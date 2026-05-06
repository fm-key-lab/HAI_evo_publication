

## plot the mean coverage and mark SNVs, indels and deletions for isolates used for phenotyping

import gzip
import glob
import pickle
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from Bio import SeqIO
import sys

plt.rcParams['font.family'] = "Helvetica"
# plt.rcParams['text.usetex'] = True 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

main_workingdir = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/')
cand_mut_table_dirs = 'case/2-candidate_mutation_table_2024_10_18/'
Ehorm_refgenome = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/P07_Ehormaechei-c1_240825/')
samples_phenotyping = {'wt\t': ['P07_2048', 'P07_2055'], 
                       'fimZ[F126L]': ['P07_2052', 'P07_2237', 'L00071_01_1', 'L00071_02_1_S3']}

chunksize = 1000
max_cov_cutoff = 10 ## maximum coverage (normalized) to be plotted, all values above will be cut to max_cov_cutoff
extract_maf_for_ambigs = True
isolate_list = list([entry for l in samples_phenotyping.values() for entry in l])
nts = np.array(['A','T','C','G'])

## from tab20 color scheme
snv_color_dict = {'A': '#d62728',
                  'T': '#9467bd',
                  'G': '#ff7f0e',
                  'C': '#17becf',
                  'N': '#c7c7c7'} ## undefined base

rename_samples = {}

def mean_over_chunks(arr, genome_length, chunk_size, max_cov_cutoff = False, normalize = False):
    remainder = genome_length % chunk_size
    if remainder != 0:
        # If not, pad the array with NaNs or another value
        arr = np.concatenate([arr, np.full((1, chunk_size - remainder), np.nan)], axis = 1)
    reshaped_arr = arr.reshape(-1, chunk_size)
    means = np.nanmean(reshaped_arr, axis=1)
    if normalize:
        means = means/normalize
    if max_cov_cutoff:
        means[means > max_cov_cutoff] = max_cov_cutoff
    return means

def genomestats(REFGENOMEFOLDER):
    # parse ref genome to extract relevant stats
    ref_genome_path=str(REFGENOMEFOLDER)+'/genome.fasta'
    fasta_file = glob.glob(ref_genome_path)
    if len(fasta_file) != 1:
        ref_genome_path=str(REFGENOMEFOLDER)+'/genome.fasta.gz'
        fasta_file_gz = glob.glob(ref_genome_path)
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    Genomelength = 0
    ChrStarts = []
    ScafNames = []
    for record in refgenome:
        ChrStarts.append(Genomelength) # chr1 starts at 0 in analysis.m
        Genomelength = Genomelength + len(record)
        ScafNames.append(record.id)
    # turn to np.arrys!
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return [ChrStarts,Genomelength,ScafNames]


## function from case module (version 24/08/02)
def chrpos2index(chrpos,chr_starts,genomelength=np.nan):
    if (np.size(chrpos,1) != 2) & (np.size(chrpos,0) == 2):
        chrpos=chrpos.T
        print('Reversed orientation of chrpos')
    elif (np.size(chrpos,1) != 2) & (np.size(chrpos,0) != 2):
        print('Wrong input format given. Provide np.array with shape (p, 2)')
        return 
    
    if len(chr_starts) == 1:
        p=chrpos[:,1]
    else:
        p=chr_starts[chrpos[:,0]-1]+chrpos[:,1]
    if not np.isnan(genomelength) and any(p>=genomelength):
        print('Invalid genome positions given with positions larger than genome length')
        return sys.exit()
    return p

################
sample_stats = {}
identified_snvs_list = []
with gzip.open(f'{main_workingdir}{cand_mut_table_dirs}/P07_Ehormaechei-c1_240825/candidate_mutation_table.pickle.gz', 'rb') as f:
    cmt = pickle.load(f)
    sampleNames = np.array(cmt[0])
    p = np.array(cmt[1])
    calls = np.array(cmt[2])
    del cmt

cov_mat_csr = scipy.sparse.load_npz(f'{main_workingdir}{cand_mut_table_dirs}/P07_Ehormaechei-c1_240825/cov_raw_sparsecsr_mat.npz')
cov_matrix = np.array(cov_mat_csr.todense())
del cov_mat_csr

[ChrStarts,genome_length,_] = genomestats(Ehorm_refgenome)

## loop over every resequenced sample
result_cov_dict = {}

for genotype, isolate_l in samples_phenotyping.items():
    for isolateid in isolate_l:
        if not isolateid in rename_samples.keys():
            sample_idx = np.where([isolateid == sample for sample in sampleNames])[0]
        else:
            sample_idx = np.where([rename_samples[isolateid] == sample for sample in sampleNames])[0]
        samplename = sampleNames[sample_idx][0]
        
        ## subset cov matrix to samples of interest 
        sample_mean = np.mean(cov_matrix[sample_idx, :])
        sample_1kb_means = mean_over_chunks(cov_matrix[sample_idx, :], genome_length, chunksize, max_cov_cutoff, sample_mean)
        result_cov_dict[(genotype, isolateid)] = [sample_1kb_means, sample_mean]


## read snv table, indel table and mge table
snp_table = pd.read_csv(f'{main_workingdir}analysis/P0007_P07_Ehormaechei-c1_240825/2026_03_26_P07_Ehormaechei-c1_240825/snp_table.csv')
snp_table_metadata_col = ['chr', 'pos', 'type', 'muts', 'locustag', 'gene', 'loc1', 'loc2', 'strand', 'product', 'nt_pos', 'aa_pos', 'nt_ref', 'nt_anc', 'nt_alt']
snp_table.columns = [col.split('__')[-1] if col not in snp_table_metadata_col else col for col in snp_table.columns ]
snp_table = snp_table[snp_table_metadata_col + isolate_list]
snp_table['genome_pos'] = chrpos2index(np.array(snp_table[['chr', 'pos']].values), ChrStarts) 
## for all ambiguous nucleotides, load the major call from the calls df
if extract_maf_for_ambigs:
    p_oi_bool = [pos in (snp_table['genome_pos'].values-1) for pos in p] ## turn to 0-based
    for isolate in isolate_list:
        isolate_idx = np.where(isolate == sampleNames)[0][0]
        if any(snp_table[isolate] == '?'):
            allele_cov = calls[isolate_idx, :4, p_oi_bool] + calls[isolate_idx, 4:, p_oi_bool]
            maf_call = nts[np.argmax(allele_cov, axis = 1)]
            hasambig_call = np.where(snp_table[isolate] == '?')[0]
            snp_table.loc[hasambig_call, isolate] = maf_call[hasambig_call]
            ## for logging print stats
            cov_of_ambig_maf = allele_cov[np.arange(allele_cov.shape[0]), np.argmax(allele_cov, axis = 1)]
            lowest_maf_ambigs = np.min(cov_of_ambig_maf / np.sum(allele_cov, axis = 1))
            lowest_cov_maf = np.min(cov_of_ambig_maf)
            print(f'{isolate}: The lowest major allele frequency was {lowest_maf_ambigs} and the lowest coverage on any position was {lowest_cov_maf}')
else:
    for isolate in isolate_list:
        snp_table.loc[snp_table[isolate] == '?', isolate] = 'N'
snp_table['genome_pos'] = snp_table['genome_pos'] / chunksize
## keep only variable positions
var_pos_bool = snp_table[isolate_list].apply(lambda row: row[row != '?'].nunique(), axis=1) >= 2 ## exclude ambiguous calls
snp_table = snp_table[var_pos_bool]


indel_table = pd.read_csv(f'{main_workingdir}analysis/P0007_P07_Ehormaechei-c1_240825/2026_03_26_P07_Ehormaechei-c1_240825/indel_table.csv')
indel_table_metadata_col = ['gene_num', 'gene_num_global', 'chr', 'pos', 'chr_boundary', 'indel_ref', 'indel_anc', 'indel_alt', 'max_indel_GL_diff_to_ref', 'indel_size_gp', 'type', 'product', 'gene', 'protein_id', 'strand', 'loc1', 'loc2', 'sequence', 'note', 'locustag', 'orthologtag', 'translation', 'indel_pos_start', 'indel_pos_end', 'alt_start', 'alt_stop', 'alt_stop_pos', 'frameshift', 'mut_translation', 'mut', 'num_isolates']
indel_table = indel_table[indel_table_metadata_col + isolate_list]
indel_table['genome_pos'] = chrpos2index(np.array(indel_table[['chr', 'pos']].values), ChrStarts)  / chunksize ## check if -1 need to be applied
## keep only variable positions
indel_table = indel_table[indel_table[isolate_list].nunique(axis=1) >= 2]

mge_table = pd.read_csv(f'{main_workingdir}analysis/P0007_P07_Ehormaechei-c1_240825/2026_03_26_P07_Ehormaechei-c1_240825/mge_analysis_ml_tree_based/mobile_genetic_elements_cov_breadth.csv')
mge_table_metadata_col = ['TYPE', 'START', 'END', 'LENGTH']
mge_table = mge_table.loc[mge_table['TYPE'] == 0, mge_table_metadata_col + isolate_list]
mge_table[['START', 'END', 'LENGTH']] = mge_table[['START', 'END', 'LENGTH']] / chunksize
        
####################
## PLOT
max_cov = max([max(arr[0]) for arr in result_cov_dict.values()])

fig, axs = plt.subplots(figsize = (6, 0.65*len(result_cov_dict)), nrows = len(result_cov_dict), sharex = True, sharey = True, gridspec_kw={'hspace': 0.15})
axs = axs.flatten()
for isolate_idx, (genotype_isolate, cov_1kb_mean) in enumerate(result_cov_dict.items()):
    genotype, isolate = genotype_isolate
    axs[isolate_idx].plot(np.arange(0, len(cov_1kb_mean[0])/1000, 0.001),cov_1kb_mean[0], color='k',linewidth=.4, zorder = 10)
    axs[isolate_idx].text(0.01, 0.95, f'{genotype}\t({isolate})', fontsize=8,horizontalalignment='left',verticalalignment='top',transform = axs[isolate_idx].transAxes)

    if not snp_table.empty:
        [axs[isolate_idx].scatter(x=snp_pos/1000, y = cov_1kb_mean[0][int(snp_pos)], color = snv_color_dict[snp_type], marker = 'o', s = 50, linewidth = 0.5, edgecolor = 'k', zorder =20) for snp_pos, snp_type in snp_table[['genome_pos', isolate]].values]
        ## mark SNVs
        [axs[isolate_idx].axvline(x=snp_pos/1000, ls = '-', linewidth = 1, color = snv_color_dict[snp_type], alpha = 0.8, zorder = 1) for snp_pos, snp_type in snp_table[['genome_pos', isolate]].values]
        [axs[isolate_idx].annotate(text = f'{gene}', xy=(snp_pos/1000, max_cov), va = 'top', ha = 'right', rotation = 90, size = 8, zorder = 25) for snp_pos, snp_type, gene in snp_table[['genome_pos', isolate, 'gene']].values]
    if not indel_table.empty:
        ## mark Indels
        [axs[isolate_idx].scatter(x=indel_pos/1000, y = cov_1kb_mean[0][indel_pos], color = 'grey', marker = 'o', s = 50, linewidth = 0.5, edgecolor = 'k', zorder =20) for indel_pos, indel_type in indel_table[['genome_pos', isolate]].values]
        [axs[isolate_idx].axvline(x=indel_pos/1000, ls = '-', linewidth = 1, color = 'grey', alpha = 0.8, zorder = 1) for indel_pos, indel_type in indel_table[['genome_pos', isolate]].values]
        [axs[isolate_idx].annotate(text = f'{gene}', xy=(indel_pos/1000, max_cov), va = 'top', ha = 'right', rotation = 90, size = 8, zorder = 25) for indel_pos, indel_type, gene in indel_table[['genome_pos', isolate, 'gene']].values]
    if not mge_table.empty:
        ## mark MGEs
        [axs[isolate_idx].axvspan(start/1000, end/1000, color='red', linewidth = 0, alpha=0.2) for start, end, cov_breadth in mge_table[['START', 'END', isolate]].values if (cov_breadth <= 0.5)]
## add grid lines
[axs[i].xaxis.set_major_locator(ticker.MultipleLocator(0.5)) for i in range(len(axs))]
[axs[i].xaxis.set_minor_locator(ticker.MultipleLocator(1)) for i in range(len(axs))]
[axs[i].yaxis.set_major_locator(ticker.MultipleLocator(4)) for i in range(len(axs))]
[axs[i].yaxis.set_minor_locator(ticker.MultipleLocator(2)) for i in range(len(axs))]
[axs[i].grid(which='major', axis = 'y', linestyle='--', linewidth=0.4, color='black', zorder = 1) for i in range(len(axs))]
## mark contig edges
[axs[i].axvline(x=x/chunksize/1000, ls = '--', linewidth = 0.5, color = 'darkgreen', alpha = 0.2, zorder = 1) for i in range(len(axs)) for x in ChrStarts]
## set labels
axs[len(axs)-1].set_xlabel('Genome position [Mb]')
[axs[i].set_ylabel('Nomalized\ncoverage') for i in range(len(axs))]
axs[len(axs)-1].set_xlim(-1/1000, (genome_length/chunksize+1)/1000)
## generate homemade legend
legend_handle = [Line2D([], [], marker='o', linestyle='None', markersize=8, markerfacecolor=color, markeredgecolor='k', label=nt) for nt, color in snv_color_dict.items()]
legend_handle += [Patch(facecolor='red', edgecolor=None, label='Deletion', alpha = 0.2)]
axs[2].legend(title = 'Allele', handles = legend_handle, loc='center left', bbox_to_anchor = (1, 0))
plt.tight_layout()
fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/fimZ_phenotyping_genome_comparisons.pdf'))
