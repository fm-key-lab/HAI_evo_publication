
## plot collectors curves for SNVs (and Indels)

import os
import glob
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = "Helvetica"
# plt.rcParams['text.usetex'] = True 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


lineage_dirs = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/P*'))
num_permutations = 100
min_isolates = 1 ## number of isolates (across all sites) required to evaluate timepoint
min_isolates_per_site = 1 ## number of isolates per individual site to consider timepoint
max_colonization_time = 0 ## max ime of observation across all lineages
num_isolates_required = 5 ## at least 5 isolates are required for the entire lineeage (across all timestamps) to be considered

## read in the isolate -> date table 
isolate_to_exact_date = pd.read_excel('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/2024_08_29_HAI_pilot_study_metadata.xlsx', sheet_name='Isolates')
isolate_to_exact_date = isolate_to_exact_date[isolate_to_exact_date['denovo_refgenome'].notna()]
isolate_to_exact_date['patient_pathogen'] = isolate_to_exact_date['denovo_refgenome'].str.split('-c').str[0]
isolate_to_exact_date['days_since_obs'] = isolate_to_exact_date.groupby('patient_pathogen')['hosp_days'].transform(lambda x: x - min(x))
isolate_to_exact_date.index = isolate_to_exact_date['ID']
isolate_to_exact_date_dict = isolate_to_exact_date['days_since_obs'].to_dict()

## get the data combined in one dict to allow panel plotting
cc_snv_dict = {}
cc_mean_dmrca_dict = {}
cc_sd_dmrca_dict = {}
cc_tp_to_days_dict = {}
cc_tp_is_inf_dict = {}

for lineage_dir in sorted(lineage_dirs):
    contains_snvs = False
    lineage_name = lineage_dir.split('/')[-1]
    patient = str(int(lineage_name.split('_')[0][1:]))
    genome = lineage_name.split('_')[2]
    
    analysis_dir = sorted(glob.glob(f'{lineage_dir}/2026_03_26*'))
    if len(analysis_dir) == 0:
        print(lineage_dir)
        continue
    else:
        analysis_dir = analysis_dir[0]

    ## simulation #SNPs identified
    # 100 times randomly choose x isolates (0<x<n) and assess mean / SD numbers SNVs identified 
    # identified versus nt_anc, based on single-iolate nt_anc (if . > discarded)
    # mean for 100sims, SD for 100 means
    if not os.path.exists(f'{analysis_dir}/snp_table.csv'):
        contains_snvs = False
        mrca_seq = []
        locations_long_names = ['Nasal', 'Oral', 'Rectal', 'Skin', 'Blood', 'Lung', 'Gastric', 'Urine']
        smpl_type_abbreviations = ['R', 'R', 'R', 'R', 'D', 'D', 'D', 'D']; ## routine = R, diagnostic = D
        ## check if some clonal isolates with only one tp are present
        goodsmpl_cmt = np.load(glob.glob(f'{analysis_dir}/*_goodpos_candindel_cmt.npz')[0], allow_pickle = True)
        is_true_sample = (goodsmpl_cmt['timepoints'] != 0)
        patient_cmt = goodsmpl_cmt['patients'][is_true_sample]
        timepoints_cmt = goodsmpl_cmt['timepoints'][is_true_sample].astype(str)
        goodsmpls_cmt = goodsmpl_cmt['sampleNames'][is_true_sample]
        dates_cmt = goodsmpl_cmt['kit_date'][is_true_sample]
        locations_cmt = goodsmpl_cmt['locations'][is_true_sample] 
        origplate_cmt = goodsmpl_cmt['original_plates'][is_true_sample]

        smpl_name = [f'{patient_cmt[idx]}_TP{timepoints_cmt[idx].zfill(2)}_{dates_cmt[idx]}_{smpl_type_abbreviations[locations_cmt[idx]]}_{locations_long_names[locations_cmt[idx]]}_{origplate_cmt[idx]}__{goodsmpls_cmt[idx]}' for idx in range(len(goodsmpls_cmt))]
        smpl_name = [sn.replace('___', '__') for sn in smpl_name]
        
        ## generate the officaial label to allow parsing below again
        df = pd.DataFrame(columns = np.concatenate([np.array(['chr', 'pos', 'type', 'muts', 'locustag', 'gene', 'loc1', 'loc2', 'strand', 'product', 'nt_pos', 'aa_pos', 'nt_ref', 'nt_anc', 'nt_alt']), 
                                                   smpl_name]))

    else:
        contains_snvs = True
        df=pd.read_csv(f'{analysis_dir}/snp_table.csv')
        ## read in MRCA sequence
        with open(f'{analysis_dir}/node_sequence_inference/ancestral_sequences.fasta') as fasta:
            for record in SeqIO.parse(fasta, 'fasta'):
                if record.id == 'NODE_0000001':
                    mrca_seq = [nt for nt in str(record.seq)]
    nt_ref_arr = np.array(df[['nt_ref']])
    nt_alt_arr = np.array(df[['nt_alt']])
    nt_anc_arr = np.array([[nt] for nt in mrca_seq])
    nt_der_arr = np.array([[alleles.replace(anc, '')] for alleles, anc in zip((nt_alt_arr + nt_ref_arr).flatten(), nt_anc_arr.flatten())])

    # collectors curve for SNVs
    np.random.seed(123)
    res_snps = {}
    mean_snp_num = {}
    # collectors curve for DMRCA
    res_mean_dmrca = {}
    res_sd_dmrca = {}
    
    ## store num isolates to later check if sufficient isolates are present to do a collectors curve (n ≥ 5)
    num_isolates_lineage = len([isolate_name for isolate_name in df.columns if 'TP' in isolate_name])
    
    unique_tps = set([isolate_name.split('TP')[1].split('_')[0] for isolate_name in df.columns if 'TP' in isolate_name])
    unique_tp_days = set([tuple([isolate_name.split('_')[1], ## timepoint
                                isolate_to_exact_date_dict[isolate_name.split('__')[-1].replace('_S2', '').replace('_S3', '')], ## isolate id to days since species observation
                                isolate_name.split('_')[3]]) ## diagnostic/routine info
                                for isolate_name in df.columns if 'TP' in isolate_name])
    tp_is_inf = {tp.replace('TP', ''): loc for tp, _, loc in sorted(unique_tp_days)}
    tp_to_days_dict = {tp.replace('TP', ''): int(days) for tp, days, _ in sorted(unique_tp_days)}

    if contains_snvs:
        max_colonization_time = max(max_colonization_time, max(tp_to_days_dict.values()))
    
    unique_tps = sorted(unique_tps) ## keep tps as tring but sort them increasingly!
    locations_of_interest = ['Body-wide']
    
    arr_mean_dmrca = np.array([])
    arr_sd_mean_dmrca = np.array([])

    for tp in unique_tps:
        df_tp = df.filter(regex=(f'_TP{tp}_')) 
        if len(df_tp.columns) < min_isolates: ## not sufficient amount of isolates at timepoint
            res_snps[tp] = {loc: np.array([]) for loc in locations_of_interest}
            mean_snp_num[tp] = {loc: np.nan for loc in locations_of_interest}
            res_mean_dmrca[tp] = {loc: np.array([]) for loc in locations_of_interest}
            res_sd_dmrca[tp] = {loc: np.array([]) for loc in locations_of_interest}
        else: ## sufficient data with mutations retained at timepoint
            res_snps[tp] = {}
            mean_snp_num[tp] = {}
            res_mean_dmrca[tp] = {}
            res_sd_dmrca[tp] = {}
            for loc in locations_of_interest:
                arr_num_snvs = np.array([])
                
                arr_mean_dmrca = np.array([])
                arr_sd_mean_dmrca = np.array([])

                df_tp_loc = df_tp.copy()
                num_isolates = df_tp_loc.shape[1] # shape represents already length of data --> (all isolates) 

                if num_isolates < min_isolates_per_site: ## need to evaluate! NOTE: utoff again as we also loop over individual sites!
                    res_snps[tp][loc] = arr_num_snvs
                    mean_snp_num[tp][loc] = np.nan
                    res_mean_dmrca[tp][loc] = arr_mean_dmrca
                    res_sd_dmrca[tp][loc] = arr_sd_mean_dmrca
                    continue
                elif df_tp_loc.empty: ## sufficient number of isolates given, but no snvs detected! --> all data are 0!
                    res_snps[tp][loc] = np.tile(0,(num_isolates,num_permutations)) ## get matrix of 
                    mean_snp_num[tp][loc] = np.mean(res_snps[tp][loc], axis = 1)[-1]
                    res_mean_dmrca[tp][loc] = np.tile(0,(num_isolates,))
                    res_sd_dmrca[tp][loc] = np.tile(0,(num_isolates,))
                    continue
                else:
                    for i in range(1,num_isolates+1): # loop through sets of isolates
                        isolate_iteration_dmrca_mean = np.array([])

                        for _ in range(num_permutations): # loop permutations
                            selection_isolates = np.random.choice(list(df_tp_loc.columns),size=i, replace=True)  # random selection of i isolates; w replacement
                            df_tp_loci = np.array(df_tp_loc[ selection_isolates ])  # dataframe for i isolates of all polarized sites 
                            
                            nt_anc_isolates = np.any((df_tp_loci == nt_anc_arr),axis=1) # any pos where any isolate has anc
                            nt_der_isolates = np.any(((df_tp_loci != nt_anc_arr) & (df_tp_loci != '?')),axis=1) # any pos where any isolate has alt; NOTE: nt_alt_arr can contain for multiallelic sites multiple different alleles! Therefore, we just check if the sequence in the isolates are not the anc seq --> mutation!
                            
                            nt_polymorphic_sites = np.all((nt_anc_isolates,nt_der_isolates),axis=0) # pos where alt and ref present, ie. polymorphic among selected isolates
                            num_variable_snv_found_in_subset = np.sum(nt_polymorphic_sites)
                            arr_num_snvs = np.append(arr_num_snvs,num_variable_snv_found_in_subset) # store # diffs

                            df_tp_loci_var = df_tp_loci[nt_polymorphic_sites,] # data only for variable positions in subset of isolates
                            dmrca_mean = np.mean(np.sum(df_tp_loci_var == nt_der_arr[nt_polymorphic_sites],axis=0)) # collect mean dMRCA for each iteration/isolate
                            isolate_iteration_dmrca_mean = np.append(isolate_iteration_dmrca_mean,dmrca_mean) # store mean dMRCAs
                        
                        arr_mean_dmrca = np.append(arr_mean_dmrca,np.mean(isolate_iteration_dmrca_mean)) # single value for each iteration
                        arr_sd_mean_dmrca = np.append(arr_sd_mean_dmrca,np.std(isolate_iteration_dmrca_mean)) # single value for each iteration based on each mean obs 

                res_snps[tp][loc] = arr_num_snvs.reshape((df_tp_loci.shape[1],num_permutations)) 
                mean_snp_num[tp][loc] = np.mean(res_snps[tp][loc], axis = 1)[-1] ## last calculation with all isolates
                res_mean_dmrca[tp][loc] = arr_mean_dmrca
                res_sd_dmrca[tp][loc] = arr_sd_mean_dmrca
    
    ## if at least one TP has data --> store, else skip
    if num_isolates_lineage >= num_isolates_required:
        cc_snv_dict[f'Patient {patient} {genome}'] = res_snps
        cc_mean_dmrca_dict[f'Patient {patient} {genome}'] = res_mean_dmrca
        cc_sd_dmrca_dict[f'Patient {patient} {genome}'] = res_sd_dmrca
        cc_tp_to_days_dict[f'Patient {patient} {genome}'] = tp_to_days_dict
        cc_tp_is_inf_dict[f'Patient {patient} {genome}'] = tp_is_inf
    
    else:
        print(f'Not storing data for {patient} {lineage_name} as it has just {num_isolates_lineage} but {num_isolates_required} are required (otherwise collectors curves do not make sense!)')
    ## get max number of days (required for the different timepoints)

    df_tp_loc_true = np.array(df.filter(regex=(f'_TP')))  # dataframe for i isolates of all polarized sites 
    if len(df_tp_loc_true) == 0:
        continue
    nt_anc_isolates_true = np.any((df_tp_loc_true == nt_anc_arr),axis=1) # any pos where any isolate has anc
    nt_der_isolates_true = np.any(((df_tp_loc_true != nt_anc_arr) & (df_tp_loc_true != '?')),axis=1) # any pos where any isolate has alt; NOTE: nt_alt_arr can contain for multiallelic sites multiple different alleles! Therefore, we just check if the sequence in the isolates are not the anc seq --> mutation!
    nt_polymorphic_sites_true = np.all((nt_anc_isolates_true,nt_der_isolates_true),axis=0) # pos where alt and ref present, ie. polymorphic among selected isolates
    num_variable_snv_found_in_subset_true = np.sum(nt_polymorphic_sites_true)
    df_tp_loci_var_true = df_tp_loc_true[nt_polymorphic_sites_true,] # data only for variable positions in subset of isolates
    dmrca_mean_true = np.mean(np.sum(df_tp_loci_var_true == nt_der_arr[nt_polymorphic_sites_true],axis=0)) # collect mean dMRCA for each iteration/isolate
    print(f'True dMRCA for {patient} {genome}: {dmrca_mean_true}')

time_to_color = mcolors.LinearSegmentedColormap.from_list("", ["white","black"])(np.linspace(0, 1, max_colonization_time+1)) #plt.get_cmap('Greys')(np.linspace(1, 0, 5)) ## 
color_to_days_sinc_obs = {days: c for days, c in zip(np.arange(max_colonization_time+1), time_to_color)}


color_inf = {'D': 'red',
             'R': 'k'}

myfontsize = 12
num_genomes = len(cc_snv_dict.keys())
ncols = 3
nrows = int((num_genomes-1)/ncols)+1

ordered_patients = ['Patient 7 Apittii-c1', 
                    'Patient 7 Bthetaiotaomicron-c1', 
                    'Patient 7 Ecoli-ST58-c5', 
                    'Patient 7 Ecoli-ST95-c3', 
                    'Patient 7 Ecoli-ST4260-c2', 
                    'Patient 7 Ehormaechei-c1', 
                    'Patient 7 Paeruginosa-c1', 
                    'Patient 10 Kmichiganensis-c1', 
                    'Patient 21 Ecoli-ST23-c2', 
                    'Patient 21 Ecoli-STnfx-c1', 
                    'Patient 21 Paeruginosa-c1', 
                    'Patient 21 Pmirabilis-c2', 
                    'Patient 21 Saureus-ST97-c3']

if len(ordered_patients) != len(cc_mean_dmrca_dict.keys()):
    ordered_patients = cc_mean_dmrca_dict.keys()

## COLLECTORS CURVES FOR dMRCA
fig, axes = plt.subplots(figsize=(4.25*ncols,2.5*nrows), ncols = ncols, nrows = nrows, sharex = True, sharey = False)
axes = axes.flatten()
## remove unecessary axes
for ax in axes[num_genomes:]:
    fig.delaxes(ax)
for ax_idx, patient in enumerate(ordered_patients):

    unique_tps_with_data = sorted([tp for tp, data in cc_mean_dmrca_dict[patient].items() if len(data['Body-wide'])!=0])
    color_list = ['#F8D9E2', '#C3B1E0', '#DCCDBE', '#F3E9B4', '#B2E0D9', '#e2b46a', 'darkgrey', 'lightgrey', 'white']
    color_to_days = {tp: color for tp, color in zip(unique_tps_with_data, color_list)}
    
    ## plot just the body wide plot 
    legend_labels = []
    legend_handle = []
    max_ylim = 0
    for tp in unique_tps_with_data:
        date = cc_tp_to_days_dict[patient][tp]
        loc_of_tp = cc_tp_is_inf_dict[patient][tp]

        legend_labels.append(f'Day {date} (n={len(cc_mean_dmrca_dict[patient][tp]["Body-wide"])})')
        legend_handle.append(Line2D([0],[0],marker='o',markersize=7,markerfacecolor=color_to_days[tp],markeredgecolor='k',markeredgewidth=0.4,color=color_to_days[tp],linestyle='-',linewidth=0.5)) #color_to_tp[tp],markeredgecolor='k',color=color_to_tp[tp],linestyle='-',linewidth=1))
        
        if cc_mean_dmrca_dict[patient][tp]['Body-wide'].size > 0:
            ## plot behind a black stroke to increase contrast
            axes[ax_idx].errorbar(range(1,cc_mean_dmrca_dict[patient][tp]['Body-wide'].shape[0]+1), 
                        cc_mean_dmrca_dict[patient][tp]['Body-wide'], 
                        cc_sd_dmrca_dict[patient][tp]['Body-wide']/2, 
                        linestyle='None', 
                        linewidth = 0.6,
                        color = 'k',
                        zorder = 5)
            axes[ax_idx].errorbar(range(1,cc_mean_dmrca_dict[patient][tp]['Body-wide'].shape[0]+1), 
                        cc_mean_dmrca_dict[patient][tp]['Body-wide'], 
                        cc_sd_dmrca_dict[patient][tp]['Body-wide']/2, 
                        linestyle='None', 
                        alpha = 0.7, 
                        marker = 'o', 
                        markersize = 4,
                        color = color_to_days[tp],
                        markeredgecolor = 'k',
                        markeredgewidth = 0.2,
                        zorder = 10)
    if axes[ax_idx].get_ylim()[1] <= 1:
        major_y_ticks = 0.2
    elif axes[ax_idx].get_ylim()[1] <= 3:
        major_y_ticks = 0.5
    elif axes[ax_idx].get_ylim()[1] > 3:
        major_y_ticks = 2
    if axes[ax_idx].get_xlim()[1] <= 20:
        major_x_ticks = 5
    else:
        major_x_ticks = 10
    if ax_idx>=(num_genomes-ncols):
        axes[ax_idx].xaxis.label.set_visible(True)
        axes[ax_idx].tick_params(labelbottom=True)
        axes[ax_idx].set_xlabel('Number of isolates')
    if ax_idx%ncols == 0:        
        axes[ax_idx].set_ylabel('Mean dMRCA',fontsize=myfontsize)
    axes[ax_idx].set_title(f'{patient.split('-c')[0]}',fontsize=myfontsize)
    axes[ax_idx].spines['right'].set_visible(False)
    axes[ax_idx].spines['top'].set_visible(False)
    axes[ax_idx].tick_params(axis='both', which='major', labelsize=myfontsize)
    axes[ax_idx].grid(which='major', axis = 'y', linestyle = '--', linewidth = 0.5, c = 'k', zorder = 0)
    axes[ax_idx].set_ylim(0, None)
    num_legend_cols = 1
    axes[ax_idx].legend(labels = legend_labels, handles = legend_handle, title = 'Days since\nfirst observation',loc='center left', bbox_to_anchor = (1, 0.5), ncols = num_legend_cols)
    ## add to all bottom axes the xlabel
    if ax_idx>=(num_genomes-ncols):
        axes[ax_idx].xaxis.label.set_visible(True)
        axes[ax_idx].tick_params(labelbottom=True)
        axes[ax_idx].set_xlabel('Number of Isolates')
    plt.tight_layout()
    fig.savefig(f'{analysis_dir}/pdf/dMRCA_SNV_BodyWide_allTPw{min_isolates}iso_{num_permutations}sim_wReplacement.pdf',dpi=300, bbox_inches = 'tight')


### COLLECTORS CURVES FOR SNVS
fig, axes = plt.subplots(figsize=(4.25*ncols,2.5*nrows), ncols = ncols, nrows = nrows, sharex = True, sharey = False)
axes = axes.flatten()
for ax in axes[num_genomes:]:
    fig.delaxes(ax)
for ax_idx, patient in enumerate(ordered_patients):

    unique_tps_with_data = [tp for tp, data in cc_snv_dict[patient].items() if len(data['Body-wide'])!=0]
    color_to_days = {tp: color for tp, color in zip(unique_tps_with_data, color_list)}
    
    ## plot just the body wide plot 
    legend_labels = []
    legend_handle = []
    max_ylim = 0
    for tp in unique_tps_with_data:
        date = cc_tp_to_days_dict[patient][tp]
        loc_of_tp = cc_tp_is_inf_dict[patient][tp]

        legend_labels.append(f'Day {date} (n={len(cc_mean_dmrca_dict[patient][tp]["Body-wide"])})')
        legend_handle.append(Line2D([0],[0],marker='o',markersize=7,markerfacecolor=color_to_days[tp],markeredgecolor='k',markeredgewidth=0.4,color=color_to_days[tp],linestyle='-',linewidth=0.5)) #color_to_tp[tp],markeredgecolor='k',color=color_to_tp[tp],linestyle='-',linewidth=1))
        
        if cc_mean_dmrca_dict[patient][tp]['Body-wide'].size > 0:
            axes[ax_idx].errorbar(range(1,cc_snv_dict[patient][tp]['Body-wide'].shape[0]+1), 
                        np.mean(cc_snv_dict[patient][tp]['Body-wide'],axis=1), 
                        np.std(cc_snv_dict[patient][tp]['Body-wide'],axis=1)/2, 
                        linestyle='None', 
                        linewidth = 0.6,
                        color = 'k',
                        zorder = 5)
            axes[ax_idx].errorbar(range(1,cc_snv_dict[patient][tp]['Body-wide'].shape[0]+1), 
                        np.mean(cc_snv_dict[patient][tp]['Body-wide'],axis=1), 
                        np.std(cc_snv_dict[patient][tp]['Body-wide'],axis=1)/2, 
                        linestyle='None', 
                        alpha = 0.7, 
                        marker = 'o', 
                        markersize = 4,
                        color = color_to_days[tp],
                        markeredgecolor = 'k',
                        markeredgewidth = 0.2,
                        zorder = 10)
    if axes[ax_idx].get_ylim()[1] <= 1:
        major_y_ticks = 0.2
    elif axes[ax_idx].get_ylim()[1] <= 3:
        major_y_ticks = 0.5
    elif axes[ax_idx].get_ylim()[1] > 3:
        major_y_ticks = 2
    if axes[ax_idx].get_xlim()[1] <= 20:
        major_x_ticks = 5
    else:
        major_x_ticks = 10
    if ax_idx>=(num_genomes-ncols):
        axes[ax_idx].xaxis.label.set_visible(True)
        axes[ax_idx].tick_params(labelbottom=True)
        axes[ax_idx].set_xlabel('Number of isolates')
    if ax_idx%ncols == 0:        
        axes[ax_idx].set_ylabel('SNVs identified',fontsize=myfontsize)
    axes[ax_idx].set_title(f'{patient.split('-c')[0]}',fontsize=myfontsize)
    axes[ax_idx].spines['right'].set_visible(False)
    axes[ax_idx].spines['top'].set_visible(False)
    axes[ax_idx].tick_params(axis='both', which='major', labelsize=myfontsize)
    axes[ax_idx].grid(which='major', axis = 'y', linestyle = '--', linewidth = 0.5, c = 'k', zorder = 0)
    axes[ax_idx].set_ylim(0, None)
    num_legend_cols = 1
    axes[ax_idx].legend(labels = legend_labels, handles = legend_handle, title = 'Days since\nfirst observation',loc='center left', bbox_to_anchor = (1, 0.5), ncols = num_legend_cols)
    ## add to all bottom axes the xlabel
    if ax_idx>=(num_genomes-ncols):
        axes[ax_idx].xaxis.label.set_visible(True)
        axes[ax_idx].tick_params(labelbottom=True)
        axes[ax_idx].set_xlabel('Number of Isolates')
    plt.tight_layout()
    fig.savefig(f'{analysis_dir}/pdf/SNV_BodyWide_allTPw{min_isolates}iso_{num_permutations}sim_wReplacement.pdf',dpi=300, bbox_inches = 'tight')