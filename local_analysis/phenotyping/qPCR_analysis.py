
########################
## Load modules 

import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np 
import scipy.stats as stats
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest
import itertools

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

########################
## If you want to save plots new, set to true
save_plt = False

########################
## Load all datasets 
## get all files starting with pattern
files = glob.glob(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/RT_qPCR/2025_05_qPCR/2025_06*Ehorm*rpoDfimZfimA_P*.csv'))
files = [f for f in files if not f.endswith('_raw.csv')]
flagella_files = glob.glob(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/qPCR/raw/2023_08_21*'))
flagella_files.extend(glob.glob(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/qPCR/raw/2023_08_24*')))

samples_to_ids_dict = {'1': 'P07_2048',
                       '2': 'P07_2055',
                       '4': 'P07_2052',
                       '5': 'P07_2237',
                       '6': 'L00071_01_1',
                       '8': 'L00071_02_1'}
map_muts = {'P07_2048': 'wt',
            'P07_2055': 'wt',
            'P07_2052': 'FimZ[F126L]',
            'P07_2237': 'FimZ[F126L]',
            'L00071_01_1': 'FimZ[F126L]',
            'L00071_02_1': 'FimZ[F126L]'}

primer_ordered = ['rpoD', 'fimZ', 'fimA', 'flhC', 'fliC']

## generate color map
genotype_col_dict = {'wt': '#c8c8c8', 
                    'FimZ[F126L]': '#f9af92'}

remove_high_gDNA_contaminations = True ## if samples with a Ct difference < 5 between sample and RTNC should be removed (True) or not (False)

## read in all files and concatenate them 
all_data_fimbriae = []
for file in files:
    if 'Primerefficiency' in file:
        continue ## do not analyse PE here!
    date = file.split('/')[-1][:10] ## get the date of each run to separate individual replicates
    plate = file.split('/')[-1].split('_')[-1].replace('.csv', '') ## get the date of each run to separate individual replicates
    growth_phase = file.split('/')[-1].split('_')[-3] ## get the date of each run to separate individual replicates
    tmp_df = pd.read_csv(file, sep = ';', decimal = ',')
    tmp_df['date'] = date
    tmp_df['plate'] = plate
    tmp_df['growth_phase'] = growth_phase
    tmp_df['qPCR_run'] = 'fimbriae'
    all_data_fimbriae.append(tmp_df)
all_data_df = pd.concat(all_data_fimbriae)

## map the genes
all_data_df['col'] = all_data_df['Well'].str[1:].astype(int)
all_data_df['row'] = all_data_df['Well'].str[0]
is_fimA_p1 = (all_data_df['plate'] == 'P1') & (all_data_df['row'].isin(['E', 'F']) | all_data_df['Well'].isin(['G1', 'G2', 'G3']))
is_fimZ_p1 = (all_data_df['plate'] == 'P1') & (all_data_df['row'].isin(['C', 'D']) | all_data_df['Well'].isin(['E1', 'E2']))
is_rpoD_p1 = (all_data_df['plate'] == 'P1') & (all_data_df['row'].isin(['A', 'B']) | all_data_df['Well'].isin(['B1']))
is_fimA_p2 = (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['C', 'G']) | all_data_df['Well'].isin(['D1', 'D2', 'D3', 'H1', 'H2', 'H3']))
is_fimZ_p2 = (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['B', 'F']) | all_data_df['Well'].isin(['C1', 'C2', 'G1', 'G2']))
is_rpoD_p2 = (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['A', 'E']) | all_data_df['Well'].isin(['B1', 'F1']))
all_data_df.loc[is_fimA_p1 | is_fimA_p2, 'Gene'] = 'fimA'
all_data_df.loc[is_fimZ_p1 | is_fimZ_p2, 'Gene'] = 'fimZ'
all_data_df.loc[is_rpoD_p1 | is_rpoD_p2, 'Gene'] = 'rpoD'

## map growth phase
is_early =  (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['A', 'B', 'C', 'D']))
is_late =  (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['E', 'F', 'G', 'H']))
all_data_df.loc[is_early & (all_data_df['growth_phase'] == 'earlyAndlatelog'), 'growth_phase'] = 'earlylog'
all_data_df.loc[is_late & (all_data_df['growth_phase'] == 'earlyAndlatelog'), 'growth_phase'] = 'latelog'

all_data_df = all_data_df[all_data_df['Exclusion'].isna()]

#################################
## read in flagella expression
all_data_flagella = []
for file in flagella_files:
    date = file.split('/')[-1][:10] ## get the date of each run to separate individual replicates
    tmp_df = pd.read_excel(file, sheet_name='Results', skiprows = 7)
    tmp_df['date'] = date
    if '-2h30-' in file: 
        tmp_df['growth_phase'] = 'earlylog'
    else:
        tmp_df['growth_phase'] = 'latelog'
    tmp_df = tmp_df[tmp_df['Target Name'].isin(['rpoD', 'flhC', 'fliC'])]
    tmp_df['qPCR_run'] = 'flagella'
    all_data_flagella.append(tmp_df) 

flagella_df = pd.concat(all_data_flagella)

## fix replicate assignment for datasets
flagella_df.loc[flagella_df['Well'].str.contains('A'), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('A'), 'Sample Name'] + '1' ## grep row A and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('B'), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('B'), 'Sample Name'] + '2' ## grep row B and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('C'), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('C'), 'Sample Name'] + '3' ## grep row C and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('D'), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('D'), 'Sample Name'] + '1' ## grep row D and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('E'), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('E'), 'Sample Name'] + '2' ## grep row E and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('F') & (flagella_df['Well'].str[1:].astype(int)%2==1), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('F') & (flagella_df['Well'].str[1:].astype(int)%2==1), 'Sample Name'] + '1' ## grep of row F all uneven columns and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('F') & (flagella_df['Well'].str[1:].astype(int)%2==0), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('F') & (flagella_df['Well'].str[1:].astype(int)%2==0), 'Sample Name'] + '2' ## grep of row F all even columns and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('G') & (flagella_df['Well'].str[1:].astype(int)%2==1), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('G') & (flagella_df['Well'].str[1:].astype(int)%2==1), 'Sample Name'] + '3' ## grep of row F all uneven columns and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('G') & (flagella_df['Well'].str[1:].astype(int)%2==0), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('G') & (flagella_df['Well'].str[1:].astype(int)%2==0), 'Sample Name'] + '1' ## grep of row F all even columns and expand with replicate id
flagella_df.loc[flagella_df['Well'].str.contains('H') & (flagella_df['Well'].str[1:].astype(int)%2==0), 'Sample Name'] = flagella_df.loc[flagella_df['Well'].str.contains('H') & (flagella_df['Well'].str[1:].astype(int)%2==0), 'Sample Name'] + '2' ## grep of row F all even columns and expand with replicate id
## match the same nameing scheme
flagella_df['Sample Name'] = flagella_df['Sample Name'].str.replace('Sample ', '').str.replace('Non template control', 'NTC')
flagella_df.loc[flagella_df['Sample Name'].str.contains('RTNC '), 'Sample Name'] = flagella_df.loc[flagella_df['Sample Name'].str.contains('RTNC '), 'Sample Name'].str.replace('RTNC ', '') + '-'

## rename the columns to allow merging of df
flagella_df = flagella_df[['Well', 'Sample Name', 'Target Name', 'Cт', 'date', 'growth_phase', 'qPCR_run']]
flagella_df.columns = ['Well', 'Sample', 'Gene', 'Ct', 'date', 'growth_phase', 'qPCR_run']

#################################

all_data_df = pd.concat([all_data_df, flagella_df]).reset_index(drop = True)

## set undetermined ct values to 40 to allow calculations on them!
all_data_df.loc[all_data_df['Ct'].isna() | (all_data_df['Ct'] == 'Undetermined'), 'Ct'] = 40
all_data_df['Ct'] = all_data_df['Ct'].astype(float)

## annotate data
all_data_df['sampleid'] = all_data_df['Sample'].str[0].map(samples_to_ids_dict)
all_data_df['replicate'] = all_data_df['Sample'].str[1]
all_data_df['sample_rep'] = all_data_df['Sample'].str[0] + '.' + all_data_df['replicate']
all_data_df.loc[all_data_df['replicate'] == 'T', 'replicate'] = ''
all_data_df['sample_type'] = np.where(all_data_df['Sample'].str.contains('-'), 'RTNC', 
                                      np.where(all_data_df['Sample'] == 'NTC', 'NTC', all_data_df['sample_rep']))

########################
## Primer efficiency
########################
efficiency_dict = {'fimZ': 103.12860110426207, 'rpoD': 106.90652565606507, 'fimA': 104.12957565472598}
r_val_dict = {'fimZ': 0.992816941116989, 'rpoD': 0.9882236015156362, 'fimA': 0.9970334188582403}
efficiency_dict_flagella = {'rpoD': 110.70544239077758, 'flhC': 105.35461938945457, 'fliC': 101.73230817989239}
r_val_dict_flagella = {'rpoD': 0.9941160965558062, 'flhC': 0.9953689773990124, 'fliC': 0.9882533183301636}

############################################################
############################################################
############################################################
## analyse data 
df_samples = all_data_df[~all_data_df['Sample'].isna()]

## extract different sample types from df
df_smpl = df_samples.loc[~df_samples['sample_type'].isin(['RTNC', 'NTC']), ['sampleid', 'sample_rep', 'Gene', 'Ct', 'growth_phase', 'qPCR_run']]

df_ctrs = df_samples[df_samples['sample_type'].isin(['RTNC', 'NTC'])]
df_ntc = df_ctrs.loc[df_ctrs['sample_type'] == 'NTC', ['Gene', 'Ct', 'growth_phase', 'qPCR_run']]
df_rtnc = df_ctrs.loc[df_ctrs['sample_type'] == 'RTNC', ['sampleid', 'sample_rep', 'Gene', 'Ct', 'growth_phase', 'qPCR_run']]
df_rtnc['sample_rep'] = df_rtnc['sample_rep'].str.replace('-', '')

## compare RTNC against NTC
df_ctr_diff = pd.merge(df_rtnc, df_ntc, how = 'left', on = ['Gene', 'growth_phase', 'qPCR_run'], suffixes=['_rtnc', '_ntc'])
df_ctr_diff['Ct_diff'] = df_ctr_diff['Ct_ntc'] - df_ctr_diff['Ct_rtnc']

## compare Samples against RTNC
df_rtnc_smpl_diff = pd.merge(df_smpl, df_rtnc, how = 'left', on = ['sampleid', 'sample_rep', 'Gene', 'growth_phase', 'qPCR_run'], suffixes=['_smpl', '_rtnc'])
df_rtnc_smpl_diff['Ct_diff'] = df_rtnc_smpl_diff['Ct_rtnc'] - df_rtnc_smpl_diff['Ct_smpl']

## remove all samples if their houskeeping gene ct difference between sample and RT- control is >5 (--> gDNA contamination is <3%!)
## cutoff based on: https://pmc.ncbi.nlm.nih.gov/articles/PMC3326333/
df_rtnc_smpl_diff_low_gDNA = df_rtnc_smpl_diff.loc[(df_rtnc_smpl_diff['Gene'] == 'rpoD') & (df_rtnc_smpl_diff['Ct_diff'] >= 5), ['sampleid', 'sample_rep', 'growth_phase', 'Ct_diff', 'qPCR_run']].drop_duplicates()
df_rtnc_smpl_diff_filt = pd.merge(df_rtnc_smpl_diff, df_rtnc_smpl_diff_low_gDNA[['sampleid', 'sample_rep', 'growth_phase', 'qPCR_run']], how = 'right', on = ['sampleid', 'sample_rep', 'growth_phase', 'qPCR_run'])
df_rtnc_smpl_diff_high_gDNA = df_rtnc_smpl_diff.loc[(df_rtnc_smpl_diff['Gene'] == 'rpoD') & (df_rtnc_smpl_diff['Ct_diff'] < 5), ['sampleid', 'sample_rep', 'growth_phase', 'Ct_diff', 'qPCR_run']].drop_duplicates()
df_rtnc_smpl_diff_rmvd = pd.merge(df_rtnc_smpl_diff, df_rtnc_smpl_diff_high_gDNA[['sampleid', 'sample_rep', 'growth_phase', 'qPCR_run']], how = 'right', on = ['sampleid', 'sample_rep', 'growth_phase', 'qPCR_run'])

print(f"Removed {len(df_rtnc_smpl_diff_high_gDNA)}/{sum(df_rtnc_smpl_diff['Gene'] == 'rpoD')} samples due to too high gDNA load:\n{df_rtnc_smpl_diff_high_gDNA.sort_values(['growth_phase', 'sample_rep', 'Ct_diff']).to_string(index = False)}")

## normalize values by removing gDNA background contamination 
if remove_high_gDNA_contaminations:
    df_rtnc_smpl_diff_corrected = df_rtnc_smpl_diff_filt.copy()
else:
    df_rtnc_smpl_diff_corrected = df_rtnc_smpl_diff.copy()

df_rtnc_smpl_diff_corrected['Ct_smpl'] = df_rtnc_smpl_diff_corrected.apply(lambda x: -np.log2( 2**( -x['Ct_smpl'] ) - 2**( -x['Ct_rtnc'] ) ), axis = 1)

########
## calculate delta ct values
dct_smpls = df_rtnc_smpl_diff_corrected.copy()
dct_smpls.loc[dct_smpls['qPCR_run'] == 'fimbriae', 'efficiency'] = dct_smpls.loc[dct_smpls['qPCR_run'] == 'fimbriae', 'Gene'].map(efficiency_dict) / 100 ## map calculated efficiencies and return to ratios
dct_smpls.loc[dct_smpls['qPCR_run'] == 'flagella', 'efficiency'] = dct_smpls.loc[dct_smpls['qPCR_run'] == 'flagella', 'Gene'].map(efficiency_dict_flagella) / 100 ## map calculated efficiencies and return to ratios

## calculate mean ct value per biological replicate and then per sample id
dct_smpls_mean_ct_smpl = df_rtnc_smpl_diff_corrected.groupby(['sample_rep', 'sampleid', 'Gene', 'growth_phase', 'qPCR_run'])['Ct_smpl'].agg('mean').reset_index()
## get mean ct value per primer and growth stage for wildtype
ref_smpls_mean_ct_smpl = dct_smpls_mean_ct_smpl.loc[dct_smpls_mean_ct_smpl['sampleid'].isin(['P07_2048', 'P07_2055']), ['Gene', 'growth_phase', 'Ct_smpl', 'qPCR_run']].groupby(['Gene', 'growth_phase', 'qPCR_run']).mean().reset_index()

dct_smpls = pd.merge(dct_smpls, ref_smpls_mean_ct_smpl, how = 'left', on = ['Gene', 'growth_phase', 'qPCR_run'], suffixes = ['', '_ref'])

dct_smpls['dct_diff'] = dct_smpls['Ct_smpl_ref'] - dct_smpls['Ct_smpl']

housekeeping_rpoD = dct_smpls[dct_smpls['Gene'] == 'rpoD'].copy()
housekeeping_rpoD.drop(['Gene'], axis = 1, inplace = True)
housekeeping_rpoD.drop_duplicates(inplace = True)

## merge data back to df
ddct_smpls = pd.merge(dct_smpls, housekeeping_rpoD, how = 'left', on = ['sample_rep', 'growth_phase','qPCR_run'], suffixes=['', '_hk_rpoD'])

## calculate fold-change per isolate and biological replicate
ddct_smpls['fc_rpoD_raw'] = ( (2*ddct_smpls['efficiency']) ** (ddct_smpls['dct_diff']) ) / ( (2*ddct_smpls['efficiency_hk_rpoD']) ** (ddct_smpls['dct_diff_hk_rpoD']) )
ddct_smpls['l2fc_rpoD_raw'] = np.log2(ddct_smpls['fc_rpoD_raw'])

## get mean per isolate and condition
ddct_smpls_grouped = ddct_smpls.groupby(['sample_rep', 'sampleid', 'Gene', 'growth_phase', 'qPCR_run'], as_index = False)[['Ct_smpl', 'Ct_rtnc', 'Ct_diff', 'efficiency', 'Ct_smpl_ref', 'dct_diff', 'Ct_smpl_hk_rpoD', 'Ct_rtnc_hk_rpoD', 'Ct_diff_hk_rpoD', 'efficiency_hk_rpoD', 'Ct_smpl_ref_hk_rpoD', 'dct_diff_hk_rpoD', 'fc_rpoD_raw']].agg('mean')
ddct_smpls_grouped['l2fc_rpoD_raw'] = np.log2(ddct_smpls_grouped['fc_rpoD_raw'])

## map different sample collections together
ddct_smpls['mut_group'] = ddct_smpls['sampleid'].map(map_muts)
ddct_smpls = ddct_smpls.sort_values('sampleid').reset_index(drop = True)

ddct_smpls_grouped['mut_group'] = ddct_smpls_grouped['sampleid'].map(map_muts)
ddct_smpls_grouped = ddct_smpls_grouped.sort_values('sampleid').reset_index(drop = True)

####################################
## group over all growth phases and locations!!
####################################
mannwhitneyu_long_name = 'mannwhitneyu test'
mannwhitneyu_short_name = ''
mannwhitneyu_func = stats.mannwhitneyu
mannwhitneyu_test = StatTest(mannwhitneyu_func, mannwhitneyu_long_name, mannwhitneyu_short_name)

sorted_order = ['rpoD', 'fimZ', 'fimA', 'flhC', 'fliC']
sorted_order_fim = ['rpoD', 'fimZ', 'fimA']
mutations = ['wt', 'FimZ[F126L]']

mutation_pairs = []
for gene in sorted_order:
    mutation_pairs.append(((gene, mutations[0]), (gene, mutations[1])))
mutation_pairs_fim = []
for gene in sorted_order_fim:
    mutation_pairs_fim.append(((gene, mutations[0]), (gene, mutations[1])))

boxplt_df = ddct_smpls_grouped.copy()
boxplt_df['x_axis'] = boxplt_df['Gene']

boxplt_df['x_axis'] = pd.Categorical(boxplt_df['x_axis'], sorted_order, ordered = True)
boxplt_df['muts'] = np.where(boxplt_df['mut_group'] == 'wt', 'wt', 'FimZ[F126L]')
boxplt_df['mut_gs'] = boxplt_df['muts'] + '\n(' + boxplt_df['growth_phase'] + ')'

boxplt_fimfli_df = boxplt_df[boxplt_df['x_axis'].isin(sorted_order)]
boxplt_fimfli_df['x_axis'] = pd.Categorical(boxplt_fimfli_df['x_axis'], sorted_order, ordered = True)

## plot same for the different growth stages
genotype_gs_col_dict = {'wt\n(earlylog)': '#c8c8c8', 
                        'wt\n(latelog)': '#989898',
                        'FimZ[F126L]\n(earlylog)': '#f9af92', 
                        'FimZ[F126L]\n(latelog)': "#f4875c"}
pval_method = 'text'

## plot as one boxplot but with different colors for the different growth stages
fig, ax = plt.subplots(figsize = (8, 6))
## plot at first white boxplots underneath (otherwise the grid shines through!)
sns.boxplot(data = boxplt_fimfli_df, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_group', 
                hue_order = genotype_col_dict.keys(), palette = ['white', 'white'],    
                linewidth = 0, fliersize = 0, zorder = 5, width = 0.8, ax = ax)
sns.boxplot(data = boxplt_fimfli_df, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_group', 
                hue_order = genotype_col_dict.keys(), palette = genotype_col_dict,    
                linewidth = 1, fliersize = 0, zorder = 10, width = 0.8, ax = ax)
for growth_phase in sorted(boxplt_fimfli_df['growth_phase'].unique()):
    stripplot_gs = boxplt_fimfli_df[boxplt_fimfli_df['growth_phase'] == growth_phase]
    stipplot_mut_gs = [mut_gs for mut_gs in genotype_gs_col_dict.keys() if mut_gs in stripplot_gs['mut_gs'].unique()]
    sns.stripplot(data = stripplot_gs, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_gs', 
                    hue_order = stipplot_mut_gs, palette = genotype_gs_col_dict, 
                    edgecolor = 'k', linewidth=0.25, jitter = 0.15, size = 8,
                    alpha = 0.8, dodge = True, zorder = 20, ax = ax)
## add significance vals
annot = Annotator(ax, mutation_pairs, data=boxplt_fimfli_df, x='x_axis', y='l2fc_rpoD_raw', order=sorted_order, hue='mut_group', hue_order=genotype_col_dict.keys())
if pval_method == 'star': 
    annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
else:
    annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
annot.apply_test()
annot.annotate(line_offset_to_group=0.04, line_offset=-1)
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
ax.set_axisbelow(True)
ax.set_ylabel('Fold change gene expression\n(log2)')
ax.set_xlabel(None)
h, l = ax.get_legend_handles_labels() ## get legend handles and labels
ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5), 
        handles = h[2:], labels = l[2:]) ## cut off the white background boxplots
plt.xticks(fontstyle = 'italic')
plt.tight_layout()
if save_plt:
    plt.savefig(os.path.expanduser(f'~/Desktop/fimZ_phenotype_validation/2025_05_qPCR/proc/plots/2026_02_xx_RTqPCR_P07_Ehorm_fimZfliC_pval{pval_method}.svg'))




boxplt_fim_df = boxplt_fimfli_df[boxplt_fimfli_df['Gene'].isin(sorted_order_fim)]

## plot for fimbriae only
fig, ax = plt.subplots(figsize = (5.5, 4))
## plot at first white boxplots underneath (otherwise the grid shines through!)
sns.boxplot(data = boxplt_fim_df, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_group', 
                hue_order = genotype_col_dict.keys(), palette = ['white', 'white'],    
                linewidth = 0, fliersize = 0, zorder = 5, width = 0.8, ax = ax)
sns.boxplot(data = boxplt_fim_df, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_group', 
                hue_order = genotype_col_dict.keys(), palette = genotype_col_dict,    
                linewidth = 1, fliersize = 0, zorder = 10, width = 0.8, ax = ax)
for growth_phase in sorted(boxplt_fim_df['growth_phase'].unique()):
    stripplot_gs = boxplt_fim_df[boxplt_fim_df['growth_phase'] == growth_phase]
    stipplot_mut_gs = [mut_gs for mut_gs in genotype_gs_col_dict.keys() if mut_gs in stripplot_gs['mut_gs'].unique()]
    sns.stripplot(data = stripplot_gs, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_gs', 
                    hue_order = stipplot_mut_gs, palette = genotype_gs_col_dict, 
                    edgecolor = 'k', linewidth=0.25, jitter = 0.15, size = 8,
                    alpha = 0.8, dodge = True, zorder = 20, ax = ax)
## add significance vals
annot = Annotator(ax, mutation_pairs_fim, data=boxplt_fim_df, x='x_axis', y='l2fc_rpoD_raw', order=sorted_order_fim, hue='mut_group', hue_order=genotype_col_dict.keys())
if pval_method == 'star': 
    annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
else:
    annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
annot.apply_test()
annot.annotate(line_offset_to_group=0.04, line_offset=-1)
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
ax.set_axisbelow(True)
ax.set_ylabel('Fold change gene expression\n(log2)')
ax.set_xlabel(None)
h, l = ax.get_legend_handles_labels() ## get legend handles and labels
ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5), 
        handles = h[2:], labels = l[2:]) ## cut off the white background boxplots
plt.xticks(fontstyle = 'italic')
plt.tight_layout()
if save_plt:
    plt.savefig(os.path.expanduser(f'~/Desktop/fimZ_phenotype_validation/2025_05_qPCR/proc/plots/2025_06_12_RTqPCR_P07_Ehorm_fimZ_pval{pval_method}.svg'))

## save the final df 
if save_plt:
    boxplt_df.to_csv(os.path.expanduser('~/Desktop/fimZ_phenotype_validation/2025_05_qPCR/proc/qPCR_analysis_values.csv'), index = False)

## calculate the impact of the citG mutation in comparison to the ones without citG (but with fimZ)
boxplt_fim_df['genotype_citG'] = np.where(boxplt_fim_df['sampleid'] == 'L00071_02_1', True, False)
for gene, gene_exp in boxplt_fim_df.groupby("Gene"):
    is_fimz = (gene_exp['mut_group'] == 'FimZ[F126L]')
    has_citG = gene_exp.loc[(gene_exp['genotype_citG'] == True) & is_fimz, 'l2fc_rpoD_raw']
    has_nocitG = gene_exp.loc[(gene_exp['genotype_citG'] == False) & is_fimz, 'l2fc_rpoD_raw']
    
    stat, p = stats.mannwhitneyu(has_citG, has_nocitG)
    print(f'For gene expression {gene} there is a p-value between citG mutants and citG wt samples (carrying all fimZ[F126L]) of: {p}')

## print median log2 upregulations
log2fold_diff = boxplt_fim_df.groupby(['Gene', 'mut_group'])['l2fc_rpoD_raw'].median().unstack().reindex(columns=['FimZ[F126L]', 'wt'])
log2fold_diff['diff'] = log2fold_diff['FimZ[F126L]'] - log2fold_diff['wt']
print(f'Identified the following log2-changes (median) between the genes:\n{log2fold_diff.to_string()}')
