
import os
import glob
import re
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

sys.path.insert(0, os.getcwd() + '/../../modules/')

import analysispy_module as apy

### functions

def read_mol_clock_data_file(file):
    with open(file, 'r') as fid:
        mol_clock_data = np.loadtxt(fid, skiprows = 1, delimiter = ',')
    if len(np.shape(mol_clock_data)) == 1: ## if 1d array --> reshape to 2d 
        mol_clock_data = np.reshape(mol_clock_data, (-1, 4))
    return mol_clock_data

def calculate_molecular_clock_rate(mol_clock_data, subj_spec, unique_tps):

    num_timepoints = len(unique_tps)

    print(f'\n\nAnalysing {subj_spec}')
    ## calculate patient specific mutation rate
    if num_timepoints > 1:
        num_samples = len(mol_clock_data)
        tps = mol_clock_data[:,1]
        num_mut_scaled = mol_clock_data[:,2]/mol_clock_data[:,3] * basescale

        regress_dict = apy.calc_regression(np.column_stack( (tps, num_mut_scaled) ),subj_spec,num_samples, num_timepoints)

        mut_rate_obs = regress_dict['slope'] * 365 ## mutation rate based on regression over all isolates 
        ts = abs(stats.t.ppf( (1-confidence) / 2, len(tps)-2)) ## (confidence interval and the degrees of freedom minus num of parameters to estimate)
        SE_mut_rate_obs = regress_dict['std_err'] * 365 
        margin_of_error = SE_mut_rate_obs * ts ## moe = std_err * t-score
        
        ci_lower = mut_rate_obs - margin_of_error
        ci_upper = mut_rate_obs + margin_of_error
        mut_rate_obs_ci95 = [ci_lower, ci_upper]
        
    else:
        print(f'{subj_spec} has only 1 tp --> skipping')
        mut_rate_obs = np.nan
        mut_rate_obs_ci95 = [np.nan, np.nan]
        num_samples = np.nan
    return [mut_rate_obs, mut_rate_obs_ci95, num_samples]

def calculate_tmrca_via_MCsimulation(num_muts_per_smpl, mut_rate, ci95_vals, n_sim=10000, seed=123, mut_rate_fit_method = 'gamma', plot_histogram = False, subj_spec = '', hist_plot_title = ''):
    ## calculate the observed point estimate for the tmrca as well as 95% CI values given Monte Carlo simulations of the mutation rate and 
    ## bootstraps of the number of observed SNVs
    ## To fit the distribution for the mutation rate, choose either a 'normal' (0-bounded), 'log-normal' or 'gamma' distribution
    ## plot histograms to evaluate model fits

    if mut_rate != mut_rate:
        print('No mutation rate given. Skipping calculation')
        return [np.nan, [np.nan, np.nan]]
    
    ## calculate patient specific tMRCA with oberved mutation rate of species
    mean_mut = np.mean(num_muts_per_smpl)
    tmrca_mutrate = mean_mut / mut_rate ## calculate observed point estimate; should be close to `np.mean(bootstrapped_tmrca)`!

    ## return just tmrca if no CI for mutationrate are given
    if ci95_vals == [0, 0]:
        print('WARNING! No confidence intervals were given')
        return [tmrca_mutrate, [np.nan, np.nan]]

    ci_diff_to_meanpoint = (ci95_vals[1]-ci95_vals[0])/2
    ci_diff_to_mean_lower = 1-(mut_rate / (ci95_vals[0]+ci_diff_to_meanpoint))
    ci_diff_to_mean_upper = 1-(mut_rate / (ci95_vals[1]-ci_diff_to_meanpoint))
    print(f'The CI value is skewed to left(negative) or right(posiive) by {ci_diff_to_mean_lower} and {ci_diff_to_mean_upper} ({ci95_vals[0]}, {mut_rate}, {ci95_vals[1]})')
    
    if (mut_rate_fit_method == 'normal') or plot_histogram:
        np.random.seed(seed)
        # fit truncated (0-bounded) normal distribution to sample from n_sim times
        sd_mut_rate = (ci95_vals[1] - ci95_vals[0]) / (2 * 1.96)
        mut_rate_dist_sampled_normal = np.random.normal(loc=mut_rate, scale=sd_mut_rate, size=n_sim)
        mut_rate_dist_sampled_truncnormal = stats.truncnorm((0 - mut_rate) / sd_mut_rate, (np.inf - mut_rate) / sd_mut_rate, loc=mut_rate, scale=sd_mut_rate).rvs(n_sim)
        if mut_rate_fit_method == 'normal':
            mut_rate_dist_sampled = mut_rate_dist_sampled_truncnormal

    if (mut_rate_fit_method == 'log-normal') or plot_histogram:
        np.random.seed(seed)
        # fit distribution to sample from n_sim times (--> mutation-rates are 0-bounded! to avoid negative values, use log-normal distribution)
        log_mut_rate_ci95_l, log_mut_rate, log_mut_rate_ci95_u = np.log([ci95_vals[0], mut_rate, ci95_vals[1]])
        sd_log_mut_rate_log = (log_mut_rate_ci95_u - log_mut_rate_ci95_l) / (2 * 1.96) ## get standard deviation to sample from distribution
        mut_rate_dist_sampled_lognormal = np.random.lognormal(mean=log_mut_rate, sigma=sd_log_mut_rate_log, size=n_sim) ## sample from mutation rate distribution n_sim times 
        if mut_rate_fit_method == 'log-normal':
            mut_rate_dist_sampled = mut_rate_dist_sampled_lognormal

    if (mut_rate_fit_method == 'gamma') or plot_histogram:
        np.random.seed(seed)
        # fit distribution to sample from n_sim times (--> mutation-rates are 0-bounded! to avoid negative values, use gamma distribution)
        ## NOTE: While this fits the CI values truely, it also performs well on credible intervals given those parameters. HOWEVER: This still does not mean that it is the true shape, but it is still rather the better guess than just using the HDP interval estimates without error propagation!!!
        gamma_shape, gamma_scale = fit_gamma_from_ci(mut_rate, ci95_vals[0], ci95_vals[1], mean_weight=1, ci_lower_weight=1, ci_upper_weight=1) ## fit gamma distribution given the mean and CI values via minimization of the difference between guessed and observed values; apply weights if one value is more important
        mut_rate_dist_sampled_gamma_fit = stats.gamma.rvs(a=gamma_shape, scale=gamma_scale, size=n_sim) ## sample from fitted gamma distribution
        if mut_rate_fit_method == 'gamma':
            mut_rate_dist_sampled = mut_rate_dist_sampled_gamma_fit
        
    if plot_histogram:
        trunc_norm_y, trunc_norm_bin_edges = np.histogram(mut_rate_dist_sampled_normal, bins=100)
        trunc_norm_x = 0.5 * (trunc_norm_bin_edges[1:] + trunc_norm_bin_edges[:-1])
        log_normal_y, log_normal_bin_edges = np.histogram(mut_rate_dist_sampled_lognormal, bins=100)
        log_normal_x = 0.5 * (log_normal_bin_edges[1:] + log_normal_bin_edges[:-1])
        gamma_fit_y, gamma_fit_bin_edges = np.histogram(mut_rate_dist_sampled_gamma_fit, bins=100)
        gamma_fit_x = 0.5 * (gamma_fit_bin_edges[1:] + gamma_fit_bin_edges[:-1])

        fig, ax = plt.subplots(figsize = (9,3.5))
        ## plot actual values
        ax.axvline(x = mut_rate, lw = 3, color = 'k', ls = '-', label = 'Observed')
        [ax.axvline(x = ci, lw = 3, color = 'k', ls = '--') for ci in ci95_vals]
        ## plot truncated normal distribution
        [ax.axvline(x = ci, lw = 2, color = '#1b9e77', ls = '--', zorder = 5) for ci in np.percentile(mut_rate_dist_sampled_truncnormal, [2.5, 97.5])]
        ax.axvline(x = np.mean(mut_rate_dist_sampled_truncnormal), lw = 2, color = '#1b9e77', ls = '-')
        ax.plot(trunc_norm_x, trunc_norm_y, color='#1b9e77', lw = 1.5, label = '0-bounded\nnormal')
        ## plot log-normal distribution
        [ax.axvline(x = ci, lw = 2, color = '#7570b3', ls = '--', zorder = 10) for ci in np.percentile(mut_rate_dist_sampled_lognormal, [2.5, 97.5])]
        ax.axvline(x = np.mean(mut_rate_dist_sampled_lognormal), lw = 2, color = '#7570b3', ls = '-')
        ax.plot(log_normal_x, log_normal_y, color='#7570b3', lw = 1.5, label = 'log-normal')
        ## plot gamma distribution
        [ax.axvline(x = ci, lw = 2, color = '#d95f02', ls = '--', zorder = 15) for ci in np.percentile(mut_rate_dist_sampled_gamma_fit, [2.5, 97.5])]
        ax.axvline(x = np.mean(mut_rate_dist_sampled_gamma_fit), lw = 2, color = '#d95f02', ls = '-')
        ax.plot(gamma_fit_x, gamma_fit_y, color='#d95f02', lw = 1.5, label = 'gamma')
        ax.set_ylim((0, None))
        ax.legend(title='Distributions with\nmean (solid lines) and\nCI95% (dashed lines)', loc = 'center left', bbox_to_anchor = (1, 0.5))
        ax.set_title(hist_plot_title)
        plt.tight_layout()
        fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_{subj_spec}_mol_clock_distribution_fit_for_MCsimulation.svg'))

    ## bootstrap n_sim times the number of mutations and get mean number of mutation 
    bootstrap_idx = np.random.randint(0, len(num_muts_per_smpl), size=(n_sim, len(num_muts_per_smpl)))
    bootstrapped_mean_mut = np.mean(num_muts_per_smpl[bootstrap_idx], axis = 1)
    bootstrapped_tmrca = bootstrapped_mean_mut / mut_rate_dist_sampled 
    tmrca_mutrate_ci95 = np.percentile(bootstrapped_tmrca, [2.5, 97.5])

    ## calculate point estimate 
    print(f'The difference between the observed tMRCA point estimate and the MC simulated point estimate is: {np.mean(bootstrapped_tmrca) - tmrca_mutrate}')
    print(f'The difference between the observed tMRCA point estimate and the median of the MC simulations is: {np.median(bootstrapped_tmrca) - tmrca_mutrate}')

    return [tmrca_mutrate, tmrca_mutrate_ci95]

def objective_function(gamma_args_tuple_log, mut_rate, ci_lower, ci_upper, mean_weight=1, ci_lower_weight=1, ci_upper_weight=1):
    ## calculate difference between guessed and observed mutation rate (mean and CIs) and 
    ## sum up (weighted) squared errors to minimize upon 
    gamma_args_tuple = np.exp(gamma_args_tuple_log)
    # validate that agruments are ≥0
    if gamma_args_tuple[0] <= 0 or gamma_args_tuple[1] <= 0:
        return np.inf

    # get mean of gamma distribution
    gamma_mean = gamma_args_tuple[0] * gamma_args_tuple[1]

    # get CI of gamma distribution
    gamma_lower_quantile = stats.gamma.ppf(0.025, a=gamma_args_tuple[0], scale=gamma_args_tuple[1])
    gamma_upper_quantile = stats.gamma.ppf(0.975, a=gamma_args_tuple[0], scale=gamma_args_tuple[1])

    ## get errors and return sum of errors with weights (optional) to apply minimization
    error_mean = ((gamma_mean - mut_rate)/mut_rate)**2 ## sum of squared error without dimensionality due to normalization
    error_lower_ci = ((gamma_lower_quantile - ci_lower)/ci_lower)**2 ## sum of squared error without dimensionality due to normalization
    error_upper_ci = ((gamma_upper_quantile - ci_upper)/ci_upper)**2 ## sum of squared error without dimensionality due to normalization

    error_value = error_mean*mean_weight + error_lower_ci*ci_lower_weight + error_upper_ci*ci_upper_weight
    return error_value

def fit_gamma_from_ci(mut_rate, ci_lower, ci_upper, mean_weight=1, ci_lower_weight=1, ci_upper_weight=1, low_val = 1e-16):
    ## fitting a gamma distribution to a mutation rate given the mean and 95% confidence intervals
    if ci_upper - ci_lower == 0:
        return None, None
    else:
        sd_mut_rate = (ci_upper - ci_lower) / (2 * 1.96) ## zsocre of 1.96 corresponds to 95% confidence interval; Note it is a good start but gamma will differ!
        var_mut_rate = sd_mut_rate**2

    ini_scale = var_mut_rate / mut_rate ## https://library.virginia.edu/data/articles/getting-started-with-gamma-regression
    ini_shape = mut_rate**2 / var_mut_rate ## https://library.virginia.edu/data/articles/getting-started-with-gamma-regression
    
    ## transform to log space that the optimizer does not need bounds which can be very fragile (optimization function returns to normal space, while the optimizer stays in log space!)
    ini_scale_log = np.log(ini_scale) 
    ini_shape_log = np.log(ini_shape)
    
    minimizer_result = minimize(objective_function,
                                (ini_shape_log, ini_scale_log),
                                args=(mut_rate, ci_lower, ci_upper, mean_weight, ci_lower_weight, ci_upper_weight),
                                method='L-BFGS-B'  # default method using bounds
                                )
    if minimizer_result.success:
        return np.exp(minimizer_result.x)
    else:
        print(f'Optimization failed: {minimizer_result.message}')
        return None, None


def get_tmrca_per_tp(mol_clock_data, unique_tps, mut_rate = np.nan, mut_rate_ci95 = [np.nan, np.nan], basescale = 1000000, req_isolate_num = 5, subj_spec = '', plot_histogram=False, mol_clock_type = '', hist_plot_title = ''):
    
    ## loop over all timepoints 
    ploted_hist = False
    for tp in unique_tps:
        tp_idx = np.where(mol_clock_data[:, 1] == tp)[0]
        ## check if tp has at least num required isolates to calculate tMRCA
        if len(tp_idx) >= req_isolate_num:
            
            ## calculate scaled number of mutations
            num_mut_scaled = mol_clock_data[tp_idx,2]/mol_clock_data[tp_idx,3] * basescale

            ## plot histogram only once!
            if ploted_hist:
                plot_histogram = False
            ploted_hist = True
            [tmrca_mean, tmrca_ci95] = calculate_tmrca_via_MCsimulation(num_mut_scaled * 365, mut_rate, mut_rate_ci95, mut_rate_fit_method = 'gamma', n_sim=10000, seed=123, plot_histogram=plot_histogram, subj_spec = subj_spec, hist_plot_title = hist_plot_title) ## note: *365 is to convert mutation rate from a yearly basis to a daily! By multiplying num_mut_scaled it is easier than implementing it to the list of CIs as well!
            
            ## get how many days since hosp admission to tp:
            d_since_hospadmission = np.unique(mol_clock_data[tp_idx, 0])[0]

            mut_rate_dict[tp + 1] = {'days_in_hosp': d_since_hospadmission,
                                    'days_to_root': np.unique(mol_clock_data[tp_idx,1])[0],
                                    'num_isolates': len(tp_idx),
                                    'mut_rate': mut_rate,
                                    'mut_rate_ci95': mut_rate_ci95,
                                    'mean_mut': np.mean(num_mut_scaled/basescale), 
                                    'ci95_mut': np.percentile(num_mut_scaled/basescale, [2.5, 97.5]), 
                                    'num_muts': f"[{','.join([str(mut) for mut in num_mut_scaled])}]",
                                    'tmrca_meanmut': tmrca_mean,
                                    'tmrca_meanmut_ci95': tmrca_ci95,
                                    'mol_clock_type': mol_clock_type
                                    }
    ## convert data to df 
    tmp_mutation_rate_df = pd.DataFrame(mut_rate_dict).T.reset_index()
    tmp_mutation_rate_df.rename(columns={'index':'timepoint'}, inplace = True)
    tmp_mutation_rate_df['subj_spec'] = subj_spec

    return tmp_mutation_rate_df

def add_metadata_to_df(tmrca_spec, inf_day_dict):

    tmrca_spec['inf_date'] = tmrca_spec['subj_spec'].map(inf_day_dict)
    tmrca_spec['subj_spec'] = tmrca_spec['subj_spec'].str.split('_').str[1:3].str.join('_')

    tmrca_spec['patient'] = tmrca_spec['subj_spec'].str.split('_').str[0]
    tmrca_spec['species'] = tmrca_spec['subj_spec'].str.split('_').str[1].str.split('-').str[0]
    tmrca_spec['species'] = tmrca_spec['species'].str[0] + '. ' + tmrca_spec['species'].str[1:]
    tmrca_spec['clade'] = tmrca_spec['subj_spec'].str.split('_').str[1].str.split('-').str[-1].str.replace('c', '').astype(int)
    tmrca_spec['spec_clade_subj_plt'] = tmrca_spec['species'] + ' clade' + tmrca_spec['clade'].astype(str) + ' (' + tmrca_spec['patient'] + ')'
    tmrca_spec['spec_subj_plt'] = tmrca_spec['species'] + ' (' + tmrca_spec['patient'] + ')'

    tmrca_spec = tmrca_spec.sort_values(['patient', 'species', 'clade'])
    return tmrca_spec

def combine_literature_inferred_clock_for_plot(tmrca_spec_df, clade_to_use_inferred_clocks = [], annotate_clock = {}, color_of_clock = {}):

    ## bring the respective plots to one column for plot
    tmrca_spec_df['tmrca_meanmut'] = tmrca_spec_df['tmrca_meanmut_lit']
    tmrca_spec_df['tmrca_meanmut_ci95'] = tmrca_spec_df['ci95_meanmut_tmrca_lit']
    tmrca_spec_df['mut_rate'] = tmrca_spec_df['mut_rate_lit']
    tmrca_spec_df['ci95_mut_rate'] = tmrca_spec_df['mut_rate_lit_ci95']
    is_ehorm_bool = (tmrca_spec_df['spec_clade_subj_plt'].isin(clade_to_use_inferred_clocks))
    tmrca_spec_df.loc[is_ehorm_bool, 'tmrca_meanmut'] = tmrca_spec_df.loc[is_ehorm_bool, 'tmrca_meanmut_pat_spec']
    tmrca_spec_df.loc[is_ehorm_bool, 'tmrca_meanmut_ci95'] = tmrca_spec_df.loc[is_ehorm_bool, 'ci95_meanmut_tmrca_pat_spec']
    tmrca_spec_df.loc[is_ehorm_bool, 'mut_rate'] = tmrca_spec_df.loc[is_ehorm_bool, 'mut_rate_obs']
    tmrca_spec_df.loc[is_ehorm_bool, 'ci95_mut_rate'] = tmrca_spec_df.loc[is_ehorm_bool, 'mut_rate_obs_ci95']
    
    ## map which species has which clock gotten // hardcoded
    tmrca_spec_inf_clades['clock_for_inference'] = tmrca_spec_inf_clades['species'].map(annotate_clock)
    tmrca_spec_inf_clades['plt_color'] = tmrca_spec_inf_clades['clock_for_inference'].map(color_of_clock)
    
    return tmrca_spec_df

def extract_listentry_pd_col(pd_col, n):
    return pd_col.apply(lambda x: x[n] if len(x) > n else None)

############################################################
############################################################
############################################################

## get genome Lengths and calculate frequency where genome length is required
genome_length_dict = {}
for ref_genome_folder in glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/*240825/')): 
    species = ref_genome_folder.split('/')[-2].split('_')[1]
    [_, genomeLength, _] = apy.genomestats(ref_genome_folder);
    genome_length_dict[species] = genomeLength


############################################################
############################################################
############################################################
## Set variables

confidence=0.95
basescale = 1000000 # rescale genome-wide molecular rate to human-readable numbers
mutation_rate_df = pd.DataFrame()

## ncbi search: ((("Pseudomonas aeruginosa") OR (Pseudomonas)) AND (("mutation rate") OR ("molecular clock") OR ("substitution rate") AND (("/year") OR ("per year")) 
## read in mutation rates from file 
mut_rates_literature = pd.read_excel('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/2025_12_09_avail_published_mol_clocks_HAI_project.xlsx', sheet_name='Sheet3')
## subset to species_of_interest
mut_rates_literature = mut_rates_literature[(mut_rates_literature['Genus'] == 'Staphylococcus') & 
                      (mut_rates_literature['species'] == 'aureus') & 
                      mut_rates_literature['Strain'].notna() & 
                      mut_rates_literature['mutation rate (credible intervals)'].notna() & 
                      (mut_rates_literature['credible interval type'] != 'IQR') & ## exclude interquartile ranges
                      (mut_rates_literature['unit'].isin(['SNV/site/year', 'SNVs/Mb/year']))
                      ].reset_index(drop = True)

mut_rates_literature.rename(columns={'mutation rate (point estimate)': 'mol_clock_pe',
                          'mutation rate (credible intervals)': 'mol_clock_ci'}, inplace=True)

mut_rates_literature['mol_clock_ci_l'] = mut_rates_literature['mol_clock_ci'].apply(lambda s: re.split(r'(?<![eE])-', s)[0] if isinstance(s, str) else s)
mut_rates_literature['mol_clock_ci_h'] = mut_rates_literature['mol_clock_ci'].apply(lambda s: re.split(r'(?<![eE])-', s)[1] if isinstance(s, str) else s)

## transform molecular clocks to same units: SNV/Mb/year
for col in ['mol_clock_pe', 'mol_clock_ci_l', 'mol_clock_ci_h']:
    mut_rates_literature[col] = mut_rates_literature[col].astype(float)
    mut_rates_literature.loc[mut_rates_literature['unit'] == 'SNV/site/year', col] = mut_rates_literature.loc[mut_rates_literature['unit'] == 'SNV/site/year', col].astype(float) * 1e6 ## transform to SNV/Mb/year
mut_rates_literature['species'] = mut_rates_literature['Genus'].str[0] + mut_rates_literature['species']

## get a dict of the tp of infections in respect to hospital admission (note: only the clades with infectious isolates will be stated!)
inf_day_dict = {'P0007_P07_Apittii-c1_240825': 50.89,
                'P0007_P07_Bthetaiotaomicron-c1_240825': 50.831,
                'P0007_P07_Ehormaechei-c1_240825': 38.688,
                'P0010_P10_Kmichiganensis-c1_240825': 3.127,
                'P0021_P21_Saureus-ST97-c3_240825': 2.011}

## plot specific variables
fontsize = 12
markersize = 25
ci_lower_idx = 0
inf_color = '#d73027'
tmrca_col = '#EAECCC'
mol_clock_col = {'Applied': '#484848',
                 'Alternative': '#c8c8c8'}

############################################################
############################################################
############################################################
## fetch files 
files = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/*/*/*_mol_clock_data_noinf.csv'))
files_inf = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/*/*/*_mol_clock_data_winf.csv'))

## subset data to certain date
files = sorted([file for file in files if '2026-03-26' in file.split('/')[-1]])
files_inf = sorted([file for file in files_inf if '2026-03-26' in file.split('/')[-1]])

############################################################
############################################################
############################################################

mutrate_noinf_obs_dict = {}
for file in files:
    req_isolate_num = 8

    subj_spec = file.split('/')[-3]
    ## read in molecular clock without infectious isolates
    mol_clock_data = read_mol_clock_data_file(file)

    ## read in molecular clock _with_ infectious isolates
    file_inf = [file for file in files_inf if subj_spec in file][0]
    mol_clock_data_inf = read_mol_clock_data_file(file_inf)

    unique_tps = set(mol_clock_data[:, 1]) ## get all unique dates since first observation of pathogen (without clinical dates!)
    
    ## remove from inf dataset the microbiome timepoints
    inf_entries_bool = ~np.isin(mol_clock_data_inf[:, 1], list(unique_tps))
    mol_clock_data_inf_only = mol_clock_data_inf[inf_entries_bool, :]
    unique_inf_tps = set(mol_clock_data_inf_only[:, 1]) ## note: it is days to the first isolate!
    ## get date if first clinical observation of pathogen
    if len(unique_inf_tps) > 0:
        first_inf_tp = min(unique_inf_tps)
    else:
        first_inf_tp = np.nan
    mutrate_noinf_obs_dict[subj_spec] = calculate_molecular_clock_rate(mol_clock_data, subj_spec, unique_tps)
    mut_rate_obs = mutrate_noinf_obs_dict[subj_spec][0]
    mut_rate_obs_ci95 = mutrate_noinf_obs_dict[subj_spec][1]
    if not subj_spec in ['P0007_P07_Ehormaechei-c1_240825']:
        print(f'Observed mutation rate for {subj_spec}: {mut_rate_obs} CI95 {mut_rate_obs_ci95}')
        mut_rate_obs = np.nan
        mut_rate_obs_ci95 = np.nan 
    
    mut_rate_dict = {}
    species_clade = subj_spec.split('_')[2]
    species = species_clade.split('-')[0]
    if species in mut_rates_literature['species'].values:
        mut_rates_lit_spec = mut_rates_literature[mut_rates_literature['species']==species]
    else:
        mut_rates_lit_spec = pd.DataFrame(columns = ['species', 'Strain', 'mol_clock_pe', 'mol_clock_ci_l', 'mol_clock_ci_h'])
    ## check if there is any timepoint within ≤48 h (needed to account for rounding errors due to the use of integers (days!); cutoff need to be checked again manually for ≤24h!) which has more than the number of isolates required
    if first_inf_tp == first_inf_tp:
        if not any([(len(np.where(mol_clock_data[:, 1] == tp)[0]) >= req_isolate_num) & (tp <= first_inf_tp+2) for tp in unique_tps]): 
            print(f'Did not find any timepoint with more than {req_isolate_num} isolates')
            continue
    else:
        ## lineage does not cause an infection
        if not any([(len(np.where(mol_clock_data[:, 1] == tp)[0]) >= req_isolate_num) for tp in unique_tps]): 
            print(f'Did not find any timepoint with more than {req_isolate_num} isolates')
            continue

    ## calculate tmrca
    tmp_mutation_rate_obs_df = get_tmrca_per_tp(mol_clock_data, unique_tps, mut_rate_obs, mut_rate_obs_ci95, basescale, req_isolate_num, subj_spec, plot_histogram=True, mol_clock_type = 'observed', hist_plot_title = subj_spec)
    tmp_mutation_rate_obs_df['Strain'] = 'observed'
    mutation_rate_df = pd.concat([mutation_rate_df, tmp_mutation_rate_obs_df])
    
    for rowidx, mut_rate_strain in mut_rates_lit_spec.iterrows():
        tmp_mutation_rate_obs_df = get_tmrca_per_tp(mol_clock_data, unique_tps, mut_rate_strain['mol_clock_pe'], [mut_rate_strain['mol_clock_ci_l'], mut_rate_strain['mol_clock_ci_h']], basescale, req_isolate_num, subj_spec, plot_histogram=True, mol_clock_type = 'literature', hist_plot_title = f'{subj_spec}\n{mut_rate_strain["paper"]}')
        tmp_mutation_rate_obs_df['Strain'] = mut_rate_strain['Strain'] + '_' + mut_rate_strain['paper'] + '_' + str(mut_rate_strain.name)
        mutation_rate_df = pd.concat([mutation_rate_df, tmp_mutation_rate_obs_df])
mutation_rate_df = mutation_rate_df.reset_index(drop = True)

## plot per species just the first TP and its tmrca with repect to the infection 
## get first tps
tmrca_spec_inf = mutation_rate_df.loc[mutation_rate_df.groupby('subj_spec')['timepoint'].transform('min').eq(mutation_rate_df['timepoint'])].reset_index(drop=True)

tmrca_spec_inf = add_metadata_to_df(tmrca_spec_inf, inf_day_dict)

####################
## plot with all clades which have enough data 
legend_info_tmrca = [Line2D([0], [0], marker = 'x', color = inf_color, linewidth = 0, markeredgewidth = 2, label = 'HAI\nonset', markersize = 8)]
legend_info_tmrca += [plt.Rectangle((0,0),1,1, facecolor = tmrca_col, edgecolor = 'w', linewidth = 0, label = 'TMRCA of pathogenic\nlineage in microbiome\n(mean & 95% CI)')]
legend_info_tmrca += [Line2D([0], [0], marker = '_', color = mol_clock_col['Applied'], markerfacecolor = mol_clock_col['Applied'], linewidth = 2, markeredgecolor = 'k', markeredgewidth=0.2, label = 'Molecular clock\n(mean & 95% CI)', markersize = 8)]

## plot with clades which contributed to inf and had enough data only
tmrca_spec_inf_clades = tmrca_spec_inf[~tmrca_spec_inf['inf_date'].isna()].copy()

## literature clock, for Ehormaechei --> inferred clock 
tmrca_spec_inf_clades = tmrca_spec_inf_clades[tmrca_spec_inf_clades['tmrca_meanmut'].notna()]

## Removal from hypermutator candidates 
candidate_hypermutators = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/candidate_hypermutators_based_on_outgroup_background.csv')
candidate_hypermutators = candidate_hypermutators.loc[candidate_hypermutators['is_candidate_hypermutator'] == True, 'genome'].values

tmrca_spec_inf_clades_nohypmut = tmrca_spec_inf_clades[~tmrca_spec_inf_clades['subj_spec'].isin(candidate_hypermutators) | (tmrca_spec_inf_clades['mol_clock_type'] == 'observed')].reset_index(drop = True)
tmrca_sel_set = {
    ('E. hormaechei clade1 (P07)', 'observed'), 
    ('S. aureus clade3 (P21)', 'ST97_Key et al 2023_5')
}

selected_tmrcas_bool = pd.Series(zip(tmrca_spec_inf_clades_nohypmut['spec_clade_subj_plt'], tmrca_spec_inf_clades_nohypmut['Strain'])).isin(tmrca_sel_set)
tmrca_spec_inf_clades_nohypmut_sel = tmrca_spec_inf_clades_nohypmut[selected_tmrcas_bool]
tmrca_spec_inf_clades_nohypmut_unsel = tmrca_spec_inf_clades_nohypmut[~selected_tmrcas_bool]
tmrca_spec_inf_clades_nohypmut_unsel['study_year'] = tmrca_spec_inf_clades_nohypmut_unsel['Strain'].str.split(' ').str[-1].str.split('_').str[0]
tmrca_spec_inf_clades_nohypmut_unsel = tmrca_spec_inf_clades_nohypmut_unsel.sort_values('study_year').reset_index(drop = True)

## Plot tMRCA as in main figure
fig, axs = plt.subplots(figsize = (6, 5), nrows = 2, sharex = True, gridspec_kw={'height_ratios': [3, 1]})
ax = axs.flatten()
ax[0].axhline(0, ls = '--', lw = 1, c = 'k', alpha = 0.5)
ax[0].annotate(text = 'Hospital admission', xy = (1.4, 0.25), fontsize=fontsize, ha = 'right')
ax[0].scatter(x= tmrca_spec_inf_clades_nohypmut_sel['species'], 
              y = tmrca_spec_inf_clades_nohypmut_sel['inf_date'], 
              linestyle='None', 
              s = 100, 
              marker='x', 
              color = inf_color, 
              zorder = 5)
ax[0].bar(x = tmrca_spec_inf_clades_nohypmut_sel['species'], 
            height = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut_ci95'], 1) - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut_ci95'], ci_lower_idx),
            bottom = tmrca_spec_inf_clades_nohypmut_sel['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut_ci95'], 1),
            color = tmrca_col, 
            width = 0.65, 
            zorder = 3)
ax[0].scatter(x= tmrca_spec_inf_clades_nohypmut_sel['species'], 
              y = tmrca_spec_inf_clades_nohypmut_sel['days_in_hosp'] - tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut'], 
              color = '#555555', 
              linestyle='None', 
              s = 3050, 
              lw = 1, 
              marker='_', 
              zorder = 7)
ax[0].set_ylim((np.nanmin(tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['tmrca_meanmut_ci95'], 1)) - 300,
            np.nanmax(tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['tmrca_meanmut_ci95'], ci_lower_idx)) + 150,
            ))
ax[0].set_ylabel('Time since hospital admission\n(in days)')
ax[0].set_yscale('symlog')
ax[0].get_yaxis().set_major_formatter(ScalarFormatter())
ax[0].legend(handles = legend_info_tmrca, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=fontsize+1)
ax[0].set_title('Acquisition time', loc='right', ha = 'right', fontsize=fontsize-2)

ax[1].scatter(x= tmrca_spec_inf_clades_nohypmut_sel['species'], 
              y = tmrca_spec_inf_clades_nohypmut_sel['mut_rate'], 
              linestyle='None', 
              s = 200, 
              marker='_',
              edgecolor = mol_clock_col['Applied'],
              lw = 2, 
              color = mol_clock_col['Applied'], 
              zorder = 15)
ax[1].vlines(x = tmrca_spec_inf_clades_nohypmut_sel['species'], 
           ymin = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['mut_rate_ci95'], 0),
           ymax = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['mut_rate_ci95'], 1),
           color = mol_clock_col['Applied'], 
           lw = 2, 
           zorder = 15)
ax[1].set_yticks(np.arange( np.ceil(max( 
                                    np.nanmax(extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['mut_rate_ci95'], 1)), 
                                    np.nanmax(tmrca_spec_inf_clades_nohypmut['mut_rate'])) 
                                    + 1))) 
ax[1].set_ylabel('Mutations/\nMb/year')
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=fontsize+1)
ax[1].set_title('Substitution rates', loc='right', ha = 'right', fontsize=fontsize-2)
ax[1].set_ylim(0, None)
[ax[idx].spines['top'].set_visible(False) for idx, a in enumerate(ax)]
[ax[idx].spines['right'].set_visible(False) for idx, a in enumerate(ax)]
plt.tick_params(labelsize=fontsize)
plt.tight_layout()
plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2026_05_15_tmrca_nonhypermut_patho_clades_min8isolates_lit_clocks_Saureus_wclocks_10kMCsim_main.svg'))




## plot the molecular clocks and their respective tMRCAs for S aureus separately all again in order of the published years (supplementary figure)
tmrca_saureus = tmrca_spec_inf_clades_nohypmut[tmrca_spec_inf_clades_nohypmut['species'] == 'S. aureus'].copy()
tmrca_saureus['pub_year'] = tmrca_saureus['Strain'].str.split(' ').str[-1].str.split('_').str[0]
tmrca_saureus['study'] = tmrca_saureus['Strain'].str.split('_').str[1:-1].str.join(' ')
tmrca_saureus['strain'] = tmrca_saureus['Strain'].str.split('_').str[0]
tmrca_saureus['study_strain'] = tmrca_saureus['study'] + '\n' + tmrca_saureus['strain']
tmrca_saureus = tmrca_saureus.sort_values(['pub_year', 'study']).reset_index(drop = True)

legend_info_tmrca = [Line2D([0], [0], marker = 'x', color = inf_color, linewidth = 0, markeredgewidth = 2, label = 'HAI\nonset', markersize = 8)]
legend_info_tmrca += [plt.Rectangle((0,0),1,1, facecolor = tmrca_col, edgecolor = 'w', linewidth = 0, label = 'TMRCA of pathogenic\nlineage in microbiome\n(mean & 95% CI)')]
legend_info_tmrca += [Line2D([0], [0], marker = '_', color = mol_clock_col['Applied'], markerfacecolor = mol_clock_col['Applied'], linewidth = 2, markeredgecolor = 'k', markeredgewidth=0.2, label = 'Molecular clock\n(mean & 95% CI)', markersize = 8)]

fig, axs = plt.subplots(figsize = (8, 6), nrows = 2, sharex = True, gridspec_kw={'height_ratios': [3, 1.25]})
ax = axs.flatten()
ax[0].axhline(0, ls = '--', lw = 1, c = 'k', alpha = 0.5)
ax[0].annotate(text = 'Hospital admission', xy = (len(tmrca_saureus['study_strain']) - 0.55, 0.25), fontsize=fontsize, ha = 'right')
ax[0].bar(x = tmrca_saureus['study_strain'], 
            height = extract_listentry_pd_col(tmrca_saureus['tmrca_meanmut_ci95'], 1) - extract_listentry_pd_col(tmrca_saureus['tmrca_meanmut_ci95'], ci_lower_idx),
            bottom = tmrca_saureus['days_in_hosp'] - extract_listentry_pd_col(tmrca_saureus['tmrca_meanmut_ci95'], 1),
            color = tmrca_col, 
            width = 0.65, 
            zorder = 3)
ax[0].scatter(x= tmrca_saureus['study_strain'], 
              y = tmrca_saureus['days_in_hosp'] - tmrca_saureus['tmrca_meanmut'], 
              color = '#555555', 
              linestyle='None', 
              s = 350, 
              lw = 1, 
              marker='_', 
              zorder = 7)
## mark HAI onset
ax[0].add_patch(mpatches.Rectangle((-0.35, np.unique(tmrca_saureus['inf_date'])[0]-0.05), len(tmrca_saureus['study_strain'])-0.3, 0.1, color=inf_color, alpha = 0.2, zorder = 2))
ax[0].annotate(text = 'Hospital admission', xy = (len(tmrca_saureus['study_strain']) - 0.25, 0.25), fontsize=fontsize, ha = 'right')

ax[0].set_ylim((np.nanmin(tmrca_saureus['days_in_hosp'] - extract_listentry_pd_col(tmrca_saureus['tmrca_meanmut_ci95'], 1)) - 300,
            np.nanmax(tmrca_saureus['days_in_hosp'] - extract_listentry_pd_col(tmrca_saureus['tmrca_meanmut_ci95'], ci_lower_idx)) + 150,
            ))
ax[0].set_ylabel('Time since hospital admission\n(in days)')
ax[0].set_yscale('symlog')
ax[0].get_yaxis().set_major_formatter(ScalarFormatter())
ax[0].legend(handles = legend_info_tmrca, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=fontsize+1)
ax[0].set_ylim(None, 101)
ax[0].set_title('Acquisition time', loc='right', ha = 'right', fontsize=fontsize-2)

ax[1].scatter(x= tmrca_saureus['study_strain'], 
              y = tmrca_saureus['mut_rate'], 
              linestyle='None', 
              s = 200, 
              marker='_',
              edgecolor = mol_clock_col['Applied'],
              lw = 1.5, 
              color = mol_clock_col['Applied'], 
              zorder = 15)
ax[1].vlines(x = tmrca_saureus['study_strain'], 
           ymin = extract_listentry_pd_col(tmrca_saureus['mut_rate_ci95'], 0),
           ymax = extract_listentry_pd_col(tmrca_saureus['mut_rate_ci95'], 1),
           color = mol_clock_col['Applied'], 
           lw = 1.5, 
           zorder = 15)
ax[1].set_yticks(np.arange( np.ceil(max( 
                                    np.nanmax(extract_listentry_pd_col(tmrca_saureus['mut_rate_ci95'], 1)), 
                                    np.nanmax(tmrca_saureus['mut_rate'])) 
                                    + 1))) 
ax[1].set_ylabel('Mutations/\nMb/year')
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=fontsize+1)
ax[1].set_title('Substitution rates', loc='right', ha = 'right', fontsize=fontsize-2)
ax[1].set_ylim(0, None)
ax[1].set_xlim(-0.75, len(tmrca_saureus['study_strain']) - 0.5)
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation = 90, ha = 'center', va = 'top')
[ax[idx].spines['top'].set_visible(False) for idx, a in enumerate(ax)]
[ax[idx].spines['right'].set_visible(False) for idx, a in enumerate(ax)]
plt.tick_params(labelsize=fontsize)
plt.tight_layout()
plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2026_05_15_tmrca_nonhypermut_patho_clades_min8isolates_lit_clocks_Saureus_wclocks_10kMCsim.svg'))




## All together within one plot (R2R)
legend_info_tmrca = [Line2D([0], [0], marker = 'x', color = inf_color, linewidth = 0, markeredgewidth = 2, label = 'HAI\nonset', markersize = 8)]
legend_info_tmrca += [plt.Rectangle((0,0),1,1, facecolor = tmrca_col, edgecolor = 'w', linewidth = 0, label = 'tMRCA of pathogenic\nlineage in microbiome\n(mean & 95% CI)')]
legend_info_tmrca += [Line2D([0], [0], marker = 'o', color = mol_clock_col['Alternative'], markerfacecolor = mol_clock_col['Alternative'], markeredgecolor = 'k', markeredgewidth=0.2, linewidth = 2, label = 'tMRCA using an\nalternative molecular clock\n(mean & 95% CI)', markersize = 8)]

legend_info_mol_clock = [Line2D([0], [0], linewidth = 0, label = 'Molecular clock\n(mean & 95% CI)')]
for label, col in mol_clock_col.items():
    legend_info_mol_clock += [Line2D([0], [0], marker = 'o', color = col, markerfacecolor = col, linewidth = 2, markeredgecolor = 'k', markeredgewidth=0.2, label = label, markersize = 8)]

fig, axs = plt.subplots(figsize = (6, 5), nrows = 2, sharex = True, gridspec_kw={'height_ratios': [3, 1]})
ax = axs.flatten()
ax[0].axhline(0, ls = '--', lw = 1, c = 'k', alpha = 0.5)
ax[0].annotate(text = 'Hospital admission', xy = (1.4, 0.25), fontsize=fontsize, ha = 'right')
ax[0].scatter(x= tmrca_spec_inf_clades_nohypmut_sel['species'], 
              y = tmrca_spec_inf_clades_nohypmut_sel['inf_date'], 
              linestyle='None', 
              s = 100, 
              marker='x', 
              color = inf_color, 
              zorder = 5)
ax[0].bar(x = tmrca_spec_inf_clades_nohypmut_sel['species'], 
            height = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut_ci95'], 1) - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut_ci95'], ci_lower_idx),
            bottom = tmrca_spec_inf_clades_nohypmut_sel['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut_ci95'], 1),
            color = tmrca_col, 
            width = 0.65, 
            zorder = 3)
ax[0].scatter(x= tmrca_spec_inf_clades_nohypmut_sel['species'], 
              y = tmrca_spec_inf_clades_nohypmut_sel['days_in_hosp'] - tmrca_spec_inf_clades_nohypmut_sel['tmrca_meanmut'], 
              color = '#555555', 
              linestyle='None', 
              s = 3050, 
              lw = 1, 
              marker='_', 
              zorder = 7)
ax[0].set_ylim((np.nanmin(tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['tmrca_meanmut_ci95'], 1)) - 300,
            np.nanmax(tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['tmrca_meanmut_ci95'], ci_lower_idx)) + 150,
            ))
ax[0].set_ylabel('Time since hospital admission\n(in days)')
ax[0].set_yscale('symlog')
ax[0].get_yaxis().set_major_formatter(ScalarFormatter())
ax[0].legend(handles = legend_info_tmrca, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=fontsize+1)
ax[0].set_title('Acquisition time', loc='right', ha = 'right', fontsize=fontsize-2)

ax[1].scatter(x= tmrca_spec_inf_clades_nohypmut_sel['species'], 
              y = tmrca_spec_inf_clades_nohypmut_sel['mut_rate'], 
              linestyle='None', 
              s = 30, 
              marker='o',
              edgecolor = 'k',
              lw = 0.2, 
              color = mol_clock_col['Applied'], 
              zorder = 15)
ax[1].vlines(x = tmrca_spec_inf_clades_nohypmut_sel['species'], 
           ymin = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['mut_rate_ci95'], 0),
           ymax = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sel['mut_rate_ci95'], 1),
           color = mol_clock_col['Applied'], 
           lw = 2, 
           zorder = 15)
ax[1].legend(handles = legend_info_mol_clock, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax[1].set_yticks(np.arange( np.ceil(max( 
                                    np.nanmax(extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['mut_rate_ci95'], 1)), 
                                    np.nanmax(tmrca_spec_inf_clades_nohypmut['mut_rate'])) 
                                    + 1))) 
ax[1].set_ylabel('Mutations/\nMb/year')
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=fontsize+1)
ax[1].set_title('Substitution rates', loc='right', ha = 'right', fontsize=fontsize-2)

## add the additional molecular clocks which have been unselected 
for spec in tmrca_spec_inf_clades_nohypmut_unsel['spec_clade_subj_plt'].unique():
    tmrca_add_mol_clocks_spec = tmrca_spec_inf_clades_nohypmut_unsel[tmrca_spec_inf_clades_nohypmut_unsel['spec_clade_subj_plt'] == spec].reset_index(drop = True)
    x_idx = np.where(tmrca_spec_inf_clades_nohypmut_sel['spec_clade_subj_plt'] == spec)[0] ## get the x axis idx for the plot 
    ## define the interspacing (-.4 to .4 with the mid zone (-.1 to .1) left empty)
    num_additional_clocks = len(tmrca_add_mol_clocks_spec)
    min_plot, max_plot, mid_free_l, mid_free_u = -0.25, 0.25, -0.05, 0.05
    xregion_left = max_plot-min_plot-(mid_free_u-mid_free_l)
    interspacing = xregion_left/num_additional_clocks
    # positions in the compressed interval without the middle gap
    compressed_x = np.linspace(min_plot, max_plot - (mid_free_u - mid_free_l), num_additional_clocks)
    # shift points that fall to the right of the gap
    x_positions = np.where(compressed_x >= mid_free_l,
                            compressed_x + (mid_free_u - mid_free_l),
                            compressed_x)
    for data_idx, datapoint in tmrca_add_mol_clocks_spec.iterrows():
        datapoints_x = x_idx+x_positions[data_idx]
        ax[0].scatter(x = datapoints_x, 
                        y = datapoint['days_in_hosp'] - datapoint['tmrca_meanmut'], 
                        linestyle='None', 
                        alpha = 0.6, 
                        marker = 'o', 
                        s = 12,
                        color = mol_clock_col['Alternative'],
                        edgecolor = 'k', 
                        lw = 0.2,
                        zorder = 10)
        ax[0].vlines(x = datapoints_x, 
                     ymin = datapoint['days_in_hosp'] - datapoint['tmrca_meanmut_ci95'][0],
                     ymax = datapoint['days_in_hosp'] - datapoint['tmrca_meanmut_ci95'][1],
                     linestyle='-', 
                     lw = 1,
                     alpha = 0.6, 
                     color = mol_clock_col['Alternative'],
                     zorder = 8)
        ax[1].scatter(x = datapoints_x, 
                        y = datapoint['mut_rate'], 
                        linestyle='None', 
                        alpha = 0.6, 
                        marker = 'o', 
                        s = 12,
                        color = mol_clock_col['Alternative'],
                        edgecolor = 'k', 
                        lw = 0.2,
                        zorder = 10)
        ax[1].vlines(x = datapoints_x, 
                     ymin = datapoint['mut_rate_ci95'][0],
                     ymax = datapoint['mut_rate_ci95'][1],
                     linestyle='-', 
                     lw = 1,
                     alpha = 0.6, 
                     color = mol_clock_col['Alternative'],
                     zorder = 8)
ax[1].set_ylim(0, None)
[ax[idx].spines['top'].set_visible(False) for idx, a in enumerate(ax)]
[ax[idx].spines['right'].set_visible(False) for idx, a in enumerate(ax)]
plt.tick_params(labelsize=fontsize)
plt.tight_layout()
plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2026_03_05_tmrca_nonhypermut_patho_clades_min8isolates_lit_clocks_excEhclock_wclocks_10kMCsim.svg'))

print(tmrca_spec_inf_clades_nohypmut_unsel[['Strain', 'mut_rate', 'mut_rate_ci95']].to_string(index = False))
