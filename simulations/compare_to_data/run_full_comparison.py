import sys
sys.path.append('../')
import pandas as pd
import os
import manifest as manifest

from simulations.compare_to_data.age_incidence_comparison import compute_inc_LL_for_all_sites, \
    plot_incidence_comparison_all_sites
from simulations.compare_to_data.age_prevalence_comparison import compute_prev_LL_for_all_sites, \
    plot_prevalence_comparison_all_sites
from simulations.compare_to_data.infectiousness_comparison import compute_infectious_LL_for_all_sites, \
    plot_infectiousness_comparison_all_sites
from simulations.compare_to_data.parasite_density_comparison import compute_parasite_density_LL_for_all_sites, \
    plot_density_comparison_all_sites
from simulations.compare_to_data.age_gamotocyte_prevalence_comparison import compute_gametocyte_prev_LL_for_all_sites, \
    plot_gametocyte_prevalence_comparison_all_sites

from simulations.compare_to_data.no_blood_comparison import compute_dead_LL_for_all_sites

def compute_LL_across_all_sites_and_metrics(numOf_param_sets = 64):
    infectious_LL = compute_infectious_LL_for_all_sites(numOf_param_sets)
    density_LL = compute_parasite_density_LL_for_all_sites(numOf_param_sets)
    prevalence_LL = compute_prev_LL_for_all_sites(numOf_param_sets)
    #gametocyte_prevalence_LL = compute_gametocyte_prev_LL_for_all_sites(numOf_param_sets)
    incidence_LL = compute_inc_LL_for_all_sites(numOf_param_sets)
    dead_LL = compute_dead_LL_for_all_sites(numOf_param_sets)


    #density_LL_w=density_LL
    #density_LL_w['ll'] = [float(val)/10 for val in density_LL['ll']]
    #print("ILL")
    #print(incidence_LL)
    combined = pd.concat([infectious_LL, incidence_LL, density_LL, prevalence_LL, dead_LL])
    
    weighting_rules = pd.read_csv('/projects/b1139/basel-hackathon-2023/reference_datasets/objective_weights.csv')
    
    b=pd.merge(combined, weighting_rules,  how='left', left_on=['site','metric'], right_on = ['site','metric'])
    b['weight'].fillna(0.1, inplace=True)
    b['ll'] = b['ll'] * b['weight']
    #print(b.to_string())
    #print(b.groupby("param_set").agg({"ll": lambda x: x.sum(skipna=False)}).reset_index().sort_values(by=['ll']).to_string())
    #print(b.groupby("param_set").agg({"weighted_LL": lambda x: x.sum(skipna=False)}).reset_index().sort_values(by=['weighted_LL']).to_string())
    
    #combined = pd.concat([density_LL])
    #print(combined.to_string())

    #fixme - Eventually, this will need to be fancier way of weighting LL across the diverse metrics/sites
    #fixme - For now, just naively add all log-likelihoods
    return b.groupby("param_set").agg({"ll": lambda x: x.sum(skipna=False)}).reset_index().sort_values(by=['ll'])

def plot_all_comparisons(param_sets_to_plot=None,plt_dir=os.path.join(manifest.simulation_output_filepath, "_plots")):
    plot_incidence_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_prevalence_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_gametocyte_prevalence_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_density_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_infectiousness_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir) 
    
def plot_all_comparisons2(param_sets_to_plot=None,plt_dir=os.path.join(manifest.simulation_output_filepath2, "_plots")):
    plot_incidence_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_prevalence_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_gametocyte_prevalence_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_density_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir)
    plot_infectiousness_comparison_all_sites(param_sets_to_plot=param_sets_to_plot,plt_dir=plt_dir) 

if __name__ == "__main__":
    #plot_all_comparisons2()
    
    ## if you are running directly from run_full_comparison you are going to probably want to 
    ## manually add a different default numOf_param_sets, for example, numOf_param_sets = 16
    ##
    print(compute_LL_across_all_sites_and_metrics())
