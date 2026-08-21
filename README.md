# Hachisuka_EEG_2026
-----

This repository includes codes and scripts for analyses and figures for the paper titled: 

Probing the neural geometry of conscious object perception in time


### *Behavioral Experiment*

Fig. 2A-C & Supplementary Fig. 2: Behavior analysis and plotting
- `~/Fig2/main_behavior_script.m`

Fig. 2D-F:
- `~/Fig2/plotsuperandbasic_byPAS.m`

Supplementary Fig. 2:
- `~/Fig2/behav_split_by_subcategory.m`

Statistics:  
- `~/Fig2/Statistics/HR_supervsbasic_30subs.jasp`
- `~/Fig2/Statistics/PASbycategoryaccuracy.jasp`
- `~/Fig2/Statistics/PASbycategoryaccuracy.jasp`
- `~/Fig2/Statistics/Main_3wayANOVA_shapirowilk.R`
- `~/Supplementary_Fig2/Statistics/percentBasicCorrect_byDetect_or_Rec.jasp`
- `~/Supplementary_Fig2/Statistics/percentSuperCorrect_byDetect_or_Rec.jasp`
- `~/Supplementary_Fig2/Statistics/Basic_BySubcategory.jasp.jasp`
- `~/Supplementary_Fig2/Statistics/Super_BySubcategory.jasp`

-----

### *Superordinate & Basic Multivariate Decoding Analysis*

Fig. 3A-C, G & Supplementary Fig. 3A-B: Timeseries decoding analysis

Subject level: 
- `~/Fig3/contentdecoding_single.m`
  
Group level: 
- `~/Fig3/contentdecoding_superbasicexemp_group.m`
  
Bootstrapping decoding onset/peak times:
- `~/Supplementary_Fig3/bootstrap_peakonset_groupfinal_main.m`

Fig. 3D-F: Searchlight decoding analysis

Subject level: 
- `~/Fig3/contentdecoding_searchlight_single.m`
  
Group level: 
- `~/Fig3/contentdecoding_searchlight_group.m`

Fig 3E-F & Supplementary Fig. 3: Timeseries cross-condition decoding analysis

Subject level: 
- `~/Supplementary_Fig3/crossdecoding_single.m`
  
Group level: 
- `~/Supplementary_Fig3/crossdecoding_group.m`
  
Bootstrapping decoding onset/peak times:
- `~/Supplementary_Fig3/bootstrap_peakonset_groupfinal_crosscond.m`

Supplementary Fig. 3: Exemplar decoding, group-level
- `~/Supplementary _Fig3/contentdecoding_exemponly_group.m`

Supplementary Fig. 4: Bayesian statistics for undetected timeseries decoding
- `~/Supplementary_Fig4/nullresult_bayesian_stats.m`

Supplementary Fig. 6: Time-generalization analysis

Subject level: 
- `~/Supplementary_Fig6/timegen_single.m`
  
Group level: 
- `~/Supplementary_Fig6/timegen_group.m`

Supplementary Fig. 7: Image position decoding analysis

Subject level: 
- `~/Supplementary_Fig7/positiondecoding_single.m`
  
Group level: 
- `~/Supplementary_Fig7/positiondecoding_group.m`

-----

### *Cross-category generalization performance (CCGP) and Parallelism Score (PS) analysis*

Fig 4B-D & Supplementary Fig. 4C-D: CCGP timeseries decoding analysis.

Subject level:
- `~/Fig4/ccgp_single.m`
  
Group level:
- `~/Fig4/ccgp_group.m`
  
Bootstrapping decoding onset/peak times:
- `~/Supplementary_Fig3/bootstrap_peakonset_groupfinal_ccgp.m`

Fig. 4F: Parallelism Score (PS) analysis & plotting

Subject-level: same script as CCGP

Group-level:
- `~/Fig4/calculatePS.m`

Supplementary Fig. 4: Bayesian statistics for CCGP and PS

- `~/Supplementary_Fig4/nullresult_bayesian_stats.m`

Supplementary Fig. 8: CCGP searchlight analysis

Subject-level: 
- `~/Supplementary_Fig8/ccgp_searchlight_single.m`
  
Group-level: 
- `~/Supplementary_Fig8/ccgp_searchlight_group.m`

-----

### *Cross-image Location Decoding Analysis*

Supplementary Fig. 5

Subject-level: 
- `~/Supplementary_Fig5/categorydecoding_crosspositions_single.m`

Group-level: 
- `~/Supplementary_Fig5/categorydecoding_crossposition_group.m`

Subject-level: 
- `~/Supplementary_Fig5/CCGPdecoding_crossposition_single.m`
- `~/Supplementary_Fig5/CCGPdecoding_crossposition_group.m`

-----

### *Distance-to-bound analysis*

Fig. 5C-E: Distance-to-bound analysis

Subject-level: 
- `~/Fig5/distancetobound_contentdecoding_bytrial.m`
  
Group-level: 
- `~/Fig5/plot_distance_to_pas_bytrial_final.m`

Bootstrapping decoding onset/peak times: 
- `~/Fig5/dist2bound_acc_bootstrap.m`

Fig. 5B-C: Plotting distance values
- `~/Fig5/plot_distance_to_pas_bytrial_final.m`

Fig. 5E: CLMM in R
- `~/Fig5/dist2bound_clmm_analysis.R`

Permutation test:
- `~/Fig5/dist2bound_run_permutation1.R`
- `~/Fig5/dist2bound_run_permutation2.R`

Group-level CLMM visualization: 
- `~/Fig5/plot_CLMM_results.m`

Supplementary Fig. 9: Distance-to-bound analysis control; predicting categorization performance

Subject level: 
- `~/Supplementary_Fig9/distancetobound_contentdecoding_bytrial_pasrating1.m`
  
Group level: 
- `~/Supplementary_Fig9/plot_distance_to_pas_bytrial_pasrating1.m`
- `~/Supplementary_Fig9/dist2bound_Lrforcorr_analysis.R

-----

### *Distance-to-bound mediation analysis*

Reformatting data for mediation analysis:
- `~/Fig6/combine_super_basic_csv_formediation.m`

Fig. 6B-D: Mediation analysis
- `~/Fig6/dist2bound_mediation.R`

-----

### *Behavioral Experiments 1 & 2 (OSOM behavioral experiments)*

Supplementary Fig. 1: Analysis & Plotting
- `~/Supplementary_Fig1and2/OSOM_behav.m`

Statistics:  
- `~/Supplementary_Fig1/Statistics/Exp1_SuperBasicHR.jasp`
- `~/Supplementary_Fig1/Statistics/Exp2_SuperBasicHR.jasp`
- `~/Supplementary_Fig1/Statistics/Exp1_PAS.jasp`
- `~/Supplementary_Fig1/Statistics/Exp2_PAS.jasp`
- `~/Supplementary_Fig1/Statistics/OSOM_shapiro_wilk_test.R`


-----

### *EEG preprocessing*

Preprocessing scripts:
- `~/Preprocessing/preprocessing_finalV2.m`
- `~/Preprocessing/interpolatebadchan.m`
- `~/Preprocessing/runICA.m`
- `~/Preprocessing/postICA_finalcleanV2.m`
- `~/Preprocessing/final_fixtrialtrigger.m`


Prepare data for multivariate decoding analysis with CoSMoMVPA:
- `~/Preprocessing/contentdecoding_preparedataset_exemplar_trialinfo.m`

ACTICAP layout:
- `~/Preprocessing/layout_acticap128chan.mat`


-----

### *toolboxes*
Fieldtrip: https://www.fieldtriptoolbox.org/

CoSMoMVPA: https://www.cosmomvpa.org/

Custom scripts: 
Replace the following scripts in `~/CoSMoMVPA/mvpa`:
- `~/custom_toolboxes/cosmo_crossvalidate.m`
- `~/custom_toolboxes/cosmo_crossvalidation_measure.m`
- `~/custom_toolboxes/cosmo_check_dataset.m`
- `~/custom_toolboxes/cosmo_searchlight.m`
  
These scripts are modified to output additional variables required for the analysis.

For CCGP and PS analysis, instead of `cosmo_classify_lda.m`, use:
- `~/custom_toolboxes/cosmo_classify_lda3.m`

For distance-to-bound analysis, replace cosmo_classify_lda.m with:
- `~/custom_toolboxes/cosmo_classify_lda4.m`

For cross position decoding analysis, also use:
- `~/custom_toolboxes/cosmo_searchlight_ccgp.m`
- `~/custom_toolboxes/cosmo_stack_ccgp.m`
- `~/custom_toolboxes/cosmo_searchlight_dist.m`
- `~/custom_toolboxes/cosmo_stack_dist.m`

Bayes Factors for timeseries decoding toolbox: https://github.com/LinaTeichmann1/BFF_repo
