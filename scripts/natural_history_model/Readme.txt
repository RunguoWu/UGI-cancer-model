This folder contains code files for the UGI-CPDI model.

functions are saved in the four file starting with "fn_".

functions in "fn_tp_recalibration.R" are not used now, as the optimisation results look good enough.

"fn_parameter_search.R" is the wrapper function to search the optimal parameters.

Currently, the master file "master_tp_optimise_by2Char.R" is used. This is for optimisation by age group over 70 or not and sex, separately.

Other optimisation master files are saved in the archive folder.

"validation_hpc_by2Char.R" is for internal validation in patients without red flag at index date.

"validation_in_high_risk.R" is for external validation in patients with imaging or 2 week wait recommended at index date.

The "job_scripts" folder has similar master files, adapted to High Performance Computer
