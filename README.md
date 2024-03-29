This repository contains files for the simulation study in the 2021 *Biometrics* article "Inverse Probability Weighted Estimators of Vaccine Effects Accommodating Partial Interference and Censoring" by Chakladar et al. [Link](https://onlinelibrary-wiley-com.libproxy.lib.unc.edu/doi/full/10.1111/biom.13459), and the full reference is below. 

Simulations were conducted in R on a Slurm-managed High Performance Computing (HPC) cluster at UNC-Chapel Hill, and the shell scripts used are specific to Slurm. Each .sh shell script performs a set of simulations under a specific set of parameters. For instance, the script m500\_ni10.sh performs simulations for m=500 groups and group size n\_i=10. Once the three aforementioned R scripts are placed on an HPC cluster, the Slurm command "sbatch m500_ni10.sh" runs a set of 1,000 simulations. All files, including shell scripts and R scripts, should be placed in the same user home directory. 

There are three R scripts that run the simulations. Before beginning, users need to **edit two global variables** in each script: 'libs' and 'user_home_directory', specifying where the user is storing (1) R libraries and (2) all R scripts and shell files on the HPC cluster. The three scripts are:
* dataset_gen.R: Generates all simulated datasets according to the parameter values passed in from SLURM and stores them as .csv files in a "base_datasets" subdirectory. Also creating a single dataset of the true estimand values and stores it in a .csv file.  
* results_gen.R: Reads in a single simulated dataset and computes the inverse probability of censoring weights (IPCW) point estimates and asymptotic standard errors (ASEs). Writes a .csv file of results to a "results" subdirectory.
* analyse_results.R: Reads in all results and writes two results tables that can be copied-and-pasted into LaTeX, one each under treatment a=1 and control a=0, to an .Rout file. 

Shell scripts corresponding to different simulation settings can be found in the different folders, and each script sets variables for the different simulation settings. In some cases, variations on the three R scripts are used to implement different scenarios. Results generation can take hours for a single simulated dataset, so the array method in Slurm is used to run a set of 1,000 simulations serially on different computing nodes. 

A LaTex version of the manuscript supplement, figures for the manuscript, and an R script 'plots.R' that generates the figures can be found in the 'supplement' directory. The directory 'R Package Tarballs' contains .tar.gz files of some R packages that are used, which may need to be loaded into the 'libs' directory on HPC clusters. 

The code is accessible at https://github.com/samrosin/interference-IPCW, where users are encouraged to submit any issues. For other questions, please email srosin@live.unc.edu. 

Article reference: Chakladar, S., Rosin, S., Hudgens, M. G., Halloran, M. E., Clemens, J. D., Ali, M., & Emch, M. E. (2021). Inverse probability weighted estimators of vaccine effects accommodating partial interference and censoring. *Biometrics* (2021)
