The simulations can in principle be run on a single computer, which here is referred to as the 'local' approach, though running a full 1,000 simulations may take at least a week to run. Using cluster computing resources where possible, where a full simulation can be run in several hours hours, is recommended. Files implementing both the local computer approach and the cluster computing approach are included. Please note that the cluster computing approach requires access to a high-performance computing cluster, and the scripts enclosed may need adaptation to a given cluster's specifications.

## Local Computer Simulation Code
Five files are included in the subdirectory local\_approach. The first four files contain only helper functions and should be run in their entirety:

* (1) sim\_dataset.R simulates the base datasets
* (2) calc\_estimands.R calculates the true values of the estimands
* (3) ase\_est.R estimates the asymptotic standard errors (ASEs) and is a helper function called by (4). 
* (4) ipcw\_eval.R estimates the IPCW point estimates and ASEs

Once files (1) through (4) are run, the user then runs the simulations using file (5), sim_local.R. Three simulation parameters must be set in this file: 

* n\_i, the number of individuals in each cluster
* m, the number of clusters
* n\_sims, the number of simulated datasets to generate

Upon running file (5) output tables of the estimator values for both a=0 and 1 are generated. 

Additionally, the parameter time.until is the time until the event of interest (days until contracting cholera, in our manuscript). This is set by default to 100, which was used throughout the manuscript simulations, but could be set to any value if different estimands are of interest. Similarly, the Bernoulli allocation probabilities of interest can be set with the parameter alphas. 

## Cluster Computing Simulation Code

First, again please note that the below instructions are most relevant for users that have access to a high-performance research computing cluster that is managed using Slurm. Users without such access would need to adapt the enclosed scripts in order to perform such a simulation. Additionally, these simulations were performed on the Longleaf cluster at the University of North Carolina at Chapel Hill, and some of the Slurm commands used in the .sh files may be specific to UNC and Longleaf's setup.

Eight files are included in the subdirectory cluster\_approach: four R scripts, and four .sh shell scripts that are each coupled to an R script and give Slurm instructions for running the R script. The example given is for the simulation with group size n_i=50, m=500 groups, and 1,000 simulations, corresponding to Web Table 7 in the manuscript, but the same approach could be taken for any combination of simulation parameters. 

Before running the files, a user must set up a home directory with the following structure:

* base_datasets
 * base\_datasets\_ni50\_m500
* truthtables
    * truths\_ni50\_m500
* results_datasets
    * results\_ni50\_m500
* output\_tables
* RLibs

RLibs is an optional subdirectory containing R packages that the cluster does not have pre-loaded.

Also, in each of the four R scripts, the user should set the value of user_home_directory to be the string corresponding to their home directory on their computing cluster. 

The user can then run the following commands from the command line to obtain final results. The user should wait for each step to fully complete (in the case of steps 1 and 3, this may take hours) before running the next step. 

* (1) sbatch --array=1-1000 dataset\_gen\_ni50\_m500.sh
* (2) sbatch truthtable\_gen\_ni50\_m500.sh
* (3) sbatch --array=1-1000 results\_gen\_ni50\_m500.sh
* (4) sbatch analyse\_ni50\_m500.sh

Upon running these commands, the file output\_tables\_ni50\_m500.Rout will appear in the output\_tables subdirectory. This .Rout file contains two results tables, one each for a=0 and a=1, which can be copied and pasted into LaTeX. 

### Further details on cluster computing
Details on what each of the four R scripts does are as follows:

* (1) dataset\_gen\_ni50\_m500.R simulates a single dataset. The array method is used to do this in a distributed fashion, so that 1,000 datasets are separately simulated on different cores and then stored in the base\_datasets\_ni50\_m500 subdirectory. They are named base\_1.csv, base\_2.csv, etc. The true estimand values are also computed and placed in subdirectory truths\_ni50\_m500 with names truth\_1.csv, truth\_2.csv, etc. 

* (2) truthtable\_gen\_ni50\_m500.R simply reads in the 1000 different datasets truth\_1.csv, truth\_2.csv, etc. and creates a single long-format dataset truthtable\_ni50\_m500 with all true estimand values for later use. 

* (3) results\_gen\_ni50\_m500.R reads in a single simulated dataset and computes the IPCW point estimates and ASEs. Like in step (1), this is done using the array method so that results can be generated near-simultaneously for, e.g., 1,000 different simulations. 

* (4) analyse\_ni50\_m500.R reads in the 1,000 different results files and outputs a file output\_tables\_ni50\_m500.Rout in the output\_tables subdirectory. This .Rout file contains two results tables, one each for a=0 and a=1, which can be copied and pasted into LaTeX. 

