# acyeh-lab/2023/TCR-simulation

Scripts are used to calculate "selective expanders" (TCRs from a single donor pool that differentially expand between 2 recipients) in Yeh AC et al, 2023.
Please read below for overview of how to run scripts

# Data preparation
Input files needed to run scripts include:
1) donor TCR sequences
2) recipient 1 TCR sequences
3) recipient 2 TCR sequences

File structure for each of the above must contain the following named columns:
1) "rearrangement" - a unique string containing CDR3 nucleotides (e.g. "CTCTGCAGCCTGGGAATCAGAACGTGCGAAGCAGAAGACTCAGCACTGTACTTGTGCTCCAGCAGTCAAAGGGGTGACACCCAGTAC")
2) "templates" - is the number of counts identified from sequencing (e.g. "11")
3) "frame_type" - is either "In" or "Out" - all frame_types set to "Out" are omitted from calculation (this column is derived from Adaptive's TCRB (v3) assay). If want to use all TCRs in file, set all templates to "frame_type" equal to "In"

# Scripts
There are 2 scripts provided (run sequentially):
1) Priors Generation Script: This is used to generate a log-log plot of TCR count vs. clonality of the donor pool, which is used as a Bayesian prior.  The output of this file named "exp_coeff_matrix.csv" will be used as input for the next Script.
2) Selective Expanders Script: This is used to calculate a p-value for each TCR specified from the input files used in the manuscript. Significant p-values represent clones that are "selectively expanded". The output of this file will have the suffix "- final probabilities.csv".

Notes:
1) The first script should run fairly quickly on a typical computing desktop (a few minutes for a file with ~1M naive TCRs)
2) The second script typically takes several hours (and potentially up to several days for pools larger clone sizes) when run on donor and recipient files with ~100-200k TCR counts.
3) The second script is  described in Supplementary Materials and Methods - "Probabilistic Simulation of Post-transplant TCR Expansion"

# Quick Start
Let's first download donor data file (see Zenodo link under samples folder). This should include the following donor file: `1901T-donor-spleen-all-supermerge.csv`.

Next, download the script "[Cluster] Publication - Priors Generation.R", change the working directory in the script to where the donor data file is located, and make sure the proper donor file is referred to:
```
setwd("C:/...") #Specify file directory
donor_file <- fread("1901T-donor-spleen-all-supermerge.csv") #Specify file name
```
The output of the above file should be named `exp_coeff_matrix.csv` and represents the coefficients for the log-log plot of the TCR clonality sizes vs. frequency. This is necessary to generate a "prior probability" for downstream Bayesian analysis.

Now, let's download 2 recipient data file (see Zenodo link under samples folder). To begin with, we can download the following files: `1901T-Gp1-1-spleen.tsv` and `1901T-Gp1-2-spleen.tsv`.  Place these two files in the same directory as `exp_coeff_matrix.csv` and '1901T-donor-spleen-all-supermerge.csv'.

Next, download the script "[Cluster] Publication - Selective Expanders.R", change the working directory in the script to where the donor data file is located, and make sure the proper donor file is referred to in the script:
```
analysis_dir <- "/..."
outputdir <- paste0(analysis_dir,"/Analysis") # Output files will be sent to "/Analysis" subdirectory
```

Next, make sure under the "File input requirements" subheading, the appropriate file are referred to.  Of note, we have to specify the donor file and two recipient files of interset to compare.
```
## File input requirements: 
## 1) Sequenced donor pool
donor_pool <- fread("1901T-donor-spleen-all-supermerge.csv")  
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type", where:
# 1) "rearrangement" is a unique string containing CDR3 nucleotides (e.g. "CTCTGCAGCCTGGGAATCAGAACGTGCGAAGCAGAAGACTCAGCACTGTACTTGTGCTCCAGCAGTCAAAGGGGTGACACCCAGTAC")
# 2) "templates" is the number of counts identified from sequencing (e.g. "11")
# 3) "frame_type" is either "In" or "Out" - we omit all frame_types that are set to "Out"

## 2) Sequenced recipient pool 1
recipient1 <- read.table(file = "1901T-Gp1-1-spleen.tsv", sep = '\t', header = TRUE)  
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type" as above

## 3) Sequenced recipient pool 2
recipient2 <- read.table(file = "1901T-Gp1-2-spleen.tsv", sep = '\t', header = TRUE)  
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type" as above

## 4) Coefficient matrix: (Generated from donor sequences using the script "[Cluster] Publication - Priors Generation")
# This file contains output for 3 variables: "a","b","x_c" as seen in Figure 1B that describe the rate of clone size decay y=a*x^(-b), with x_c being the x-intercept
coefficient_matrix <- "exp_coeff_matrix.csv"
```

Next, we have several script parameters that can be modified.  The cutoffs below reflect what is used in Yeh AC, et al, Figure 3D (B6 to B6D2F1 twin transplant systems):
```
## Numerical Variables (adjust below):
## 1) Donor pool size: - number of donor T-cells donated to each recipient (adjust for each experiment / recipient pair)
D1 <- 2000000
D2 <- 2000000

## 2) Combined recipient TCR clone count threshold - default set to 10 (see Figure S4)
TCR_cutoff <- 10

## 3) Donor TCR count probability threshold - default cutoff set to 0.001 (see Figure S4)
numerical_cutoff <- 1e-3

## 4) Negative binomial PMF variance - default set to 2 (see Figure S5)
PMF <- 2
```

We can also modify the intermediate output steps below.  Can set both to "FALSE" to streamline disk space usage; in this case, only the final output file with all p-values for each TCR clonotype will be returned.
```
## Debugging output files (assumes large amounts of disk space)
output_priors <- FALSE # This will output .csv files of p(N), p(B), p(B|N) and p(N|B)
output_clone_files <- TRUE # This will output .csv file of background clone probabilities
```

Now we can run the script "[Cluster] Publication - Selective Expanders.R" in the directory.
The final output file should be labeled "Exp1 - final probabilities.csv", which contains p-values for each TCR clonotype that represents the probability that the particular clonotype expanded similarly between the two recipients given.

An example of the output file is as follows:
|experiment |donor_count |donor_total |recipient1_count	|recipient1_total	|recipient2_count	|recipient2_total	|rearrangement	|final_prob |
|-----------|------------|------------|-----------------|-----------------|-----------------|-----------------|---------------|-----------|
|Exp1	|2	|652099	|9739	|171490	|129	|103692	|ACTCTGAAGATCCAGAGCACGCAACCCCAGGACTCAGCGGTGTATCTTTGTGCAAGCAGCTTAGCTGGGGGGGACTATGAACAGTAC	|1.04E-28
|Exp1	|0	|652099	|1	|171490	|216	|103692	|TCTCTCACTGTGACATCTGCCCAGAAGAACGAGATGGCCGTTTTTCTCTGTGCCAGCACCGTACTGGGGGGGCGTCATGAACAGTAC	|2.17E-16
|Exp1	|0	|652099	|1	|171490	|136	|103692	|TCTCTCACTGTGACATCTGCCCAGAAGAACGAGATGGCCGTTTTTCTCTGTGCCAGCAGTCGGGACTGGGGGATATATGAACAGTAC	|1.51E-11
|Exp1	|26	|652099	|0	|171490	|61	|103692	|CTACTTTTACATATATCTGCCGTGGATCCAGAAGACTCAGCTGTCTATTTTTGTGCCAGCAGCCAAGACTGGGGGTATGAACAGTAC	|1.05E-10

The "count" columns reflect clonotype count, "total" columns represent total number of TCRs from given sample, and "final_prob" columns gives the final calculated p-values for each TCR clonotype.

