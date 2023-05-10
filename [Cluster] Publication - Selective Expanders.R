library(plyr)
library(data.table)

## This script calculates probability values for finding selective expanders as described in "Materials and Methods" section: 
## "Probabilistic Simulation of Post-transplant TCR Expansion" in Yeh Ac, et al.
## Note before running this script, need to generate donor TCR log-log plot coefficients for priors estimate (Bayesian step, part 3) - 
## This is run in a separate script and loaded in file inpu requirements part 4) below

## Experiment name:
experiment <- "run_08_17_22_Exp1_Gp4_5" # Used for output file, change as needed

## Set directories:
## Place all required files in the analysis_dir directory set below:
analysis_dir <- "/fh/fast/hill_g/Albert/TCR Sequences"
outputdir <- paste0(analysis_dir,"/Analysis_081722_Exp1_Gp4_5") # Output files will be sent to "/Analysis" subdirectory

setwd(analysis_dir)
message(paste0("Creating directory ", outputdir))
dir.create(outputdir)

## File input requirements: 
## 1) Sequenced donor pool
donor_pool <- fread("1901T-B6_45_1-donor-spleen-all.csv")  
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type", where:
# 1) "rearrangement" is a unique string containing CDR3 nucleotides (e.g. "CTCTGCAGCCTGGGAATCAGAACGTGCGAAGCAGAAGACTCAGCACTGTACTTGTGCTCCAGCAGTCAAAGGGGTGACACCCAGTAC")
# 2) "templates" is the number of counts identified from sequencing (e.g. "11")
# 3) "frame_type" is either "In" or "Out" - we omit all frame_types that are set to "Out"

## 2) Sequenced recipient pool 1
recipient1 <- read.table(file = "1901T-Day2-Gp0-spleen.tsv", sep = '\t', header = TRUE)  
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type" as above

## 3) Sequenced recipient pool 2
recipient2 <- read.table(file = "1901T-Day2-Gp1-spleen.tsv", sep = '\t', header = TRUE)  
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type" as above

## 4) Coefficient matrix: (Generated from donor sequences using the script "[Cluster] Publication - Priors Generation")
# This file contains output for 3 variables: "a","b","x_c" as seen in Figure 1B that describe the rate of clone size decay y=a*x^(-b), with x_c being the x-intercept
coefficient_matrix <- "1901T-B6_45_1-donor-spleen-all_coeff_matrix.csv"
coefficients <- as.numeric(fread(coefficient_matrix))
# For example, if a=0.75, b=3.03, x_c=82.0 (as in Fig 1b), then coefficients <- c(0.75,3.03,82.0)
#coefficients <- c(0.75,3.03,82.0) # Delete this if line if importing coefficient variables


## Numerical Variables (adjust below):
## 1) Donor pool size: - number of donor T-cells donated to each recipient (adjust for each experiment / recipient pair)
D1 <- 5000000
D2 <- 5000000

## 2) Combined recipient TCR clone count threshold - default set to 10 (see Figure S4)
TCR_cutoff <- 10

## 3) Donor TCR count probability threshold - default cutoff set to 0.001 (see Figure S5)
numerical_cutoff <- 1e-3

## 4) Negative binomial PMF variance - default set to 2 (see Figure S6)
PMF <- 2


## Debugging output files (assumes large amounts of disk space)
output_priors <- FALSE # This will output .csv files of p(N), p(B), p(B|N) and p(N|B)
output_clone_files <- FALSE # This will output .csv file of background clone probabilities


#---------------------#  
#----- FUNCTIONS -----#  
#---------------------#
message("Loading Functions")

merge_counts <- function(list1, list2)
  # Merges in-frame templates only - combines template frequencies
{
  list1_inframe_rearrangements <- as.character(list1$rearrangement[list1$frame_type=="In"])
  list1_inframe_template_counts <- as.numeric(list1$templates[list1$frame_type=="In"])
  
  list2_inframe_rearrangements <- as.character(list2$rearrangement[list2$frame_type=="In"])
  list2_inframe_template_counts <- as.numeric(list2$templates[list2$frame_type=="In"])
  
  intersect_list <- intersect(list1_inframe_rearrangements,list2_inframe_rearrangements)
  
  intersect_list_freq <- NULL
  
  if(length(intersect_list)==0)
  { 
    message("no intersections detected")
    rearrangement <- c(list1_inframe_rearrangements, list2_inframe_rearrangements)
    templates <- c(list1_inframe_template_counts, list2_inframe_template_counts)
  }
  if(length(intersect_list)!=0)
  {
    for(i in 1:length(intersect_list))
    {
      temp_count <- list1_inframe_template_counts[list1_inframe_rearrangements==intersect_list[i]] + 
        list2_inframe_template_counts[list2_inframe_rearrangements==intersect_list[i]]
      intersect_list_freq <- c(intersect_list_freq,temp_count)
    }
    
    rearrangement <- c(list1_inframe_rearrangements[!(list1_inframe_rearrangements%in%intersect_list)],
                       list2_inframe_rearrangements[!(list2_inframe_rearrangements%in%intersect_list)],
                       intersect_list)
    
    templates <- c(list1_inframe_template_counts[!(list1_inframe_rearrangements%in%intersect_list)],
                   list2_inframe_template_counts[!(list2_inframe_rearrangements%in%intersect_list)],
                   intersect_list_freq)
  }
  
  frame_type <- rep("In",length(rearrangement))
  return(data.frame(rearrangement,templates,frame_type))
}

merge_counts_2 <- function(list1, list2)
  # Merges in-frame templates only - this also gives template frequency in from each list, unlike "merge_counts" function above
{
  list1_inframe_rearrangements <- as.character(list1$rearrangement[list1$frame_type=="In"])
  list1_inframe_template_counts <- as.numeric(list1$templates[list1$frame_type=="In"])
  
  list2_inframe_rearrangements <- as.character(list2$rearrangement[list2$frame_type=="In"])
  list2_inframe_template_counts <- as.numeric(list2$templates[list2$frame_type=="In"])
  
  intersect_list <- intersect(list1_inframe_rearrangements,list2_inframe_rearrangements)
  
  intersect_list_freq_1 <- NULL
  intersect_list_freq_2 <- NULL
  
  if(length(intersect_list)==0)
  { 
    message("no intersections detected")
    rearrangement <- c(list1_inframe_rearrangements, list2_inframe_rearrangements)
    templates <- c(list1_inframe_template_counts, list2_inframe_template_counts)
    templates_1 <- c(list1_inframe_template_counts, rep(0,length(list2_inframe_rearrangements)))
    templates_2 <- c(rep(0,length(list1_inframe_rearrangements)),list2_inframe_template_counts)
  }
  if(length(intersect_list)!=0)
  {
    for(i in 1:length(intersect_list))
    {
      intersect_list_freq_1 <- c(intersect_list_freq_1,
                                 list1_inframe_template_counts[list1_inframe_rearrangements==intersect_list[i]])
      
      intersect_list_freq_2 <- c(intersect_list_freq_2,
                                 list2_inframe_template_counts[list2_inframe_rearrangements==intersect_list[i]])
    }
    
    rearrangement <- c(list1_inframe_rearrangements[!(list1_inframe_rearrangements%in%intersect_list)],
                       list2_inframe_rearrangements[!(list2_inframe_rearrangements%in%intersect_list)],
                       intersect_list)
    
    unique_list_1 <- length(list1_inframe_template_counts[!(list1_inframe_rearrangements%in%intersect_list)])
    unique_list_2 <- length(list2_inframe_template_counts[!(list2_inframe_rearrangements%in%intersect_list)])
    
    
    templates_1 <- c(list1_inframe_template_counts[!(list1_inframe_rearrangements%in%intersect_list)],
                     rep(0,unique_list_2),
                     intersect_list_freq_1)
    
    templates_2 <- c(rep(0,unique_list_1),
                     list2_inframe_template_counts[!(list2_inframe_rearrangements%in%intersect_list)],
                     intersect_list_freq_2)
    templates <- templates_1+templates_2
  }
  
  frame_type <- rep("In",length(rearrangement))
  return(data.frame(rearrangement,templates,templates_1,templates_2))
}

merger_donor_recip <- function(donor_list,recip_list)
  # Takes input from merge_counts_2 and adds template count from donor list for each clone
{
  index<-length(recip_list$rearrangement)
  donor_frequency<- NULL
  for(i in 1:index)
  {
    message(paste0("Merging ",i," of ", index, " TCRs"))
    sequence <- as.character(recip_list$rearrangement[i])
    sequence_freq <- 0
    if(sequence%in%donor_list$rearrangement)
    #If the recipient clone is found in the donor pool, update count, otherwise default is set to 0
    {
      sequence_freq <- donor_list$templates[as.character(donor_list$rearrangement)==sequence]
    }
    
    donor_frequency <- c(donor_frequency,sequence_freq)
  }
  
  rearrangement <- recip_list$rearrangement
  recip_templates <- recip_list$templates
  recip_templates_1 <- recip_list$templates_1
  recip_templates_2 <- recip_list$templates_2
  return(data.frame(rearrangement,donor_frequency,recip_templates,recip_templates_1,recip_templates_2))
}

calc_nbinom_parameters <- function(mean, variance)
  # Outputs the parameters n (# of successes desired) and p (prob success) from the negative binomial distribution
  # Given fixed mean and variance
  # note mean = n*(1-p)/p, and variance = n*(1-p)/(p^2) 
  # thus variance = mean/p, and p = mean/variance
  # returns n, p
{
  p <- mean/variance
  n <- mean*p/(1-p)
  return(c(n,p))
}



#---- Probability calculations----#
{
  
  # Each Part listed below roughly corresponds to the numbering used in Supplemental Materials and Methods - "Probabilistic Simulation of Post-transplant TCR Expansion"
  
  # Part 1: Definitions and Assumptions
  {
    message("Part 1: Loading definitions and assumptions")
    message(Sys.time())
    message(paste0("Donor pool size 1 = ",D1))
    message(paste0("Donor pool size 2 = ",D2))
    # "B" = background count (from donor sequencing results of the TCR)
    # "C1" and "C2" = counts in two recipient samples.  By definition, at least one of these will be nonzero, and C1+C2 >= 10 (or specified)
    # "T_0", "T_1R",  "T_2R", "T_3R" are the sizes of the background pool and recipient pools (i.e. total sequenced)
    # "D1" and "D2" is the total number of cells transplanted to corresponding recipients.
  
    T_0 <- sum(count(donor_pool$templates[donor_pool$frame_type=="In"])[,1]*
      count(donor_pool$templates[donor_pool$frame_type=="In"])[,2])
    
    #Now load recipient pool counts
    message("Defining recipient pool counts")
    {
      T1 <- sum(count(recipient1$templates[recipient1$frame_type=="In"])[,1]*
                count(recipient1$templates[recipient1$frame_type=="In"])[,2])
      T2 <- sum(count(recipient2$templates[recipient2$frame_type=="In"])[,1]*
                count(recipient2$templates[recipient2$frame_type=="In"])[,2])     
    }

    message(paste0("T_0: ", T_0))
    message(paste0("T1: ", T1))
    message(paste0("T2: ", T2))
  
    total_pool_count <- T_0+D1+D2
    message(paste0("Total_pool_count: ", total_pool_count))
    
    
    # We must loop from 1 to the max clonality count+1 (to account for "B+1"), which in this case is:
    max_donor_clonality_count <- max(count(donor_pool$templates[donor_pool$frame_type=="In"])[,1])
    
    message(paste0("Max_donor_clonality_count: ", max_donor_clonality_count))
    
    upper_end = floor(1.5*max_donor_clonality_count* total_pool_count / T_0)
    # Upper end represents the maximum estimated clone size of the total TCR pool (2*D+sequenced pool) based on the maximum seen in the sequenced pool
    # This variable is used in Part 3 below to define an upper-bound for calculations, as likelihood decreases exponentially. 
    # It is also denoted as the variable "L" in Supplemental Materials and Methods - "Probabilistic Simulation of Post-transplant TCR Expansion", part 3. 
    message(paste0("Upper_End: ", upper_end)) 
  }
  
  # Part 2: Select all the TCRs of interest between each pair of recipients
  {
    message("Part 2: Selecting TCRs")
    message(Sys.time())

    # we will now create a matrix of each specific TCR including the Final TCR counts between each pair
    # col 1 = TCR sequence, col 2 = frequency in group 1, col 3 = frequency in group 2
    # Also apply the TCR_cutoff and eliminate pairs with total TCR less than the cutoff point
    # Also Now we have lists of TCR sequences for each group that occur in at least 10 (or specified) instances
    # Now we expand that list so we also get the count from the initial donor sample
    
    Group_Merge <- merge_counts_2(recipient1, recipient2)
    Group_Merge_freq10 <- Group_Merge[Group_Merge[,2]>=TCR_cutoff,]
    final_donor_recip_data <- merger_donor_recip(donor_pool,Group_Merge_freq10)
    
    setwd(outputdir) #Change to output directory as defined as a global variable
    fwrite(final_donor_recip_data,"final_donor_recip_data.csv") # For recordkeeping
  }
  
  # Part 3: Estimate probability distribution over total clones in the merged pool (T0 + 2*D) either sequenced or transplanted
  {

    ## What we want is a probability distribution of the total number of clone cells for any given TCR.
    # This starts at B+1, given "B" = background count and our assumption that we only look at samples with at least 1 recipient cell
    # Truncate after awhile as probabilities decay rapidly
    # Estimated by hypergeometric distribution, to model background pool seqeucing, and a "prior" on clone size with shape set by background sequencing data
    # This prior has the effect of shifting the clone size distribution smaller than we would get by using the flat prior, thus making
    # p-values more conservative, since the founder pools are smaller and subject to fluctuation.
    # Could also use negative binomial distribution here (see below) to account for noise (over-dispersion) in the sequencing counts.
    
    # Bayesian statistics: we want p(N | B) over possible background counts, 
    # where B = sequenced pool
    # where N = actual starting clone size = merged pool (T0 + 2*D), starting at B+1 as at least one recipient must have had a TCR
    # Note p(N|B) = p(N) * p(B|N) / p(B).  Note p(N) is the prior based on the power law slope, which we estimate

    ## part 3.1: Calculating our PRIOR, p(N) 
    {
      message("Part 3-1: Generating prior probability p(N) for bayesian analysis")
      message(Sys.time())
      message("Coefficients preloaded") #Generated from another script
      message(paste0("a=",coefficients[1]))
      message(paste0("b=",coefficients[2]))
      message(paste0("x_c=",coefficients[3]))    
      
      message(paste0("Calculating prior vector from 1 to ",upper_end))
      p_N <- NULL
      for(i in 1:upper_end) #Note that p_N starts at 1, as only look at clones with at least 1 recipient (so p_N must be > 0)
      {
        p_N <- c(p_N,coefficients[1]*(i^-coefficients[2]))
      }
      if(output_priors)
      {
        write.csv(p_N,"p(N).csv",row.names=paste0("N=",1:upper_end))
      }
      
      message(paste0("sum(p_N):",sum(p_N))) #This should be as close to 1 as possible
    }
  
    ## Part 3.2: Calculating p(B|N)
    {
      message("Part 3-2: Calculating p(B|N)")
      message(Sys.time())

      # dhyper(x,m,n,k) gives probability mass function of observing "x" successes
      # m = number of success states in the population
      # n = number of failure states in the population
      # k = number of draws
      # x = number of observed successes
      # m <- will vary - this is the number of "true success" states in the population 
      
      (N <- total_pool_count) # Total pool including sequenced pool + 2* donor pool (for each mouse)
      (k <- T_0) # Number of donor TCRs sequenced (a.k.a. "draws")
      
      p_BgivenN <- NULL # We will store results p(B|N) in this matrix
      for(m in 1:upper_end) # This is setting the N for p(B|N).  Limit is "upper_end", arbitrarily defined - should refine for final paper
      {
        message(paste0("Calculating p(B|N): ", m," of ", upper_end, " done"))
        prob <- rep(0,upper_end+1) #We make default value 0 so the matrix is even
        #prob_adj <- rep(0,upper_end+1)
        
        for(x in 0:m) 
          # Must start at 0, as we can have 0 show up in the sequenced pool even if the total pool can't be 0.
          # Note that x cannot be more than m, so prob will be 0 in that case, but we need to output the 0, hence setting the default value to 0 above
        {
          value <- dhyper(x,m,N-m,k) 
          prob[x+1] <- value 
        }
        p_BgivenN<-cbind(p_BgivenN,prob) 
      }
      
      colnames(p_BgivenN)<-paste0("N=",1:upper_end)
      if(output_priors)
      {
        write.csv(p_BgivenN,"p(B|N).csv",
                  col.names = TRUE,
                  row.names = paste0("B=",0:upper_end))
      }
      message(paste0("p(B|N) for B=0, N=1: ",p_BgivenN[1,1]))
    }
    
    ## Part 3.3: Calculating p(B)
    {
      message("Part 3-3: Calculating p(B)")
      message(Sys.time())
      
      # Now calculate p(B) for each B
      p_B <- rep(0,upper_end+1)
      
      for(i in 0:upper_end) # Loops to the maximum TCR clone size observed - for each "B"
      {
        for(j in 1:upper_end) #upper bound set by the # of .csv files we put out
        {
          temp <- p_BgivenN[i+1,j]*p_N[j] # representing p(B|N) * p(N) for each N
          p_B[i+1]<-p_B[i+1] + temp
        }
      }
      if(output_priors)
      {
        write.csv(p_B, "p(B).csv", row.names=paste0("B=",0:upper_end))
      }
    }
    
    ## Part 3.4: Calculating p(N|B) = p(N)
    {
      # We need to put together Bayes' Theorem: p(N|B) = p(N) * p(B|N) / p(B)
      p_NgivenB <- matrix(0,nrow=upper_end,ncol=max_donor_clonality_count+1) # Set the max background count to the actual max background count; +1 as you can have "0"
      
      for(N in 1:upper_end)
      {
        for(B in 0:max_donor_clonality_count){
          p_NgivenB[N,B+1] <- p_N[N] * (p_BgivenN[B+1,N] / p_B[B+1])
        }
        
      }
      colnames(p_NgivenB)<-paste0("B=",0:max_donor_clonality_count)
      if(output_priors)
      {
        write.csv(p_NgivenB,"p(N|B).csv",
                  col.names = TRUE,
                  row.names = paste0("N=",1:upper_end))
      }
    }
  }
  
  # Part 4/5: Calculate an expansion rate (Exp_R) for each observed clone in a given pair of mouse
  # Exp_R = F_final / F_start
  # Where f_start = starting clone frequency in both mice =
  # (S-B)/(2*D) where S = starting clone size from the merged donor pool, B = donor sequencing count, D = # of cells donated to a mouse
  # S-B represents the clone size that is transplanted to the mouse
  # And f_final = final clone frequency in the donor population = 
  # (C1+C2)/(T1+T2) where C1,C2 represent TCR count in the recipient sample and T1,T2 represent the # of T-cells sequenced
  {
    message("Part 4: Calculating expansion rate")
    message(Sys.time())
    
    pair <- final_donor_recip_data # We have to iterate through each pair of mice
    
    all_names <- as.character(pair$rearrangement)
    all_S <- list(as.character(pair$rearrangement))
    all_prob_S <- list(as.character(pair$rearrangement))
    all_S_minus_B <- list(as.character(pair$rearrangement))
    all_recip_1_count <- list(as.character(pair$rearrangement))
    all_recip_2_count <- list(as.character(pair$rearrangement))
    all_exp_r <- list(as.character(pair$rearrangement))
    
    for(i in 1:length(all_names))
    {
      message(paste0("Calculating expansion rates: ",i," of ", length(pair$rearrangement)))
      
      # Here we calculate f_final = (C1+C2)/(T1+T2) and f_start = (S-B)/(2*D)
      name<-as.character(pair$rearrangement[i])
      
      #----- calculating f_final -----# 
      final_clone_size <- pair$recip_templates[i] # This is C1+C2
      f_final <- final_clone_size / (T1+T2)
      #-------------------------------#
      
      B<-pair$donor_frequency[i]
      S_vector<-p_NgivenB[,B+1]
        
      
      ##### This part calculates the cutoff to set below which the probabilities is deemed too small to affect our overall calculations; 
      # cutoff set via global variable
      index_cutoff <- max(which(S_vector>numerical_cutoff)) 
      # Issue here is if there is no cutoff assigned (nothing > numerical cutoff), we should just take the max value (peak)
      if(index_cutoff=="-Inf")
      {
        index_cutoff <- which.max(S_vector)
      }
      message(paste0("index cutoff: ", index_cutoff))
      ##### ----- #####
      
      S <- 1:index_cutoff
      S_minus_B <- S-B # a vector
      prob_S <- S_vector[1:index_cutoff] 
      
      ##Need to make this adjustment so that prob_S is shifted to only be available for positive "S_minus_B's"
      ##We will then normalize prob_S so that the sum = 1
      {
        prob_S[match(0,S_minus_B)]<-0 #If no cells are donated, than we have to set that probability to 0 as we only take T-cells that have at least 1 recipeint clone
        prob_S <- prob_S/sum(prob_S) 
      }
      
      #Note that prob_S is a probability distribution
      f_start <- (S-B)/(D1+D2)
      exp_r <- f_final/f_start
   
      recip_total <- final_clone_size
      recip_1_count <- pair$recip_templates_1[i]
      recip_2_count <- pair$recip_templates_2[i]

      ##### Output file here #####
      temp<-data.frame(S,prob_S,S_minus_B,
                       recip_1_count,
                       recip_2_count,
                       exp_r)
      # The resulting output is a  file for each TCR clone that appears in either of the two recipient pools of interest with 3 columns:
      # 1) Each possible total starting clone size "S"
      # 2) The probabiltiy of that starting clone size
      # 3) S-B (1)-(2)
      # 4) Sequenced recipient pool 1 clone size 
      # 5) Sequenced reicpient pool 2 clone size
      # 6) The calculated expansion rate with that clone size (Exp_R)
      ######################################
      
      all_S[[i]] <- S
      all_prob_S[[i]] <- prob_S
      all_S_minus_B[[i]] <- S_minus_B
      all_recip_1_count[[i]] <- recip_1_count
      all_recip_2_count[[i]] <- recip_2_count
      all_exp_r[[i]] <- exp_r

      if(output_clone_files)
      {
        fwrite(temp,paste0(experiment,"_clone_files - ",name,".csv")) #Write out for debugging purposes
      }
    }
  }
  
  # Part 6A: Iterating over all possible splits in the donor pool, perfomring binomial splits and calculate respective probabilities
  # In this part, we loop over all possible splits "A1" and "A2" of transplanted cells, with A1+A2 = (S-B).  We can use the
  # binomial distribution to assign probability of each split for each S-B. Then, for each asplit, we compute the final frequency
  # f1 = exp_r*A1/D and f2 = exp_r*A2/D using the same shared expansion rate (exp_r) and starting frequencies
  # Now assume C1 < C2 (choose mouse with smaller sequenced count).  Then we calculate the probability of seeing of seeing only C1 counts assuming
    # that we sequenced T1 and the actual frequency at the end is f1
  {
    message("Part 5: Performing binomial splits")
    message(Sys.time())
    
    for(m in 1:length(all_names)) # Loop through all TCRs
    {
      gc()
      message(paste0("Calculating binomial donor splits: ", m," out of ",length(all_names)," clones"))
      setwd(outputdir)
      
      dir.create(all_names[m]) #Make folders without the ".csv" file extension)

      #Now, for each "S_minus_B", we split "A1" and "A2", and we compute the final frequencies f1 = exp_r*A1/D, f2 = exp_r*A2/D
      start_index <- min(which(all_S_minus_B[[m]]>0)) #We can't have negative values or 0 for # of T-cells transplanted, so set boundary to >0
      if(start_index=="Inf"){
        start_index <- 1 # have to start somewhere
      }
      
      C1<-as.numeric(all_recip_1_count[[m]])
      C2<-as.numeric(all_recip_2_count[[m]])
      
      p_values <- NULL
      
      for(i in start_index:length(all_S_minus_B[[m]]))
      {
        temp <- NULL
        binominal_split_prob <- dbinom(0:all_S_minus_B[[m]][i],all_S_minus_B[[m]][i],0.5) #binomial distribution to calculate prob of A1/A2 split

        for(j in 0:(all_S_minus_B[[m]][i])) # There are N+1 ways to split this (e.g. if N=3, we have (0,3), (1,2), (2,1), (3,0))
        {

          rate<-all_exp_r[[m]][i]
          A1<-j
          A2<-all_S_minus_B[[m]][i]-j
          prob_A1A2_split<-binominal_split_prob[j+1]  #Have to add 1 because j starts at 0
          f1<-rate*A1/D1 # frequency of desired TCR for mouse 1 based on assumed rate and binomial split
          f2<-rate*A2/D2 # frequency of desired TCR for mouse 2 based on assumed rate and binomial split
          p_values <- NULL
          
          #----- Now calculate respective p-values for seeing C1 counts assuming C1/T1 < C2/T2 (pick the smaller mouse) -----#
          # Note that for this, one could use 
          # (1) hypergeometric distribution if we have the full donor pool size
          # (2) poisson distribution if sequencing were perfect, or 
          # (3) negative binomial distribution, which is like the poisson but has a larger variance 
          # (to account for extra noise in the read counts). 
          
          # Need to change variance here by setting Mu
          if(C1/T1 <= C2/T2) # We want to calculate on the one with less clones proportionally
          {
            if(A2==0) # If the bigger value is 0, not possible
            {
              p_values <- 0
            }else if(A1==0){
              if((C1*C2)>0){
                p_values <- 0  #If there is at least 1 count in the lower numbered recipient but it got "0" clones. 
              }else{
                p_values <- 1  #If there are no counts in the lower numbered recipient and it got "0" clones
              }
            }else{
              expected_mean <- f1*T1
              temp2<-calc_nbinom_parameters(expected_mean, expected_mean*PMF) #note if expected_mean = 0, this is undefined, hence above
              p_values <- pnbinom(C1,temp2[1],temp2[2])
            }

          }else{
            if(A1==0) # If the bigger value is 0
            {
              p_values <- 0
            }else if(A2==0){
              if((C1*C2)>0){
                p_values <- 0  #If there is at least 1 count in the lower numbered recipient but it got "0" clones
              }else{
                p_values <- 1  #If there are no counts in the lower numbered recipient and it got "0" clones
              }
            }else{
              expected_mean <- f2*T2
              temp2<-calc_nbinom_parameters(expected_mean, expected_mean*PMF) #note if expected_mean = 0, this is undefined, hence above
              p_values <- pnbinom(C2,temp2[1],temp2[2]) 
            }
          }
          
          temp<-rbind(temp,data.frame(all_S_minus_B[[m]][i],
                                      all_prob_S[[m]][i],
                                      A1, A2, prob_A1A2_split,
                                      rate, f1, f2, C1, C2, T1, T2,
                                      p_values))
        }

        setwd(paste0(outputdir,"/",all_names[m]))
        
        colnames(temp)<-c("S_minus_B","prob_S","A1","A2","prob_A1A2_split","rate","f1","f2", "C1", "C2", "T1", "T2", "p_values")
        
        #####  Output file here #####
        fwrite(temp,paste0("[",all_S_minus_B[[m]][i],"].csv"))      
        #############################
      }
      
    }
  }
  
  # Part 6B: Calculating final p-value
  # This part sums up all the probabilites and spits them in a csv file listing each clone
  {
    message("Part 6: Summing up probabilities")
    message(Sys.time())
    setwd(outputdir) 

    # This loop then calculates the final p_value for each clone, and outputs everything in one large file
    # Multiples each S_B by each binomial split probaiblity and by each final p-value
    # Note that number of files in the "files" directory is less than the number of rows in "test" file as we only calculate
    # The binomial splits for all probabilities > set global variable "numerical cutoff"
    # Note have to multiple final probability by 2 (1-sided vs. 2-sided)
    
    final_data_matrix <- NULL
    
    for(i in 1:length(all_names))
    {
      message(paste0("Summing final p-values per clone: ",i," in ", length(all_names)))
      setwd(paste0(outputdir,"/",all_names[i]))
      message(paste0("Clone: ",all_names[i]))
      
      files <- list.files()
      
      
      vector_S_minus_B <- NULL
      for(j in 1:length(files))
      {
        temp2<-read.csv(files[j])
        temp_prob<-sum(temp2$prob_S[1]*temp2$prob_A1A2_split*temp2$p_values)
        vector_S_minus_B[j]<-temp_prob
      }

      final_prob<- sum(vector_S_minus_B) *2  # Note have to multiple final probability by 2 (1-sided vs. 2-sided thing)
      
      rearrangement <- all_names[i]
      recipient1_count <- all_recip_1_count[[i]][1]
      recipient2_count <- all_recip_2_count[[i]][1]
      recipient1_total <- T1
      recipient2_total <- T2
      donor_count <- as.numeric(all_S[[i]][1])-as.numeric(all_S_minus_B[[i]][1])
      donor_total <- T_0
      
      final_data_matrix<-rbind(final_data_matrix,data.frame(experiment, 
                                                            donor_count, donor_total, 
                                                            recipient1_count, recipient1_total, 
                                                            recipient2_count, recipient2_total, 
                                                            rearrangement, 
                                                            final_prob))
    }
    #####  Final Output file here #####
    setwd(outputdir) 
    fwrite(final_data_matrix,paste0(experiment," - final probabilities.csv"))    
    ###################################
  }

  message(Sys.time())
  message("Script Complete")   
}

