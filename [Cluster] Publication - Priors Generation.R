library(data.table)
library(plyr)


#----- Updated 08/15/22 -----#
# This script generates a prior matrix for donor TCR probabilities to be used for part 3 of Selective Expander calculations
# This as a standalone script is because it only needs to be run once for all recipient pair comparisons for any given recipient a it only depends on the donor counts
# Will need to save this file in the same folder as the TCR file Library


#----- Choose file to generate log-log plot (Priors) on -----#
setwd("C:/Users/ayeh/Desktop/Projects/[Project] Mouse TCR-seq/Adaptive Submission/[Data] 1901T") #Loading sequences
donor_file <- fread("1901T-B6_45_1-donor-spleen-all.csv")
# File structure must contain the following named columns: 1) "rearrangement", 2) "templates", 3) "frame_type", where:
# 1) "rearrangement" is a unique string containing CDR3 nucleotides (e.g. "CTCTGCAGCCTGGGAATCAGAACGTGCGAAGCAGAAGACTCAGCACTGTACTTGTGCTCCAGCAGTCAAAGGGGTGACACCCAGTAC")
# 2) "templates" is the number of counts identified from sequencing (e.g. "11")
# 3) "frame_type" is either "In" or "Out" - we omit all frame_types that are set to "Out"

#---------------------#  
#----- FUNCTIONS -----#  
#---------------------#
message("Loading Functions")

generate_power_law_variables <- function(countMatrix, filename, plot=FALSE)
  # This function fits best line for power law plots based on the JI paper - PMID: 12682227 (Naumov et al 2003 JI)
  
  # Monte Carlo simulation from JI paper - note that frequency (max = 1) is used on the y axis
  # power law: y = a * x ^ -b, where x is rank (clonality), and y is the frequency for a given rank
  # in otherwords: log(y) = log(a) - b*log(x)
  
  # x_c = e^((log N + log a)/b), where x_c separates the two portions of the log-log curve
  # where x_c is the rank at which the frequency is equal to the reciprocal of the number of different clonotypes (in a frequency axis that is not normalized, this is equal to a frequency of 1)
  # thus log(y) = log(a) - b* ((log(N)+log(a))/b) = log(a) - log(N) - log(a) = -log(N)
  # or y = 1/N
  
# When xi <= x_c, then b = sum(log(a/yi)/logxi)/M), where M = number of ranks
# when xi > x_c, then b = log(a*N)/log(x), whiich means log(y) = log(a) - log(a*N) = log(a) - log(a) - log(N) = - log(N), where N is the number of unique clonotypes, or  y = 1/N
{
  
  #--- Sample data for troubleshooting ---#
  #countMatrix <- count(Sample_1901T_B6Spleen3$templates[Sample_1901T_B6Spleen3$frame_type=="In"])
  #countMatrix
  
  N <- sum(countMatrix[,2]) # N = number of unique clonotypes
  ranks <- countMatrix[,1]
  #-- Now changes frequencies to % from absolute values --#
  freqs <- countMatrix[,2] / N
  #-----#
  
  a <- NULL
  b <- NULL
  x_c <- NULL
  y1 <- freqs[1]
  
  
  #-- Initial estimate --#
  a0 <- y1
  b0 <- 1-a0
  
  # Now calculate residuals for each a_j, b_j, x_cj
  residuals_vector <- NULL
  a_vector <- NULL
  b_vector <- NULL
  x_c_vector <- NULL
  for(j in -100:100) #These are iterations for monte carlo experiment - j allows for random walk for the "a0" parameter with increments of 0.0001
  {
    
    # Goal is to find a residual for each a_j and corresponding b_j and x_cj
    a_j <- a0 + (0.001 * j)
    if(a_j > 1)
    {
      a_j <- 1 # Cap at 1 as you can't get more than that as a fraction
    }
    
    # Now calculate b_j
    b_j <- 0
    bs <- NULL
    for(i in 2:length(ranks))
    {
      # For each rank, assign new value for "b" for both the first and second part of the curve
      # Note that we do not assign a "b" value for the first ran, as log(1) = 0, and we get a divide by "0" error
      #-- Assign "b" for first part of curve --#
      b <- log(a_j/freqs[i])/log(ranks[i]) # Assigns b value for 1st part of curve
      
      #message(paste("b=",b))
      bs <- c(bs,b)
      
    }
    b_j <- sum(bs) / (length(ranks)-1)
    x_cj <- exp((log(N) + log(a_j))/b_j) # Assigns x_cj for each walk
    
    
    for(k in -250:250) # Now vary the b_j by increments to determine best fit
    {
      b_js <- b_j + (0.005 * k)
      # Calculate residual
      sum_residual_squared <- 0
      for(i in 1:length(ranks))
      {
        if(ranks[i] <= x_cj)
        {    
          expected_freq <- a_j / (ranks[i]^b_js)
          residual_squared <- (freqs[i] - expected_freq)^2
          
          # Use the log-scale to calculate residual as to not oveemphasize the larger numbers
          log_expected_freq <- log(expected_freq*N,N)
          log_residual_squared <- (log((freqs[i]*N),N) - (log_expected_freq))^2
          
          
          #sum_residual_squared <- sum_residual_squared + residual_squared
          sum_residual_squared <- sum_residual_squared + log_residual_squared
        }
      }
      residuals_vector <- c(residuals_vector,sum_residual_squared)
      a_vector <- c(a_vector,a_j)
      b_vector <- c(b_vector,b_js)
      x_c_vector <- c(x_c_vector,x_cj)
      
    }
  }
  

  
  min_residual <- min(residuals_vector)
  min_index <- match(min_residual,residuals_vector)
  a <- a_vector[min_index]
  b <- b_vector[min_index]
  x_c <- x_c_vector[min_index]
  
  if(plot==TRUE)
  {
    pdf(filename,width=10,height=10)
    par(mfrow=c(2,2),mar=c(8.1,5.1,2.1,1.1))
    
    plot(freqs~countMatrix[,1], log="xy", type='b', col="blue", pch=19, cex.lab=1.5, cex.main=1.5,
         ylab="Frequency (% of Total)", xlab="Clonality Count",main="Clonal Distribution")
    
    
    lines(1:floor(x_c),
          power_law(1:floor(x_c),a,b),
          col="red")
    
    text(max(ranks)*0.4,max(freqs)*0.9,paste("a = ", round(a,3)))
    text(max(ranks)*0.4,max(freqs)*0.5,paste("b = ", round(b,3)))
    text(max(ranks)*0.4,max(freqs)*0.25,paste("x_c = ", round(x_c,3)))
    
    
    plot(freqs~countMatrix[,1], type='b', col="blue", pch=19, cex.lab=1.5, cex.main=1.5,
         ylab="Frequency", xlab="Clonality Count",main="Clonal Distribution")
    
    lines(1:floor(x_c),
          power_law(1:floor(x_c),a,b),
          col="red")
    
    text(max(ranks)*0.85,max(freqs)*0.95,paste("a = ", round(a,3)))
    text(max(ranks)*0.85,max(freqs)*0.9,paste("b = ", round(b,3)))
    text(max(ranks)*0.85,max(freqs)*0.85,paste("x_c = ", round(x_c,3)))
    
    dev.off() 
  }
  return(c(a,b,x_c))
}

# Power law - this returns y for y = a/(x^b)
power_law <- function(x, a, b)
{
  y <- a / (x^b)
  return (y)
}


message("Generating p(N) coefficients for bayesian analysis")
# coefficients is defined by power law: y = a*x^-b, where x is rank (clonality), and y is the frequency for a given rank
coefficients <- generate_power_law_variables(count(donor_file$templates[donor_file$frame_type=="In"]), 
                                             paste0("priors_log-log_plot.pdf"), plot=TRUE)


coefficients_matrix <- matrix(coefficients,nrow=1,ncol=3)
colnames(coefficients_matrix) <- c("a","b","x_c")
write.csv(coefficients_matrix,paste0("exp_coeff_matrix.csv"),row.names=FALSE)
