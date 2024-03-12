## ----setup, include = FALSE---------------------------------------------------------------------------
library(knitr)
library(formatR)
library(tidyverse)
library(gifski)
knitr::opts_chunk$set(tidy = TRUE, tidy.opts = list(width.cutoff = 72), 
                      cache =TRUE, 
                      echo = TRUE)


## ----warning=FALSE, message=FALSE---------------------------------------------------------------------
#BiocManager::install("LEA")
#devtools::install_github("tdhock/directlabels")

library(dartRverse)

library(ggplot2)
library(data.table)

library(leaflet.minicharts)
library(LEA)


## ----warning=FALSE, message=FALSE, include=FALSE------------------------------------------------------
# Load the data
glsim <- gl.load("./data/sess2_glsim.rds")


## ----eval = FALSE-------------------------------------------------------------------------------------
## # Load the data
## glsim <- gl.load("./data/sess2_glsim.rds")


## -----------------------------------------------------------------------------------------------------
#check the population size and number of populations
table(pop(glsim))
#number of loci
nLoc(glsim)



## ----include=FALSE------------------------------------------------------------------------------------
glw <- gl.load("./data/sess2_glDes.rds")
#check the population size and number of populations
table(pop(glw)) # G stands for trapped in a grid, T for transect
#number of loci
nLoc(glw)

glw <- gl.filter.monomorphs(glw, verbose = 5)


## ----eval = FALSE-------------------------------------------------------------------------------------
## glw <- gl.load("./data/sess2_glDes.rds")
## #check the population size and number of populations
## table(pop(glw)) # G stands for trapped in a grid, T for transect
## #number of loci
## nLoc(glw)
## 
## glw <- gl.filter.monomorphs(glw, verbose = 5)


## -----------------------------------------------------------------------------------------------------
# Run EMIBD9 to detect related individuals. Output saved to save time
#ibd9 <- gl.run.EMIBD9(glw, Inbreed = FALSE, 
 #                     emibd9.path = "E:/Software/Genetic software/EMIBD9")
#saveRDS(ibd9, file="Gen Diversity/ibd9.rds")

# Load EMIBD9 results
ibd9 <- readRDS("./data/sess2_ibd9.rds")

# Do some manipulation to have everything in one table
ibd9Tab <- ibd9[[2]]
# Kick out self comparisons
ibd9Tab <- ibd9Tab[ibd9Tab$Indiv1 != ibd9Tab$Indiv2,  c(1, 2, 21)]
# Add trapping des for each individuals
Des1 <- pop(glw)[as.numeric(ibd9Tab$Indiv1)]
Des2 <- pop(glw)[as.numeric(ibd9Tab$Indiv2)]
# Flag pairs trapping within G, T and in between (BW)
PDes <- ifelse((Des1 == "G" & Des2 == "G"), yes = "G", 
       no = ifelse((Des1 == "T" & Des2 == "T"), yes = "T", "BW"))
# Combine together
ibd9DT <- data.table(cbind(ibd9Tab, Des1, Des2, PDes))
# Compute the mean relatedness
ibd9DT[, mean(as.numeric(`r(1,2)`)), by=PDes]



## -----------------------------------------------------------------------------------------------------
ibd9DT[, sum(`r(1,2)`>=0.125), by=PDes] # HS or more
ibd9DT[, sum(`r(1,2)`>=0.25), by=PDes] # FS or PO

# Recall that the number of pairwise comparisons within 'G' is
nG <- sum(pop(glw) == "G")
nG*(nG-1)/2


## -----------------------------------------------------------------------------------------------------
hwe <- gl.report.hwe(glw, multi_comp = TRUE)



## -----------------------------------------------------------------------------------------------------
het <- gl.report.heterozygosity(glw)


## -----------------------------------------------------------------------------------------------------
gl.report.heterozygosity(glsim)
pcaglsim <- gl.pcoa(glsim)
gl.pcoa.plot(glPca = pcaglsim, x = glsim)

hwe <- gl.report.hwe(glsim, multi_comp = TRUE)



## -----------------------------------------------------------------------------------------------------
# gl <- dartR::gl.read.dart("Report_DUp20-4995_1_moreOrders_SNP_mapping_2.csv",
#                           ind.metafile = "Uperoleia_metadata.csv")

load('./data/session_2.RData') # data named gl



## -----------------------------------------------------------------------------------------------------
gl.map.interactive(gl)


## -----------------------------------------------------------------------------------------------------
#Get rid of really poorly sequenced loci
#But don’t cut hard
gl.report.callrate(gl)
gl.1 <- gl.filter.callrate(gl, method = "loc", threshold = 0.8)

#Very low filter – this is only to get rid of your really bad individuals
gl.report.callrate(gl.1, method = "ind")
gl.2 <- gl.filter.callrate(gl.1, method = "ind", threshold = 0.25)

#Always run this after removing individuals – removes loci that are no longer variable
gl.3 <- gl.filter.monomorphs(gl.2)

#Get rid of unreliable loci
gl.report.reproducibility(gl.3)
gl.4 <- gl.filter.reproducibility(gl.3)

#Get rid of low and super high read depth loci
#do twice so you can zoom in
gl.report.rdepth(gl.4)
gl.5 <- gl.filter.rdepth(gl.4, lower = 0, upper = 25)
gl.clean <- gl.filter.rdepth(gl.5, lower = 8, upper = 17)

nInd(gl.clean)
nLoc(gl.clean)

#look at the data to see if you see any obvious issues and redo if you do.
plot(gl.clean)

rm(gl.1, gl.2, gl.3, gl.4, gl.5)





## ----output = FALSE-----------------------------------------------------------------------------------
gl.report.callrate(gl.clean)
gl.1 <- gl.filter.callrate(gl.clean, method = "loc", threshold = 0.98)

#Remove minor alleles
#I usually set up the threshold so it is just 
# removing singletons to improve computation time
gl.report.maf(gl.1)
gl.2 <- gl.filter.maf(gl.1,  threshold = 1/(2*nInd(gl.1)))

#check that the data looks fairly clean
#this starts to show some obvious population banding
plot(gl.2)

#remove secondary SNPs on the same fragment
#Always do this as the last loci filter so that you’ve cut for quality 
# before you cut because there are two SNPs
gl.3 <- gl.filter.secondaries(gl.2)

#Filter on individuals. You can usually be a bit flexible at this point.
#individuals look a whole lot better now
#make note of any idnviduals with a low call rate. Keep them in for now
#but if they act weird in the analysis, you may want to consider removing
gl.report.callrate(gl.3, method = "ind")
gl.4 <- gl.filter.callrate(gl.3, method = "ind", threshold = .9)
#Always run this after removing individuals
gl.structure <- gl.filter.monomorphs(gl.4)

#this is your cleaned dataset for a population structure analysis              
plot(gl.structure)

nInd(gl.structure)
nLoc(gl.structure)



## ----output = FALSE-----------------------------------------------------------------------------------
#takes quite a while to run (about 15 minutes)
#in case you want to run it uncomment the following lines and comment the readRDS line

# LEA requires the genotype style file
gl2faststructure(gl.structure, outfile = "gl_structure.fstr",
                outpath = './data/')
struct2geno("./data/gl_structure.fstr", ploidy = 2, FORMAT = 2)
###this hates any loci with all heterozygotes

snmf.Sy.K1_8.10 = snmf("./data/gl_structure.fstr.geno", K = 1:8,
                      entropy = T, ploidy = 2, project="new", repetitions = 10)





## -----------------------------------------------------------------------------------------------------
plot(snmf.Sy.K1_8.10)
k <- 2 #chose best based on lowest cross entropy in graph
ce = cross.entropy(snmf.Sy.K1_8.10, K = k)
best <- which.min(ce)
par (mfrow = c(1,1))
barplot(t(Q(snmf.Sy.K1_8.10, K = k, run = best)), col = 1:k)


#Do a PCoA plot. I hate these but they are also good for first pass visualisation.
pc <- gl.pcoa(gl.structure)
gl.pcoa.plot(pc, gl.structure)

rm(ce, gl.1, gl.2, gl.3, gl.4, snmf.Sy.K1_8.10, best, k, pc)




## -----------------------------------------------------------------------------------------------------
get_tajima_D <- function(x){
  require(dartRverse) # possibly not needed for a function in an R package?
  
  # Find allele frequencies (p1 and p2) for every locus in every population
  allele_freqs <- dartR::gl.percent.freq(x)
  #detach(dartR)
  names(allele_freqs)[names(allele_freqs) == "frequency"] <- "p1"
  allele_freqs$p1 <- allele_freqs$p1 / 100
  allele_freqs$p2 <- 1 - allele_freqs$p1
  
  # Get the names of all the populations
  pops <- unique(allele_freqs$popn)
  
  #split each population
  allele_freqs_by_pop <- split(allele_freqs, allele_freqs$popn)
  
  # Internal function to calculate pi
  calc_pi <- function(allele_freqs) {
    n = allele_freqs$nobs * 2  # vector of n values
    pi_sqr <- allele_freqs$p1 ^ 2 + allele_freqs$p2 ^ 2
    h = (n / (n - 1)) * (1 - pi_sqr) # vector of values of h
    sum(h) # return pi, which is the sum of h across loci
  }
  
  get_tajima_D_for_one_pop <- function(allele_freqs_by_pop) {
    pi <- calc_pi(allele_freqs_by_pop)
    
# Calculate number of segregating sites, ignoring missing data (missing data will not appear in teh allele freq calcualtions)
    #S <- sum(!(allele_freqs_by_pop$p1 == 0 | allele_freqs_by_pop$p1 == 1))
    S <- sum(allele_freqs_by_pop$p1 >0 & allele_freqs_by_pop$p1 <1)
    if(S == 0) {
      warning("No segregating sites")
      data.frame(pi = NaN, 
                 S = NaN, 
                 D = NaN, 
                 Pval.normal = NaN, 
                 Pval.beta = NaN)
    }
    
    n <- mean(allele_freqs_by_pop$nobs * 2 )
    
    tmp <- 1:(n - 1)
    a1 <- sum(1/tmp)
    a2 <- sum(1/tmp^2)
    b1 <- (n + 1)/(3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    
    
    #calculate D and do beta testing
    D <- (pi - S/a1) / sqrt(e1 * S + e2 * S * (S - 1))
    Dmin <- (2/n - 1/a1)/sqrt(e2)
    Dmax <- ((n/(2*(n - 1))) - 1/a1)/sqrt(e2)
    tmp1 <- 1 + Dmin * Dmax
    tmp2 <- Dmax - Dmin
    a <- -tmp1 * Dmax/tmp2
    b <- tmp1 * Dmin/tmp2
  
    
    data.frame(pi = pi, 
               S = S, 
               D = D)
  }
  
  output <- do.call("rbind", lapply(allele_freqs_by_pop, 
                                    get_tajima_D_for_one_pop))
  data.frame(population = rownames(output), output, row.names = NULL)
}



## -----------------------------------------------------------------------------------------------------

#This function is written to calculate Tajima's D with a fair amount of missing data
#so we are going to filter lightly here
gl.report.callrate(gl.clean)
gl.1 <- gl.filter.callrate(gl.clean, method = "loc", threshold = 0.9)

#In this first round, we are going to actually remove singletons (our rare alleles) to see what happens
gl.report.maf(gl.1)
nLoc(gl.1)
gl.2 <- gl.filter.maf(gl.1,  threshold = 0.05)
nLoc(gl.2)

#check that the data looks fairly clean
#this starts to show some obvious banding which are the two populations
plot(gl.2)

#we are also going to remove secondary SNPs on the same fragment in this first round
gl.3 <- gl.filter.secondaries(gl.2)

#Filter on individuals. You can usually be a bit flexible at this point.
#make note of any individuals with a low call rate. Keep them in for now
#but if they act weird in the analysis, you may want to consider removing
gl.report.callrate(gl.3, method = "ind")
gl.4 <- gl.filter.callrate(gl.3, method = "ind", threshold = .9)
#Always run this after removing individuals
gl.D_withfiltering <- gl.filter.monomorphs(gl.4)


#calculate tajima's D with removing singletons and secondaries

D_w_filtering <- get_tajima_D(gl.D_withfiltering)
D_w_filtering

rm(gl.1, gl.2, gl.3, gl.4)



## -----------------------------------------------------------------------------------------------------
#This function is written to calculate Tajima's D with a fair amount of missing data
#we are going to filter lightly here 
gl.report.callrate(gl.clean)
gl.1 <- gl.filter.callrate(gl.clean, method = "loc", threshold = 0.9)


#check that the data looks fairly clean
#this starts ot show some obvious population banding
plot(gl.1)

#Filter on individuals. You can usually be a bit flexible at this point.
#make note of any idnviduals with a low call rate. Keep them in for now
#but if they act weird in the analysis, you may want to consider removing
gl.report.callrate(gl.1, method = "ind")
gl.2 <- gl.filter.callrate(gl.1, method = "ind", threshold = .9)
#Always run this after removing individuals
gl.D_withOutfiltering <- gl.filter.monomorphs(gl.2)

rm(gl.1, gl.2)

D_wOUT_filtering <- get_tajima_D(gl.D_withOutfiltering)



## -----------------------------------------------------------------------------------------------------
#This function is written to calculate Tajima's D with a fair amount of missing data
# we are going to filter lightly here 
gl.report.callrate(gl.clean)
gl.1 <- gl.filter.callrate(gl.clean, method = "loc", threshold = 0.9)


#check that the data looks fairly clean
#this starts ot show some obvious population banding
plot(gl.1)

#filter to loci with a lower read depth, so that we are really confident
#that our base calls are correct
gl.report.rdepth(gl.1)
gl.2 <- gl.filter.rdepth(gl.1, lower = 12, upper = 17)


#Filter on individuals. You can usually be a bit flexible at this point.
#make note of any idnviduals with a low call rate. Keep them in for now
#but if they act weird in the analysis, you may want to consider removing
gl.report.callrate(gl.2, method = "ind")
gl.3 <- gl.filter.callrate(gl.2, method = "ind", threshold = .9)
#Always run this after removing individuals
gl.D_withOutfilteringRDepth <- gl.filter.monomorphs(gl.3)

rm(gl.1, gl.2, gl.3)

#calculate tajima's D with removing singletons and secondaries

D_wOUT_filtering_Rdepth <- get_tajima_D(gl.D_withOutfilteringRDepth)



## -----------------------------------------------------------------------------------------------------
D_w_filtering
D_wOUT_filtering
D_wOUT_filtering_Rdepth

