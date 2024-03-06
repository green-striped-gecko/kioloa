library(dartRverse)
library(ggplot2)

#######how do estimates of heterozygosity change as the number of individuals changes#####
#######in your SNP calling?

#This uses data from Litoria rubella, a very abundant and widespread frog species
#we are using data following https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13947

#let's focus on the Kimberley, where we have lots of samples, and we expect the 
#diversity to be really high because there are millions of them everywhere

#SNPs were called on different numbers of individuals, in increments of 5 using
#the Stacks pipeline. 

#we are only going to do a simple filter for call rates because we already
#filtered for allele depth in Stacks. SNPs were called independent of other
#populations, which we will get to later.

load('./data/Session3_data.RData')

#create a list of the kimberley genlights
kimberley_names <- ls(pattern = "^Kimberley")
#put all the genlights into a mega list
kimberley <- mget(kimberley_names)

#we're going to do this in a loop for speed, applying the same filters
# Iterate over the names of the kimberley list
for(name in names(kimberley)){
  # Extract the genlight object from the kimberley list using its name
  genlight_object <- kimberley[[name]]
  # Apply the filter call rate function
  # Assuming gl.filter.callrate is a function that operates on a genlight object
  filtered_object <- gl.filter.callrate(genlight_object, threshold = 0.7, mono.rm = TRUE)
  # Assign the filtered object back to the environment with a new name
  assign(paste0(name, "_0.7"), filtered_object)
}

####now calculate heterozygosity 

# List all object names in the environment
all_names <- ls()

# Use grep() to match names that start with "Kimberley" and end with "0.7", .+ indicates any characters in between
kimberley_filtered <- grep("^Kimberley.+0\\.7$", all_names, value = TRUE)#put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_filtered)

#Initialize an empty data frame
heterozygosity_reports_df <- data.frame()

# Iterate over the kimberley list to apply gl.report.heterozygosity and bind the results
for(name in names(kimberley)) {
  # Apply the function
  report <- gl.report.heterozygosity(kimberley[[name]])
  
  # Add 'ObjectName' as the first column of the report
  report <- cbind(ObjectName = name, report)
  
  # Bind this report to the main data frame
  heterozygosity_reports_df <- bind_rows(heterozygosity_reports_df, report)
}

# heterozygosity_reports_df now contains all the reports with an additional column for object names
View(heterozygosity_reports_df)


# Example using ggplot2 to plot the data
library(ggplot2)
kimberley_Ho_0.7callrate <- ggplot(heterozygosity_reports_df, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity by Sample number", x = "Sample Number", y = "Observed Heterozygosity (Ho)")

kimberley_Ho_0.7callrate

#As you can see, different numbers samples can substantially 
#change your heterozygosity estimate.


######lets redo this, but see how changing your call rate filter impacts estimates#####

# List all object names in the environment
all_names <- ls()

# Use grep() to match names that start with "Kimberley" 
kimberley_names <- grep("^Kimberley.*\\.vcf$", all_names, value = TRUE) #put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_names)

#we're going to do this in a loop for speed, applying the same filters
# Iterate over the names of the kimberley list
for(name in names(kimberley)){
  # Extract the genlight object from the kimberley list using its name
  genlight_object <- kimberley[[name]]
  # Apply the filter call rate function
  # Assuming gl.filter.callrate is a function that operates on a genlight object
  filtered_object <- gl.filter.callrate(genlight_object, threshold = .95, mono.rm = TRUE)
  # Assign the filtered object back to the environment with a new name
  assign(paste0(name, "_0.95"), filtered_object)
}

####now calculate heterozygosity 

# List all object names in the environment
all_names <- ls()

# Use grep() to match names that start with "Kimberley" and end with "0.95", .+ indicates any characters in between
kimberley_filtered <- grep("^Kimberley.+0\\.95$", all_names, value = TRUE)#put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_filtered)


#Initialize an empty data frame
heterozygosity_reports_df_0.95 <- data.frame()

# Iterate over the kimberley list to apply gl.report.heterozygosity and bind the results
for(name in names(kimberley)) {
  # Apply the function
  report <- gl.report.heterozygosity(kimberley[[name]])
  
  # Add 'ObjectName' as the first column of the report
  report <- cbind(ObjectName = name, report)
  
  # Bind this report to the main data frame
  heterozygosity_reports_df_0.95 <- bind_rows(heterozygosity_reports_df_0.95, report)
}

# heterozygosity_reports_df now contains all the reports with an additional column for object names
View(heterozygosity_reports_df_0.95)


# Example using ggplot2 to plot the data
library(ggplot2)
kimberley_Ho_0.95callrate <- ggplot(heterozygosity_reports_df_0.95, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity by Sample number 0.95 Call Rate", x = "Sample Number", y = "Observed Heterozygosity (Ho)")

kimberley_Ho_0.95callrate

par(mfrow = c(2,1))
kimberley_Ho_0.7callrate +kimberley_Ho_0.95callrate

#Higher call rate filters can reduce your slightly reduce your heterozygosity estimate


#######does subsetting to have the same sample numbers fix the issue? ########
#######compare heterozygosity when SNPs are called on 5 individuals, versus called
#######on 40 individuals and then is subset to 5 individuals


#If you called SNPs on more individuals than you wanted to make them equal, can you just remove
#individuals without re-calling SNPS? Let's test.


# Use grep() to match names that start with "Kimberley" 
kimberley_names <- grep("^Kimberley.*\\.vcf$", all_names, value = TRUE) #put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_names)

#Now lets subsample the datasets down to the same five individuals
#so that the only difference is our SNP calling
#who are the individuals
inds <- indNames(Kimberley_n_05.vcf_0.7)

#Initialize an empty data frame
heterozygosity_when_subsampling <- data.frame()

for(name in names(kimberley)) {
  # Access the genlight object from your list
  genlight_object <- kimberley[[name]]
  
  # Subset the individuals
  x <- gl.keep.ind(genlight_object, ind.list = inds, mono.rm = TRUE)
  
  #filter on call rate
  
  x <- gl.filter.callrate(x, threshold = .7)
  
  # Apply the function
  report <- gl.report.heterozygosity(x)
  
  # Add 'ObjectName' as the first column of the report
  report$ObjectName <- name # Adjusting to add column without cbind to maintain data frame classes
  
  # Bind this report to the main data frame
  heterozygosity_when_subsampling <- bind_rows(heterozygosity_when_subsampling, report)
}

kimberley_Ho_subsampling <- ggplot(heterozygosity_when_subsampling, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity when subsampling", x = "SNP Calling", y = "Observed Heterozygosity (Ho)")

kimberley_Ho_subsampling

#there are some minor differences but it's not too bad


####Is this true for all populations?

#create a list of the southeast genlights
Southeast_names <- ls(pattern = "^SouthEast")
#put all the genlights into a mega list
southeast <- mget(Southeast_names)


#If you called SNPs on more individuals than you wanted to make them equal, can you just remove
#individuals without re-calling SNPS? Let's test.

#Now lets subsample the datasets down to the same five individuals
#so that the only difference is our SNP calling
#who are the individuals
inds <- indNames(SouthEast_n_05.vcf)

#Initialize an empty data frame
heterozygosity_when_subsampling_southeast <- data.frame()

for(name in names(southeast)) {
  # Access the genlight object from your list
  genlight_object <- southeast[[name]]
  
  # Subset the individuals
  x <- gl.keep.ind(genlight_object, ind.list = inds, mono.rm = TRUE)
  
  #filter on call rate
  
  x <- gl.filter.callrate(x, threshold = .7)
  
  # Apply the function
  report <- gl.report.heterozygosity(x)
  
  # Add 'ObjectName' as the first column of the report
  report$ObjectName <- name # Adjusting to add column without cbind to maintain data frame classes
  
  # Bind this report to the main data frame
  heterozygosity_when_subsampling <- bind_rows(heterozygosity_when_subsampling, report)
}

southeast_Ho_subsampling <- ggplot(heterozygosity_when_subsampling, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity when subsampling", x = "SNP Calling", y = "Observed Heterozygosity (Ho)")

southeast_Ho_subsampling


#create a list of the central australian genlights
central_names <- ls(pattern = "^Central")
#put all the genlights into a mega list
central <- mget(central_names)


#If you called SNPs on more individuals than you wanted to make them equal, can you just remove
#individuals without re-calling SNPS? Let's test.

#Now lets subsample the datasets down to the same five individuals
#so that the only difference is our SNP calling
#who are the individuals
inds <- indNames(Central_n_05.vcf)

#Initialize an empty data frame
heterozygosity_when_subsampling_central <- data.frame()

for(name in names(central)) {
  # Access the genlight object from your list
  genlight_object <- central[[name]]
  
  # Subset the individuals
  x <- gl.keep.ind(genlight_object, ind.list = inds, mono.rm = TRUE)
  
  #filter on call rate
  
  x <- gl.filter.callrate(x, threshold = .7)
  
  # Apply the function
  report <- gl.report.heterozygosity(x)
  
  # Add 'ObjectName' as the first column of the report
  report$ObjectName <- name # Adjusting to add column without cbind to maintain data frame classes
  
  # Bind this report to the main data frame
  heterozygosity_when_subsampling <- bind_rows(heterozygosity_when_subsampling, report)
}

central_Ho_subsampling <- ggplot(heterozygosity_when_subsampling, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity when subsampling 5 individuals", x = "SNP Calling", y = "Observed Heterozygosity (Ho)")

central_Ho_subsampling


###so we're not seeing too much impact of subsampling if SNPs are called 
###separately on populations.

# BUT BUT BUT BUT BUT BUT BUT
# The kimberley population has the highest actual heterozygosity (see the paper)
# the difference here is that this population lost a lot of variable sites 
# because 7% of the variable sites have more than two alleles
# So a comparison between Central & SE and Kimberley Rubella would find the highest in the
# wrong population.

# This is very problematic and is an ongoing bioinformatic issue.

#Therefore, calling SNPs separately on populations and keeping sample sizes 
#the same does not resolve the issues with heterozygosity based on SNPs.


#### Calling populations together - Does this resolve the issue? #########


###Calling SNPs together but with different sample sizes ##########




#######Calling SNPs together with different sample sizes for each population then #########
#######subsetting to equal sample sizes ########




###Calling SNPs together, but with the same sample sizes ###########




### SNP-based heterozygosity is bad and cannot be used to compare different
### Populations, and there are no clear workarounds that will give you
### reliable answers.

### You need to report autosomal heterozygosity at a minimum, which means
### doing your own bioinformatics from raw data at this time.