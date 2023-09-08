#!/usr/bin/Rscript

# This script prepares a dataset of TPC predictor variables for each 
# thermal performance dataset. These variables can then used by the 
# ctree_fitting.R script to train multi-output conditional inference 
# trees for predicting the AICc weights of all models.

library(matrixStats)

#####################
# F U N C T I O N S #
#####################

# This function reads all the fits for each dataset and generates 
# an AICc weights matrix.
get_AICc_table <- function()
{
	
	# Read the fits of the dataset with an id of 264 to use as a template. 
	# This dataset includes fits with the maximum number of model parameters.
	load('../Results/model_fits/264.Rda')

	# Prepare the template matrix of AICc weights.
	template_AICc_weights <- AICc_weights
	template_AICc_weights[1,] <- 0
	rm(AICc_weights)

	# Get all the datasets for which we have fits and go through each 
	# dataset.
	model_fits <- list.files('../Results/model_fits/', pattern = '*.Rda')

	# Prepare the final AICc weights matrix.
	AICc_table <- matrix(nrow = length(model_fits), ncol = 83)
	rownames(AICc_table) <- 1:nrow(AICc_table)
	colnames(AICc_table) <- colnames(template_AICc_weights)

	# For each dataset ...
	for ( i in 1:length(model_fits) )
	{
		
		# ... load the corresponding fits.
		load(paste('../Results/model_fits/', model_fits[i], sep = ''))
		
		# Store all AICc weights for this dataset.
		current_AICc_weights <- template_AICc_weights
		
		for ( j in 1:length(colnames(current_AICc_weights)) )
		{
			if ( colnames(current_AICc_weights)[j] %in% colnames(AICc_weights) )
			{
				current_AICc_weights[
					1, colnames(current_AICc_weights)[j]
				] <- AICc_weights[1, colnames(current_AICc_weights)[j]]
			} else
			{
				next
			}
		}
			
		AICc_table[i,] <- current_AICc_weights
		rownames(AICc_table)[i] <- gsub('.Rda', '', model_fits[i])
	}
	
	# Return the final matrix of AICc weights.
	return(AICc_table)
}

# This function extracts the p-value of a linear model.
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

####################
# M A I N  C O D E #
####################

# Read all the thermal performance datasets.
entire_dataset <- read.csv('../Data/thermal_performance_datasets.csv', stringsAsFactors = FALSE)

# Get the AICc weights matrix and write it to an output file.
AICc_table <- get_AICc_table()

write.csv(AICc_table, file = '../Results/AICc_table.csv')

# Prepare vectors to store the predictors.
broad_trait_group <- rep(NA, nrow(AICc_table))
specific_trait_group <- rep(NA, nrow(AICc_table))
kingdom <- rep(NA, nrow(AICc_table))
phylum <- rep(NA, nrow(AICc_table))
n_data_points <- rep(NA, nrow(AICc_table))
n_data_points_rise <- rep(NA, nrow(AICc_table))
n_data_points_fall <- rep(NA, nrow(AICc_table))
n_unique_temps <- rep(NA, nrow(AICc_table))
n_unique_temps_rise <- rep(NA, nrow(AICc_table))
n_unique_temps_fall <- rep(NA, nrow(AICc_table))
median_distance_between_temps <- rep(NA, nrow(AICc_table))
median_distance_between_temps_rise <- rep(NA, nrow(AICc_table))
median_distance_between_temps_fall <- rep(NA, nrow(AICc_table))
median_distance_between_trait_vals <- rep(NA, nrow(AICc_table))
median_distance_between_trait_vals_rise <- rep(NA, nrow(AICc_table))
median_distance_between_trait_vals_fall <- rep(NA, nrow(AICc_table))
sampled_TPC_width <- rep(NA, nrow(AICc_table))
sampled_rise_width <- rep(NA, nrow(AICc_table))
sampled_fall_width <- rep(NA, nrow(AICc_table))
min_temp <- rep(NA, nrow(AICc_table))
min_value_rise_peak <- rep(NA, nrow(AICc_table))
min_value_fall_peak <- rep(NA, nrow(AICc_table))
max_temp <- rep(NA, nrow(AICc_table))
T_pk_value <- rep(NA, nrow(AICc_table))
skew_scalar <- rep(NA, nrow(AICc_table))
rise_cor <- rep(NA, nrow(AICc_table))
fall_cor <- rep(NA, nrow(AICc_table))
rise_exponent <- rep(NA, nrow(AICc_table))
fall_exponent <- rep(NA, nrow(AICc_table))

# Initialise a counter to track the current row number in the AICc 
# weights matrix.
counter <- 1

# For each row name (which stands for a different dataset) ...
for ( i in row.names(AICc_table) )
{
	cat("Now at ", i, "...\n", sep = "")

	# ... extract the corresponding dataset.
	matching_subset_in_db <- entire_dataset[
		entire_dataset$id == i,
	]
	
	# Get the taxonomy and trait information.
	kingdom[counter] <- matching_subset_in_db$kingdom[1]
	phylum[counter] <- matching_subset_in_db$phylum[1]
	
	broad_trait_group[counter] <- matching_subset_in_db$broad_trait_group[1]
	specific_trait_group[counter] <- matching_subset_in_db$specific_trait_group[1]

	# Load the corresponding fits.
	load(paste('../Results/model_fits/', i, '.Rda', sep = ''))

	# Calculate the number of data points, the number of distinct temperatures,
	# the median distance among temperatures, the median distance among 
	# trait values, the width of the TPC, and the minimum and maximum temperatures.
	n_data_points[counter] <- nrow(dataset)
	n_unique_temps[counter] <- length(unique(dataset$temp))
	
	median_distance_between_temps[counter] <- median(diff(sort(unique(dataset$temp))))

	dataset_averaged_by_temperature <- aggregate(
		trait_value ~ temp, dataset, FUN = median
	)

	median_distance_between_trait_vals[counter] <- median(
		abs(
			diff(
				dataset_averaged_by_temperature$trait_value[
					order(dataset_averaged_by_temperature$temp)
				]
			)
		)
	)/(max(dataset$trait_value) - min(dataset$trait_value))
	
	sampled_TPC_width[counter] <- max(dataset$temp) - min(dataset$temp)
	min_temp[counter] <- min(dataset$temp)
	max_temp[counter] <- max(dataset$temp)

	# Fit a quadratic curve to check if the data cover the entire TPC 
	# range, only its rise or fall, or if this cannot be objectively 
	# determined.
	quadratic_fit <- lm(trait_value ~ poly(temp, 2, raw = TRUE), data = dataset)
	
	a <- quadratic_fit$coefficients[3]
	b <- quadratic_fit$coefficients[2]
	g <- quadratic_fit$coefficients[1]
	
	# If these conditions are true, the optimum temperature (T_pk) can be 
	# estimated from the data.
	if ( a < 0 && lmp(quadratic_fit) <= 0.05 )
	{
		
		# Calculate T_pk.
		T_pk <- -b/(2 * a)
		
		# Make sure that T_pk lies within the sampled temperatures.
		if ( T_pk < max_temp[counter] && T_pk > min_temp[counter] )
		{
			T_pk_value[counter] <- T_pk
		}
		
		# Calculate the proportion of data points at temperatures below 
		# the optimum.
		current_prop_rise <- length(dataset$temp[dataset$temp <= T_pk])/nrow(dataset)
		
		# If the proportion of data points at the rise of the curve 
		# are between 0 and 1, then both the rise and the fall of the 
		# curve were sampled.
		if ( current_prop_rise > 0 && current_prop_rise < 1 )
		{

			# For the rise and fall of the TPC, separately, calculate:
			#
			# - the median distance among temperatures,
			#
			# - the median distance among trait values,
			# 
			# - the ratio of the minimum to the maximum trait values,
			#
			# - the skew scalar (if the skew-normal model was able to 
			#   produce a fit),
			#
			# - the number of unique temperatures,
			#
			# - the correlation between temperatures and trait values,
			#
			# - the temperature range width,
			#
			# - the exponent of temperature with which trait values 
			#   increase (or decrease).
			median_distance_between_temps_rise[counter] <- median(diff(sort(unique(dataset$temp[dataset$temp <= T_pk]))))
			
			dataset_averaged_by_temperature_rise <- dataset_averaged_by_temperature[
				dataset_averaged_by_temperature$temp <= T_pk,
			]
			
			median_distance_between_trait_vals_rise[counter] <- median(
				abs(
					diff(
						dataset_averaged_by_temperature_rise$trait_value[
							order(dataset_averaged_by_temperature_rise$temp)
						]
					)
				)
			)/(max(dataset$trait_value) - min(dataset$trait_value))
									
			min_value_rise_peak[counter] <- min(dataset$trait_value[dataset$temp <= T_pk])/max(dataset$trait_value)
			min_value_fall_peak[counter] <- min(dataset$trait_value[dataset$temp > T_pk])/max(dataset$trait_value)
			
			if ( !is.null(fitted_models$skew_normal_4_pars) )
			{
				skew_scalar[counter] <- -coefficients(fitted_models$skew_normal_4_pars)[3]
			}
			
			n_data_points_rise[counter] <- nrow(dataset[dataset$temp < T_pk,])
			n_data_points_fall[counter] <- nrow(dataset[dataset$temp > T_pk,])
			
			n_unique_temps_rise[counter] <- length(unique(dataset$temp[dataset$temp < T_pk]))
			n_unique_temps_fall[counter] <- length(unique(dataset$temp[dataset$temp > T_pk]))
						
			rise_cor[counter] <- cor(
				dataset$trait_value[dataset$temp < T_pk], 
				dataset$temp[dataset$temp < T_pk]
			)
			
			fall_cor[counter] <- cor(
				dataset$trait_value[dataset$temp > T_pk], 
				dataset$temp[dataset$temp > T_pk]
			)

			dataset_rise <- dataset[dataset$temp < T_pk & dataset$temp > 0,]
			sampled_rise_width[counter] <- max(dataset_rise$temp) - min(dataset_rise$temp)
			
			if (nrow(dataset_rise) >= 2 && length(unique(dataset_rise$temp)) >= 2 )
			{
				exponent_fit <- lm(log(trait_value) ~ log(temp), data = dataset_rise)
				rise_exponent[counter] <- coefficients(exponent_fit)[2]
			}

			median_distance_between_temps_fall[counter] <- median(diff(sort(unique(dataset$temp[dataset$temp > T_pk]))))
			
			dataset_averaged_by_temperature_fall <- dataset_averaged_by_temperature[
				dataset_averaged_by_temperature$temp > T_pk,
			]
			
			median_distance_between_trait_vals_fall[counter] <- median(
				abs(
					diff(
						dataset_averaged_by_temperature_fall$trait_value[
							order(dataset_averaged_by_temperature_fall$temp)
						]
					)
				)
			)/(max(dataset$trait_value) - min(dataset$trait_value))
						
			dataset_fall <- dataset[dataset$temp > T_pk,]
			sampled_fall_width[counter] <- max(dataset_fall$temp) - min(dataset_fall$temp)

			if ( nrow(dataset_fall) >= 2 && length(unique(dataset_fall$temp)) >= 2 )
			{
				exponent_fit <- lm(log(trait_value) ~ log(temp), data = dataset_fall)
				fall_exponent[counter] <- coefficients(exponent_fit)[2]				
			}
			
		} else if ( current_prop_rise == 0 )
		{
			
			# If the proportion of data points at the rise of the curve 
			# is 0, only the fall of the curve was sampled. Therefore, 
			# calculate only the applicable variables for the fall.
			median_distance_between_temps_fall[counter] <- median(diff(sort(unique(dataset$temp))))
			
			median_distance_between_trait_vals_fall[counter] <- median(
				abs(
					diff(
						dataset_averaged_by_temperature$trait_value[
							order(dataset_averaged_by_temperature$temp)
						]
					)
				)
			)/(max(dataset$trait_value) - min(dataset$trait_value))
			
			sampled_rise_width[counter] <- 0			
			sampled_fall_width[counter] <- max(dataset$temp) - min(dataset$temp)
			
			exponent_fit <- lm(log(trait_value) ~ log(temp), data = dataset)
			fall_exponent[counter] <- coefficients(exponent_fit)[2]
			
			fall_cor[counter] <- cor(dataset$trait_value, dataset$temp)
			
			n_data_points_rise[counter] <- 0
			n_data_points_fall[counter] <- nrow(dataset)
			
			n_unique_temps_rise[counter] <- 0
			n_unique_temps_fall[counter] <- length(unique(dataset$temp))
		} else
		{
			
			# If the proportion of data points at the rise of the curve 
			# is 1, only the rise of the curve was sampled. Therefore, 
			# calculate only the applicable variables for the rise. 
			median_distance_between_temps_rise[counter] <- median(diff(sort(unique(dataset$temp))))
			
			median_distance_between_trait_vals_rise[counter] <- median(
				abs(
					diff(
						dataset_averaged_by_temperature$trait_value[
							order(dataset_averaged_by_temperature$temp)
						]
					)
				)
			)/(max(dataset$trait_value) - min(dataset$trait_value))
			
			sampled_rise_width[counter] <- max(dataset$temp) - min(dataset$temp)
			sampled_fall_width[counter] <- 0

			n_data_points_rise[counter] <- nrow(dataset)
			n_data_points_fall[counter] <- 0
			
			n_unique_temps_rise[counter] <- length(unique(dataset$temp))
			n_unique_temps_fall[counter] <- 0
			
			rise_cor[counter] <- cor(dataset$trait_value, dataset$temp)
						
			if ( nrow(dataset[dataset$temp > 0,]) >= 2 && length(unique(dataset$temp[dataset$temp > 0])) >= 2 )
			{
				exponent_fit <- lm(log(trait_value) ~ log(temp), data = dataset[dataset$temp > 0,])
				rise_exponent[counter] <- coefficients(exponent_fit)[2]
			}			
		}
	} else
	{
		
		# If the optimum cannot be objectively estimated, calculate 
		# the correlation between trait values and temperatures.
		current_cor <- cor.test(dataset$trait_value, dataset$temp)
		
		# If the p_value of the correlation is 0.05 or lower ...
		if ( current_cor$p.value <= 0.05 )
		{
			
			# ... then check if the fall or rise was sampled.
			if ( current_cor$estimate < 0 )
			{
				
				# If the correlation estimate is negative, the fall 
				# was sampled. Therefore, calculate only the applicable 
				# variables for the fall. 				
				median_distance_between_temps_fall[counter] <- median(diff(sort(unique(dataset$temp))))
				
				median_distance_between_trait_vals_fall[counter] <- median(
					abs(
						diff(
							dataset_averaged_by_temperature$trait_value[
								order(dataset_averaged_by_temperature$temp)
							]
						)
					)
				)/(max(dataset$trait_value) - min(dataset$trait_value))

				sampled_rise_width[counter] <- 0
				sampled_fall_width[counter] <- max(dataset$temp) - min(dataset$temp)
				
				n_data_points_rise[counter] <- 0
				n_data_points_fall[counter] <- nrow(dataset)
			
				n_unique_temps_rise[counter] <- 0
				n_unique_temps_fall[counter] <- length(unique(dataset$temp))
				
				fall_cor[counter] <- cor(dataset$trait_value, dataset$temp)
				
				if ( nrow(dataset[dataset$temp > 0,]) >= 2 && length(unique(dataset$temp[dataset$temp > 0])) >= 2 )
				{
					exponent_fit <- lm(log(trait_value) ~ log(temp), data = dataset[dataset$temp > 0,])
					fall_exponent[counter] <- coefficients(exponent_fit)[2]
				}
			} else
			{
				
				# If the correlation estimate is positive, the rise 
				# was sampled. Therefore, calculate only the applicable 
				# variables for the rise.
				median_distance_between_temps_rise[counter] <- median(diff(sort(unique(dataset$temp))))
				
				median_distance_between_trait_vals_rise[counter] <- median(
					abs(
						diff(
							dataset_averaged_by_temperature$trait_value[
								order(dataset_averaged_by_temperature$temp)
							]
						)
					)
				)/(max(dataset$trait_value) - min(dataset$trait_value))
								
				sampled_rise_width[counter] <- max(dataset$temp) - min(dataset$temp)
				sampled_fall_width[counter] <- 0
				
				n_data_points_rise[counter] <- nrow(dataset)
				n_data_points_fall[counter] <- 0
			
				n_unique_temps_rise[counter] <- length(unique(dataset$temp))
				n_unique_temps_fall[counter] <- 0
				
				rise_cor[counter] <- cor(
					dataset$trait_value, dataset$temp
				)
								
				if ( nrow(dataset[dataset$temp > 0,]) >= 2 && length(unique(dataset$temp[dataset$temp > 0])) >= 2 )
				{
					exponent_fit <- lm(log(trait_value) ~ log(temp), data = dataset[dataset$temp > 0,])
					rise_exponent[counter] <- coefficients(exponent_fit)[2]
				}
			}
		}
	}
	
	counter <- counter + 1
}

# Store all the predictor variables to a data frame.
dataset_of_predictors <- data.frame(
	kingdom = kingdom,
	phylum = phylum,
	broad_trait_group = broad_trait_group,
	specific_trait_group = specific_trait_group,
	n_data_points = n_data_points,
	n_data_points_rise = n_data_points_rise,
	n_data_points_fall = n_data_points_fall,
	n_unique_temps = n_unique_temps,
	n_unique_temps_rise = n_unique_temps_rise,
	n_unique_temps_fall = n_unique_temps_fall,
	median_distance_between_temps = median_distance_between_temps,
	median_distance_between_temps_rise = median_distance_between_temps_rise,
	median_distance_between_temps_fall = median_distance_between_temps_fall,
	median_distance_between_trait_vals = median_distance_between_trait_vals,
	median_distance_between_trait_vals_rise = median_distance_between_trait_vals_rise,
	median_distance_between_trait_vals_fall = median_distance_between_trait_vals_fall,
	sampled_TPC_width = sampled_TPC_width,
	sampled_rise_width = sampled_rise_width,
	sampled_fall_width = sampled_fall_width,
	min_temp = min_temp,
	min_value_rise_peak = min_value_rise_peak,
	min_value_fall_peak = min_value_fall_peak,
	max_temp = max_temp,
	T_pk_value = T_pk_value,
	skew_scalar = skew_scalar,
	rise_cor = rise_cor,
	fall_cor = fall_cor,
	rise_exponent = rise_exponent,
	fall_exponent = fall_exponent
)

# Specify the ids of the thermal performance datasets as row names.
row.names(dataset_of_predictors) <- row.names(AICc_table)

# Write the predictor variables to an output file.
write.csv(
	dataset_of_predictors, file = '../Results/TPC_predictors_dataset.csv'
)
