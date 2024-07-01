#!/usr/bin/env Rscript

# This script handles the fitting of TPC models to thermal performance 
# datasets, as well as the calculation of R-squared, AICc weights, and 
# dataset-specific distance matrices among model fits.
#
# The script requires a command-line argument to determine if models
# will be fitted to trait or enzyme data:
#
# ./TPC_fitting_pipeline.R traits
#
# OR
#
# ./TPC_fitting_pipeline.R enzymes

library(doParallel)
library(minpack.lm)
library(Matrix)
library(MuMIn)
library(mvMORPH)
library(nls.multstart)
library(pracma)
library(puma)
library(stringr)

#####################
# F U N C T I O N S #
#####################

# Read all the TPC fitting functions from the other script.
source('TPC_fitting_functions.R')

# This function takes all the model fits for a given thermal performance 
# dataset and calculates a distance matrix among them.
calculate_distance_matrix_across_model_fits <- function(
		fitted_models, model_list, dataset, R_squared_vals 
	)
{
	
	# Prepare the distance matrix.
	dist_matrix <- matrix(nrow = length(fitted_models), ncol = length(fitted_models))
	
	rownames(dist_matrix) <- names(model_list)
	colnames(dist_matrix) <- names(model_list)
	
	# Set the diagonals (i.e., the distance between a model fit and 
	# itself) to zero.
	diag(dist_matrix) <- 0
	
	# Initialise the maximum height across all fitted curves to zero.
	max_predicted_height <- 0
	
	# For each model ...
	for ( i in 1:(length(fitted_models) - 1) )
	{
		
		# ... make sure that a fit with an R-squared of at least 0.5 exists.
		if ( is.null(fitted_models[[i]]) || (R_squared_vals[i] < 0.5) )
		{
			next
		} else
		{
			
			# ... go through all other model fits with an R-squared of 
			# at least 0.5.
			for ( j in (i+1):length(fitted_models) )
			{
				if ( is.null(fitted_models[[j]]) || (R_squared_vals[j] < 0.5) )
				{
					next
				} else
				{
					
					# Calculate the height of the fitted curves at 
					# all temperatures from the minimum one to the 
					# maximum one in the dataset, with a step size 
					# of 0.1°C.
					temps_for_prediction <- seq(
						min(dataset$temp), max(dataset$temp), 0.1
					)
					
					pred_1 <- exp(
						predict(
							fitted_models[[i]], newdata = data.frame(
								temp = temps_for_prediction
							)
						)
					)
					
					pred_2 <- exp(
						predict(
							fitted_models[[j]], newdata = data.frame(
								temp = temps_for_prediction
							)
						)
					)

					# Make sure that there are no infinite heights in 
					# both fits that are being compared.
					if ( any(is.infinite(pred_1)) || any(is.infinite(pred_2)) )
					{
						next
					} else
					{
						
						# Check if the maximum height of the two fits is 
						# higher than what we already have observed.
						max_predicted_height <- max(
							c(
								max_predicted_height,
								pred_1,
								pred_2
							)
						)
						
						# Calculate the Euclidean distance between the 
						# two fits.
						dist_matrix[i,j] <- sqrt(sum((pred_1 - pred_2)^2))
					}
				}
			}
		}
	}
	
	# Return the distance matrix, with maximum normalisation based on 
	# the maximum observed TPC height.
	return(forceSymmetric(dist_matrix / max_predicted_height, uplo = 'U'))
}

# This function fits each of the 83 models to a thermal performance dataset.
fit_models <- function(dataset, ID, data_type)
{
	fitting_status <- NA
	
	# Determine the maximum number of parameters that can be used 
	# (one fewer than the number of distinct temperatures).
	max_pars_start <- length(unique(dataset$temp)) - 1
	if ( max_pars_start > 7 )
	{
		max_pars_start <- 7
	}
	
	# Try to fit models with the maximum number of parameters. If 
	# there is a model that returns an AICc of -infinity, then 
	# decrease the maximum number of parameters and try again, until 
	# you get to models with 3 parameters.
	for ( max_pars in max_pars_start:3 )
	{
		
		# Identify the models that can be fitted.
		model_list <- get_model_list(max_pars)
		
		# Start 8 parallel threads.
		cl <- makeCluster(8)
		registerDoParallel(cl)
		
		# Fit each model.
		fitted_models <- foreach( 
			x = model_list,
			.final = function(x) setNames(x, names(model_list)),
			.packages = c('minpack.lm', 'nls.multstart', 'pracma')
		) %dopar% 
		{
			set.seed(1337)
			x(dataset)
		}
		
		# Stop the threads.
		stopCluster(cl)
		
		# Get ready to store the AICc weights and R-squared values.
		AICc_weights <- matrix(nrow = 1, ncol = length(fitted_models))
		colnames(AICc_weights) <- names(fitted_models)
	
		AICc_vals <- rep(NA, length(fitted_models))
		R_squared_vals <- rep(NA, length(fitted_models))
		
		# For each model ...
		for ( i in 1:length(fitted_models))
		{
			
			# ... if a fit could not be obtained, set AICc to a huge 
			# number.
			if ( is.null(fitted_models[[i]]) )
			{
				AICc_vals[i] <- 99999999
			} else
			{
				
				# Otherwise, calculate the R-squared.
				rss <- sum((dataset$trait_value - exp(fitted(fitted_models[[i]])))^2)
				tss <- sum((dataset$trait_value - mean(dataset$trait_value))^2)
				
				r_squared <- 1 - (rss/tss)
				R_squared_vals[i] <- r_squared
				
				# If the R-squared is at least 0.5, then calculate the 
				# AICc value for the fit.
				if ( r_squared >= 0.5 )
				{
					AICc_vals[i] <- AICc(fitted_models[[i]])
				} else
				{
					AICc_vals[i] <- 99999999
				}
			}
		}
		
		# If any of the AICc values was -infinity, go to the beginning 
		# and try to fit models with fewer parameters.
		if ( any(is.infinite(AICc_vals)) )
		{
			fitting_status <- '-Inf AICc for at least one fit'
			next
		} else if ( all(AICc_vals == 99999999) )
		{
			
			# If no model could produce an R-squared of at least 0.5, 
			# stop trying to fit models to this dataset.
			fitting_status <- 'no R2 >= 0.5'
			break
		} else
		{
			
			# Otherwise, return the maximum number of model parameters, 
			# the AICc weights, calculate the distance matrix among 
			# fits, and store all the objects to an .Rda file.
			fitting_status <- paste(
				'fitted models with up to ', max_pars, 
				' parameters', sep = ''
			)
			
			AICc_weights[1,] <- aicw(AICc_vals)$aicweights
			
			dist_matrix <- calculate_distance_matrix_across_model_fits(
				fitted_models, model_list, dataset, R_squared_vals
			)

			if ( data_type == 'traits' )
			{
				save(
					dataset, fitted_models, AICc_weights, R_squared_vals,
					dist_matrix, file = paste(
						'../Results/model_fits/', ID, '.Rda', sep = ''
					)
				)
			} else if ( data_type == 'enzymes' )
			{
				save(
					dataset, fitted_models, AICc_weights, R_squared_vals,
					dist_matrix, file = paste(
						'../Results/enzyme_model_fits/', ID, '.Rda', sep = ''
					)
				)
			}
			break
		}
	}
	
	# Return the result of the fitting process.
	return(fitting_status)
}

# This function returns a list of models that have no more than 
# X parameters (the "max_pars" argument).
get_model_list <- function(max_pars)
{
	model_list <- list(
		Johnson_Lewin_4_pars = fit_Johnson_Lewin_4_pars,
		
		extended_Johnson_Lewin_5_pars = fit_extended_Johnson_Lewin_5_pars,
		
		simplified_Johnson_Lewin_4_pars = fit_simplified_Johnson_Lewin_4_pars,
		
		simplified_extended_Johnson_Lewin_5_pars = fit_simplified_extended_Johnson_Lewin_5_pars,
		
		Sharpe_Schoolfield_6_pars = fit_Sharpe_Schoolfield_6_pars,
		
		extended_Sharpe_Schoolfield_7_pars = fit_extended_Sharpe_Schoolfield_7_pars,

		simplified_Sharpe_Schoolfield_6_pars = fit_simplified_Sharpe_Schoolfield_6_pars,
		
		simplified_extended_Sharpe_Schoolfield_7_pars = fit_simplified_extended_Sharpe_Schoolfield_7_pars,
				
		Gaussian_3_pars = fit_Gaussian_3_pars,
		
		double_Gaussian_4_pars = fit_double_Gaussian_4_pars,
		
		modified_Gaussian_4_pars = fit_modified_Gaussian_4_pars,
		
		exponentially_modified_Gaussian_5_pars = fit_exponentially_modified_Gaussian_5_pars,

		Gaussian_Gompertz_5_pars = fit_Gaussian_Gompertz_5_pars,
		
		skew_normal_4_pars = fit_skew_normal_4_pars,
		
		Weibull_4_pars = fit_Weibull_4_pars,

		Logan_I_4_pars = fit_Logan_I_4_pars,
		
		Logan_II_5_pars = fit_Logan_II_5_pars,
		
		Logan_III_4_pars = fit_Logan_III_4_pars,
		
		Lactin_I_3_pars = fit_Lactin_I_3_pars,
		
		Lactin_II_4_pars = fit_Lactin_II_4_pars,
				
		Asbury_Angilletta_6_pars = fit_Asbury_Angilletta_6_pars,
		
		simplified_Asbury_Angilletta_4_pars = fit_simplified_Asbury_Angilletta_4_pars,
		
		Thomas_I_4_pars = fit_Thomas_I_4_pars,
		
		Thomas_II_5_pars = fit_Thomas_II_5_pars,
		
		Ratkowsky_4_pars = fit_Ratkowsky_4_pars,
		
		second_order_polynomial_3_pars = fit_2nd_order_polynomial_3_pars,
				
		third_order_polynomial_4_pars = fit_3rd_order_polynomial_4_pars,
		
		fourth_order_polynomial_5_pars = fit_4th_order_polynomial_5_pars,
		
		fifth_order_polynomial_6_pars = fit_5th_order_polynomial_6_pars,
		
		Ashrafi_I_3_pars = fit_Ashrafi_I_3_pars,
		
		Ashrafi_II_3_pars = fit_Ashrafi_II_3_pars,
		
		Ashrafi_III_3_pars = fit_Ashrafi_III_3_pars,
		
		Ashrafi_IV_4_pars = fit_Ashrafi_IV_4_pars,
		
		Ashrafi_V_4_pars = fit_Ashrafi_V_4_pars,
		
		Enzyme_Assisted_Arrhenius_5_pars = fit_Enzyme_Assisted_Arrhenius_5_pars,
		
		Eubank_3_pars = fit_Eubank_3_pars,
		
		Stinner_4_pars = fit_Stinner_4_pars,
		
		Briere_I_3_pars = fit_Briere_I_3_pars,
		
		simplified_Briere_I_3_pars = fit_simplified_Briere_I_3_pars,
		
		Briere_II_4_pars = fit_Briere_II_4_pars,
		
		simplified_Briere_II_4_pars = fit_simplified_Briere_II_4_pars,
		
		extended_Briere_5_pars = fit_extended_Briere_5_pars,
		
		simplified_extended_Briere_5_pars = fit_simplified_extended_Briere_5_pars,
		
		Janisch_I_3_pars = fit_Janisch_I_3_pars,
		
		Janisch_II_4_pars = fit_Janisch_II_4_pars,
		
		simplified_beta_type_3_pars = fit_simplified_beta_type_3_pars,
		
		Ritchie_4_pars = fit_Ritchie_4_pars,
		
		Thornton_Lessem_6_pars = fit_Thornton_Lessem_6_pars,
		
		Huey_Stevenson_5_pars = fit_Huey_Stevenson_5_pars,
		
		Stevenson_6_pars = fit_Stevenson_6_pars,
		
		Wang_Lan_Ding_7_pars = fit_Wang_Lan_Ding_7_pars,
		
		ONeill_4_pars = fit_ONeill_4_pars,
						
		modified_Deutsch_4_pars = fit_modified_Deutsch_4_pars,
				
		Johnk_5_pars = fit_Johnk_5_pars,
		
		Rezende_Bozinovic_4_pars = fit_Rezende_Bozinovic_4_pars,
				
		Rice_Clock_6_pars = fit_Rice_Clock_6_pars,
		
		Bilinear_4_pars = fit_Bilinear_4_pars,
		
		modified_Bilinear_6_pars = fit_modified_Bilinear_6_pars,
		
		Dent_like_5_pars = fit_Dent_like_5_pars,
		
		modified_Dent_like_7_pars = fit_modified_Dent_like_7_pars,
		
		Yan_Hunt_4_pars = fit_Yan_Hunt_4_pars,
		
		Cardinal_Temperature_4_pars = fit_Cardinal_Temperature_4_pars,
		
		Cardinal_Temperature_with_inflection_4_pars = fit_Cardinal_Temperature_with_inflection_4_pars,
		
		Ross_Ratkowsky_5_pars = fit_Ross_Ratkowsky_5_pars,
		
		Hobbs_4_pars = fit_Hobbs_4_pars,
				
		Mitchell_Angilletta_3_pars = fit_Mitchell_Angilletta_3_pars,
		
		Wang_Engel_4_pars = fit_Wang_Engel_4_pars,
		
		Regniere_6_pars = fit_Regniere_6_pars,
		
		Analytis_Allahyari_5_pars = fit_Analytis_Allahyari_5_pars,
		
		Analytis_Kontodimas_3_pars = fit_Analytis_Kontodimas_3_pars,
		
		linear_logistic_5_pars = fit_linear_logistic_5_pars,

		Tomlinson_Menz_4_pars = fit_Tomlinson_Menz_4_pars,
		
		Tomlinson_Phillips_3_pars = fit_Tomlinson_Phillips_3_pars,
		
		Finstad_Jonsson_4_pars = fit_Finstad_Jonsson_4_pars,
		
		Ruiz_4_pars = fit_Ruiz_4_pars,
		
		Boatman_5_pars = fit_Boatman_5_pars,
		
		Hinshelwood_4_pars = fit_Hinshelwood_4_pars,
		
		Kumaraswamy_5_pars = fit_Kumaraswamy_5_pars,
		
		Vant_Hoff_4_pars = fit_Vant_Hoff_4_pars,
		
		Warren_Dreyer_3_pars = fit_Warren_Dreyer_3_pars,
		
		Atkin_3_pars = fit_Atkin_3_pars,
		
		Newbery_4_pars = fit_Newbery_4_pars,
		
		Taylor_Sexton_3_pars = fit_Taylor_Sexton_3_pars
	)
	
	# Remove models that have more parameters than the maximum 
	# requested.
	if ( max_pars < 7 )
	{
		for ( i in (max_pars + 1):7 )
		{
			model_list <- within(
				model_list, rm(
					list = grep(
						paste('_', i, '_pars', sep = ''), 
						names(model_list), value = TRUE
					)
				)
			)
		}
	}
	
	return(model_list)
}

# This function takes all the thermal performance datasets and 
# separates them by id.
split_data <- function(dataset)
{
	
	# Prepare a list to store each dataset separately.
	split_dataset <- list()
	
	# Rename the 'temperature' column as 'temp' for simplicity.
	names(dataset)[names(dataset) == 'temperature'] <- 'temp'
	
	# Remove rows with missing or non-positive trait measurements.
	dataset <- dataset[!is.na(dataset$trait_value) & dataset$trait_value > 0,]

	# Get a vector of all the unique ids.
	unique_IDs <- unique(dataset$id)
	
	# For each id ...
	for ( i in unique_IDs )
	{
		
		# ... isolate the corresponding dataset.
		temp_dataset <- dataset[dataset$id == i,]
		
		# If the trait was expressed as time or duration, compute the 
		# inverse to get the rate.
		if (
			grepl('[^\\/]time', temp_dataset$trait[1], ignore.case = TRUE) ||
			grepl('duration$', temp_dataset$trait[1], ignore.case = TRUE)
		)
		{
			temp_dataset$trait_value <- 1/temp_dataset$trait_value
		}
		
		# Store the dataset to the list.
		split_dataset[[as.character(i)]] <- temp_dataset[,c('id', 'trait_value', 'temp')]
	}
	
	# Return the list of datasets.
	return(split_dataset)
}

####################
# M A I N  C O D E #
####################

# Check if there is a command-line argument of 'traits' or 'enzymes'.
# Otherwise, exit.
args <- commandArgs(TRUE)
if ( length(args) == 0 )
{
	stop("Usage: ./TPC_fitting_pipeline.R traits OR ./TPC_fitting_pipeline.R enzymes")
} else
{
	data_type <- args[1]
}

# Read all the thermal performance datasets.
if ( data_type == 'traits' )
{
	all_data <- read.csv('../Data/thermal_performance_datasets.csv', stringsAsFactors = FALSE)
} else if ( data_type == 'enzymes' )
{
	all_data <- read.csv('../Data/enzyme_thermal_performance_datasets.csv', stringsAsFactors = FALSE)
}

# Separate the datasets by id.
split_dataset <- split_data(all_data)

# Remove the original thermal performance datasets object to save memory.
rm(all_data)
gc()

# Prepare a data frame to store the fitting results (e.g., whether 
# only 3-parameter models could be fitted) to a particular dataset.
results_summary <- data.frame(
	ID = rep(NA, length(split_dataset)),
	outcome = rep(NA, length(split_dataset))
)

# Iterate over the thermal performance datasets and try to fit models.
for ( i in 1:length(split_dataset) )
{
	cat("Now at TPC ", i, "/", length(split_dataset), " ...\n", sep = '')
	
	results_summary$ID[i] <- split_dataset[[i]]$id[1]

	results_summary$outcome[i] <- fit_models(
		split_dataset[[i]], split_dataset[[i]]$id[1],
		data_type
	)
	gc()
}

# Write the fitting results to an output file.
if ( data_type == 'traits' )
{
	write.csv(
		results_summary, file = '../Results/results_summary.csv',
		row.names = FALSE
)
} else if ( data_type == 'enzymes' )
{
	write.csv(
		results_summary, file = '../Results/enzyme_results_summary.csv',
		row.names = FALSE
	)
}
