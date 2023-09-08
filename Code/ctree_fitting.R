#!/usr/bin/env Rscript

# This script fits multi-output conditional inference trees to predict 
# AICc weights of all models based on 29 predictor variables.

library(partykit)

#####################
# F U N C T I O N S #
#####################

# This function calculates the R-squared value of a conditional inference 
# tree against training or testing data.
get_R_squared <- function(current_model, data, training_or_testing)
{
	if ( training_or_testing == 'training' )
	{
		return(
			summary(
				lm(
					unlist(
						as.vector(predict(current_model))
					) ~ unlist(
						as.vector(data[,1:83])
					)
				)
			)$r.squared
		)
	} else if ( training_or_testing == 'testing' )
	{
		return(
			summary(
				lm(
					unlist(
						as.vector(predict(current_model, newdata = data))
					) ~ unlist(
						as.vector(data[,1:83])
					)
				)
			)$r.squared
		)
	}
}

####################
# M A I N  C O D E #
####################

# Read the table of AICc weights of the 83 models for all datasets.
AICc_table <- read.csv('../Results/AICc_table.csv', row.names = 1)

# Read the dataset of predictor variables.
predictors <- read.csv('../Results/TPC_predictors_dataset.csv', row.names = 1)
predictors$kingdom <- as.factor(predictors$kingdom)
predictors$phylum <- as.factor(predictors$phylum)
predictors$broad_trait_group <- as.factor(predictors$broad_trait_group)
predictors$specific_trait_group <- as.factor(predictors$specific_trait_group)

# Split the data into training and testing subsets (0.8% and 0.2% of the 
# data, respectively).
set.seed(1337)
training_IDs <- sample(row.names(AICc_table), round(0.8 * nrow(AICc_table)))

# Add the predictor columns to the AICc weights table.
all_data <- cbind(AICc_table, predictors)
training_data <- all_data[training_IDs,]
testing_data <- all_data[!row.names(all_data) %in% training_IDs,]

# Fit 15 different conditional inference trees.

# Taxonomy predictors
set.seed(1)
ctree_taxonomy <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ kingdom + phylum'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Trait predictors
set.seed(2)
ctree_trait <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ broad_trait_group + specific_trait_group'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# TPC shape predictors
set.seed(3)
ctree_TPC_shape <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Sampling intensity predictors
set.seed(4)
ctree_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# TPC shape + sampling intensity predictors
set.seed(5)
ctree_TPC_shape_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Trait + TPC shape predictors
set.seed(6)
ctree_trait_TPC_shape <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ broad_trait_group + specific_trait_group + skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Trait + Sampling intensity predictors
set.seed(7)
ctree_trait_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ broad_trait_group + specific_trait_group + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Taxonomy + Sampling intensity predictors
set.seed(8)
ctree_taxonomy_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ kingdom + phylum + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Taxonomy + TPC shape predictors
set.seed(9)
ctree_taxonomy_TPC_shape <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ kingdom + phylum + skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Taxonomy + Trait predictors
set.seed(10)
ctree_taxonomy_trait <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ kingdom + phylum + broad_trait_group + specific_trait_group'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Taxonomy + Trait + TPC shape predictors
set.seed(11)
ctree_taxonomy_trait_TPC_shape <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
                        '~ kingdom + phylum + broad_trait_group + specific_trait_group + skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Trait + TPC shape + Sampling intensity predictors
set.seed(12)
ctree_trait_TPC_shape_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ broad_trait_group + specific_trait_group + skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Taxonomy + TPC_shape + Sampling intensity predictors
set.seed(13)
ctree_taxonomy_TPC_shape_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
                        '~ kingdom + phylum + skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Taxonomy + Trait + Sampling intensity predictors
set.seed(14)
ctree_taxonomy_trait_sampling_intensity <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ kingdom + phylum + broad_trait_group + specific_trait_group + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# All possible predictors
set.seed(15)
ctree_all_pred <- ctree(
	eval(
		parse(
			text = paste(
				paste(colnames(AICc_table), collapse = ' + '),
				'~ kingdom + phylum + broad_trait_group + specific_trait_group + skew_scalar + min_temp + max_temp + T_pk_value + rise_exponent + fall_exponent + min_value_rise_peak + min_value_fall_peak + n_data_points + n_data_points_rise + n_data_points_fall + n_unique_temps + n_unique_temps_rise + n_unique_temps_fall + rise_cor + fall_cor + sampled_TPC_width + sampled_rise_width + sampled_fall_width + median_distance_between_temps + median_distance_between_temps_rise + median_distance_between_temps_fall + median_distance_between_trait_vals + median_distance_between_trait_vals_rise + median_distance_between_trait_vals_fall'
			)
		)
	),
	data = training_data,
	control = ctree_control(
		maxdepth = 4,
		cores = 8
	)
)

# Calculate and store the R-squared values for all conditional inference 
# trees against the training dataset.
results_training_data <- data.frame(
	Model = c(
		'Taxonomy', 'Trait', 'TPC shape', 'Sampling intensity', 
		'TPC shape + Sampling intensity', 'Trait + TPC shape',
		'Trait + Sampling intensity', 'Taxonomy + Sampling intensity',
		'Taxonomy + TPC shape', 'Taxonomy + Trait',
		'Taxonomy + Trait + TPC shape', 
		'Trait + TPC shape + Sampling intensity',
		'Taxonomy + TPC shape + Sampling intensity',
		'Taxonomy + Trait + Sampling intensity',
		'Taxonomy + Trait + TPC shape + Sampling intensity'
	),
	R_squared = c(
		get_R_squared(ctree_taxonomy, training_data, 'training'),
		get_R_squared(ctree_trait, training_data, 'training'),
		get_R_squared(ctree_TPC_shape, training_data, 'training'),
		get_R_squared(ctree_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_TPC_shape_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_trait_TPC_shape, training_data, 'training'),
		get_R_squared(ctree_trait_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_taxonomy_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_taxonomy_TPC_shape, training_data, 'training'),
		get_R_squared(ctree_taxonomy_trait, training_data, 'training'),
		get_R_squared(ctree_taxonomy_trait_TPC_shape, training_data, 'training'),
		get_R_squared(ctree_trait_TPC_shape_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_taxonomy_TPC_shape_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_taxonomy_trait_sampling_intensity, training_data, 'training'),
		get_R_squared(ctree_all_pred, training_data, 'training')
	)
)

# Order the table according to the R-squared column.
results_training_data <- results_training_data[rev(order(results_training_data$R_squared)),]

# Write the results to an output file.
write.csv(results_training_data, file = '../Results/ctree_R_squared_training_data.csv', row.names = FALSE)

# Calculate the R-squared value of the best model against the testing 
# data, and store it to an output file.
sink(file = '../Results/R_squared_of_best_ctree_against_testing_data.txt')

cat(get_R_squared(ctree_trait_TPC_shape_sampling_intensity, testing_data, 'testing'))

sink()
