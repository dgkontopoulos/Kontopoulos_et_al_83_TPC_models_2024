#!/usr/bin/env Rscript

# This script calculates a median matrix of Euclidean distances across 
# TPC models, and infers from this a dendrogram based on hierarchical 
# clustering with average linkage (the UPGMA algorithm).

library(ape)
library(phytools)

# Read the fits of the dataset with an id of 264 to use as a template. 
# This dataset includes fits with the maximum number of model parameters.
load('../Results/model_fits/264.Rda')

# Prepare the template matrix.
template_distance_matrix <- dist_matrix
template_distance_matrix[!is.na(template_distance_matrix)] <- NA
rm(dist_matrix)

# Create a list to store the dataset-specific distance matrices.
distance_matrices <- list()

# Get all the datasets for which we have fits and go through each 
# dataset.
model_fits <- list.files('../Results/model_fits/', pattern = '*.Rda')

for ( i in 1:length(model_fits) )
{
	cat('Now at ', model_fits[i], '...\n', sep = '')
	
	# Load the fits for the dataset.
	load(paste('../Results/model_fits/', model_fits[i], sep = ''))
	
	# Store all measured pairwise distances that were inferred from this 
	# dataset.
	current_dist_matrix <- template_distance_matrix
	
	for ( j in 1:length(rownames(current_dist_matrix)) )
	{
		if ( rownames(current_dist_matrix)[j] %in% rownames(dist_matrix) )
		{
			for ( k in 1:length(colnames(current_dist_matrix)) )
			{
				if ( colnames(current_dist_matrix)[k] %in% colnames(dist_matrix) )
				{
					current_dist_matrix[
						rownames(current_dist_matrix)[j], 
						colnames(current_dist_matrix)[k]
					] <- dist_matrix[
						rownames(current_dist_matrix)[j], 
						colnames(current_dist_matrix)[k]
					]
				} else
				{
					next
				}
			}
		} else
		{
			next
		}
	}
	
	distance_matrices[[i]] <- current_dist_matrix
}

# Calculate the median of each cell across all the distance matrices, 
# where possible.
averaged_distance_matrix <- apply(
	array(
		do.call(cbind, distance_matrices), 
		dim = c( dim(distance_matrices[[1]]), length(distance_matrices) )
	),
	c(1,2), median, na.rm = TRUE
)

colnames(averaged_distance_matrix) <- colnames(distance_matrices[[1]])
rownames(averaged_distance_matrix) <- rownames(distance_matrices[[1]])

# Infer a dendrogram using the UPGMA algorithm, save it to an output file, 
# and also save the median distance matrix to an .Rda file.
tree <- as.phylo(hclust(as.dist(averaged_distance_matrix), method = 'average'))

write.tree(tree, file = '../Results/tree_of_models.nwk')

save(averaged_distance_matrix, file = '../Results/averaged_distance_matrix.Rda')
