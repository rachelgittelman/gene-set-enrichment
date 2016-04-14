## This code will build an empirical null model to determine the significance of overlap
## between two genes sets. It accomplishes the same things as a hypergeometric test.
## by Rachel Gittelman

####################################################################################################
## you will need to enter in the following variables in R:

## background_set: the list of all possible genes that could have had methylated cysteines
## foreground_set: the list of genes that did have methylated cysteines
## overlap_set: the second set of genes that had methylated cysteines that you want to determine if
                there is significant overlap with.
####################################################################################################
## here is some quick example data:

background_set <- seq(1,1000)
foreground_set <- seq(80,180)
overlap_set <- seq(1,100)

####################################################################################################
## Now construct the empirical null and determine where the actual overlap falls within that
## distribution. "pval.result" will be the final result, and "empirical_null" will be the 
## list of null overlaps. Right now the number of permutations is set to 1000, but you can set it to 
## 10000 if you'd like.

n_per <- 1000
actual_overlap <- length(intersect(foreground_set, overlap_set))
n_to_sample <- length(foreground_set)

sample_genes <- function(background_set, n_to_sample, overlap_set)
{
	permutation <- sample(background_set, size=n_to_sample, replace=F)
	n_overlap <- length(intersect(permutation,overlap_set))
	return(n_overlap)
}

empirical.pval <- function(x, dist) {sum(x < dist)/length(dist)}

empirical_null <- sort(unlist(replicate(n_per, sample_genes(background_set, n_to_sample, overlap_set))))

pval_result <- empirical.pval(actual_overlap, empirical_null)

####################################################################################################
## Now you can do the exact same thing with the hypergeometric test:

phyper(actual_overlap,len(overlap_set),1-len(overlap_set),n_to_sample, lower.tail=F)

