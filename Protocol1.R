
# sourcing the protocols required functions
# Packeges required:
# ape, parallel, Rfast, Rfast2
source('scripts/mitoP.R')

## To create the mitoP "full_data"
#1 reading the individual alignments
full_data_loci <- lapply(X = dir(path = 'main_datasets/mitoP_original_seqIDs/',
                            pattern = '\\.fast$',
                            full.names = TRUE),
  function(x) {
    x <- as.matrix(read.FASTA(file = x,
                              type = 'AA'));
    x[sort(labels(x)),]
    })

#2 Remove sequence unique identifiers prior to concatenation
full_data_loci <- lapply(full_data_loci, function(z) {
  updateLabel(x = z,
              old = rownames(z),
              new = sub("^((?:[^_]+_){3}).*", "\\1", rownames(z)))
})

#3 Concatenate
mitoP_full_data <- catAAbin(full_data_loci)

# Write 'mitoP_full_data.fasta' file:
write.FASTA(x = mitoP_full_data,
            file = "main_datasets/mitoP_full_data.fasta")

## To Generate the mitoP "full_data + P1 mask"
  #1) create the subsamples set: 
    # For reproducibility, the subsets and their inferred ML trees used in the
    #  paper are saved in the object mitoP_P1_subsets.RData in the path 
    #  (intermediate/P1/subsets)
mitoP_P1_subsets <- subsamG(x = full_data_loci, y = 14, z = 2000) # This will generate 1976 subsamples

## To reproduce the analyses of the paper:
load('intermediate_data/P1/P1_subsets_trees.RData')

P1_results <- txerr(tr = P1_subsets_trees[[4]],
                    subg = P1_subsets_trees,
                    rm = 0.10,
                    rgn = 0.5,
                    rtx = 0.5,
                    outrgx = "^ub") 

# Write 'mitoP_P1_mask.fasta' file:
write.FASTA(x = P1_results[[1]],
            file = "main_datasets/mitoP_P1_data.fasta")

# Write updated single gene alignment
dir.create( p1out <- 'P1_output')
gaps <- charToRaw(c('-?xX'))
invisible(lapply(seq_along(P1_results[[2]]), function(z) {
  z1 <- P1_results[[2]][[z]]
  z2 <- apply(z1, 1, function(z3) !all(z3 %in% gaps) )
  z1 <- z1[z2, ]
  write.FASTA(x = z1, file = sprintf('P1_output/p1mitoP%02d.fast', z))
}))
