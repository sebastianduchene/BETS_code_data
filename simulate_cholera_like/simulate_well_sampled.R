library(NELSI)
library(phangorn)
i <- 1 # Change for simulation replicates using a for loop
ntaxa <- 50
  # For the ultrametric tree
  tree <- rcoal(ntaxa)
  tree$edge.length <- tree$edge.length * 100 # Convert branch lengths for a more realistic population size
  root_height <- max(allnode.times(tree))
  root_height
  # sampling times are drawn from an exponential distro because the time-tree is ultrametric  
  sampling_times <- abs( round(sort( c(0,  rexp(ntaxa-1, 10 / root_height))) - 2009, 2) )
  sampling_times
  
  tree$tip.label <- paste(tree$tip.label, sampling_times, sep = "_")
  
  tree$root.edge 
  plot.tree.lines(tree, plot.new = T, line.type = "l", rotation.angle = 3*pi/2)
  
  # Simulate sampling times
  clock_rate <-  7e-7 * 2948589 / 5000 # Similar site patterns as empirical data
  clock_rate
  
  phylogram <- tree
  phylogram$edge.length <- phylogram$edge.length * clock_rate
  
  aln <- as.DNAbin(simSeq(phylogram, l = 5000))
  
  print(root_height)
  print(aln)
  print('variable sites')
  print(length(seg.sites(aln)))
  print(max(sampling_times) - root_height)
  cat("root height is: ", root_height, "\n")
  cat("Proportion of sampling times is: ", diff(range(sampling_times)) / root_height, "\n")
  cat("variable sites", length(seg.sites(aln)), "\n")
  write.dna(aln, file = paste0("ultrametric_sample_clock_rate0.0004_sim", i, ".fasta"), format = "fasta", 
            nbcol = -1, colsep = "")
  write.tree(tree, file = paste0("ultrametric_sample_clock_rate0.0004_sim", i, ".tree"))
