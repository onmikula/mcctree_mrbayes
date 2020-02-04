### TOOLS
library(ape)
library(phangorn)
source("func/comm_anc_heights.R")
source("func/mp2list.R")
source("func/probClades.R")
source("func/writeNexus.R")


### GENUS NAME
genus <- "Fukomys"

# output file
output <- paste("trees", paste(genus, "trees", sep="."), sep="/")
# burn-in as a fraction of posterior sample
burnin <- 0.2

# combining trees
treefiles <- list.files("mb", pattern=paste(genus, "run", sep="."), full.names=TRUE)
treefiles <- treefiles[grepl("t$", treefiles)]
trees <- lapply(lapply(treefiles, ape::read.nexus), function(x) x[-seq(ceiling(burnin * length(x)))])
trees <- unlist(lapply(trees, mp2list), recursive=FALSE)
writeNexus(trees, output)

# rooting & discarding of outgroups
outgroup_file <- paste("mb", paste(strsplit(genus, "_")[[1]][1], "outgroups.txt", sep="_"), sep="/")
source(outgroup_file)
trees <- ape::read.nexus(output)
trees <- lapply(trees, ape::root, outgroup=outgroups, resolve.root=TRUE)
trees <- lapply(trees, ape::drop.tip, tip=outgroups)
writeNexus(trees, gsub("\\.trees$", "_rooted.trees", output))

# MCC tree calculation
mcct <- phangorn::maxCladeCred(trees, tree=TRUE, part=NULL, rooted=TRUE)
mcct <- comm_anc_heights(mcct, trees)
mcct$node.label <- formatC(probClades(mcct, trees), format="f", digits=4)
mcct$node.label <- paste("[&posterior=", mcct$node.label, "]", sep="")

# export in nexus and newick format (newick is for mptp)
writeNexus(mcct, gsub("\\.trees$", "_mcct.tre", output))
nwk <- mcct[names(mcct) != "node.label"]
class(nwk) <- "phylo"
ape::write.tree(nwk, gsub("\\.trees$", "_mcct.nwk", output))

