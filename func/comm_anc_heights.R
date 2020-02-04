comm_anc_heights <- function(phy, trees) {
	R <- Ntip(phy) + 1
	tips <- phy$tip.label
	torders <- lapply(lapply(trees, "[[", "tip.label"), order)
	theights <- matrix(, length(tips), length(trees))
	for (i in seq_along(trees)) {
		tpaths <- ape::nodepath(trees[[i]])
		tbranches <- lapply(tpaths, function(x) match(x, trees[[i]]$edge[,2]))
		theights[,i] <- sapply(tbranches, function(x) sum(trees[[i]]$edge.length[x], na.rm=TRUE))[torders[[i]]]
	}
	TH <- rowMeans(theights)
	TH <- TH[match(tips, sort(tips))]
	clades <- lapply(prop.part(phy), function(x) tips[x])
	nheights <- matrix(, length(clades), length(trees))
	for (i in seq_along(trees)) {
		anc <- sapply(clades, function(x) getMRCA(trees[[i]],x))
		npaths <- lapply(anc, function(a) ape::nodepath(trees[[i]], from=R, to=a))
		nbranches <- lapply(npaths, function(x) match(x, trees[[i]]$edge[,2]))
		nheights[,i] <- sapply(nbranches, function(x) sum(trees[[i]]$edge.length[x], na.rm=TRUE))
	}
	NH <- rowMeans(nheights)
	H <- c(TH, NH)
	phy$edge.length <- as.numeric(diff(t(matrix(H[phy$edge], Nedge(phy), 2))))
	return(phy)
}

