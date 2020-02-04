probClades <- function(phy, psample, complete=FALSE, species=NULL) {
	clades <- ape::prop.part(phy)
	clades <- lapply(clades, function(x) attr(clades, "labels")[x])
	clades <- lapply(clades, sort)
	part <- ape::prop.part(psample)
	pppart <- attr(part, "number") / length(psample)
	part <- lapply(part, function(x) attr(part, "labels")[x])
	part <- lapply(part, sort)
	ord <- sapply(seq_along(clades),  function(i) which(sapply(part, identical, y=clades[[i]])))
	pp <- pppart[ord]
	if (isTRUE(complete)) {
		if (!is.null(species)) {
			binarize <- function(x, maxn) as.numeric(seq(maxn) %in% x)
			part <- lapply(lapply(part, function(x) phy$tip.label[x]), match, table=species)
			part <- apply(sapply(part, binarize, maxn=length(species)), 2, paste, collapse="")
		}
		attr(pp, "part") <- data.frame(part=part, pp=ppart, stringsAsFactors=FALSE)
	}
	return(pp)
}



