# X - input file with phylogenetic trees or DNA alignment
#	trees: object of class "phylo", "multiPhylo" or a list of "phylo" objects
#	DNA: character matrix or object of class "DNAbin" 
# file - name of output file

writeNexus <- function(X, file, translate=TRUE, tree.names=TRUE, part=NULL) {
	cls <- class(X)
	if (cls[1] == "list") {
		cls <- class(X[[1]])
	}
	if (cls %in% c("phylo", "multiPhylo")) {
		if (class(X) == "phylo") {
			taxa <- X$tip.label
			trees <- write.tree(X, file="")
		} else {
			taxa <- X[[1]]$tip.label
			trees <- unlist(lapply(X, write.tree, file=""))
		}
		n <- length(taxa)
		header <- c("Begin taxa;", paste("\tDimensions ntax=", n, ";", sep=""), "\t\tTaxlabels", paste("\t\t", taxa, sep=""), "\t;", "End;")
		if (isTRUE(translate)) {
			transl <- paste(paste("\t\t", seq(n), sep=""), paste(taxa, ",", sep=""))
			transl[length(transl)] <- sub("\\,$", "", transl[length(transl)])
			transl <- c("\tTranslate", transl, ";")
			for (i in seq(n)) {
				trees <- gsub(paste(taxa[i], ":", sep=""), paste(i, ":", sep=""), trees)
			}
		} else {
			transl <- NULL
		}
		if (isTRUE(tree.names)) {
			if(class(X) %in% c("list","multiPhylo") & !is.null(names(X))) {
				nam <- names(X)
			} else {
				nam <- paste("no", seq_along(trees), sep="_")
			}
			trees <- paste("tree", nam, "=", trees)
		}
		nexus <- c("#NEXUS", "\n", header, "\n", "Begin trees;", transl, trees, "End;")
	} else if (cls %in% c("matrix", "DNAbin")) {
		if (cls == "DNAbin") {
			X <- as.matrix(X)
		}
		dims <- dim(X)
		header <- c("#NEXUS\n", "BEGIN DATA;",
			paste("DIMENSIONS", "  ", "NTAX=", dims[1], " ", "NCHAR=", dims[2], ";", sep=""),
			"FORMAT DATATYPE=DNA GAP=~ MISSING=-;", "MATRIX\n")
		sequences <- c(rownames(X), unname(apply(X, 1, paste, collapse="")))
		sequences <- sequences[rep(seq(dims[1]), each=2) + (seq(2 * dims[1]) %% 2 == 0) * dims[1]]
		nexus <- c(header, sequences, ";\n", "END;")
		if (!is.null(part)) {
			p <- apply(part, 1, function(x) paste(x, collapse=" - "))
			p <- paste(names(p), p, sep=" = ")
			p <- paste(paste(paste("\t", "CHARSET", sep=""), p), ";", sep="")
			p <- c("BEGIN assumptions;", p, "END;")
			nexus <- c(nexus, "\n", p)
		}
	}
	writeLines(nexus, con=file, sep="\n", useBytes=FALSE)
}

