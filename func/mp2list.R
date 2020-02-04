mp2list <- function(mp) {
	return(lapply(seq_along(mp), function(i) mp[[i]]))
}
