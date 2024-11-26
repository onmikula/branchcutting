#' Roots the tree and drops the outgroup.
#' 
#' The function performs outgroup-rooting of the tree and then drops the outgroup(s) out.
#' 
#' @param phy an object of class `phylo` or `multiPhylo`.
#' @param outgroup a character vector with outgroup names or their unique identifier(s).
#' @returns A rooted tree without outgroup(s).
#' @export

rootanddrop <- function(phy, outgroup="outgroup") {
	rad <- function(phy, outgroups) {
		phy <- ape::root(phy, outgroup=outgroups, resolve.root=TRUE, edgelabel=TRUE)
		phy <- ape::drop.tip(phy, tip=outgroups)
		phy <- ape::read.tree(text=ape::write.tree(ape::rotateConstr(phy, constraint=phy$tip.label)))
		return(phy)
	}
	if (!inherits(phy, c("phylo", "multiPhylo"))) {
		stop("'phy' must be of class 'phylo' or 'multiPhylo'")
	} else if (inherits(phy, "phylo")) {
		tips <- phy$tip.label
	} else {
		tips <- phy[[1]]$tip.label
	}
	outgroups <- unlist(lapply(outgroup, grep, x=tips, value=TRUE))
	outgroups <- sort(unique(c(outgroups, unlist(lapply(outgroup, grep, x=tips, value=TRUE, fixed=TRUE)))))
	if (length(outgroups) == 0) {
		stop("there are no tree tips conforming the outgroup definition")
	}
	if (inherits(phy, "phylo")) {
		phy <- rad(phy, outgroups)
	} else if (inherits(phy, "multiPhylo") | inherits(phy[[1]], "phylo")) {
		phy <- lapply(phy, rad, outgroups=outgroups)
		class(phy) <- "multiPhylo"
	}
	return(phy)
}
