#' Roots the tree and drops the outgroup.
#' 
#' The function performs outgroup-rooting of the tree and then drops the outgroup(s) out.
#' 
#' @param phy an object of class `phylo` or `multiPhylo`.
#' @param outgroup a character vector with outgroup names or their unique identifier(s).
#' @returns A rooted tree without outgroup(s).
#' @examples
#' rootanddrop(phy=aethomys_tree, outgroup="outgroup")
#' @export

rootanddrop <- function(phy, outgroup="outgroup") {
	rad <- function(phy, outgroup) {
		outgroups <- unlist(lapply(outgroup, grep, x=phy$tip.label, value=TRUE))
		outgroups <- sort(unique(c(outgroups, unlist(lapply(outgroup, grep, x=phy$tip.label, value=TRUE, fixed=TRUE)))))
		phy <- ape::root(phy, outgroup=outgroups, resolve.root=TRUE)
		phy <- ape::drop.tip(phy, tip=outgroups)
		ntip <- ape::Ntip(phy)
		if (!is.null(phy$node.label)) {
			root <- ntip + 1
			empty <- phy$edge[phy$edge[,1] == root, 2]
			empty <- empty[empty > ntip]
			empty <- empty[nchar(phy$node.label[empty - ntip]) == 0]
			if (length(empty) == 1) {
				swap <- c(root, empty) - ntip
				phy$node.label[swap] <- phy$node.label[rev(swap)]
			}
		}
		return(phy)
	}
	if (inherits(phy, "phylo")) {
		phy <- rad(phy, outgroup)
	} else if (inherits(phy, "multiPhylo") | inherits(phy[[1]], "phylo")) {
		phy <- lapply(phy, rad, outgroup=outgroup)
		class(phy) <- "multiPhylo"
	} else {
		stop("'phy' must be of class 'phylo' or 'multiPhylo'")
	}
	return(phy)
}
