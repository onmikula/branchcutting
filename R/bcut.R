#' Branchcutting on a single tree.
#' 
#' @description
#' `bcut` performs branchcutting on a single tree.
#' 
#' @param phy An object of class `phylo` representing phylogenetic tree with branch lengths.
#' @param heterosp A two-column matrix with rows containing tree tips constrained to be heterospecific.
#' @param consp A two-column matrix with rows containing tree tips constrained to be conspecific.
#' @param rescale Logical, whether to rescale branch lengths (default is TRUE).
#' @param bw A character string indicating bandwidth function or a number giving the bandwidth;
#'   defaults to `mdif`, a bandwith calculation rule designed to the problem in hand.
#' @param kernel A character string indicating function used in the kernel density estimation,
#'   defaults to `epanechnikov`, but see [stats::density()] for other options.
#' @param return.bcut Logical, whether to return `bcut` object (the full output)
#'   or just a list of delimited species.
#' @returns An object of class `bcut` or, if `return.bcut=FALSE`,
#'   a list of character vectors with tree tips assigned to the same species.
#'   
#'   The `bcut` object is a list with the following components:
#' * `species` a data frame with columns `ID` (tree tip labels)
#'   and `Species` (classification to species labelled BCUT1, BCUT2, ...)
#' * `phy` the input tree (as a `phylo` object)
#' * `score` a numeric vector with importance scores of the tree branches
#'   (listed in the same order as in `phy$edge`)
#' * `density` an object of class `density` returned by [stats::density()]
#'   in the outlier detection step of the branchcutting algorithm
#' * `bpoint` the breakpoint, an importance score value separating outliers from the rest
#' * `backbone` a numeric vector with indices of branches belonging to the interspecific backbone
#' * `criteria` a data frame with as many rows as the number of tree branches and four components
#'   with logical vectors. The first one (`backbone`) indicates which branches belong to the backbone
#'   and the others (`outliers`, `heterosp` and `consp`) upon which criteria they were included.
#' * `nspecies` the number of delimited species 
#' * `rescaled` the tree with rescaled branch length (if the rescaling was applied)
#' @examples
#' bcspecies <- bcut(dendromus)
#' @export
bcut <- function(phy, heterosp=NULL, consp=NULL, rescale=TRUE, bw="mdif", kernel="epanechnikov", return.bcut=TRUE) {

	if (isTRUE(rescale)) {
		rphy <- rescale_branchlengths(phy)
		score <- calculate_importance(rphy)
	} else {
		rphy <- NULL
		score <- calculate_importance(phy)
	}
	
	ntip <- length(phy$tip.label)
	nedge <- nrow(phy$edge)
	if (!is.null(heterosp)) {
		maxdeepanc <- sapply(seq(nrow(heterosp)), function(i) ape::getMRCA(phy, tip=heterosp[i,]))
		hetconstrain <- find_backbone(phy, nodes=maxdeepanc)
	} else {
		hetconstrain <- rep(FALSE, nedge)
	}
	if (!is.null(consp)) {
		conconstrain <- find_crowns(phy, tips=consp)
	} else {
		conconstrain <- rep(FALSE, nedge)
	}
	if (any(hetconstrain & conconstrain)) {
		stop("Heterospecific and conspecific constraints are in conflict.")
	}
	
	if (bw == "mdif") {
		bw <- mean(diff(sort(score)))
	}
	outliers <- find_outliers(score, bw=bw, kernel=kernel, backbone=hetconstrain)
	dens <- outliers$density
	bpoint <- outliers$bpoint
	outliers <- outliers$outliers
	retained <- sort(which((outliers | hetconstrain) & !conconstrain))
	if (length(retained) > 0) {
		retnodes <- unique(sort(phy$edge[retained,]))
		backbone <- find_backbone(phy, nodes=retnodes)
		ancestors <- setdiff(phy$edge[backbone, 2], phy$edge[backbone, 1])
		singletons <- ancestors <= ntip
		species <- vector("list", length(ancestors))		
		species[singletons] <- phy$tip.label[ancestors[singletons]]
		species[!singletons] <- lapply(ancestors[!singletons], function(i) ape::extract.clade(phy, node=i)$tip.label)
	} else {
		species <- list(phy$tip.label)
	}

	if (isTRUE(return.bcut)) {
		species <- data.frame(ID=unlist(species), Species=rep(seq_along(species), sapply(species, length)), stringsAsFactors=FALSE)
		species$Species <- paste0("BCUT", formatC(species$Species, format="d", flag=0, width=nchar(max(species$Species))))
		criteria <- data.frame(backbone=backbone, outliers=outliers, heterosp=hetconstrain, consp=conconstrain)
		output <- list(species=species, phy=phy, score=score, density=dens, bpoint=bpoint, backbone=which(backbone), criteria=criteria, nspecies=length(unique(species[,2])), rescaled=rphy)
		class(output) <- "bcut"
	} else {
		species <- lapply(species, sort)
		output <- species
	}
	return(output)
	
}


#' Rescaling of branch lengths.
#' 
#' Rescales branch lengths according to a proper Laplacian graph transform.
#' 
#' @param phy An object of class `phylo`.
#' @returns An object of class `phylo` with rescaled branch lengths.
#' @export
rescale_branchlengths <- function(phy) {
	nodes <- lapply(seq(max(phy$edge)), function(i) which(phy$edge[,1] == i))
	deg <- sapply(nodes, function(x) sum(phy$edge.length[x]))
	if (ape::is.rooted(phy)) {
		phy$edge.length <- phy$edge.length * deg[phy$edge[,1]]		
	} else {
		phy$edge.length <- phy$edge.length * deg[phy$edge[,1]] * deg[phy$edge[,2]]		
	}
	return(phy)
}


#' Calculation of importance scores.
#' 
#' Calculates importance score for every branch in the tree.
#' 
#' @param phy An object of class `phylo`.
#' @returns A numeric vector with importance scores of the tree branches.
#' @export
calculate_importance <- function(phy) {
	edge <- phy$edge
	ntip <- min(edge[,1]) - 1
	noff <- as.numeric(edge[,2] <= ntip)
	last <- which(noff == 1)
	while(any(noff == 0)) {
		term <- last[duplicated(edge[last,1])]
		term <- which(edge[,1] %in% edge[term,1])
		kept <- setdiff(last, term)
		last <- which(edge[,2] %in% edge[term,1])
		new <- tapply(noff[term], edge[term,1], sum)
		noff[last] <- new[order(edge[last,2])]
		last <- c(last, kept)
	}
	npath <- noff * (ntip - noff)
	score <- npath * phy$edge.length / choose(ntip, 2)
	return(score)
}


#' Detection of outlying branches.
#' 
#' Finds branches with outlying (very high) importance scores.
#'
#' @param score A numeric vector with importance scores of the tree branches.
#' @param bw The bandwidth parameter of the kernel density function.
#' @param kernel A character string indicating the kernel density function.
#' @param backbone A logical vector indicating branches constrained to be in the interspecific backbone.
#' @returns A list with the components `outliers` (logical vector indicating outlying branches),
#'   `density` (an object with class `density` returned by the [stats::density()]) and
#'   `bpoint` (the breakpoint, an importance score value separating outliers from the rest).
#' @import ape
#' @export
find_outliers <- function(score, bw, kernel, backbone) {
	dens <- suppressWarnings(stats::density(score[!backbone], bw=bw, kernel=kernel, cut=0))
	zero <- sqrt(.Machine$double.eps)
	nonzeros <- which(dens$y <= zero)
	if (length(nonzeros) > 0) {
		bpoint <- dens$x[min(nonzeros)]
	} else {
		bpoint <- max(score)
	}
	outliers <- score > bpoint
	result <- list(outliers=outliers, density=dens, bpoint=bpoint)
	return(result)	
}


#' Finding tree backbone.
#' 
#' Finds basal / central backbone of a tree anchored by specified nodes.
#' 
#' @param phy An object of class `phylo`.
#' @param nodes A numeric vector with node numbers.
#' @returns A logical vector indicating whether the tree branches are in the backbone.
#' @export
find_backbone <- function(phy, nodes) {	
	if (ape::is.rooted(phy)) {
		edge <- phy$edge
		last <- edge[,2] %in% nodes
		backbone <- last
		while(any(last)) {
			backbone <- backbone | last
			last <- edge[,2] %in% edge[last,1]
		}
		backbone <- backbone | edge[,1] %in% edge[backbone,1]  & !backbone
	} else {
		edge <- apply(apply(phy$edge, 1, sort), 2, paste, collapse="-")
		backbone <- rep(FALSE, length(edge))
		terminal <- sort(nodes)
		while (length(terminal) > 1) {
			path <- ape::nodepath(phy, terminal[1], terminal[length(terminal)])
			terminal <- setdiff(terminal, path)
			path <- cbind(path[-1], path[-length(path)])
			path <- apply(apply(path, 1, sort), 2, paste, collapse="-")
			backbone <- backbone | edge %in% path
		}
		if (length(terminal) == 1) {
			path <- ape::nodepath(phy, terminal, setdiff(nodes, terminal)[1])
			path <- cbind(path[-1], path[-length(path)])
			path <- apply(apply(path, 1, sort), 2, paste, collapse="-")
			backbone <- backbone | edge %in% path			
		}
		bbnodes <- unique(as.numeric(phy$edge[backbone,]))
		backbone <- phy$edge[,1] %in% bbnodes | phy$edge[,2] %in% bbnodes
	}
	return(backbone)
}


#' Finding tree crowns.
#' 
#' Finds tree crowns, meaning a set of subtrees consisting of the most recent
#' common ancestors of the specified tips and all their offspring.
#' 
#' @param phy An object of class `phylo`.
#' @param tips A two-column matrix with pairs of tips (indicated by names or numbers).
#' @returns A logical vector indicating whether the tree branches are in any crown.
#' @export
find_crowns <- function(phy, tips) {
	if (mode(tips) != "numeric") {
		tips[,] <- match(tips, phy$tip.label)
		mode(tips) <- "numeric"
	}
	nodes <- sapply(seq(nrow(tips)), function(i) ape::getMRCA(phy, tip=tips[i,]))
	edge <- phy$edge
	nedge <- nrow(edge)
	ntip <- length(phy$tip.label)
	crowns <- rep(FALSE, nedge)
	for (i in seq_along(nodes)) {
		path1 <- ape::nodepath(phy, nodes[i], tips[i,1])
		path2 <- ape::nodepath(phy, nodes[i], tips[i,2])
		paths <- rbind(cbind(path1[-length(path1)], path1[-1]), cbind(path2[-length(path2)], path2[-1]))
		crowns <- crowns | duplicated(rbind(edge, paths), fromLast=TRUE)[1:nedge]
		offshoots <- !crowns & edge[,1] %in% paths
		while (any(offshoots)) {
			crowns <- crowns | offshoots
			offshoots <- !crowns & edge[,1] %in% edge[offshoots,2]
		}
	}
	return(crowns)
}



#' @export
summary.bcut <- function(bc) {
	cat("Object of class 'bcut':\n")
	cat(paste(paste("\tNo. of species:", bc$nspecies), "\n", sep=""))
	cat(paste(paste("\tNo. of tree tips:", length(bc$phy$tip.label)), "\n", sep=""))
	cat(paste(paste("\tNo. of retained branches:", length(bc$backbone)), "\n", sep=""))
}

							
#' @export
print.bcut <- function(bc) {
	k <- bc$nspecies
	nn <- table(bc$species[,2])

	cat(paste(paste("Object of class 'bcut' with", k, "delimited species"), ":\n", sep=""))
	for (i in seq(k)) {
		six <- head(bc$species[,1][bc$species[,2] == names(nn)[i]])
		if (nn[i] > 6)
			six <- c(six, "...")
		species <- paste(paste("\t", paste("Species", i), ":", sep=""), paste(six, collapse=", "))
		cat(paste(species, "\n", sep=""))
	}
}
