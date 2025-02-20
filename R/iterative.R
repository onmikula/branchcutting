#' Iterative branchcutting.
#' 
#' Iterative application of the branchcutting species delimitation to the same tree
#' using delimitation from a previous run as heterospecific constraint to the next one.
#' 
#' @param phy phylogenetic tree with branch lengths supplied as an object of class `phylo`,
#'   a character string in newick format or a name of file with the tree in newick or nexus format.
#' @param heterosp as in [branchcutting], initial heterospecific constraint.
#' @param consp as in [branchcutting], initial conspecific constraint.
#' @param rescale logical, whether to rescale branch lengths (default is TRUE).
#' @param bw a character string indicating bandwidth function or a number giving the bandwidth;
#'   defaults to `mdif`, a bandwith calculation rule designed to the problem in hand.
#' @param kernel a character string indicating function used in the kernel density estimation,
#'   defaults to `epanechnikov`, but see [stats::density()] for other options.
#' @param outgroup a character vector with outgroup names or their unique identifier(s).
#' @returns a list with two components: `iterations` containing a list of `bcut` objects
#'   (results from subsequent iterations) and `species`, a data frame with the final delimitation.
#' @export
iterative <- function(phy, heterosp=NULL, consp=NULL, rescale=TRUE, bw="mdif", kernel="epanechnikov", outgroup=NULL) {
	
	if (is.character(phy) & length(phy) == 1) {
		newick <- grepl("^\\(.+;$", phy)
		if (isTRUE(newick)) {
			phy <- ape::read.tree(text=phy)
		} else {
			first <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(phy, n=10))
			first <- first[nchar(first) > 0][1]
			newick <- substr(first, 1, 1) == "("
			if (isTRUE(newick)) {
				phy <- ape::read.tree(phy)
			} else {
				phy <- ape::read.nexus(phy)
			}
		}
	} else {
		phy <- ape::as.phylo(phy)
	}
	if (!ape::is.rooted(phy)) {
		stop("Heterospecific constraints and hence iterative branchcutting cannot be applied to unrooted tree.")
	}
	if (!is.null(outgroup)) {
		phy <- rootanddrop(phy, outgroup=outgroup)
	}
	
	tips <- phy$tip.label
	if (!is.null(heterosp)) {
		if (length(heterosp) == 1) {
			heterosp <- read.delim(heterosp, header=FALSE)
		} else if (is.vector(heterosp) & length(heterosp) == 2) {
			heterosp <- matrix(heterosp, 1, 2)
		}
		heterosp <- as.matrix(heterosp)[,1:2,drop=FALSE]
		hetmiss <- matrix(!heterosp %in% tips, nrow(heterosp), 2)
		hetmiss <- which(apply(hetmiss, 1, any))
		if (length(hetmiss) > 0) {
			heterosp <- heterosp[-hetmiss,,drop=FALSE]
			hetmiss <- paste(hetmiss, collapse=", ")
			warning(paste("heterospecific constraint(s) no.", hetmiss, "missing among tree tips and removed"))
			if (nrow(heterosp) == 0) {
				heterosp <- NULL
			}
		}
	}
	if (!is.null(consp)) {
		if (length(consp) == 1) {
			consp <- read.delim(consp, header=FALSE)
		} else if (is.vector(consp) & length(consp) == 2) {
			consp <- matrix(consp, 1, 2)
		}
		consp <- as.matrix(consp)[,1:2,drop=FALSE]
		conmiss <- matrix(!consp %in% tips, nrow(consp), 2)
		conmiss <- which(apply(conmiss, 1, any))
		if (length(conmiss) > 0) {
			consp <- consp[-conmiss,,drop=FALSE]
			conmiss <- paste(conmiss, collapse=", ")
			warning(paste("conspecific constraint(s) no.", conmiss, "missing among tree tips and removed"))
			if (nrow(consp) == 0) {
				consp <- NULL
			}
		}
	}

	i <- 1
	iterations <- list(bcut(phy=phy, heterosp=heterosp, consp=consp, rescale=rescale, bw=bw, kernel=kernel, return.bcut=TRUE))
	convergence <- FALSE
	while (convergence == FALSE) {
		heterosp <- prepare_constraints(phy=phy, spdelim=iterations[[i]]$species, type="heterospecific")
		iterations[[i+1]] <- bcut(phy=phy, heterosp=heterosp, consp=NULL, rescale=rescale, bw=bw, kernel=kernel, return.bcut=TRUE)
		current <- iterations[[i]]$species[order(iterations[[i]]$species[,1]),2]
		current <- match(current, unique(current))
		candidate <- iterations[[i+1]]$species[order(iterations[[i+1]]$species[,1]),2]
		candidate <- match(candidate, unique(candidate))
		i <- i + 1
		convergence <- all(current == candidate)
	}
	iter <- list(iterations=iterations[-i], species=iterations[[i-1]]$species)
	
	return(iter)
	
}



#' Preparing constraints.
#' 
#' Prepares constraints for species delimitation using branchcutting algorithm.
#' 
#' @param phy phylogenetic tree with branch lengths supplied as an object of class `phylo`,
#'   a character string in newick format or a name of file with the tree in newick or nexus format.
#' @param spdelim an object of class `bcut` or a data frame or matrix
#'   with tree tips in the first column and their classification into species in the second
#'   or a name of tab-delimited file containing such data.
#' @param type character string specifies the type of constraint, either "heterospecific" (default)
#'   or "conspecific" or their unambiguous abbreviation.
#' @returns a two-column matrix specifying the constraints.
#' @export
prepare_constraints <- function(phy, spdelim, type="heterospecific") {

	if (inherits(phy, "bcut")) {
		phy <- phy$phy
	} else if (is.character(phy) & length(phy) == 1) {
		newick <- grepl("^\\(.+;$", phy)
		if (isTRUE(newick)) {
			phy <- ape::read.tree(text=phy)
		} else {
			first <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(phy, n=10))
			first <- first[nchar(first) > 0][1]
			newick <- substr(first, 1, 1) == "("
			if (isTRUE(newick)) {
				phy <- ape::read.tree(phy)
			} else {
				phy <- ape::read.nexus(phy)
			}
		}
	} else {
		phy <- ape::as.phylo(phy)
	}
	
	if (inherits(spdelim, "bcut")) {
		spdelim <- bcut$species
	} else if (is.character(spdelim) & length(spdelim) == 1) {
		spdelim <- read.delim(spdelim)[,1:2]
	} else {
		spdelim <- as.data.frame(spdelim[,1:2])
	}

	type <- grep(tolower(type), c("heterospecific", "conspecific"), value=TRUE)[1]
	tips <- intersect(phy$tip.label, spdelim[,1])
	spdelim <- spdelim[match(tips, spdelim[,1]),]
	if (type == "heterospecific") {
		spdelim <- spdelim[match(unique(spdelim[,2]),spdelim[,2]),1]
		constr <- t(combn(spdelim, 2))
		attr(constr, "type") <- "heterospecific"
	}
	if (type == "conspecific") {
		spdelim[,1] <- match(spdelim[,1], phy$tip.label)
		spdelim <- split(spdelim[,1], spdelim[,2])
		spdelim <- spdelim[sapply(spdelim, length) > 1]
		for (i in seq_along(spdelim)) {
			mrca <- ape::getMRCA(phy, spdelim[[i]])	
			paths <- lapply(spdelim[[i]], function(to) ape::nodepath(phy, from=mrca, to=to))
			basal <- sapply(paths, "[[", 2)
			spdelim[[i]] <- sapply(lapply(paths[match(unique(basal), basal)], rev), "[[", 1)
			spdelim[[i]] <- phy$tip.label[spdelim[[i]]]
		}
		constr <- as.matrix(do.call(rbind, spdelim))
		attr(constr, "type") <- "conspecific"
	}
	
	return(constr)
	
}
