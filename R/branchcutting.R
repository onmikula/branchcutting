#' Branchcutting species delimitation.
#' 
#' The function `branchcutting` delimits species on a single tree and/or
#' a posterior sample of trees using branch-cutting algorithm.
#' 
#' @param phy phylogenetic tree with branch lengths supplied as an object of class `phylo`,
#'   a character string in newick format or a name of file with the tree in newick or nexus format.
#' @param psample a posterior sample of trees as a list of `phylo` obbjects, `multiPhylo` object,
#'   or a name of file with the trees in newick or nexus format.
#' @param heterosp a two-column matrix with rows containing tree tips constrained to be heterospecific
#'   or a name of tab-delimited file containing such matrix.
#' @param consp a two-column matrix with rows containing tree tips constrained to be conspecific.
#'   or a name of tab-delimited file containing such matrix.
#' @param rescale logical, whether to rescale branch lengths (default is TRUE).
#' @param bw a character string indicating bandwidth function or a number giving the bandwidth;
#'   defaults to `mdif`, a bandwith calculation rule designed to the problem in hand.
#' @param kernel a character string indicating function used in the kernel density estimation,
#'   defaults to `epanechnikov`, but see [stats::density()] for other options.
#' @param outgroup a character vector with outgroup names or their unique identifier(s).
#' @returns an object of class `bcut` or, if only `psample` is supplied,
#'   a list with components `consdelim` (data frame with the consensus species delimitation) and
#'   `psample` (a list containing information about delimited species and their relative supports).	
#' @examples
#' hetspmat <- matrix(c("TZ30480", "TZ30767"), 1, 2)
#' conspmat <- matrix(c("RS1629", "TZ30767"), 1, 2)
#' bcspecies <- branchcutting(phy=aethomys_tree, psample=aethomys_psample, heterosp=hetspmat, consp=conspmat, outgroup="outgroup")
#' @export

branchcutting <- function(phy, psample, heterosp=NULL, consp=NULL, rescale=TRUE, bw="mdif", kernel="epanechnikov", outgroup=NULL) {

	if (!missing(phy)) {
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
	} else {
		phy <- NULL
	}
	if (!missing(psample)) {
		if (is.character(psample) & length(psample) == 1) {
			first <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(psample, n=10))
			first <- first[nchar(first) > 0][1]
			newick <- substr(first, 1, 1) == "("
			if (isTRUE(newick)) {
				psample <- ape::read.tree(psample)
			} else {
				psample <- ape::read.nexus(psample)
			}
		}
	} else if (is.list(phy) & !inherits(phy, "phylo")) {
		psample <- phy
		phy <- NULL
		warning("cannot work with 'phy' as it is not in 'phylo' format")
	} else {
		psample <- NULL
	}
	if (inherits(phy, "phylo")) {
		if (is.null(phy$edge.length)) {
			phy <- NULL
			warning("cannot work with 'phy' as it has no branch lengths")	
		}
	}
	if (inherits(psample[[1]], "phylo")) {
		if (is.null(psample[[1]]$edge.length)) {
			psample <- NULL
			warning("cannot work with 'psample' as it has no branch lengths")
		}
	}
	if (is.null(phy) & is.null(psample)) {
		stop("at least one of 'phy' and 'psample' must be a valid input")
	}
	
	if (!is.null(outgroup)) {
		if (!is.null(phy)) {
			phy <- rootanddrop(phy, outgroup=outgroup)		
		}
		if (!is.null(psample)) {
			psample <- rootanddrop(psample, outgroup=outgroup)		
		}
	}

	if (!is.null(heterosp)) {
		if (!ape::is.rooted(phy)) {
			heterosp <- NULL
			warning("Heterospecific constraints cannot be applied to unrooted tree.")
		} else {
			if (length(heterosp) == 1) {
				heterosp <- read.delim(heterosp, header=FALSE)
			} else if (is.vector(heterosp) & length(heterosp) == 2) {
				heterosp <- matrix(heterosp, 1, 2)
			}
			heterosp <- as.matrix(heterosp)[,1:2,drop=FALSE]
		}
	}
	if (!is.null(consp)) {
		if (length(consp) == 1) {
			consp <- read.delim(consp, header=FALSE)
		} else if (is.vector(consp) & length(consp) == 2) {
			consp <- matrix(consp, 1, 2)
		}
		consp <- as.matrix(consp)[,1:2,drop=FALSE]
	}

	if (!is.null(heterosp) | !is.null(consp)) {
		if (!is.null(phy)) {
			tips <- phy$tip.label
		} else if (!is.null(psample)) {
			tips <- psample[[1]]$tip.label
		}
		if (!is.null(heterosp)) {
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
	}

	if (!is.null(phy)) {
		bc <- bcut(phy=phy, heterosp=heterosp, consp=consp, rescale=rescale, bw=bw, kernel=kernel, return.bcut=TRUE)		
	}
	
	if (!is.null(psample)) {
		ps <- lapply(psample, bcut, heterosp=heterosp, consp=consp, rescale=rescale, bw=bw, kernel=kernel, return.bcut=FALSE)
		ps <- lapply(unlist(ps, recursive=FALSE), sort)
		rs <- table(sapply(ps, paste, collapse="-")) / length(psample)
		rs <- rs[order(rs, names(rs), method="radix", decreasing=c(TRUE, FALSE))]
		species <- strsplit(names(rs), "-")
		names(species) <- paste0("sp", formatC(seq_along(species), format="d", flag=0, width=nchar(length(species))))
		tips <- psample[[1]]$tip.label
		mat <- matrix(0, length(species), length(species))
		for (i in 1:(length(species) - 1)) {
			for (j in (i+1):length(species)) {
				mat[i,j] <- mat[i,j] <- as.numeric(length(intersect(species[[i]], species[[j]])) == 0)
			}
		}
		consdelim <- 1
		complete <- length(setdiff(tips, unlist(species[consdelim]))) == 0
		while (complete == FALSE) {		
			consdelim <- c(consdelim, min(which(colSums(mat[consdelim,,drop=FALSE] == 1) == length(consdelim))))
			complete <- length(setdiff(tips, unlist(species[consdelim]))) == 0
		}
		spdf <- data.frame(Species=names(species), Support=as.numeric(rs), Consensus=seq_along(species) %in% consdelim, row.names=NULL)
		map <- matrix(0, length(species), length(tips), dimnames=list(NULL, tips))
		for (i in seq_along(species)) {
			map[i,species[[i]]] <- 1
		}
		psample <- list(support=spdf, classification=map)
		if (!is.null(phy)) {
			splist <- lapply(split(bc$species[,1], bc$species[,2]), sort)
			bc$support <- data.frame(Species=names(splist), Support=as.numeric(rs[match(splist, species)]))
			bc$psample <- psample
		}
	}
	
	if (!is.null(phy)) {
		output <- bc
		class(output) <- "bcut"
	} else {
		consdelim <- data.frame(ID=unlist(species[spdf$Consensus]), Species=rep(names(species[spdf$Consensus]), sapply(species[spdf$Consensus], length)), stringsAsFactors=FALSE, row.names=NULL)
		output <- list(consdelim=consdelim, psample=psample)
	}
	return(output)
	
}
