#' Branchcutting species delimitation.
#' 
#' The function `branchcutting` delimits species on a single tree and/or
#' a posterior sample of trees using branch-cutting algorithm.
#' 
#' @param phy Phylogenetic tree with branch lengths supplied as an object of class `phylo`,
#'   a character string in newick format or a name of file with the tree in newick or nexus format.
#' @param psample A posterior sample of trees as a list of `phylo` obbjects, `multiPhylo` object,
#'   or a name of file with the trees in newick or nexus format.
#' @param heterosp A two-column matrix with rows containing tree tips constrained to be heterospecific
#'   or a name of tab-delimited file without header containing such matrix.
#' @param consp A two-column matrix with rows containing tree tips constrained to be conspecific.
#'   or a name of tab-delimited file without header containing such matrix.
#' @param rescale Logical, whether to rescale branch lengths (default is TRUE).
#' @param bw A character string indicating bandwidth function or a number giving the bandwidth;
#'   defaults to `mdif`, a bandwith calculation rule designed to the problem in hand.
#' @param kernel A character string indicating function used in the kernel density estimation,
#'   defaults to `epanechnikov`, but see [stats::density()] for other options.
#' @param outgroup A character vector with outgroup names or their unique identifier(s).
#' @returns An object of class `bcut` or, if only `psample` is supplied,
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
				first <- Filter(function(x) nchar(x) > 0, sub("^[[:blank:]]*", "", readLines(phy, n=10)))[1]
				newick <- substr(first, 1, 1) == "("
				if (isTRUE(newick)) {
					phy <- ape::read.tree(phy)
				} else {
					phy <- ape::read.nexus(phy)
				}
			}
		}
	} else {
		phy <- NULL
	}
	if (is.list(phy) & !inherits(phy, "phylo") & missing(psample)) {
		psample <- phy
		phy <- NULL
	} else if (is.list(phy) & !inherits(phy, "phylo") & !missing(psample)) {
		phy <- NULL
		warning("cannot work with 'phy' as it is not in 'phylo' format")
	} else if (missing(psample)) {
		psample <- NULL
	}
	if (inherits(phy, "phylo") & is.null(phy$edge.length)) {
		phy <- NULL
		warning("cannot work with 'phy' as it has no branch lengths")
	}
	if (inherits(psample[[1]], "phylo") & is.null(psample[[1]]$edge.length)) {
		psample <- NULL
		warning("cannot work with 'psample' as it has no branch lengths")
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
		if (length(heterosp) == 1) {
			heterosp <- read.delim(heterosp, header=FALSE)
		}
		heterosp <- as.matrix(heterosp)[,1:2,drop=FALSE]
	}
	if (!is.null(consp)) {
		if (length(consp) == 1) {
			consp <- read.delim(consp, header=FALSE)
		}
		consp <- as.matrix(consp)[,1:2,drop=FALSE]
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
		spdf <- data.frame(Species=names(species), RS=as.numeric(rs), Consensus=seq_along(species) %in% consdelim, row.names=NULL)
		map <- matrix(0, length(species), length(tips), dimnames=list(NULL, tips))
		for (i in seq_along(species)) {
			map[i,species[[i]]] <- 1
		}
		psample <- list(rs=spdf, classification=map)
		if (!is.null(phy)) {
			splist <- lapply(split(bc$species[,1], bc$species[,2]), sort)
			bc$rs <- data.frame(Species=names(splist), RS=as.numeric(rs[match(splist, species)]))
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
