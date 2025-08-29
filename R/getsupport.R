#' Branchcutting support of species. 
#' 
#' Extracts support for specified species from branchcutting results.
#' 
#' @param species Species for which support is to be assessed. Either a two-column matrix or data frame
#'   (1st column tip labels, 2nd column species) or a (list of) vector(s) with tip labels
#'   defining the species.
#' @param bc An object of class `bcut`.
#' @details The input `bcut` object must include component `psample`.
#' @export
getsupport <- function(species, bc) {

	if (!is.null(bc$psample)) {
		psample <- bc$psample
	} else if (all(c("support", "classification", "consensus") %in% names(bc)) ) {
		psample <- bc
	} else {
		stop("The input 'bcut' object must contain the component 'psample'.")
	}

	if (isTRUE(ncol(species) >= 2)) {
		species <- split(as.character(species[,1]), as.character(species[,2]))
	} else if (is.list(species)) {
		if (is.null(names(species))) {
			names(species) <- paste0("querysp", formatC(seq_along(species), format="d", flag=0, width=nchar(length(species))))
		}
	} else if (is.character(species)) {
		species <- list(querysp=species)
	}
	tips <- colnames(psample$classification)
	species <- lapply(species, function(sp, tips) as.numeric(tips %in% sp), tips=tips)

	k <- length(species)
	support <- setNames(numeric(k), names(species))
	for (i in seq(k)) {
		dupl <- duplicated(rbind(species[[i]], psample$classification))
		if (any(dupl)) {
			sp <- rownames(psample$classification)[which(dupl) - 1]
			support[i] <- as.numeric(psample$support$Support[psample$support$Species == sp])
		}
	}
	return(support)

}


#' Consensus species delimitation. 
#' 
#' Extracts onsensus species delimitation from a branchcutting analysis. 
#' 
#' @param bc An object of class `bcut`.
#' @details The input `bcut` object must include component `psample`.
#' @export
getconsensus <- function(bc) {

	if (inherits(bc, "bcut")) {
		consensus <- bc$psample$consensus
	} else if (all(c("support", "classification", "consensus") %in% names(bc)) ) {
		consensus <- bc$consensus
	} else {
		stop("The input 'bcut' object must contain the component 'psample'.")
	}
	return(consensus)
	
}


