#' Exporting species delimitations.
#' 
#' Generic function, intended for exporting species delimitations into plain text tables and/or data frames.
#'
#' @export
export <- function(x) {
  UseMethod("export")
}



#' Export of branchcutting delimitation.
#' 
#' `export` method for class `bcut`.
#' 
#' @param bc An object of class `bcut`.
#' @param file A name for the tab-delimited output file.
#' @param return A logical, whether to return result as a data frame.
#' @param taxonomy A data frame with an existing classification, informing about species names.
#' @param col.names Preferred column names (1st must refer to tree tips, 2nd to the species).
#' @export
export.bcut <- function(bc, file=NULL, return=TRUE, taxonomy=NULL, col.names=c("ID", "Species")) {
	species <- bc$species
	if (!is.null(taxonomy)) {
		taxonomy <- taxonomy[match(species[,1], taxonomy[,1]), 2]
		tab <- table(species[,2], taxonomy)
		rows <- rowSums(tab != 0) == 1
		cols <- colSums(tab != 0) == 1
		rows <- names(which(rowSums(tab[rows,cols,drop=FALSE] != 0) == 1))
		cols <- names(which(colSums(tab[rows,cols,drop=FALSE] != 0) == 1))
		tab <- tab[rows,cols,drop=FALSE]
		map <- setNames(colnames(tab)[apply(tab > 0, 1, which)], rownames(tab))
		named <- species[,2] %in% names(map)
		species[named,2] <- map[species[named,2]]
	}
	species <- setNames(species, col.names)
	if (!is.null(file)) {
		write.table(species, file, row.names=FALSE, quote=FALSE, sep="\t")
	}
	if (isTRUE(return)) {
		attr(species, "nspecies") <- bc$nspecies
		return(species)
	}
}
