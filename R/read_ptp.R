#' Reads files produced by PTP / mPTP.
#' 
#' The function reads files with results of PTP / mPTP species delimitation.
#' 
#' @param file Character, name of the file.
#' @returns An object of class `spdelim`.
#' @export
read_ptp <- function(file) {
	lin <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(file))
	lin <- lin[nchar(lin) > 0]
	cmd <- grep("command", lin, ignore.case=TRUE, value=TRUE)
	method <- c("ptp", "mptp")[grepl("--multi", cmd) + 1]
	start <- grep("^Species", lin) + 1
	end <- c(start[-1] - 2, length(lin))
	sp <- lapply(Map(":", start, end), function(ii) lin[ii])
	nsp <- length(sp)
	species <- data.frame(ID=unlist(sp), Species=rep(seq(nsp), sapply(sp, length)), stringsAsFactors=FALSE)
	output <- list(species=species, nspecies=nsp, method=method)
	class(output) <- "spdelim"
	return(output)
}


#' @export
print.spdelim <- function(x) {
	k <- x$nspecies
	nn <- table(x$species[,2])
	m <- paste0("'", x$method, "'")

	cat(paste0(paste("Object of class 'spdelim' with", k, "species delimited by", m, "method"), ":\n"))
	for (i in seq(k)) {
		six <- head(x$species[,1][x$species[,2] == names(nn)[i]])
		if (nn[i] > 6)
			six <- c(six, "...")
		species <- paste(paste("\t", paste("Species", i), ":", sep=""), paste(six, collapse=", "))
		cat(paste(species, "\n", sep=""))
	}
}
