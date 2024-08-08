#' Graphical display of `branchcutting` results.
#' 
#' `plot` method for class "bcut".
#' 
#' @param bc an object of class "bcut".
#' @param which a type of plot to be created, options are "tree" (default), "species", "score" and "density",
#'   unambiguous abbreviations are accepted.
#' @param edge.width an object of class "bcut".
#' @param show.tip.label logical, whether to show the tip labels (as in `plot.phylo`).
#' @param tip.cex factor scaling tip labels when `show.tip.label=TRUE`.
#' @param pt.cex factor scaling points when `which="species"`.
#' @param show.legend not implemented yet.
#' @param ... Other graphical parameters.
#' @export

plot.bcut <- function(bc, which="tree", edge.width=2, show.tip.label=FALSE, tip.cex=0.7, pt.cex=2, show.legend=FALSE, ...) {

	draw_phylogram <- function(phy, xy, lwd, col, ...) {
		lwd <- rep_len(lwd, nrow(phy$edge))
		col <- rep_len(col, nrow(phy$edge))
		for (i in seq(nrow(phy$edge))) {
			edg <- xy[phy$edge[i,],]
			lines(matrix(edg[c(1,2,4,4)],2,2), lwd=lwd[i], col=col[i], ...)
			lines(matrix(edg[c(1,1,3,4)],2,2), lwd=lwd[i], col=col[i], ...)
		}
	}
	draw_unrooted <- function(phy, xy, lwd, col, ...) {
		lwd <- rep_len(lwd, nrow(phy$edge))
		col <- rep_len(col, nrow(phy$edge))
		for (i in seq(nrow(phy$edge))) {
			edg <- xy[phy$edge[i,],]
			lines(edg[,1], edg[,2], lwd=lwd[i], col=col[i], ...)
		}
	}

	which <- grep(which, c("tree", "species", "score", "density"), value=TRUE)
	if (which == "tree") {
		col <- ifelse(bc$criteria$backbone, 3, 1)
		col <- ifelse(bc$criteria$outliers, 2, col)
		col <- ifelse(!bc$criteria$outliers & bc$criteria$heterosp, 4, col)
		col <- ifelse(bc$criteria$outliers & bc$criteria$consp, 8, col)
		type <- ifelse(ape::is.rooted(bc$phy), "phylogram", "unrooted")
		draw_tree <- list(phylogram=draw_phylogram, unrooted=draw_unrooted)[[type]]
		if (type == "phylogram") {
			xy <- ape::plotPhyloCoor(bc$phy, type=type, direction="rightwards", use.edge.length=TRUE)
		}
		if (type == "unrooted") {
			nb.sp <- ape::node.depth(bc$phy)
			xy <- ape::unrooted.xy(Ntip=length(bc$phy$tip.label), Nnode=bc$phy$Nnode, edge=bc$phy$edge, edge.length=bc$phy$edge.length, nb.sp=nb.sp, rotate.tree=0)$M
		}
		xlim <- range(xy[,1])
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
			par(mai=c(0.22, 0.22, 0.22, 0.22))
		}
		if (isTRUE(show.tip.label)) {
			ntip <- length(bc$phy$tip.label)
			offset <- 0.25 * mean(bc$phy$edge.length[bc$phy$edge[,2] <= ntip])
			xlim[2] <- xlim[2] + offset
			strw <- tip.cex * max(strwidth(bc$phy$tip.label, units="inches", font=3))
			strw <- strw * 1.08 * diff(xlim) / (par("fin")[2] - sum(par("mai")[c(2,4)]))
			xlim[2] <- xlim[2] + strw
		}
		plot(xy, type="n", bty="n", axes=FALSE, ann=FALSE, xlim=xlim)
		draw_tree(bc$phy, xy, lwd=edge.width, col=col)
		if (isTRUE(show.tip.label)) {
			xytips <- xy[seq(ntip),]
			xytips[,1] <- xytips[,1] + offset
			text(xytips, labels=bc$phy$tip.label, cex=tip.cex, font=3, adj=c(0,0.25))
		}
	}

	if (which == "species") {
		draw_ellipse <- function(size, w, h) {
			rs <- seq(0, 2 * pi, len=200)
			pts <- scale(cbind(0.5 * w * cos(rs), 0.5 * h * sin(rs)), scale=FALSE)
			return(pts %*% diag(rep(size / w, 2)))
		}
		species <- split(bc$species[,1], bc$species[,2])
		single <- sapply(species, length) == 1
		anc <- bc$ancestors[names(species)]
		type <- c("unrooted","cladogram")[ape::is.rooted(bc$phy) + 1]
		if (type == "cladogram") {
			xy <- ape::plotPhyloCoor(bc$phy, type=type, direction="rightwards", edge.width=edge.width, use.edge.length=TRUE)
		}
		if (type == "unrooted") {
			nb.sp <- ape::node.depth(bc$phy)
			xy <- ape::unrooted.xy(Ntip=length(bc$phy$tip.label), Nnode=bc$phy$Nnode, edge=bc$phy$edge, edge.length=bc$phy$edge.length, nb.sp=nb.sp, rotate.tree=0)$M
		}
		offset <- mean(bc$phy$edge.length[bc$phy$edge[,2] <= length(bc$phy$tip.label)])
		circle <- draw_ellipse(size=offset, w=diff(range(xy[,1])), h=diff(range(xy[,2])))
		hulls <- vector("list", bc$nspecies)
		for (i in seq(bc$nspecies)) {
			tips <- match(bc$species[bc$species[,2] == names(species)[i],1], bc$phy$tip.label)
			nodes <- unique(unlist(lapply(tips, function(x) ape::nodepath(bc$phy, x, anc[i]))))
			nodes <- nodes[grDevices::chull(xy[nodes,])]
			hulls[[i]] <- lapply(nodes, function(n, circle) circle + rep(1, 200) %*% t(xy[n,]), circle=circle)
			hulls[[i]] <- do.call(rbind, hulls[[i]])
			hulls[[i]] <- hulls[[i]][grDevices::chull(hulls[[i]]),]
		}
		xylim <- do.call(rbind, lapply(c(list(xy), hulls), function(x) apply(x, 2, range)))
		xylim <- apply(xylim, 2, range)
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
			par(mai=c(0.22, 0.22, 0.22, 0.22))
		}
		plot(xy, type="n", bty="n", axes=FALSE, ann=FALSE, xlim=xylim[,1], ylim=xylim[,2])
		for (i in seq(bc$nspecies)) {
			polygon(hulls[[i]], col="grey80", border="transparent")
		}
		for (i in seq(nrow(bc$phy$edge))) {
			lines(xy[bc$phy$edge[i,],], lwd=edge.width, col="black")
		}
		if (any(single)) {
			off <- match(unlist(species[single]), bc$phy$tip.label)
			points(xy[off,,drop=FALSE], pch=16, cex=pt.cex, col="black")
		}
		if (any(!single)) {
			points(xy[anc[!single],,drop=FALSE], pch=16, cex=pt.cex, col="black")
		}
	}
	
	if (which == "score") {
		col <- ifelse(bc$criteria$backbone, 3, 1)
		col <- ifelse(bc$criteria$outliers, 2, col)
		col <- ifelse(!bc$criteria$outliers & bc$criteria$heterosp, 4, col)
		col <- ifelse(bc$criteria$outliers & bc$criteria$consp, 8, col)
		col <- col[order(bc$score)]
		score <- sort(bc$score)
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
		}
		par(mai=c(1.02, 0.92, 0.42, 0.42))
		plot(seq_along(score), score, type="n", xlab="Rank", ylab="Importance score", main="", cex.axis=1.2, cex.lab=1.3, ...)
		points(which(col == 1), score[col == 1], pch=16, col=1, ...)
		points(which(col == 3), score[col == 3], pch=16, col=3, ...)
		points(which(col == 4), score[col == 4], pch=16, col=4, ...)
		points(which(col == 8), score[col == 8], pch=16, col=8, ...)
		points(which(col == 2), score[col == 2], pch=16, col=2, ...)
	}
	
	if (which == "density") {
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
		}
		par(mai=c(1.02, 0.92, 0.42, 0.42))
		plot(bc$density, xlab="Importance score", main="", cex.axis=1.2, cex.lab=1.3, lwd=2, ...)
		if (any(bc$criteria$outliers)) {
			abline(v=bc$bpoint, col="red", lwd=2, lty=3, ...)
		}
	}
}
