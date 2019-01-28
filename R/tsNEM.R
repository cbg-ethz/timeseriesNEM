#' Internal function for making a row in the canonical ILP constraint matrix.
#'
#' @param edgeCoordinates Node indices of an edge in the adjacency matrix of the
#'   transitively reduced signalling graph.
#' @param dimension The total number of nodes in the signalling graph.
#' @return A row of length \code{dimension} for the canonical ILP constraint
#'   matrix.
rowMaker <- function(edgeCoordinates, dimension) {
  result <- rep(0, dimension)
  result[edgeCoordinates] <- c(1, -1)
  return(result)
}

#' Internal function for running the ILP.
#'
#' @param NamedEgeneProfile Effect profile as a named vector.
#' @param nemObject A \code{nem} object.
#' @param constraintMatrix The corresponding canonical ILP constraint matrix.
#' @return The ILP solution of the signalling states of the nodes in the graph.
#' @import lpSolve
runILP <- function(NamedEgeneProfile, nemObject, constraintMatrix) {
  # BEWARE: nem `inference = "score"' returns nested nemObject$mappos$mappos.
  obj.func.coefs <- sapply(nemObject[["control"]][["Sgenes"]], function(Sgene) {
    sum(NamedEgeneProfile[nemObject[["mappos"]][[Sgene]]], na.rm = TRUE)
  })
  ilp <- lp(direction = "max",
            objective.in = obj.func.coefs,
            const.mat = constraintMatrix,
            const.dir = rep("<=", nrow(constraintMatrix)),
            const.rhs = rep(0, nrow(constraintMatrix)),
            all.int = TRUE,
            all.bin = TRUE)
  activity <- ilp$solution
  names(activity) <- names(ilp$objective)
  return(activity)
}

#' Map an observational E-gene profile onto a NEM.
#'
#' @param nemObject A \code{nem} object.
#' @param observationalLogDensities A matrix (or vector) of effect profiles
#'   recorded as log densities. The (row) names must correspond to the E-gene
#'   names in \code{nemObject}. The column names should identify the
#'   observational profiles, for instance with time points.
#' @return A matrix of signalling states for the nodes (rows) in the various
#'   time points (columns).
#' @export
tsNEM <- function(nemObject, observationalLogDensities) {
  if(is.vector(observationalLogDensities)) {
    observationalLogDensities <- as.matrix(observationalLogDensities)
  }
  # Construct constraint matrix for the ILP (canonical form). We need one row
  # per edge even in the transitive reduction, since there may be branches.
  adjacency <- as(nemObject[["graph"]], "matrix")
  reduction <- transitive.reduction(adjacency)
  edges <- which(reduction == 1, arr.ind = TRUE)
  edges <- edges[order(edges[,1]),]
  constraint.matrix <- t(apply(edges, 1, rowMaker, nrow(reduction)))
  colnames(constraint.matrix) <- nemObject[["control"]][["Sgenes"]]
  # Run ILP on all effect profiles.
  apply(observationalLogDensities, 2, runILP, nemObject = nemObject,
        constraintMatrix = constraint.matrix)
}
