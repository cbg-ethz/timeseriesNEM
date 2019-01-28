require(tidyverse)

#' Wrapper function for tsNEM to return tidy output.
#' 
#' @param nemPath File path to nem object.
#' @param obsLogDen See \code{tsNEM()}.
#' @return A tibble with node, contrast, and binary activity call.
tidy_tsNEM <- function(nemPath, obsLogDen) {
  tsNEM(nemObject = readRDS(nemPath), observationalLogDensities = obsLogDen) %>%
    as.data.frame() %>%
    rownames_to_column("node") %>%
    as_tibble() %>%
    gather(contrast, activity, -node)
}

#' Tidy adjacency matrix.
#' 
#' @param nemPath File path to nem object.
#' @param reduce Boolean whether to transitively reduce the graph.
#' @return Tidy adjacency matrix.
tidy_adjacency <- function(nemPath, reduce = FALSE) {
  adjacency <- as(readRDS(nemPath)[["graph"]], "matrix")
  if(reduce)
    adjacency <- transitive.reduction(adjacency)
  as.data.frame(adjacency) %>%
    rownames_to_column(var = "parent") %>%
    as_tibble() %>%
    gather(key = "child", value = "directed edge", -parent)
}

#' Tidy support from list of objects returned by tidy_adjacency().
#' 
#' @param list.of.tidy.adjacencies A list of \code{tidy_adjacency()} objects.
#' @return The average of each directed edge in the list of adjacency matrices.
tidy_support <- function(list.of.tidy.adjacencies) {
  bind_rows(list.of.tidy.adjacencies) %>%
    group_by(child, parent) %>%
    summarise(support = mean(`directed edge`)) %>%
    ungroup()
}

#' Remove "miR-" and "-[35]p" from the start and end of microRNA names.
#' 
#' @param string A (vector of) string(s).
#' @return Abbreviated microRNA names.
abbrv_miR <- function(string) {
  string %>% sub("^miR-", "", .) %>% sub("-[35]p$", "", .)
}
