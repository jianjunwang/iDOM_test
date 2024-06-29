
#' @title Calculate molecular taxonomic diversity indices
#' @description Calculate various diversity indices for molecular data.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions and columns represent individual molecules. Default: mol.data.
#' @param type A character vector specifying the diversity indices to calculate. Options: "Richness", "Shannon", "Simpson", "Invsimpson", "PielouEven", "Chao1". Default: c("Richness", "Shannon", "Simpson", "Invsimpson", "PielouEven", "Chao1").
#' @return A data frame containing calculated diversity indices for each community.
#' @seealso 
#'  \code{\link[vegan]{decostand}}, \code{\link[vegan]{diversity}}, \code{\link[vegan]{specpool}}
#' @rdname chemoDiv
#' @export 
#' @importFrom vegan decostand specnumber diversity estimateR
commTD <- function(mol.data, type = c("Richness", "Shannon", "Simpson", "Invsimpson", "PielouEven", "Chao1")) {
  
  norm.mol.data = vegan::decostand(mol.data, method = "total")
  
  # Validate the input types
  valid_types <- c("Richness", "Shannon", "Simpson", "Invsimpson", "PielouEven", "Chao1")
  if (!all(type %in% valid_types)) {
    stop("Invalid type provided. Valid types are: Richness, Shannon, Simpson, Invsimpson, PielouEven, and Chao1.")
  }
  
  # Initialize the result data frame
  molTD <- as.data.frame(matrix(data = NA, nrow = nrow(mol.data), ncol = length(type)))
  colnames(molTD) <- type
  
  # Define a list of functions for dynamic call
  diversity_functions <- list(
    Richness = function(data) vegan::specnumber(data),
    Shannon = function(data) vegan::diversity(data, index = "shannon"),
    Simpson = function(data) vegan::diversity(data, index = "simpson"),
    Invsimpson = function(data) vegan::diversity(data, index = "invsimpson"),
    PielouEven = function(data) vegan::diversity(data, index = "shannon") / log(vegan::specnumber(data)),
    Chao1 = function(data) {
      norm_mol_counts_100k = round(norm.mol.data * 100000)
      return(vegan::estimateR(norm_mol_counts_100k, method = "chao")[2,])
    }
  )
  
  # Iterate over the selected types and calculate the indices
  for (t in type) {
    if (t %in% names(diversity_functions)) {
      molTD[[t]] <- diversity_functions[[t]](norm.mol.data)
    }
  }
  
  return(molTD)
}
