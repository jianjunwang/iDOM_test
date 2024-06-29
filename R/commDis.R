
#' @title Calculate dissimilarity between molecular community
#' @description Calculate dissimilarity metrics between molecular communities based on given method.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions and columns represent individual molecules. Default: mol.data
#' @param method Dissimilarity method to use. Options: "jaccard", "braycurtis", "unweighted_unifrac", "weighted_unifrac". Default: "braycurtis".
#' @param mol.dendrogram The dendrogram object representing the relationships among molecules. Default: NULL
#' @return Dissimilarity matrix
#' @seealso 
#'  \code{\link[phyloseq]{phyloseq}}, \code{\link[phyloseq]{UniFrac}}
#' @rdname commDis
#' @export 
#' @importFrom phyloseq phyloseq UniFrac
commDis <- function(mol.data, method = "braycurtis", mol.dendrogram = NULL) {
  # Load required libraries inside the function
  library(vegan)
  library(phyloseq)
  
  # Check if mol.data is a data frame
  if (!is.data.frame(mol.data)) {
    stop("Input 'mol.data' must be a data frame.")
  }
  
  # Check method parameter
  if (!method %in% c("jaccard", "braycurtis", "unweighted_unifrac", "weighted_unifrac")) {
    stop("Unsupported dissimilarity method. Supported methods: 'jaccard', 'braycurtis', 'unweighted_unifrac', 'weighted_unifrac'.")
  }
  
  # Calculate dissimilarity based on chosen method
  if (method == "jaccard") {
    dis <- vegdist(mol.data, method = "jaccard")
  } else if (method == "braycurtis") {
    dis <- vegdist(mol.data, method = "bray")
  } else if (method == "unweighted_unifrac" || method == "weighted_unifrac") {
    # Check if mol.dendrogram is provided
    if (is.null(mol.dendrogram)) {
      stop("UniFrac dissimilarity calculation requires a dendrogram 'mol.dendrogram'.")
    }
    
    # Convert mol.data to phyloseq object
    ps <- phyloseq::phyloseq(
      otu_table(as.matrix(mol.data), taxa_are_rows = F),
      phy_tree(mol.dendrogram)
    )
    
    # Calculate UniFrac dissimilarity
    if (method == "unweighted_unifrac") {
      dis <- phyloseq::UniFrac(ps, weighted = FALSE)
    } else {  # method == "weighted_unifrac"
      dis <- phyloseq::UniFrac(ps, weighted = TRUE)
    }
    
  } else {
    stop("Unsupported dissimilarity method.")
  }
  
  return(dis)
}
