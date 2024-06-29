
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mol.data Sample/DOM matrix with samples in the rows and DOM molecules in the columns (compositional/abundant data).
#' @param dendrogram PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION, Default: c("DD", "MPD", "MNTD")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[picante]{pd}}, \code{\link[picante]{mpd}}, \code{\link[picante]{mntd}}
#' @rdname chemoDD
#' @export 
#' @importFrom picante pd mpd mntd

commDD <- function(mol.data, dendrogram, type = c("DD", "MPD", "MNTD")){
  library(picante)
  
  norm.mol.data = vegan::decostand(mol.data, method = "total")
  
  # Validate the input types
  valid_types <- c("DD", "MPD", "MNTD")
  if (!all(type %in% valid_types)) {
    stop("Invalid type provided. Valid types are: DD, MPD, MNTD.")
  }
  
  # Initialize the result data frame
  mol_DD <- as.data.frame(matrix(data = NA, nrow = nrow(mol.data), ncol = length(type)))
  colnames(mol_DD) <- type
  
  # Define a list of functions for dynamic call
  diversity_functions <- list(
    DD = function(data, dendrogram){
      combined <- match.phylo.comm(dendrogram, data)
      
      dendrogram.matched <- combined$phy
      data.matched <- combined$comm
      
      return(picante::pd(data.matched, dendrogram.matched)[,"PD"])
    },
    
    MPD = function(data, dendrogram){
      combined <- match.phylo.comm(dendrogram, data)
      
      dendrogram.matched <- combined$phy
      data.matched <- combined$comm
      
      dendrogram.dist <- cophenetic(dendrogram.matched)
      return(picante::mpd(data.matched, dendrogram.dist, abundance.weighted = T))
    },
    
    MNTD = function(data, dendrogram){
      combined <- match.phylo.comm(dendrogram, data)
      
      dendrogram.matched <- combined$phy
      data.matched <- combined$comm
      
      dendrogram.dist <- cophenetic(dendrogram.matched)
      return(picante::mntd(data.matched, dendrogram.dist, abundance.weighted = T)) 
    }
  )
  
  # Iterate over the selected types and calculate the indices
  for (t in type) {
    if (t %in% names(diversity_functions)) {
      mol_DD[[t]] <- diversity_functions[[t]](norm.mol.data, dendrogram)
    }
  }
  
  return(mol_DD)
  
}
