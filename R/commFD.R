#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mol.data PARAM_DESCRIPTION
#' @param mol.trait PARAM_DESCRIPTION
#' @param trait.name PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname chemoFD
#' @export 

commFD <- function(mol.data, mol.trait, trait.names){
  library(fundiversity)
  
  norm.mol.data = vegan::decostand(mol.data, method = "total")
  
  fd.raoq = fd_raoq(sp_com = norm.mol.data, traits = mol.trait[, trait.names, drop = F])
  
  return(fd.raoq)
}
