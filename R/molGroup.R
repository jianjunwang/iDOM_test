
#' @title Molecular Grouping Analysis
#' @description This function groups molecular data based on specific traits.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions and columns represent individual molecules. Default: mol.data
#' @param mol.trait Data frame containing molecular trait data such as mass or other properties associated with each molecule. Default: mol.trait
#' @param group The grouping criteria, Default: 'RA'
#' @param Trans_threshold The threshold for transformations, Default: c(1, 10)
#' @return A data frame with grouped molecules.
#' @rdname molGroup
#' @export 
molGroup <- function(mol.data, mol.trait, group = "RA", Trans_threshold = c(1,10)) {
  
  library(dplyr)
  
  # Checking row names consistency
  if(!identical(colnames(mol.data), rownames(mol.trait))){
    stop("Mismatch in row names between mol.data and mol.trait")
  }
  
  if (group == "RA") {
    # Reactivity: labile: H/C ≥ 1.5; recalcitrant: H/C < 1.5
    mol.HC = mol.trait %>% 
      select(HtoC_ratio) %>% 
      mutate(peak = rownames(.)) %>% 
      relocate(peak, .before = HtoC_ratio)
    
    if(exists("Transformation.results")){
      peak.profile <- Transformation.results[[2]]
    } else {
      Transformation.results <- molTrans(mol.data = mol.data, mol.trait = mol.trait, error.term = 0.000010, type = "Dataset")
      peak.profile <- Transformation.results[[2]]
    }
    
    mol.trait.group <- peak.profile %>% 
      left_join(mol.HC, by = "peak") %>% 
      rename(HtoC_ratio = colnames(.)[4]) %>% 
      mutate(group = case_when(HtoC_ratio >= 1.5 & num.trans.involved.in > Trans_threshold[2] ~"Labile_Active",
                               HtoC_ratio >= 1.5 & num.trans.involved.in <= Trans_threshold[1] ~"Labile_Inactive",
                               HtoC_ratio < 1.5 & num.trans.involved.in > Trans_threshold[2] ~"Recalcitrant_Active",
                               HtoC_ratio < 1.5 & num.trans.involved.in <= Trans_threshold[1] ~"Recalcitrant_Inactive"))
    
  }else if(group == "El_comp"){
    mol.trait.group <- mol.trait %>% 
      select(El_comp) %>% 
      mutate(peak = rownames(.)) %>% 
      relocate(peak, .before = El_comp) %>% 
      mutate(group = El_comp)
      
  }else if(group == "Class")
    mol.trait.group <- mol.trait %>% 
      select(Class) %>% 
      mutate(peak = rownames(.)) %>% 
      relocate(peak, .before = Class) %>% 
      mutate(group = Class)
  
  return(mol.trait.group)
}
