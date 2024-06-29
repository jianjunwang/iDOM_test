
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mol.data PARAM_DESCRIPTION
#' @param uncharacterized.mol.data PARAM_DESCRIPTION
#' @param mol PARAM_DESCRIPTION
#' @param occu.rate PARAM_DESCRIPTION, Default: 0.5
#' @param bootstrap PARAM_DESCRIPTION, Default: 100
#' @param Network.size PARAM_DESCRIPTION, Default: 200
#' @param sparcc.R PARAM_DESCRIPTION, Default: sparcc.r.sig
#' @param sparcc.R.threshold PARAM_DESCRIPTION, Default: 0.3
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[igraph]{graph_from_adjacency_matrix}}, \code{\link[igraph]{V}}, \code{\link[igraph]{degree}}
#'  \code{\link[bipartite]{specieslevel}}
#'  \code{\link[dplyr]{rename}}
#' @rdname iDME
#' @export 
#' @importFrom igraph graph_from_adjacency_matrix V degree
#' @importFrom bipartite specieslevel
#' @importFrom dplyr rename
iDME <- function(mol.data, uncharacterized.mol.data, mol, occu.rate = 0.5, bootstrap = 100, Network.size = 200, sparcc.R = sparcc.r.sig, sparcc.R.threshold = 0.3){

  # Checking row names consistency
  if(identical(x = rownames(mol.data), y = rownames(uncharacterized.mol.data)) == FALSE){
    stop("Something is incorrect in your row names")
  }else{
    mol.data.all = cbind(mol.data, uncharacterized.mol.data)
    mol.data.all.RA = decostand(mol.data.all, method = "total")
    mol.data.all.round = round(100000 * mol.data.all.RA, 0)
  }

  mol.data.all.occu <- mol.data.all.round[, colSums(mol.data.all.round > 0) >= nrow(mol.data.all.round) * occu.rate];dim(mol.data.all.occu)

  dim(mol.data.all.occu)

  identical(colnames(mol.data.all.occu), rownames(sparcc.R))
  identical(colnames(mol.data.all.occu), colnames(sparcc.R))

  sparcc.R.thresh <- sparcc.R %>%
    mutate_all(~ifelse(abs(.) >= 0.3, abs(.), 0))

  sparcc.R.thresh[1:10,1:10]

  ii = 1; jj = 1
  for (ii in 1:nrow(mol.data.all.occu)) {
    mol.data.sample = mol.data.all.occu[ii, , drop = F];dim(mol.data.sample)
    mol.data.sample = mol.data.sample[, colSums(mol.data.sample) > 0, drop = F]; dim(mol.data.sample)

    mol.known.occu <- mol.data.sample[, colnames(mol.data.sample) %in% colnames(mol.data), drop = F];dim(mol.known.occu)
    mol.unknown.occu <- mol.data.sample[, colnames(mol.data.sample) %in% colnames(uncharacterized.mol.data), drop = F];dim(mol.unknown.occu)

    for (jj in 1:bootstrap) {
      # Network nodes
      Known.seq <- sample(ncol(mol.known.occu), Network.size, replace = FALSE)
      K1.seq <- sample(Known.seq, Network.size/2, replace = FALSE)
      K2.seq <- Known.seq[!Known.seq %in% K1.seq]

      Unknown.seq <- sample(ncol(mol.unknown.occu), Network.size/2, replace = FALSE)

      #
      Net.K1 = mol.known.occu[, K1.seq, drop = F]; dim(Net.K1)
      Net.K2 = mol.known.occu[, K2.seq, drop = F]; dim(Net.K2)
      Net.D = mol.unknown.occu[, Unknown.seq, drop = F]; dim(Net.D)

      # Known networks
      Known_Net.all = sparcc.R.thresh[c(colnames(Net.K1), colnames(Net.K2)), c(colnames(Net.K1), colnames(Net.K2))];dim(Known_Net.all)
      Known.all <- igraph::graph_from_adjacency_matrix(as.matrix(Known_Net.all), mode = "undirected", weighted = TRUE, diag = FALSE)
      Known.all.degree <- data.frame(
        nodes_id = igraph::V(Known.all)$name,Degree = igraph::degree(Known.all)
      );rownames(Known.all.degree) = NULL

      Net_Known.Inner.K1 <- Known_Net.all[colnames(Net.K1), colnames(Net.K1)];dim(Net_Known.Inner.K1)
      Net_Known.Inner.K2 <- Known_Net.all[colnames(Net.K2), colnames(Net.K2)];dim(Net_Known.Inner.K2)

      Known.Inner.K1 <- igraph::graph_from_adjacency_matrix(as.matrix(Net_Known.Inner.K1), mode = "undirected", weighted = TRUE, diag = FALSE)
      Known.Inner.K2 <- igraph::graph_from_adjacency_matrix(as.matrix(Net_Known.Inner.K2), mode = "undirected", weighted = TRUE, diag = FALSE)

      Known.Inner.degree <- data.frame(
        nodes_id = c(igraph::V(Known.Inner.K1)$name,igraph::V(Known.Inner.K2)$name),Degree = c(igraph::degree(Known.Inner.K1),igraph::degree(Known.Inner.K2))
      );rownames(Known.Inner.degree) = NULL

      Net_Known.Inter <- Known_Net.all[colnames(Net.K2), colnames(Net.K1)];dim(Net_Known.Inter)
      Known.Inter.higher <- bipartite::specieslevel(Net_Known.Inter, index="degree", level = "higher");dim(Known.Inter.higher)
      Known.Inter.lower <- bipartite::specieslevel(Net_Known.Inter, index="degree", level = "lower");dim(Known.Inter.lower)
      Known.Inter.degree = data.frame(nodes_id = c(rownames(Known.Inter.higher),rownames(Known.Inter.lower)), Degree = c(Known.Inter.higher$degree,Known.Inter.lower$degree))

      if (identical(Known.Inter.degree$nodes_id,Known.all.degree$nodes_id) == F) {
        Known.Inter.degree.0 = data.frame(nodes_id = Known.all.degree$nodes_id[!Known.all.degree$nodes_id %in% Known.Inter.degree$nodes_id],Degree = 0)
        Known.Inter.Degree = rbind(Known.Inter.degree,Known.Inter.degree.0)
      }else{
        Known.Inter.Degree = Known.Inter.degree
      }

      Net_Known_degree = Known.all.degree %>%
        full_join(Known.Inner.degree,by = "nodes_id") %>%
        full_join(Known.Inter.Degree,by = "nodes_id") %>%
        dplyr::rename(KK_all.degree = Degree.x, KK_inner.degree = Degree.y, KK_inter.degree = Degree) %>%
        select(KK_all.degree, KK_inner.degree, KK_inter.degree) %>%
        summarise_all(mean, na.rm = TRUE)
        # mutate(class = case_when(nodes_id %in% colnames(Net.K1) ~ "K1",
        #                          nodes_id %in% colnames(Net.K2) ~ "K2"))

      # DK_Networks -------------------------------------------------------------
      DK_Net.all = sparcc.R.thresh[c(colnames(Net.K1), colnames(Net.D)), c(colnames(Net.K1), colnames(Net.D))];dim(DK_Net.all)

      DK.all <- igraph::graph_from_adjacency_matrix(as.matrix(DK_Net.all), mode = "undirected", weighted = TRUE, diag = FALSE)
      DK.all.degree <- data.frame(
        nodes_id = igraph::V(DK.all)$name,Degree = igraph::degree(DK.all)
      );rownames(DK.all) = NULL

      Net_DK.Inner.K1 <- DK_Net.all[colnames(Net.K1), colnames(Net.K1)];dim(Net_DK.Inner.K1)
      Net_DK.Inner.D <- DK_Net.all[colnames(Net.D), colnames(Net.D)];dim(Net_DK.Inner.D)

      DK.Inner.K1 <- igraph::graph_from_adjacency_matrix(as.matrix(Net_DK.Inner.K1), mode = "undirected", weighted = TRUE, diag = FALSE)
      DK.Inner.D <- igraph::graph_from_adjacency_matrix(as.matrix(Net_DK.Inner.D), mode = "undirected", weighted = TRUE, diag = FALSE)

      DK.Inner.degree <- data.frame(
        nodes_id = c(igraph::V(DK.Inner.K1)$name,igraph::V(DK.Inner.D)$name),Degree = c(igraph::degree(DK.Inner.K1),igraph::degree(DK.Inner.D))
      );rownames(DK.Inner.degree) = NULL

      # Inter
      Net_DK.Inter <- DK_Net.all[colnames(Net.D), colnames(Net.K1)];dim(Net_DK.Inter)
      DK.Inter.higher <- bipartite::specieslevel(Net_DK.Inter, index="degree", level = "higher");dim(DK.Inter.higher)
      DK.Inter.lower <- bipartite::specieslevel(Net_DK.Inter, index="degree", level = "lower");dim(DK.Inter.lower)
      DK.Inter.degree = data.frame(nodes_id = c(rownames(DK.Inter.higher),rownames(DK.Inter.lower)), Degree = c(DK.Inter.higher$degree,DK.Inter.lower$degree))

      if (identical(DK.Inter.degree$nodes_id, DK.all.degree$nodes_id) == F) {
        DK.Inter.degree.0 = data.frame(nodes_id = DK.all.degree$nodes_id[!DK.all.degree$nodes_id %in% DK.Inter.degree$nodes_id],Degree = 0)
        DK.Inter.Degree = rbind(DK.Inter.degree, DK.Inter.degree.0)
      }else{
        DK.Inter.Degree = DK.Inter.degree
      }

      Net_DK_degree = DK.all.degree %>%
        full_join(DK.Inner.degree,by = "nodes_id") %>%
        full_join(DK.Inter.Degree,by = "nodes_id") %>%
        dplyr::rename(DK_all.degree = Degree.x,DK_inner.degree = Degree.y,DK_inter.degree = Degree) %>%
        select(DK_all.degree, DK_inner.degree, DK_inter.degree) %>%
        # mutate(class = case_when(nodes_id %in% colnames(Net.K1) ~ "K1",
        #                          nodes_id %in% colnames(Net.D) ~ "D")) %>%
        summarise_all(mean, na.rm = TRUE)

      test = cbind(Net_Known_degree, Net_DK_degree, bootstrap = jj, sample = rownames(mol.data.sample))

      if (ii == 1 & jj == 1) {
        iDME = test
      }else{
        iDME = rbind(iDME, test)
      }
      print(ii)
    }
    ii = 1
    for (ii in 1:length(samples)) {
      iDME.sample = iDME %>% filter(sample %in% samples[ii])
      
      iDME.test = ((mean(iDME.sample$DK_all.degree)/mean(iDME.sample$KK_all.degree))-1)*100
      iDME.intra = (((mean(iDME.sample$DK_inner.degree)-mean(iDME.sample$KK_inner.degree))/mean(iDME.sample$KK_all.degree)))*100
      iDME.inter = (((mean(iDME.sample$DK_inter.degree)-mean(iDME.sample$KK_inter.degree))/mean(iDME.sample$KK_all.degree)))*100
      
      iDME.sample.test = data.frame(iDME = iDME.test, iDME.intra, iDME.inter, sample = unique(iDME.sample$sample))
      
      if (ii == 1) {
        iDME.Sample = iDME.sample.test
      }else{
        iDME.Sample = rbind(iDME.Sample, iDME.sample.test)
      }
    }
    
    return(iDME.Sample)
  }
}
