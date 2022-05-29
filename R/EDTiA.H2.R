



#' @title The specialization index H2' of DOM-microbe associations
#' @description This function calculates the network-level specialization of all interacting trophic levels in DOM-microbe bipartite networks, including full, negative and positive networks.
#' @param Comm.Microbe PARAM_DESCRIPTION
#' @param Comm.DOM PARAM_DESCRIPTION
#' @param occurrence.threshold PARAM_DESCRIPTION, Default: 0.5
#' @param threshold.r PARAM_DESCRIPTION, Default: 0.3
#' @param N PARAM_DESCRIPTION, Default: 100
#' @param Null.model PARAM_DESCRIPTION, Default: 'swap.web'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[SpiecEasi]{sparcc}}
#'  \code{\link[bipartite]{nullmodel}}, \code{\link[bipartite]{networklevel}}
#' @rdname DOM.H2
#' @export 
#' @importFrom SpiecEasi sparcc
#' @importFrom bipartite nullmodel networklevel



EDTiA.H2 <- function(Comm.Microbe,Comm.DOM,occurrence.threshold = 0.5,threshold.r = 0.3,N = 100,Null.model = "swap.web") {
  library(vegan)
  
  Comm.Microbe.total = 10000 * decostand(Comm.Microbe, method = "total")
  Comm.Microbe.total = round(Comm.Microbe.total, 0)  
  colnames(Comm.Microbe.total) = paste("Bac", colnames(Comm.Microbe.total), sep="_")
  
  Comm.DOM.total = 10000 * decostand(Comm.DOM, method = "total")
  Comm.DOM.total = round(Comm.DOM.total, 0)
  colnames(Comm.DOM.total) = paste("DOM", colnames(Comm.DOM.total), sep="_")
  
  ####
  Comm.Microbe.pa = decostand(Comm.Microbe, method = "pa")
  Comm.Microbe.keep = Comm.Microbe.total[, colSums(Comm.Microbe.pa) >= nrow(Comm.Microbe.pa) * occurrence.threshold]
  
  Comm.DOM.pa = decostand(Comm.DOM, method = "pa")
  Comm.DOM.keep = Comm.DOM.total[, colSums(Comm.DOM.pa) >= nrow(Comm.DOM.pa) * occurrence.threshold]
  
  if (identical(rownames(Comm.Microbe.keep),rownames(Comm.DOM.keep))) {
    Comm.Bac.DOM = cbind(Comm.Microbe.keep, Comm.DOM.keep)
  }
  
  ####
  sparcc.Bac.DOM <- SpiecEasi::sparcc(Comm.Bac.DOM, iter = 20, inner_iter = 10, th = 0.1)
  sparcc.cor.Bac.DOM = sparcc.Bac.DOM$Cor
  
  adj.cor.Bac.DOM = ifelse(abs(sparcc.cor.Bac.DOM) > threshold.r, sparcc.cor.Bac.DOM, 0)
  
  colnames(adj.cor.Bac.DOM) = colnames(Comm.Bac.DOM)
  rownames(adj.cor.Bac.DOM) = colnames(Comm.Bac.DOM)
  
  sel.adj.cor.Bac.DOM = adj.cor.Bac.DOM[grep("Bac_",rownames(adj.cor.Bac.DOM)),
                                        grep("DOM_",colnames(adj.cor.Bac.DOM))]
  
  cor.Bac.DOM.int = round(10000 * sel.adj.cor.Bac.DOM, 0)
  
  classes <- c("Full","Positive","Nagative")
  
  ii = 1;jj = 1;i = 1
  for (ii in 1:length(classes)) {
    class = classes[ii]
    
    if (class == "Full") {
      adj.cor.int = ifelse(abs(cor.Bac.DOM.int) > 0, abs(cor.Bac.DOM.int), 0)
    }else if(class == "Positive"){
      adj.cor.int = ifelse(cor.Bac.DOM.int > 0, cor.Bac.DOM.int, 0) 
    }else{
      adj.cor.int = ifelse(cor.Bac.DOM.int < 0, abs(cor.Bac.DOM.int), 0) 
    }
    
    nulls <- bipartite::nullmodel(adj.cor.int, N = N, method = Null.model)
    
    weighted.indices = c("H2")
    
    for (jj in 1:nrow(Comm.Bac.DOM)) {
      Comm.Bac.DOM.site = Comm.Bac.DOM[jj,,drop = F];dim(Comm.Bac.DOM.site)
      Comm.Bac.DOM.site = Comm.Bac.DOM.site[,colSums(Comm.Bac.DOM.site) > 0];dim(Comm.Bac.DOM.site)
      
      site.Bac <- Comm.Bac.DOM.site[,grep("Bac_",colnames(Comm.Bac.DOM.site))];dim(site.Bac)
      site.DOM <- Comm.Bac.DOM.site[,grep("DOM_",colnames(Comm.Bac.DOM.site))];dim(site.DOM)
      
      adj.cor.int.site = adj.cor.int[colnames(site.Bac),colnames(site.DOM)];dim(adj.cor.int.site)
      
      nulls.site = list()
      
      for (i in 1:N) {
        colnames(nulls[[i]]) = colnames(adj.cor.int)
        rownames(nulls[[i]]) = rownames(adj.cor.int)
        
        nulls.site[[i]] <- nulls[[i]][colnames(site.Bac),colnames(site.DOM)];dim(nulls.site[[i]])
        # print(dim(nulls.site[[i]]))
      }
      
      index.obs = bipartite::networklevel(adj.cor.int.site, weighted = T, index = weighted.indices) 
      index.null = unlist(sapply(nulls.site, bipartite::networklevel, weighted = T, index = weighted.indices))
      
      index.rand.mean = mean(index.null)
      index.rand.sd = sd(index.null)
      index.ses = (index.obs - index.rand.mean) / index.rand.sd
      praw = sum(index.null > index.obs) / length(index.null)
      index.obs.p = ifelse(praw > 0.5, 1-praw, praw)
      index = data.frame(Index = names(index.obs), Observed = index.obs, Standardised = index.ses, P.value = index.obs.p,
                         Network.type = class, 
                         Site = rownames(Comm.Bac.DOM.site))
      rownames(index) = NULL
      
      if (ii == 1 & jj == 1) {
        Index = index
      }else{
        Index = rbind(Index,index)
      }
      # print(jj)
    }
  }
  return(Index)
}
