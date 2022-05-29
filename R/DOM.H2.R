



#' @title The specialization index H2' of DOM-microbe associations
#' @description This function calculates the network-level specialization of all interacting trophic levels in DOM-microbe bipartite networks, including full, negative and positive networks.
#' @param Comm.Microbe Sample/Microbe matrix with samples in the rows and Microbial species in the columns (compositional/abundant data).
#' @param Comm.DOM Sample/DOM matrix with samples in the rows and DOM molecules in the columns (compositional/abundant data).
#' @param occurrence.threshold The threshold of retaining bacterial species or DOM molecules observed in more than half of the total samples. Default: 0.5.
#' @param threshold.r The threshold of including correlations between DOM molecules and bacterial species in bipartite networks. Defaulted: 0.3.
#' @param N Number of null models to be generated; defaults to 100 (more might be better, less probably not).
#' @param Null.model Null model type. Can be given as an integer or name: 1/"r2dtable", 2/"swap.web", 3/"vaznull", 4/"shuffle.web"; allows for partial match of names.
#' @return Returns a data.frame, which contains standardised H2', Observed H2', P-value and Network type. Standardised: Standardised H2'. It is standardised by using a null modelling approach. 
#' Observed: Observed H2'. It ranges between 0 (complete generalization) and 1 (complete specialization).p.value: The significance of differences between observed and random H2'.Network.type: Network type. Can be given as a name: Full (full network), 
#' Negative (negative network), Positive (positive network).
#' @details H2' is a network-level property to describe how much the two trophic 
#' levels are interacting with each other in a bipartite network. For example,
#' H2' is used to quantify the specialization of DOM-microbe associations at a 
#' network level. Specifically, elevated H2' values convey that there is a high 
#' degree of specialization between DOM and microbes. By contrast, lower H2' 
#' values reflect a more generalized bipartite network where different DOM 
#' molecules can be used by a large range of bacterial taxa.
#' @references Hu A, Choi M, Tanentzap AJ, Liu J, Jang K-S, Lennon JT, Liu Y, Soininen J, Lu X, Zhang Y, Shen J, Wang J. 
#'   Quantifying microbial associations of dissolved organic matter under global change. \emph{bioRxiv}, 2021.
#' @references Hu A, Choi M, Tanentzap AJ, Liu J, Jang K-S, Lennon JT, Liu Y, Soininen J, Lu X, Zhang Y, Shen J, Wang J. 
#'   Quantifying the associations between dissolved organic matter and microbes
#'   under global change. \emph{Nature communications}, 2022.
#' @examples 
#' \dontrun{
#' # Example data of a Microbial compositional table (50 samples by 100 bacterial species)
#' Microbial.data
#' # Example data of a DOM compositional table (50 samples by 100 DOM molecules)
#' DOM.data
#' # Calculation of H2' index
#' DOM.H2(Comm.Microbe = Microbial.data,
#'        Comm.DOM = DOM.data,
#'        occurrence.threshold = 0.5,
#'        threshold.r = 0.3,
#'        N = 100,
#'        Null.model = "swap.web")
#' }
#' @seealso 
#'  \code{\link[SpiecEasi]{sparcc}}
#'  \code{\link[bipartite]{nullmodel}}, \code{\link[bipartite]{networklevel}}
#' @rdname DOM.H2
#' @export 
#' @importFrom SpiecEasi sparcc
#' @importFrom bipartite nullmodel networklevel


DOM.H2 <- function(Comm.Microbe,Comm.DOM,occurrence.threshold = 0.5,threshold.r = 0.3,N = 100,Null.model = "swap.web") {
  library(vegan)
  library(bipartite)
  
  Comm.Microbe.total = 10000 * decostand(Comm.Microbe, method = "total")
  Comm.Microbe.total = round(Comm.Microbe.total, 0)  
  colnames(Comm.Microbe.total) = paste("Bac", colnames(Comm.Microbe.total), sep="_")
  
  Comm.DOM.total = 10000 * decostand(Comm.DOM, method = "total")
  Comm.DOM.total = round(Comm.DOM.total, 0)
  colnames(Comm.DOM.total) = paste("DOM", colnames(Comm.DOM.total), sep="_")
  
  ####
  Comm.Microbe.keep = Comm.Microbe.total[, colSums(Comm.Microbe.total > 0) >= nrow(Comm.Microbe.total) * occurrence.threshold]
  Comm.DOM.keep = Comm.DOM.total[, colSums(Comm.DOM.total > 0) >= nrow(Comm.DOM.total) * occurrence.threshold]
  
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
                         Sample = rownames(Comm.Bac.DOM.site))
      rownames(index) = NULL
      
      if (ii == 1 & jj == 1) {
        Index = index
      }else{
        Index = rbind(Index,index)
      }
      
      if (jj %% 5 == 0) {
        print(paste0("Sample:",jj,"; ","Null model type:",classes[ii]))
      }
    }
  }
  return(Index)
}





