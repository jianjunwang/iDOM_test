
qp.ratios <- list()
qp.results <- list()


commProc <- function(mol.data, dendrogram, group){
  # Initialize lists to store results
  qp.ratios <- list()
  qp.results <- list()

  if (group == "RA") {
    
    mol.trait.group = molGroup
    
    for (ii in 1:4) {
      mol.group = mol.trait.group %>% 
        filter(group %in% groups[jj])
      
      mol.data.group = mol.data[,colnames(mol.data) %in% mol.group$peak];dim(mol.data.group)
      
      phylo = picante::match.phylo.comm(phy = dendrogram, comm = mol.data.group)
      tree = phylo$phy; comm = phylo$comm
      
      pd = cophenetic(tree)
      qp = qpen(comm = comm, pd = pd, rand.time = 1000, nworker = 10, ab.weight = F)
      
      # Store results
      qp.ratios[[length(qp.ratios) + 1]] <- data.frame(ratio = qp$ratio, group = groups[ii])
      qp.results[[length(qp.results) + 1]] <- data.frame(result = qp$result, group = groups[ii])
    }
  }else{
    phylo = picante::match.phylo.comm(phy = dendrogram, comm = mol.data)
    tree = phylo$phy; comm = phylo$comm
    
    pd = cophenetic(tree)
    qp = qpen(comm = comm, pd = pd, rand.time = 1000, nworker = 10, ab.weight = F)
    
    # Store results
    qp.ratios <- data.frame(ratio = qp$ratio)
    qp.results <- data.frame(result = qp$result)
  }
}