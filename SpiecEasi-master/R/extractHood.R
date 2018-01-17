#' Spiec-Easi function to extract neighborhood surrounding nodes of interest.
#' Input is Spiec-Easi network and node or nodes of interest as vector such as c(1:5)
#' This code was generalized from the code used to make figure 4 in Tipton, Müller, et. al.
#' written by Laura Tipton; 1/5/18
#' @export

extract_hood <- function(selntwrk, node){
  refit <- selntwrk$refit
  assoc <- which(selntwrk$refit[,node[1]]==1)
  assoc <- c(node, assoc)
  if(length(node) > 1){
    for(nod in 2:length(node)){
      tempassoc <- which(selntwrk$refit[,node[nod]]==1)
      assoc <- unique(c(assoc, tempassoc))
    }
  }
  refit2 <- refit[assoc,assoc]
  colnames(refit2) = rownames(refit2) <- assoc
  return(refit2)
  rm(assoc)
}