# Code taken from huva package published on github (Bonaguro et al, 2022)

plot_binned_gene_mod <- function(goi, huva_experiment) {
  
  if (class(huva_experiment)!= "huva_experiment") {
    print("use huva_experiment class object for reliable results")
  }
  
  container <- list()
  
  for (i in names(huva_experiment)) {
    
    for (j in names(huva_experiment[[i]][["data"]])) {
      
      if (sum(goi %in% rownames(huva_experiment[[i]][["data"]][[j]]))>=1) {
        
        tmp <- as.data.frame(huva_experiment[[i]][["data"]][[j]])
        tmp <- t(tmp[goi,])
        
        #q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        
        tmp <- merge(huva_experiment[[i]][["anno"]][[j]][, c("Row.names", "group")], tmp, by.x = "Row.names", by.y = "row.names")
        tmp <- suppressMessages(melt(tmp, value.name = "expression", variable.name = "gene"))
        
        tmp <- ggplot(tmp, aes(x = gene, y = expression, fill = group))+
          geom_boxplot() + xlab("") + ylab("expression") + ggtitle(paste(unlist(strsplit(j, "_"))[-1], collapse = "_")) + 
          stat_compare_means(method = "t.test",  label = "p.signif", bracket.size = 1) +
          scale_fill_manual(values = c("#6b1414", "#005596")) + 
          ylab("vst expression")+
          theme_bw() + 
          theme(panel.grid = element_blank())
        
        container[[j]] <- tmp
        
      }
      
    }
    
  }
  
  return(container)
  
}