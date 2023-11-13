# Code taken from huva package published on github (Bonaguro et al, 2022)

# Modification turns direction of comparison originally used in huva package (low vs high -> now: high vs low) 

run_huva_experiment_rev <- function(data = datasets, gene, quantiles, gs_list, summ = T, datasets_list = NULL, adjust.method = "none") {
  
  if (class(data)!= "huva_dataset") {
    error("Use huva_dataset class object to run the huva_experiment function.")
  }
  
  container <- list()
  
  print(paste("Binning on ", gene, " expression", sep = ""))
  
  if (is.null(datasets_list)==F) {
    # This is to allow the selection of the datasets to use in the analysis
    data <- data[datasets_list]
  }
  
  for (i in names(data)) {
    
    for (j in names(data[[i]][["data"]])) {
      
      if (gene %in% rownames(data[[i]][["data"]][[j]])) {
        
        expr <- as.data.frame(data[[i]][["data"]][[j]][gene,])
        colnames(expr) <- c("expression")
        
        # Calculation of the percentiles
        
        expr$group <- ifelse(expr$expression<= quantile(expr$expression, c(1-quantiles,quantiles), na.rm = T)[2], "low",
                             ifelse(expr$expression >= quantile(expr$expression, c(1-quantiles,quantiles), na.rm = T)[1],
                                    "high", "none"))
        
        #q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        q2 <- paste("DE", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        q3 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        q4 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        
        anno_tmp <- merge(x= expr[expr$group != "none",], y= data[[i]][["anno"]][[j]], by="row.names")
        
        data_tmp <- data[[i]][["data"]][[j]][, anno_tmp$Row.names]
        
        sample <- as.factor(anno_tmp$group)
        design.mat <- model.matrix(~0+sample)
        colnames(design.mat) <- levels(sample)
        
        contrast.matrix <- makeContrasts(Diff = high - low, levels = design.mat)
        
        fit <- lmFit (data_tmp, design.mat)
        fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(fit)
        
        DE_table <- topTable(fit, coef = "Diff", p.value = 1, adjust.method = adjust.method, lfc = log2(1), number = 100000)
        rank <- fit$coefficients[order(fit$coefficients[,1],decreasing = T),]
        
        gse <- suppressWarnings(fgsea(pathways = gs_list,
                                      stats = rank,
                                      minSize = 1,
                                      maxSize = Inf,
                                      nperm = 1000))
        
        container[[i]][["anno"]][[paste(i,j, sep = "_")]] <- anno_tmp
        container[[i]][["data"]][[paste(i,j, sep = "_")]] <- data_tmp
        container[[i]][["DE_genes"]][[paste(i,j, sep = "_")]] <- DE_table
        container[[i]][["Rank_genelist"]][[paste(i,j, sep = "_")]] <- rank
        container[[i]][["gsea"]][[paste(i,j, sep = "_")]] <- gse
        
        if (summ==T) {
          
          container[["summary"]][["Rank"]][[paste(i,j, sep = "_")]] <- rank
          container[["summary"]][["gsea"]][[paste(i,j, sep = "_")]] <- gse
          
          # Add the summary of the anno here
          tmp <- anno_tmp
          
          cont_anno <- list()
          
          for (n in colnames(tmp)[-c(1,2,3)]) {
            
            if (is.numeric(tmp[[n]])==T) {
              list <- tmp[,c("group", n)]
              list <- t.test(tmp[tmp$group=="high",][[n]], tmp[tmp$group=="low",][[n]], paired = F, var.equal = T, )
            }
            if (is.numeric(tmp[[n]])==F) {
              list <- tmp[,c("group", n)]
              list <- table(list)
              list <- prop.table(list, margin = 1)*100
            }
            
            cont_anno[[n]] <- list
          }
          
          container[["summary"]][["anno"]][[paste(i,j, sep = "_")]] <- cont_anno
          
        }
        
        if (length(names(data[[i]][["metadata"]]))>0) {
          
          for (k in names(data[[i]][["metadata"]])) {
            
            tmp_metadata <- merge(expr[expr$group != "none",], data[[i]][["metadata"]][[k]], by = "row.names")
            
            container[[i]][["metadata"]][[paste(i, j, k, sep = "_")]] <- tmp_metadata
            
            if (summ==T) {
              
              tmp_metadata$expression <- NULL
              tmp_metadata <- suppressMessages(melt(tmp_metadata))
              tmp_metadata2 <- summarySE(tmp_metadata, groupvars = c("group", "variable"), measurevar = "value", na.rm = T)[,c(1,2,4)]
              tmp_metadata2 <- merge(tmp_metadata2[tmp_metadata2$group=="high",], tmp_metadata2[tmp_metadata2$group=="low",], by= "variable")
              tmp_metadata2 <- data.frame(variable = tmp_metadata2$variable, high_mean = tmp_metadata2$value.x, low_mean = tmp_metadata2$value.y, fc_high_low = tmp_metadata2$value.x/tmp_metadata2$value.y)
              
              # Calculate the pvalue
              tmp_metadata_pval <- compare_means(value~group, data = tmp_metadata, method = "t.test", paired = F, group.by = "variable", var.equal = FALSE, p.adjust.method = "none")
              
              # Merging with the sign value
              tmp_metadata <- merge(tmp_metadata2, tmp_metadata_pval[,c(1,5,8)], by = "variable")
              
              container[["summary"]][["metadata"]][[paste(i, j, k, sep = "_")]] <- tmp_metadata
              
            }
            
          }
          
        }
        
      }
      
      else {
        print(paste(gene, "is not present in", j, sep = " "))
      }
      
    }
    
  }
  
  class(container) <- "huva_experiment"
  
  return(container)
}
