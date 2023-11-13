plot_cluster_heatmap_mean <- function(matrix_groupmean,
                                      stats,
                                      col_order = NULL, 
                                      row_order = NULL, 
                                      cluster_columns = T,
                                      cluster_rows = T,
                                      k = 0, 
                                      return_HM = T, 
                                      file_name = "module_heatmap_mean.pdf"){
  
  # filter for included clusters (non-white)
  c_df <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], cluster_included == "yes")
  mat_heatmap <- NULL
  
  # merge result of t-test 
  c_df <- left_join(c_df, stats[,c("color", "p.adj")])
  
  
  if(!base::is.null(row_order)){
    for (c in row_order){
      #get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      # Mean expression per module
      c_GFCs <- dplyr::filter(matrix_groupmean, Gene %in% genes)
      c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))], 2, base::mean)
      
      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
      
    }
    base::rownames(mat_heatmap) <- row_order
    
  }else{
    for (c in base::unique(c_df$color)){
      #get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      
      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(matrix_groupmean, Gene %in% genes)
      
      if(base::is.vector(c_GFCs)){
        c_GFC_means <- cGFCs
      }else{
        c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))] %>% 
                                     base::as.data.frame(), 2, base::mean)
      }
      
      
      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
      
    }
    base::rownames(mat_heatmap) <- c_df$color
  }
  
  
  base::colnames(mat_heatmap) <- base::colnames(matrix_groupmean)[1:(base::ncol(matrix_groupmean)-1)]
  
  if(!base::is.null(col_order)){
    mat_heatmap <- mat_heatmap %>% base::as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order] %>% base::as.matrix()
  }
  
  
  if(!base::is.null(row_order)){
    cluster_colors <- base::factor(row_order)
    base::names(cluster_colors) <- row_order
    c_df <- c_df[base::match(row_order, c_df$color),]
  }else{
    cluster_colors <- base::factor(c_df$color)
    base::names(cluster_colors) <- c_df$color
    row_order <- base::unique(c_df$color)
  }
  
  
  # row annotation
  col_fun <- colorRamp2(breaks = c(0, 0.01, 0.05, 1), 
                        colors = c("#3F007D", "#807DBA", "#BCBDDC", "#FCFBFD"))
  ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                                simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                          p.adj = ComplexHeatmap::anno_simple(c_df$p.adj, col = col_fun,
                                                                               width = grid::unit(1.5, "cm"), 
                                                                               gp = grid::gpar(col = "black"), 
                                                                               pch = format(c_df$p.adj, scientific = FALSE, trim = TRUE), 
                                                                               pt_size = unit(1, "snpc")*0.3,
                                                                               pt_gp = gpar(col = "white")),
                                          genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),
                                          # gene_nums = ComplexHeatmap::anno_text(c_df$gene_no, width = grid::unit(1.5, "cm"), gp = grid::gpar(fontsize = 10)),
                                          which = "row",
                                          width = grid::unit(4.5, "cm"),
                                          annotation_name_side = "top",
                                          gap = grid::unit(2, "mm"),
                                          annotation_name_rot = 0,
                                          annotation_name_gp = grid::gpar(fontsize = 8))
  
  lgd1 <- Legend(col_fun = col_fun, title = "p.adj", at = c(0, 0.01, 0.05, 1), labels = format(c(0, 0.01, 0.05, 1))) 
  
  # col annotation
  ca <- ComplexHeatmap::HeatmapAnnotation(group = ComplexHeatmap::anno_simple(gsub('[0-9]+', '', col_order), col = col_group,                                          
                                                                              simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                          which = "column",
                                          height = grid::unit(.5, "cm"),
                                          annotation_name_side = "left",
                                          gap = grid::unit(2, "mm"),
                                          annotation_name_rot = 0,
                                          annotation_name_gp = grid::gpar(fontsize = 8))
  
  lgd2 <- Legend(labels = c("NEO", "AD"), title = "group", legend_gp = gpar(fill = col_group))
  
  lgd_list <- packLegend(lgd1, lgd2)
  
  all_conditions <- NULL
  
  for(setnum in 1:base::length(hcobject[["layers"]])){
    all_conditions <- base::c(all_conditions, base::as.character(dplyr::pull(hcobject[["data"]][[base::paste0("set", setnum, "_anno")]], hcobject[["global_settings"]][["voi"]])))
  }
  all_conditions <- base::table(all_conditions) %>%
    base::as.data.frame() %>%
    dplyr::filter(., all_conditions %in% base::colnames(mat_heatmap))
  all_conditions <- all_conditions[base::match(base::colnames(mat_heatmap), base::as.character(all_conditions$all_conditions)),]
  all_conditions <- base::paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")
  
  
  Cairo::Cairo(file = paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", file_name),
               width = 50,
               height = 30,
               pointsize=11,
               dpi=300,
               type = "pdf",
               units = "in")
  
  hm <- ComplexHeatmap::Heatmap(mat_heatmap,
                                right_annotation = ha,
                                top_annotation = ca,
                                col = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(51),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                width = ncol(mat_heatmap)*unit(10, "mm"), 
                                height = nrow(mat_heatmap)*unit(10, "mm"),
                                column_names_rot = 90,
                                row_names_gp = grid::gpar(fontsize = 8),
                                column_names_gp = grid::gpar(fontsize = 8),
                                rect_gp = grid::gpar(col = "black"),
                                heatmap_legend_param = list(title = "scaled mean expr.", legend_height = grid::unit(3, "cm")), column_km = k)
  
  
  hm_w_lgd <- ComplexHeatmap::draw(hm, 
                                   annotation_legend_list = lgd_list, 
                                   merge_legends = T,
                                   padding = grid::unit(c(2, 2, 2, 30), "mm"))
  
  grDevices::dev.off()
  
  print(hm_w_lgd)
  if(return_HM){
    hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]] <<- hm_w_lgd
  }
  
}
