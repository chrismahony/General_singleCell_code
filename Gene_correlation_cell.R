#function definition

matrix_to_expression_df<- function(x, obj){
        df<- x %>%
                as.matrix() %>% 
                as.data.frame() %>%
                tibble::rownames_to_column(var= "gene") %>%
                tidyr::pivot_longer(cols = -1, names_to = "cell", values_to = "expression") %>%
                tidyr::pivot_wider(names_from = "gene", values_from = expression) %>%
                left_join(obj@meta.data %>% 
                                  tibble::rownames_to_column(var = "cell"))
        return(df)
}


get_expression_data<- function(obj, assay = "RNA", slot = "data", 
                               genes = NULL, cells = NULL){
        if (is.null(genes) & !is.null(cells)){
                df<- GetAssayData(obj, assay = assay, slot = slot)[, cells, drop = FALSE] %>%
                        matrix_to_expression_df(obj = obj)
        } else if (!is.null(genes) & is.null(cells)){
                df <- GetAssayData(obj, assay = assay, slot = slot)[genes, , drop = FALSE] %>%
                        matrix_to_expression_df(obj = obj)
        } else if (is.null(genes & is.null(cells))){
                df <- GetAssayData(obj, assay = assay, slot = slot)[, , drop = FALSE] %>%
                        matrix_to_expression_df(obj = obj)
        } else {
                df<- GetAssayData(obj, assay = assay, slot = slot)[genes, cells, drop = FALSE] %>%
                        matrix_to_expression_df(obj = obj)
        }
        return(df)
}


#plotting

genes<- c("gene1", "gene2")

expression_data<- get_expression_data(obj, genes = genes)



ggplot(expression_data, aes(x= gene1, y = gene2)) + 
        #geom_smooth(method="lm") +
        geom_point(size = 0.8) +
        facet_wrap(~seurat_annotations) +
        ggpubr::stat_cor(method = "pearson")
