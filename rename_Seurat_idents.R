obj$new_ident <- obj@meta.data[["integrated_snn_res.0.06"]]
Idents(obj) <- 'new_ident'
Levels(obj)
current.sample.ids <- c("Linning", "0","1","2", "3", "4", "5","6","7")
new.sample.ids <- c("Linning", "0","1","2", "3", "4", "5","0","0")

obj@meta.data[["new_ident"]] <- plyr::mapvalues(x = obj@meta.data[["new_ident"]], from = current.sample.ids, to = new.sample.ids)
