data_nfkb<-readRDS(file="data-raw/data_nfkb.rds")

devtools::use_data(data_nfkb,overwrite=TRUE)