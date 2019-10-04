data_nfkb<-readRDS(file="data-raw/data_nfkb.rds")

data_nfkb<-data_nfkb[,c("signal","response_0","response_3",
                        "response_21","response_90","response_120")]
use_data(data_nfkb,overwrite=TRUE)
