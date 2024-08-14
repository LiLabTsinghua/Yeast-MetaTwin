library('classyfireR')

setwd('D:/LHH/YMDB/ymdb/data')

data <- read.csv("ymdb.csv", encoding="UTF-8")
data$INCHIKEY -> key

Classification_List <- purrr::map(key, get_classification)
Classification_List = Classification_List[!sapply(Classification_List,is.null)]
l <- lengths(list(lengths(Classification_List)))
list1 <- list()
for ( i in 1:l) {
  list1[i] = Classification_List[[i]]@classification$Classification[3]
}
list1 = list1[!sapply(list1,is.null)]
write.csv(list1,"YMDB_class.csv", row.names = FALSE)
list1 <- list()
for ( i in 1:l) {
  list1[i] = Classification_List[[i]]@classification$Classification[2]
}
list1 = list1[!sapply(list1,is.null)]
write.csv(list1,"YMDB_super_class.csv", row.names = FALSE)
list1 <- list()
for ( i in 1:l) {
  list1[i] = Classification_List[[i]]@meta$inchikey
}
list1 = list1[!sapply(list1,is.null)]
write.csv(list1,"YMDB_inchikey.csv", row.names = FALSE)

