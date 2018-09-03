setwd("C:/Users/kevinautin/Google Drive/Thesis/")
Sys.setenv(LANG = "en")


#data.txt file
data = read.table("data.txt", header=TRUE)
colnames(data)[1] = "GeneID"
data = data[, -2]

#regulators file
data_regu = read.csv("TF_list.csv", header = TRUE, sep=",")
data_regu = as.data.frame(data_regu[, 1])
colnames(data_regu)[1] = "GeneID"

#Isolate the regulators and the targets from the data
regu = merge(data, data_regu, by="GeneID")
library(dplyr)
targ_all = anti_join(data, data_regu, by="GeneID")

# 1000 targets with highest variance
targ_all_t = t(targ_all)
colnames(targ_all_t) <- targ_all_t[1,]
targ_all_t <- targ_all_t[-1,] 
variances = diag(var(targ_all_t))
var_1000 = sort(variances, decreasing = TRUE)[1:1000]
var_1000 = as.data.frame(names(var_1000))
colnames(var_1000)[1] = "GeneID"
targ_1000 = merge(targ_all, var_1000, by="GeneID")


#reshape data.frame
regu = t(regu)
colnames(regu) <- paste0("X", regu[1,])
regu <- regu[-1,] 
targ_1000 = t(targ_1000)
colnames(targ_1000) = paste0("Y", targ_1000[1,])
targ_1000 = targ_1000[-1,]

#rescale data
regu = scale(regu)
targ = scale(targ_1000)

write.table(targ, "targ.txt")
write.table(regu, "regu.txt")