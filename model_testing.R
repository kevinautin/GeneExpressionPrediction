setwd("C:/Users/kevinautin/Google Drive/Thesis/Code/GeneExpressionPrediction")
#setwd("C:/Users/user/Thesis/")
Sys.setenv(LANG = "en")
source("my_func.R")
library(caret)

targ = read.table("targ.txt")
regu = read.table("regu.txt")

set.seed(5)
folds = sample(rep(1:5, length.out = nrow(regu)), size = nrow(regu), replace = F)



#Lasso
tune_grid = expand.grid(alpha=1,
				lambda=10^seq(0,-3, length=5))
preds_file = "preds_lasso.txt"
true_file = "true_lasso.txt"
method = "glmnet"
sequence = seq(1, 1000)
testing(method, folds, tune_grid, sequence,  targ, regu, preds_file, true_file)
test_error(preds_file, true_file)
stats = read.table("stats.txt")
summary(stats)


#Elastic net
tune_grid = expand.grid(alpha=10^seq(0,-3, length=5),
				lambda=10^seq(0,-3, length=5))
preds_file = "preds_elanet.txt"
true_file = "true_elanet.txt"
method = "glmnet"
sequence = seq(1, 1000)
testing(method, folds, tune_grid, sequence,  targ, regu, preds_file, true_file)
test_error(preds_file, true_file)
stats = read.table("stats.txt")
summary(stats)
compareModels(preds_file, true_file, "preds_lasso.txt", "true_lasso.txt")
