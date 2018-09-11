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
compareModels(preds_file, true_file, "preds_lasso.txt", "true_lasso.txt")




#Elanet
tune_grid = expand.grid(alpha=10^seq(0,-3, length=8), 
				lambda=10^seq(0,-3, length=8)) 
preds_file = "preds_elanet.txt"
true_file = "true_elanet.txt"
method = "glmnet"
sequence = seq(1:1000)
testing(method, folds, tune_grid, sequence, targ, regu, preds_file, true_file)
test_error(preds_file, true_file)
stats = read.table("stats.txt")
summary(stats)
compareModels(preds_file, true_file, "preds_lasso.txt", "true_lasso.txt")



testing_comb = function(method, folds, tune_grid, sequence, targ, regu, preds_file, true_file){

	preds_all_genes = c()
	true_all_genes = c()
	stats = c()

	for(index_targ in sequence){
		display(paste0("Target Gene: ", index_targ))

		CV <- lapply(1:5, function(x){
			control <- trainControl(method="cv", number=5, search="grid")
			model_en <- train(regu[folds != x,], targ[folds != x, index_targ], 
				method=method, tuneGrid=tune_grid, trControl=control)
			ind_best = rownames(model_en$bestTune)
			rmse_elanet = model_en$results[ind_best,3]


			rfe_control <- rfeControl(functions=rfFuncs, method="cv", 
				number=5, 
				returnResamp = "all", 
				saveDetails=TRUE)
			model_rf <- rfe(regu[folds != x,], targ[folds != x, index_targ],
				sizes=c(300),
				rfeControl=rfe_control)
			display(model_rf$results)
			#ind_best = rownames(model_rf$fit$bestTune)
			rmse_rf = min(model_rf$results[,2])

			use_elanet = rmse_elanet < rmse_rf
			display(use_elanet)

			if(use_elanet){model = model_en}
			else{model = model_rf}
					
			preds <- predict(model, regu[folds == x,])
			return(data.frame(preds, real = targ[folds == x, index_targ]))
		})
		listcv = do.call(rbind, CV)
		preds = listcv[[1]]
		true = listcv[[2]]
		preds_all_genes = cbind(preds_all_genes, preds)
		true_all_genes = cbind(true_all_genes, true)
	}

	write.table(preds_all_genes, preds_file)
	write.table(true_all_genes, true_file)
}
