rmse <- function(actual, predicted){
	sqrt(mean((actual - predicted)^2))
}

meanRmse <- function(preds, true){
	total_rmse = 0
	for(i in 1:ncol(preds)){
		total_rmse = total_rmse + rmse(preds[,i], true[,i])
	}
	mean_rmse = total_rmse/ncol(preds)
	return(mean_rmse)
}

display = function(text){
	print(text)
	Sys.sleep(0.0001)
	flush.console()
}

test_error = function(preds_file, true_file){
	preds = read.table(preds_file)
	true = read.table(true_file)
	meanRmse(preds, true)
}

testing = function(method, folds, tune_grid, sequence, targ, regu, preds_file, true_file){

	preds_all_genes = c()
	true_all_genes = c()
	stats = c()

	for(index_targ in sequence){
		display(paste0("Target Gene: ", index_targ))

		CV <- lapply(1:5, function(x){
			control <- trainControl(method="cv", number=5, search="grid")
			model <- train(regu[folds != x,], targ[folds != x, index_targ], 
				method=method, tuneGrid=tune_grid, trControl=control)
			display(model)
			stats <<- rbind(stats, model$bestTune)		
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
	write.table(stats, "stats.txt")
}

compareModels = function(preds_file_1, true_file_1, preds_file_2, true_file_2){
	preds_1 = read.table(preds_file_1)
	true_1 = read.table(true_file_1)
	preds_2 = read.table(preds_file_2)
	true_2 = read.table(true_file_2)
	print(dim(preds_1))
	rmse_1 = c()
	rmse_2 = c()
	for(i in 1:ncol(preds_1)){
		rmse_1 = cbind(rmse_1, rmse(preds_1[,i], true_1[,i]))
		rmse_2 = cbind(rmse_2, rmse(preds_2[,i], true_2[,i]))
	}
	display(t.test(rmse_1, rmse_2, paired=TRUE))
	plot(rmse_2~rmse_1, main="RMSE", xlab = preds_file_1, ylab = preds_file_2, xlim=c(0.2, 0.8), ylim=c(0.2, 0.8))
	lines(c(0,1),c(0,1))
	print(paste("Model 1 better than model 2:", mean(rmse_1 < rmse_2)))
	print(paste("Bad predictions model l:", mean(rmse_1 > 0.4)))
	print(paste("Bad predictions model 2:", mean(rmse_2 > 0.4)))
	print(paste("SD model l:", sd(rmse_1)))
	print(paste("SD model 2:", sd(rmse_2)))
}