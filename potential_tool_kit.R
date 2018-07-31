# an initial look at some potential tools to find weeds in fields
# general ploting data manipulation
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
# classifiers
library(extraTrees)
library(gamsel)
library(randomForest)
library(gbm)

setwd('/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/')

################################################################################
# make a simple data generator that makes the categorical data, and categories change
# as simple linear function 10 predictors, option to stratify so some fileds will
# be all high and others all low
make_cat_data  = function(n_samp = 100, n_field = 1, noise = 2.5, strat = FALSE){

	# make the dummy X variables
	X = matrix(runif(n_samp * n_field * 10, -2, 2), ncol = 10)
	colnames(X) = paste0('x', 1:10)

	# estiomates on linear predictors
	B = c(-2, 0.5, 0.2, 0.2, 0.6, -1.7, 0.2)
	# function to priduce the liner predictor
	f = function(x){
		return(x[1] * B[1] + x[2] * B[2] + x[3] * B[3] + x[4] * B[4] +
			x[3] * x[4] * B[5] + x[5] * B[6] + (x[5] ^ 2) * B[7])
	}

	# make the linear predictors
	lp = apply(X, MARGIN = 1, FUN = f)
	Xlp = as.data.frame(cbind(X, lp))

	# assinge field IDs
	if(strat){

		Xlp = Xlp[order(Xlp$lp), ]
		Xlp$field = rep(paste0('f', 1:n_field), each = n_samp)

	}else{

		Xlp$field = rep(paste0('f', 1:n_field), each = n_samp)

	}

	#probit cum norm to get probs from linear to categorical
	# set cut points
	C = quantile(Xlp$lp, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))

	# apply to cumulative norm
	Xlp = mutate(Xlp, p0 = pnorm(C[2], mean = lp, sd = noise),
		pl = pnorm(C[3], mean = lp, sd = noise) - p0,
		pm = pnorm(C[4], mean = lp, sd = noise) - p0 - pl,
		ph = pnorm(C[5], mean = lp, sd = noise) - p0 - pl - pm,
		pv = pnorm(C[6], mean = lp, sd = noise) - p0 - pl - pm - ph)

	# make the random draws to assing the categories using multinomial
	Xlp$ob_den_state = apply(select(Xlp, p0:pv), MARGIN = 1, FUN = function(x){
		ob = rmultinom(n = 1, size = 1, prob = x)
		DS = c('0','L','M','H','V')
		return(DS[ob == 1])
	})

	return(Xlp)

}

################################################################################
# make a confustion matrix dataframe,
# tol controls where the cutoff for labeling predictions Uncertian
make_CM_df = function(obs, pred, group_ID, tol = 0.2){

	pred_cat = colnames(pred)

	# get the classification thresholds
	up_t = (1 / length(pred_cat)) * (1 + tol)

	# find the most proabal class
	most_prob = apply(pred, MARGIN = 1, FUN = function(x){
		max_prob = max(x)
		pred_class = ifelse(max_prob < up_t, 'Un', pred_cat[which(x == max_prob)])
		return(pred_class)
	})

	# put it all together into a dataframe
	df = data.frame(group_ID = group_ID, obs_class = obs, pred_class = most_prob) %>%
		group_by(group_ID, obs_class, pred_class) %>%
		summarise(counts = n())

	return(df)

}

################################################################################
# generate a dataframe to play with, start with a random one, no stratification
rand_df = make_cat_data(n_samp = 200, n_field = 10, noise = 2.5, strat = FALSE)

# check I can recover the underlying model of the linear predictor
m1 = lm(lp ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x3:x4 + I(x5^2), data = rand_df)
summary(m1)
# as expected recovers everything perfectly as expected given no noise term on the
# lp, interestingly even with no noise at all three predcitors that have no effect
# have p-values between 0.05 and 0.1

# plot to see the relationship between the lp and the observed category for each
# field
ggplot(rand_df, aes(x = ob_den_state, y = lp)) + geom_boxplot() +
 facet_wrap(~field, ncol = 4)

ggplot(rand_df, aes(x = ob_den_state)) + geom_bar() +
 facet_wrap(~field, ncol = 4)

# create a stratified dataset
strat_df = make_cat_data(n_samp = 200, n_field = 10, noise = 2.5, strat = TRUE)

ggplot(strat_df, aes(x = ob_den_state, y = lp)) + geom_boxplot() +
 facet_wrap(~field, ncol = 4)

ggplot(strat_df, aes(x = ob_den_state)) + geom_bar() +
 facet_wrap(~field, ncol = 4)

# split the data into test sets and training sets take the classic
# 80:20 split, but split on fields as I want to predict to new
# fields
train_rand = filter(rand_df, !field %in% c('f1', 'f4'))
test_rand = filter(rand_df, field %in% c('f1', 'f4'))

train_strat = filter(strat_df, !field %in% c('f1', 'f4'))
test_strat = filter(strat_df, field %in% c('f1', 'f4'))

# make some binary versions as a wider range of tools exists for this class of data
train_rand_bin = mutate(train_rand,
		ob_bin = ifelse(ob_den_state %in% c('0', 'L'), 'L', 'H'))
test_rand_bin = mutate(test_rand,
		ob_bin = ifelse(ob_den_state %in% c('0', 'L'), 'L', 'H'))

train_strat_bin = mutate(train_strat,
		ob_bin = ifelse(ob_den_state %in% c('0', 'L'), 'L', 'H'))
test_strat_bin = mutate(test_strat,
		ob_bin = ifelse(ob_den_state %in% c('0', 'L'), 'L', 'H'))

###############################################################################
# Look at classifying using extraTrees, with posible multi-task learing
# exTree wants the input and output data as matricies and vectors so make these
exT_MT_x = as.matrix(select(train_rand_bin, x1:x10))
exT_MT_y = as.factor(train_rand_bin$ob_bin)
exT_MT_task = as.numeric(sapply(strsplit(train_rand_bin$field, split = 'f'),
	FUN = function(x) x[2]))

ET_MT_rand = extraTrees(x = exT_MT_x, y = exT_MT_y, ntree = 5000, nodesize = 5,
		numRandomCuts = 4, evenCuts = TRUE, numThreads = 2, tasks = exT_MT_task,
		probOfTaskCuts = 0.25)

test_exT_MT_x = as.matrix(select(test_rand_bin, x1:x10))
test_exT_MT_task = as.numeric(sapply(strsplit(test_rand_bin$field, split = 'f'),
	FUN = function(x) x[2]))

ET_MT_rand_OoS = predict(ET_MT_rand, newdata = test_exT_MT_x, probability = TRUE,
	newtasks = test_exT_MT_task)

ET_MT_conmat = make_CM_df(obs = test_rand_bin$ob_bin, pred = ET_MT_rand_OoS,
		group_ID = test_rand_bin$field, tol = 0.1) %>%
	mutate(classifier = 'exTree_MT_rand')

#plot the confustion matrix
ET_MT_perf = ggplot(ET_MT_conmat, aes(x = obs_class, y = pred_class)) +
	geom_tile(aes(fill = counts)) +
	geom_text(aes(label = counts)) +
	facet_grid(classifier ~ group_ID)

ET_MT_perf

# no multi-task learning
ET_rand = extraTrees(x = exT_MT_x, y = exT_MT_y, ntree = 5000, nodesize = 5,
		numRandomCuts = 4, evenCuts = TRUE, numThreads = 2)

ET_rand_OoS = predict(ET_rand, newdata = test_exT_MT_x, probability = TRUE)

ET_conmat = make_CM_df(obs = test_rand_bin$ob_bin, pred = ET_rand_OoS,
		group_ID = test_rand_bin$field, tol = 0.1) %>%
	mutate(classifier = 'exTree_rand')

ET_perf = ggplot(ET_conmat, aes(x = obs_class, y = pred_class)) +
	geom_tile(aes(fill = counts)) +
	geom_text(aes(label = counts)) +
	facet_grid(classifier ~ group_ID)

ET_perf

################################################################################
# Try a Lasso GAM
gam_x = as.matrix(select(train_rand_bin, x1:x10))
gam_y = ifelse(train_rand_bin$ob_bin == 'H', 1, 0)
GAM_rand = cv.gamsel(x = gam_x, y = gam_y, lambda = 26 * exp((-0.15 * 0:50)),
	family = 'binomial', degrees = rep(5, dim(gam_x)[2]), dfs = rep(5, dim(gam_x)[2]),
	type.measure =  "class", nfolds = 5, parallel = FALSE)

# find the best lambda
op_lam = GAM_rand$lambda[GAM_rand$cvm == min(GAM_rand$cvm)][1]

GAM_rand_fit = gamsel(x = gam_x, y = gam_y, lambda = op_lam, family = 'binomial',
	degrees = rep(10, dim(gam_x)[2]), dfs = rep(5, dim(gam_x)[2]))

test_GAM = as.matrix(select(test_rand_bin, x1:x10))
pred_OoS = predict(GAM_rand_fit, newdata = test_GAM, type = 'response')
pred_OoS = cbind(pred_OoS, 1 - pred_OoS)
colnames(pred_OoS) = c('H', 'L')

GAM_conmat = make_CM_df(obs = test_rand_bin$ob_bin, pred = pred_OoS,
		group_ID = test_rand_bin$field, tol = 0.1) %>%
	mutate(classifier = 'lasso_GAM')

################################################################################
# Random forest
RF_y = as.factor(train_rand_bin$ob_bin)
RF_x = as.matrix(select(train_rand_bin, x1:x10))

RF_rand = randomForest(x = RF_x, y = RF_y, ntree = 500, replace = TRUE,
	strata = RF_y, sampsize = c(100, 100), nodesize = 10,
	data = train_rand_bin)

RF_test_x = as.matrix(select(test_rand_bin, x1:x10))
pred_OoS = predict(RF_rand, newdata = RF_test_x, type = 'prob')

RF_conmat = make_CM_df(obs = test_rand_bin$ob_bin, pred = pred_OoS,
		group_ID = test_rand_bin$field, tol = 0.1) %>%
	mutate(classifier = 'rand_forest')

################################################################################
# generlaized boosted regression tree
BRT_y = ifelse(train_rand_bin$ob_bin == 'H', 1, 0)
BRT_x = as.matrix(select(train_rand_bin, x1:x10))

BRT_rand_cv = gbm(BRT_y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
	n.trees = 10000, interaction.depth = 5, shrinkage = 0.005, cv.folds = 5,
	data = train_rand_bin)

op_tree = gbm.perf(BRT_rand_cv)
BRT_rand = gbm.fit(x = BRT_x, y = BRT_y, n.trees = op_tree, shrinkage = 0.005)

BRT_test_x = as.matrix(select(test_rand_bin, x1:x10))
pred_OoS = predict(BRT_rand, newdata = BRT_test_x, n.trees = op_tree,
	type = 'response')

pred_OoS = cbind(pred_OoS, 1 - pred_OoS)
colnames(pred_OoS) = c('H', 'L')

BRT_conmat = make_CM_df(obs = test_rand_bin$ob_bin, pred = pred_OoS,
		group_ID = test_rand_bin$field, tol = 0.1) %>%
	mutate(classifier = 'BRT')

################################################################################
# Over all comparison plots
df_all_bin = bind_rows(ET_MT_conmat, ET_conmat, GAM_conmat, RF_conmat, BRT_conmat)
plt_all_bin = ggplot(df_all_bin, aes(x = obs_class, y = pred_class)) +
	geom_tile(aes(fill = counts)) +
	geom_text(aes(label = counts)) +
	scale_fill_viridis() +
	facet_grid(classifier ~ group_ID)

pdf(file = 'sim_data_rand_confustion_mat.pdf', width = 8, height = 12)
	plt_all_bin
dev.off()

################################################################################
# Now try it with stratified data
################################################################################
exT_MT_x = as.matrix(select(train_strat_bin, x1:x10))
exT_MT_y = as.factor(train_strat_bin$ob_bin)
exT_MT_task = as.numeric(sapply(strsplit(train_strat_bin$field, split = 'f'),
	FUN = function(x) x[2]))

ET_MT_strat = extraTrees(x = exT_MT_x, y = exT_MT_y, ntree = 5000, nodesize = 5,
		numRandomCuts = 4, evenCuts = TRUE, numThreads = 2, tasks = exT_MT_task,
		probOfTaskCuts = 0.25)

test_exT_MT_x = as.matrix(select(test_strat_bin, x1:x10))
test_exT_MT_task = as.numeric(sapply(strsplit(test_strat_bin$field, split = 'f'),
	FUN = function(x) x[2]))

ET_MT_strat_OoS = predict(ET_MT_strat, newdata = test_exT_MT_x, probability = TRUE,
	newtasks = test_exT_MT_task)

ET_MT_conmat = make_CM_df(obs = test_strat_bin$ob_bin, pred = ET_MT_strat_OoS,
		group_ID = test_strat_bin$field, tol = 0.1) %>%
	mutate(classifier = 'exTree_MT_strat')
################################################################################
exT_x = as.matrix(select(train_strat_bin, x1:x10))
exT_y = as.factor(train_strat_bin$ob_bin)
ET_strat = extraTrees(x = exT_x, y = exT_y, ntree = 5000, nodesize = 5,
		numRandomCuts = 4, evenCuts = TRUE, numThreads = 2)

test_exT_x = as.matrix(select(test_strat_bin, x1:x10))

ET_strat_OoS = predict(ET_strat, newdata = test_exT_x, probability = TRUE)

ET_conmat = make_CM_df(obs = test_strat_bin$ob_bin, pred = ET_strat_OoS,
		group_ID = test_strat_bin$field, tol = 0.1) %>%
	mutate(classifier = 'exTree_strat')

################################################################################
gam_x = as.matrix(select(train_strat_bin, x1:x10))
gam_y = ifelse(train_strat_bin$ob_bin == 'H', 1, 0)
GAM_strat = cv.gamsel(x = gam_x, y = gam_y, lambda = 26 * exp((-0.15 * 0:50)),
	family = 'binomial', degrees = rep(5, dim(gam_x)[2]), dfs = rep(5, dim(gam_x)[2]),
	type.measure =  "class", nfolds = 5, parallel = FALSE)

# find the best lambda
op_lam = GAM_strat$lambda[GAM_strat$cvm == min(GAM_strat$cvm)][1]

GAM_strat_fit = gamsel(x = gam_x, y = gam_y, lambda = op_lam, family = 'binomial',
	degrees = rep(10, dim(gam_x)[2]), dfs = rep(5, dim(gam_x)[2]))

test_GAM = as.matrix(select(test_strat_bin, x1:x10))
pred_OoS = predict(GAM_strat_fit, newdata = test_GAM, type = 'response')
pred_OoS = cbind(pred_OoS, 1 - pred_OoS)
colnames(pred_OoS) = c('H', 'L')

GAM_conmat = make_CM_df(obs = test_strat_bin$ob_bin, pred = pred_OoS,
		group_ID = test_strat_bin$field, tol = 0.1) %>%
	mutate(classifier = 'lasso_GAM')

################################################################################
RF_y = as.factor(train_strat_bin$ob_bin)
RF_x = as.matrix(select(train_strat_bin, x1:x10))

RF_strat = randomForest(x = RF_x, y = RF_y, ntree = 500, replace = TRUE,
	strata = RF_y, sampsize = c(50, 50), nodesize = 20,
	data = train_strat_bin)

RF_test_x = as.matrix(select(test_strat_bin, x1:x10))
pred_OoS = predict(RF_strat, newdata = RF_test_x, type = 'prob')

RF_conmat = make_CM_df(obs = test_strat_bin$ob_bin, pred = pred_OoS,
		group_ID = test_strat_bin$field, tol = 0.1) %>%
	mutate(classifier = 'rand_forest')

################################################################################
BRT_y = ifelse(train_strat_bin$ob_bin == 'H', 1, 0)
BRT_x = as.matrix(select(train_strat_bin, x1:x10))

BRT_strat_cv = gbm(BRT_y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
	n.trees = 10000, interaction.depth = 5, shrinkage = 0.005, cv.folds = 5,
	data = train_strat_bin)

op_tree = gbm.perf(BRT_strat_cv)
BRT_strat = gbm.fit(x = BRT_x, y = BRT_y, n.trees = op_tree, shrinkage = 0.005)

BRT_test_x = as.matrix(select(test_strat_bin, x1:x10))
pred_OoS = predict(BRT_strat, newdata = BRT_test_x, n.trees = op_tree,
	type = 'response')

pred_OoS = cbind(pred_OoS, 1 - pred_OoS)
colnames(pred_OoS) = c('H', 'L')

BRT_conmat = make_CM_df(obs = test_strat_bin$ob_bin, pred = pred_OoS,
		group_ID = test_strat_bin$field, tol = 0.1) %>%
	mutate(classifier = 'BRT')


################################################################################
# overall comparision
df_all_bin = bind_rows(ET_conmat, GAM_conmat, RF_conmat, BRT_conmat, ET_MT_conmat)
plt_all_bin = ggplot(df_all_bin, aes(x = obs_class, y = pred_class)) +
	geom_tile(aes(fill = counts)) +
	geom_text(aes(label = counts)) +
	scale_fill_viridis() +
	facet_grid(classifier ~ group_ID)

pdf(file = 'sim_data_stritified_confustion_mat.pdf', width = 8, height = 12)
	plt_all_bin
dev.off()

# As expected none of the methods do that well on stratified data lasso_GAM
# works the best, but is still not that good, random forest is the worst by far.
# Could try differnt cut offs for classification to try and address the unbalanced data
