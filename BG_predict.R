# Goal is to try and classify how much black grass is in a given cell, given
# some hyper-spectral imagery, and see how well that prediction works in and
# out of training date.

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(viridis)
# classifiers
library(extraTrees)
library(gamsel)
library(randomForest)
library(gbm)

veg_ind_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/Data/sample_data/UoS/hyperspectral/vegetationindices/'
den_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/Data/dens_data/den_shape_2018/'
data_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/Data/'
plot_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/outputs/plots'

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

# Turn a Factor vector, FV, into a matrix for categorical observation, use factor
# to maintaine desired order, which is needed for RPSS
vect2cat_mat = function(FV){

	fact_levels = levels(FV)

	m = t(sapply(FV, FUN = function(x){

			cat_rep = rep(0, length(fact_levels))
			cat_rep[fact_levels == x] = 1
			return(cat_rep)
		}))

	colnames(m) = fact_levels

	return(m)

}

# Calculate the Ranked Probability Skill Score (absolute), as the absolute has no biase,
# See Muller et al. 2005, Journal of Climate 18:1513-1523 and
# Weigel et al. 2007, Monthly Weather Review 135:118-124
# RPSS compairs a given forcast to a reference forcast, in this case we assume
# random classification as the reference (i.e. p_i = 1/K for calss i)
# pred = N x K matrix (K is number of calsses)
# obs = N x K matrix, elements are [0, 1], each row sums to 1
RPSS = function(pred, obs){

	# turn the pred and obs matricies into cumulative prob
	pred_cum = t(apply(t(pred), 2, cumsum))
	obs_cum = t(apply(t(obs), 2, cumsum))

	# expected Ranked Prob Skill from prediction
	E_RPS = mean(apply(abs(pred_cum - obs_cum), 1, sum))

	# make the reference, if the classifery is just picking random classes the
	# results will average to be the propotion in each class
	rand_guess = cumsum(colSums(obs) / dim(obs)[1])
	E_refS = mean(apply(t(abs(rand_guess - t(obs_cum))), 1, sum))

	return(1 - (E_RPS / E_refS))

}
# bit of random test data for the RPSS
test_obs = sample(c('l', 'm', 'h', 'v'), size = 1000, replace = TRUE)
test_df = data.frame(test_obs = test_obs,
		X1 = sample.int(n = 10, size = length(test_obs), replace = TRUE),
		X2 = sample.int(n = 10, size = length(test_obs), replace = TRUE),
		X3 = sample.int(n = 10, size = length(test_obs), replace = TRUE))

RF_y = factor(test_df$test_obs, levels = c('l', 'm', 'h', 'v')) # order will matter for RPSS
RF_x = as.matrix(select(test_df, X1:X3))
test_RF = randomForest(x = RF_x, y = RF_y, ntree = 500, replace = TRUE,
	strata = RF_y, sampsize = 20, nodesize = 10)
RF_pred_x = data.frame(X1 = sample.int(n = 10, size = length(test_obs), replace = TRUE),
		X2 = sample.int(n = 10, size = length(test_obs), replace = TRUE),
		X3 = sample.int(n = 10, size = length(test_obs), replace = TRUE))
test_pred = predict(test_RF, newdata = as.matrix(RF_pred_x), type = 'prob')
# test
test_RPSS = RPSS(pred = test_pred, obs = vect2cat_mat(RF_y)) # note use of RF_y, need ordered factor
# Think about how I might randomize this to get an approiate null distribution.
# Assume the null model the predictions are random with respect to the observations
# draw the predictions from  a multinomial with prob equal to fraction of each
# class in observations
RPSS_null = function(obs, n_sim){

	# realtive frequency of each class in the data
	rel_freq = colSums(obs) / dim(obs)[1]
	rpss_null = numeric(n_sim)

	# simulate random guessing by the classifier
	for(i in 1:n_sim){

		rand_pred = t(rmultinom(n = dim(obs)[1], size = 1, prob = rel_freq))
		rpss_null[i] = RPSS(pred = rand_pred, obs = obs)

		if(i %% 100 == 0) print(paste0('sim: ', i))
	}

	return(rpss_null)

}

RPSS_null_dist = RPSS_null(obs = vect2cat_mat(RF_y), n_sim = 500)
hist(RPSS_null_dist)

################################################################################
# get the data
setwd(data_loc)
# Start with the PCA data, already proccesed, so each grid square is summaries
# by the frist 4 moments of the distribution of each PCA dimention in each cell
PCA_dat = na.omit(read.csv('hyper_spec_dist_sum_PCA.csv', header = TRUE, stringsAsFactor = FALSE))

# break the data in to training and test sets
# take one of the strips and see if the models can predict to new strips,
# both those grid cells in the training set (but in differnt strips), and
# compleatly new ones
PCA_st5 = filter(PCA_dat, strip == 'St5')
PCA_dat = filter(PCA_dat, strip != 'St5')

# sample the remaining data to make a test set
PCA_test = sample_n(PCA_dat, size = 17)
PCA_train = filter(PCA_dat, !(uID %in% PCA_test$uID))

################################################################################
# fit a random forest model, although, not really enough data to make it work.
RF_y = as.factor(PCA_train$den_st)
RF_x = as.matrix(select(PCA_train, PC1_mean:PC7_kur))
RF_x_test = as.matrix(select(PCA_test, PC1_mean:PC7_kur))
RF_y_test = vect2cat_mat(as.factor(PCA_test$den_st))[, c('l', 'm', 'h')] # note reorder

# Tune the classifyer a bit, not much data so test to a small amount of out of sample data
tune_table = data.frame(n_trees = rep(c(250, 500, 1000, 2000), each = 4),
	samp_size = rep(c(10, 20, 50, 100), times = 4),
	OoS_skill_test = NA)

for(n in seq_along(tune_table$n_trees)){

	RF = randomForest(x = RF_x, y = RF_y, ntree = tune_table$n_trees[n], replace = TRUE,
		strata = RF_y, sampsize = tune_table$samp_size[n], nodesize = 10)

		pred_OoS = predict(RF, newdata = RF_x_test, type = 'prob')
		tune_table$OoS_skill_test[n] =  RPSS(pred = pred_OoS[, c('l', 'm', 'h')],
			obs = RF_y_test)

		print(paste0('ntrees = ', tune_table$n_trees[n], ' | sampsize = ', tune_table$samp_size[n]))
}

# looks like sample size is the most important for predictive skill, even 250 trees does fine
# test on an unseen strip (same field, different flight)
RF_tuned = randomForest(x = RF_x, y = RF_y, ntree = 250, replace = TRUE,
	strata = RF_y, sampsize = 100, nodesize = 10)

RF_st5_x = as.matrix(select(PCA_st5, PC1_mean:PC7_kur))
pred_OoStrip = predict(RF_tuned, newdata = RF_st5_x, type = 'prob')

RF_conmat = make_CM_df(obs = PCA_st5$den_st, pred = pred_OoStrip,
		group_ID = PCA_st5$strip, tol = 0.1) %>%
	ungroup() %>%
	mutate(obs_class = factor(obs_class, levels = c('l', 'm', 'h', 'Un')),
		pred_class = factor(pred_class, levels = c('l', 'm', 'h', 'Un')))


plt_con_mat = ggplot(RF_conmat, aes(x = obs_class, y = pred_class)) +
	geom_tile(aes(fill = counts)) +
	geom_text(aes(label = counts)) +
	scale_fill_viridis()

OoS_skill = RPSS(pred = pred_OoStrip[, c('l', 'm', 'h')],
	obs = vect2cat_mat(as.factor(PCA_st5$den_st))[, c('l', 'm', 'h')]) # reorder the coloumns to make cumsum work

OoS_skill_null = RPSS_null(obs = vect2cat_mat(as.factor(PCA_st5$den_st))[, c('l', 'm', 'h')], n_sim = 1000)
hist(OoS_skill_null)
lines(x = c(OoS_skill, OoS_skill), y = c(0, 300), col = 'red')
# Even with this very modest data this works better than I thought it would, need
# The skill to out of sample prediction probably better than chance, but there is
# a marginal prob our apprent skill is just chance due to small test set.
# Need the full data set to really test it, this will have 2 effects, it will expand the
# test set shrinking the confidence limits on the skill scores, second it will
# give the classifer more data to train on and hopefully improve its performance.

# map these predictions to see how the errors realte to each other and the underlying field.


#
