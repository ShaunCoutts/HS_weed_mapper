# script to get the hyper spectral data and vegetation indices, also map them
# and explore the images
# install.packages('moments', repos = "https://cloud.r-project.org")

library(rgdal)
library(ggplot2)
library(hexbin)
library(scales)
library(viridis)
library(GISTools)
library(dplyr)
library(tidyr)
library(ggmap)
library(RgoogleMaps)
library(dismo)
library(rgeos)
library(raster)
library(fastcluster)
library(caret)
library(e1071)
library(moments)

veg_ind_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/Data/sample_data/UoS/hyperspectral/vegetationindices/'
den_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/Data/dens_data/den_shape_2018/'
data_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/Data/'
plot_loc = '/Users/shauncoutts/Dropbox/projects/hyper_spec_weed_mapper/outputs/plots'

################################################################################
# Black grass density of fields, I want a polygone of each 20m by 20m grid, with a
# density state and field ID
sites = unique(list.files(den_loc, pattern = "\\.shp$")) %>% #_2016",pattern = "\\.shp$"))
	as.vector() %>%
	as.character() %>%
	(function(x) gsub('.shp', '', x))

# get the list of strip names
strips = unique(list.files(veg_ind_loc, pattern = "\\.img$"))

# use the naming convention of using the Pcode to remove the sites already done
proc_fields =  unique(list.files(data_loc, pattern = "^den_HS_st"))

site_names = gsub(" _combi", '', sites)
done_names = unique(substr(proc_fields, start = 12, stop = 26))
# sites = sites[!site_names %in% done_names]

################################################################################
# make a list of dataframes to hold the meta data each cell joined with extracted HS data
# loop over each strip to see which fields match which strips, note the extraction
# takes a couple of hours per field
for(st in seq_along(strips)){

	setwd(veg_ind_loc)
	link_veg = brick(strips[st])

	# Loop through each site
	for(i in 1:length(sites)){

		###load merged shapefile
		setwd(den_loc)
	    shape = readOGR(dsn = ".",sites[i])  #dsn = the directory (without a trailing backslash); layer = the shapefile name without the .shp.
		# change the projection to match the hyper-spec dataset
		shape_proj = spTransform(shape, crs(link_veg))
		# now remove the outer layer of grid squares
		shape_crop = subset(shape_proj, !is.na(DnstySt))
		shape_crop = subset(shape_crop, apply(gTouches(shape_crop, byid = TRUE,
			returnDense = TRUE), MARGIN = 1, FUN = sum) == 8)

		# check if the field polygons overlap the strip
		if(!is.null(raster::intersect(extent(link_veg), extent(shape_crop)))){

			HS_list = list()
			counter = 1

			for(j in 1:dim(shape_crop)[1]){

				print(paste0(st, ':', i, ':', j))

				# get the bit of the density data I want
				single_grid = shape_crop[j, ]
				single_grid@data$uID = paste0('st', st, ':', single_grid$farm, ':',
					single_grid$Pcode, ':', single_grid$quadcnv)

				HS_extract = raster::extract(link_veg, single_grid) # check that 20 / xres(link_veg_S1), is roughly equal to sqrt(dim(HS_extract[[1]])[1])

				#take out all the NAs and join to the density state data
				HS_list[[counter]] = filter(as_data_frame(HS_extract[[1]]), !is.na(Red.Green.Ratio.Image)) %>%
					mutate(uID = single_grid$uID) %>%
					inner_join(single_grid@data, by = 'uID')

				#TODO, in here I need to also add the coordinates for each pixel so I can match up the clustered data to
				# check what various clusters and metrics mean in terms of the images, I need to locate those pixels in the
				# images.

				counter = counter + 1

			}

			# save the resulting data
			setwd(data_loc)
			Pcode_ID = as.character(unique(na.omit(sapply(HS_list, FUN = function(x) x$Pcode[1])))[1])
			write.csv(bind_rows(HS_list),
				file = paste0('den_HS_st', st, '_', Pcode_ID, '.csv'))

		}
	}
}

################################################################################
# read in and explore the data a bit
df_names = unique(list.files(data_loc, pattern = "^den_HS_st"))

setwd(data_loc)
HS_den_df = read.csv(df_names[1], stringsAsFactor = FALSE)

# make a summary for each cell
HS_sum = group_by(HS_den_df, field, Pcode, quad) %>%
	summarise(crop = paste(unique(crop), sep = ':'), X = mean(x), Y = mean(y),
		den_state = unique(DnstySt),
		ARI_1_mean = mean(Anthocyanin.Reflectance.Index),
		ARI_2_mean = mean(Anthocyanin.Reflectance.Index.2..Modified.Anthocyanin.Reflectance.Index.),
		CRI_mean = mean(Carotenoid.Reflectance.Index.1),
		red_edge_NDVI_mean = mean(Modified.Red.Edge.NDVI),
		SAV_mean = mean(Modified.Soil.Adjusted.Vegetation.Index.2),
		inf_red_mean = mean(Normalized.Difference.Infrared.Index),
		lign_mean = mean(Normalized.Difference.Lignin.Index, rm.na = TRUE),
		veg_ind_mean = mean(Normalized.Difference.Vegetation.Index),
		H2O_ind_mean = mean(Normalized.Difference.Water.Index),
		PCR_mean = mean(Photochemical.Reflectance.Index),
		red_edge_mean = mean(Red.Edge.Position.Index..linear.interpolation.),
		RG_mean = mean(Red.Green.Ratio.Image),
		SIP_mean = mean(Structure.Insensitive.Pigment.Index.1),
		chlor_abs_mean = mean(Transformed.Chlorophyll.Absorption.in.Reflectance.Index),
		chlor_abs_soil_mean = mean(Transformed.Chlorophyll.Absorption.in.Reflectance.Index.Optimised.Soil.Adjusted.Vegetation.Index),
		tri_veg = mean(Triangular.Vegetation.Index),
		volgA = mean(Vogelmann.Red.Edge.A.or.1),
		volgB = mean(Vogelmann.Red.Edge.B.or.2),
		volgC = mean(Vogelmann.Red.Edge.C.or.3),
		water_band = mean(Water.Band.Index)) %>% # resacel so max value = 1
	mutate(ARI_1_mean = rescale(ARI_1_mean), ARI_2_mean = rescale(ARI_2_mean),
		CRI_mean = rescale(CRI_mean), red_edge_NDVI_mean = rescale(red_edge_NDVI_mean),
		SAV_mean = rescale(SAV_mean), inf_red_mean = rescale(inf_red_mean),
		lign_mean = rescale(lign_mean), veg_ind_mean = rescale(veg_ind_mean),
		H2O_ind_mean = rescale(H2O_ind_mean), PCR_mean = rescale(PCR_mean),
		red_edge_mean = rescale(red_edge_mean), RG_mean = rescale(RG_mean),
		SIP_mean = rescale(SIP_mean), chlor_abs_mean = rescale(chlor_abs_mean),
		chlor_abs_soil_mean = rescale(chlor_abs_soil_mean), tri_veg = rescale(tri_veg),
		volgA = rescale(volgA), volgB = rescale(volgB), volgC = rescale(volgC),
		water_band = rescale(water_band)) %>%
	gather(key = 'index', value = 'value', ARI_1_mean:water_band)

# map each index as a color and show density class
sum_plt = ggplot(HS_sum, aes(x = X, y = Y)) +
	geom_tile(aes(fill = value), color = 'white') +
	geom_text(aes(label = den_state), color = grey(0.7), size = 2) +
	scale_fill_viridis() +
	facet_wrap(~index, nrow = 5, scales = 'free')

################################################################################
# make some summaries of the varaibles to see how they look, and how they relate to
# each other. There are going to be lots of rows, far more than I need and it will
# make it unwieldly to work with, sample from within each grid cell and strip.

proc_fields =  unique(list.files(data_loc, pattern = "^den_HS_st"))

setwd(data_loc)
HS_samp_list = list()

# read in a strip, take a sample, then read in another
for(i in seq_along(proc_fields)){

	print(i)
	HS_den_df = read.csv(proc_fields[i])

	HS_samp_list[[i]] = group_by(HS_den_df, field, Pcode, quad) %>%
		sample_frac(0.05, replace = TRUE) %>%
		mutate(strip = paste0('St', i))

}

HS_samp = ungroup(bind_rows(HS_samp_list)) %>% as.data.frame()

################################################################################
# make a paris plot with 2D histograms for each pair of predictors
veg_inds = names(HS_samp)[2:21]
setwd(plot_loc)
pdf(file = 'pairs_plots_veg_inds.pdf', width = 6, height = 8)

	for(i in 1:(length(veg_inds) - 1)){

		print(paste0(i, '/', (length(veg_inds) - 1)))

		for(j in (i + 1):length(veg_inds)){

			plt = ggplot(HS_samp, aes_string(x = veg_inds[i], y = veg_inds[j])) +
				geom_hex() +
				scale_fill_viridis() +
				facet_wrap(~strip, ncol = 2)

			print(plt)

		}
	}

dev.off()

# immediatly can see that a couple of the indexes are very strongly corrleated with
# each other Volgmann A, B and C all very corrlated with each other, only need one
# struct.insene.pigment.1 and norm.dif.veg.index have correlation near 1, only need one

HS_samp_redu = dplyr::select(HS_samp, Anthocyanin.Reflectance.Index:Normalized.Difference.Infrared.Index,
	Normalized.Difference.Vegetation.Index:Red.Green.Ratio.Image,
	Transformed.Chlorophyll.Absorption.in.Reflectance.Index:Vogelmann.Red.Edge.A.or.1,
	Water.Band.Index:strip)

# look at the reduced set to make sure they are none that are very very correlated
veg_inds = names(HS_samp_redu)[1:16]
setwd(plot_loc)
pdf(file = 'pairs_plots_veg_inds_redu.pdf', width = 6, height = 8)

	for(i in 1:(length(veg_inds) - 1)){

		print(paste0(i, '/', (length(veg_inds) - 1)))

		for(j in (i + 1):length(veg_inds)){

			plt = ggplot(HS_samp, aes_string(x = veg_inds[i], y = veg_inds[j])) +
				geom_hex() +
				scale_fill_viridis() +
				facet_wrap(~strip, ncol = 2)

			print(plt)

		}
	}

dev.off()

# No multi-modal distributions, which is good, becuase they would be really annoying to
# deal with

######################################################################################
# Still lots of correlation, try some clustering tools/unsupervised learners, then PCA
# do this by strip to see if there is consistency between strips, if the strip
# clusters are not consistent there is probably not much hope for cleaning the data
# in this way

strip_ID = unique(HS_samp_redu$strip)

setwd(plot_loc)
pdf(file = 'hier_clust_by_strip.pdf', width = 20, height = 20)

	for(i in strip_ID){

		print(i)

		samp_st = filter(HS_samp_redu, strip == i)
		preds_st = data.matrix(dplyr::select(samp_st, Anthocyanin.Reflectance.Index:Water.Band.Index))
		cl_st = hclust.vector(preds_st, method = 'ward', metric = 'euclidean')

		plot(cl_st)

	}

dev.off()

# They do cluster out in a similar way, try k-means clustering to see how the
# clusters break down across predictors. Try all the strips together, try and break
# the data into the two main clusters

# first, center and scale the preictors the dataset
HS_samp_pred = dplyr::select(HS_samp_redu, Anthocyanin.Reflectance.Index:Water.Band.Index) %>%
	base::scale() %>%
	as_data_frame()

clust_2 = kmeans(HS_samp_pred, 2, nstart = 10)

# This clustering didn't work very well, the heiricial clustering is pretty slow
# and hard to work out what the clusters mean. This does not seem like a good way to
# filter the data, might be better to just keep the whole thing and PCA it all
prePro_HS = preProcess(HS_samp_pred, method = c('center', 'scale', 'BoxCox', 'pca'),
	thresh = 0.95)

HS_PCA = predict(prePro_HS, HS_samp_pred)

################################################################################
# Do this to for the full data, rather than just the samples, As the idea is to
# take different summaries of the distribution within each grid cell
proc_fields =  unique(list.files(data_loc, pattern = "^den_HS_st"))

setwd(data_loc)
HS_full_list = list()

# read in a strip, take a sample, then read in another
for(i in seq_along(proc_fields)){

	print(i)
	HS_den_df = read.csv(proc_fields[i])

	HS_full_list[[i]] = mutate(HS_den_df, strip = paste0('St', i))

}

HS_full = ungroup(bind_rows(HS_full_list)) %>% as.data.frame()

write.csv(HS_full, file = 'hyper_spec_full_extracted.csv')

################################################################################
# PCA reduced the dimentionally down to 8 dimentions
setwd(data_loc)
HS_full = read.csv('hyper_spec_full_extracted.csv', stringsAsFactor = FALSE)

HS_full_redu = dplyr::select(HS_full, Anthocyanin.Reflectance.Index:Normalized.Difference.Infrared.Index,
		Normalized.Difference.Vegetation.Index:Red.Green.Ratio.Image,
		Transformed.Chlorophyll.Absorption.in.Reflectance.Index:Vogelmann.Red.Edge.A.or.1,
		Water.Band.Index:strip) %>%
	filter(Modified.Red.Edge.NDVI > 0.4) # filters out the bare earth like tractor tracks.

HS_full_pred = dplyr::select(HS_full_redu, Anthocyanin.Reflectance.Index:Water.Band.Index) %>%
	as_data_frame()

prePro_HS = preProcess(HS_full_pred, method = c('center', 'scale', 'pca'),
	thresh = 0.95)

HS_PCA = predict(prePro_HS, HS_full_pred) # there a 7 dimetions needed to capture 95% of variance

HS_full_PCA = cbind(dplyr::select(HS_full_redu, uID:strip), HS_PCA)

setwd(data_loc)
write.csv(HS_full_PCA, file = 'hyper_spec_full_PCA.csv')

################################################################################
# These files are very large to work with, and I need to down scale them anyway, to a
# distribution then take summaries of that distribution, mean, varaince, skewness
# kurtosis, start with the PCA data as it will be nicer to work with
setwd(data_loc)
HS_full_PCA = read.csv('hyper_spec_full_PCA.csv', header = TRUE, stringsAsFactor = FALSE)

# Take the four summaries for each PCA dimention to up scale the hyper-spec data
DS_PCA = group_by(HS_full_PCA, strip, Pcode, quad, uID) %>%
 	summarise(den_st = unique(DnstySt),
		PC1_mean = mean(PC1), PC1_var = var(PC1), PC1_sk = skewness(PC1), PC1_kur = kurtosis(PC1),
		PC2_mean = mean(PC2), PC2_var = var(PC2), PC2_sk = skewness(PC2), PC2_kur = kurtosis(PC2),
		PC3_mean = mean(PC3), PC3_var = var(PC3), PC3_sk = skewness(PC3), PC3_kur = kurtosis(PC3),
		PC4_mean = mean(PC4), PC4_var = var(PC4), PC4_sk = skewness(PC4), PC4_kur = kurtosis(PC4),
		PC5_mean = mean(PC5), PC5_var = var(PC5), PC5_sk = skewness(PC5), PC5_kur = kurtosis(PC5),
		PC6_mean = mean(PC6), PC6_var = var(PC6), PC6_sk = skewness(PC6), PC6_kur = kurtosis(PC6),
		PC7_mean = mean(PC7), PC7_var = var(PC7), PC7_sk = skewness(PC7), PC7_kur = kurtosis(PC7))

setwd(data_loc)
write.csv(DS_PCA, file = 'hyper_spec_dist_sum_PCA.csv')







################################################################################
#see if aieral imagry lines up for the one field we know we have
setwd(veg_ind_loc)
link_veg = brick(strips[5])
plot(link_veg[[20]])

setwd(den_loc)
shape = readOGR(dsn = ".",sites[2])  #dsn = the directory (without a trailing backslash); layer = the shapefile name without the .shp.
# change the projection to match the hyper-spec dataset
shape_proj = spTransform(shape, crs(link_veg))
# now remove the outer layer of grid squares
shape_crop = subset(shape_proj, !is.na(DnstySt))
shape_crop = subset(shape_crop, apply(gTouches(shape_crop, byid = TRUE,
	returnDense = TRUE), MARGIN = 1, FUN = sum) == 8)

plot(shape_crop, add = TRUE)
