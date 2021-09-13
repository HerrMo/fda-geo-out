### Qualitative analysis of real world data

# source(here::here("R/setup.R")) # uncomment setup.R for stand-alone usage

## ecg200 data

set.seed(446)
ecg200_tr <- read.table("~/data/ECG200/ECG200_TRAIN.txt")
ecg200_te <- read.table("~/data/ECG200/ECG200_TEST.txt")

lbls <- c(ecg200_tr$V1, ecg200_te$V1)
ecg200 <- rbind(ecg200_tr, ecg200_te)[ , -1]

out_ecg200 <- ecg200
d_l2_out_ecg200 <- dists(ecg200)
emb_out_ecg200 <- embed(d_l2_out_ecg200, "mds", k = 5)


# contamination
# set.seed(446)
# smpl <- sample(which(lbls == -1), 0.1 * length(lbls == 1))
# out_ecg200 <- rbind(ecg200[lbls == 1, ], ecg200[smpl, ])
#
# d_l2_out_ecg200 <- dists(out_ecg200)
# emb_out_ecg200 <- embed(d_l2_out_ecg200, "mds", k = 5)

##-----------------------------------------------
# Embedding dimensions for the following datasets
dims <- 5
##-----------------------------------------------

### octane
set.seed(946)

data(Octanespectrum)
data(Octanevalues)

oct <- t(Octanespectrum$y)
oct_lbls <- Octanevalues
d_l2_oct <- dists(oct)
oct_emb <- embed(d_l2_oct, "mds", k = dims)

mPts <- 0.75 * nrow(extract_points(oct_emb))
oct_lof <- lof(oct_emb, minPts = mPts)

plt_emb_oct <- plot_emb_temp(oct_emb,
                             col = oct_lof,
                             pch = ifelse(frank(oct_emb, ndim = dims) > 12, 1, 3))
plt_funs_oct <- plot_funs(t(t(oct) - colMeans(oct))[rank(-oct_lof),], col = oct_lof[rank(-oct_lof)])


### spanish weather
set.seed(946)
data(aemet)
weather <- aemet$temp$data

d_l2_weather <- dists(weather)
emb_weather <- embed(d_l2_weather, "mds", k = dims)

mPts <- 0.75 * nrow(extract_points(emb_weather))
plt_emb_we <- plot_emb_temp(emb_weather,
                            color = lof(emb_weather, minPts = mPts),
                            pch = ifelse(frank(emb_weather, ndim = dims) > 12, 1, 3))
plt_funs_we <- plot_funs(weather, col = lof(emb_weather, minPts = mPts))


### tecator
set.seed(534)

data(tecator)
tec <- tecator$absorp.fdata$data

plot_funs(tec, args = tecator$absorp.fdata$argvals)
d_l2_tec <- dists(tec)
emb_tec <- embed(d_l2_tec, "mds", k = dims)

mPts <- 0.75 * nrow(extract_points(emb_tec))
plt_emb_tec <- plot_emb_temp(emb_tec,
                             color = lof(emb_tec, minPts = mPts),
                             pch = ifelse(frank(emb_tec, ndim = dims) > 12, 1, 3))
plt_funs_tec <- plot_funs(t(t(tec)-colMeans(tec)), col = lof(emb_tec, minPts =  mPts))


### wine data
set.seed(297)
wine_tr <- read.table("~/data/Wine/Wine_TRAIN.txt")
wine_te <- read.table("~/data/Wine/Wine_TEST.txt")

wine_lbls <- c(wine_tr$V1, wine_te$V1)
wine <- rbind(wine_tr, wine_te)[ , -1]

d_l2_wine <- dists(wine)
emb_wine <- embed(d_l2_wine, "mds", k = dims)

mPts <- 0.75 * nrow(extract_points(emb_wine))
wine_lof <- lof(emb_wine, minPts = mPts)

plt_emb_wine <- plot_emb_temp(emb_wine,
                              col = wine_lof,
                              pch = ifelse(frank(emb_wine, ndim = dims) > 12, 1, 3))
plt_funs_wine <- plot_funs(t(t(wine) - colMeans(wine))[rank(wine_lof),], col = wine_lof[rank(wine_lof)])


