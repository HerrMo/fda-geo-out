library(fda.usc) # tecator data, weather data
library(fds) # octane data

## ecg200 data --------------------------------------------------------------------------------------------------------

set.seed(446)
ecg200_tr <- read.table("data/raw/ECG200_TRAIN.txt")
ecg200_te <- read.table("data/raw/ECG200_TEST.txt")

lbls <- c(ecg200_tr$V1, ecg200_te$V1)
ecg200 <- rbind(ecg200_tr, ecg200_te)[ , -1]

out_ecg200 <- ecg200
save(out_ecg200, file = "data/processed/ECG.RData")

# contamination
#set.seed(446)
# smpl <- sample(which(lbls == -1), 0.1 * length(lbls == 1))
#
# out_ecg200 <- rbind(ecg200[lbls == 1, ], ecg200[smpl, ])
#
# save(out_ecg200, file = "data/processed/ECG_cont.RData" )


### octane -------------------------------------------------------------------------------------------------------------
set.seed(946)

data(Octanespectrum)
data(Octanevalues)

oct <- t(Octanespectrum$y)
oct_lbls <- Octanevalues

save(oct, file = "data/processed/Octane.RData")

### spanish weather ----------------------------------------------------------------------------------------------------
set.seed(946)
data(aemet)
weather <- aemet$temp$data

save(weather, file = "data/processed/Weather.RData")

### tecator ------------------------------------------------------------------------------------------------------------
set.seed(534)
data(tecator)
tec <- tecator$absorp.fdata$data

save(tec, file = "data/processed/Tecator.RData")

### wine data ----------------------------------------------------------------------------------------------------------
set.seed(297)
wine_tr <- read.table("data/raw/Wine_TRAIN.txt")
wine_te <- read.table("data/raw/Wine_TEST.txt")

wine_lbls <- c(wine_tr$V1, wine_te$V1)
wine <- rbind(wine_tr, wine_te)[ , -1]

save(wine, file = "data/processed/Wine.RData")


