# script for cleaning demographics columns (excludes all testing and cancer history columns)
## ---------------------------------------------------------------------------------------
# recode gender
datw <-
  dat %>%
  mutate(Gender = recode(Gender, 'F' = 0, 'M' = 1),
         Gender = na_if(Gender, "U"))

## ---------------------------------------------------------------------------------------
# drop missing in FatherID, MotherID
datw <- datw[!(is.na(datw$FatherID) & is.na(datw$MotherID)), ]  

## ---------------------------------------------------------------------------------------
len <- dim(datw)[1]
datw$Twins <- rep(0, len)
i_list <- which(datw$isProband == "Identical Twin")
j = 1
for(i in i_list){
  # find their proband twin, give twin pair unique number
  if(datw$Twins[i] == 0){
    datw$Twins[which((datw$isProband == "Proband" | datw$isProband == "Identical Twin") & 
                       datw$PedigreeID == datw$PedigreeID[i] & 
                       datw$FatherID == datw$FatherID[i] & 
                       datw$MotherID == datw$MotherID[i])] <- j
    datw$Twins[i] <- j
    
    j = j + 1
  }
}

# manually fix twin pair 11
datw[which(datw$Twins == 11 & datw$year.of.DOB != 1979), ]$Twins <- 0

## ---------------------------------------------------------------------------------------
# recode
datw$isProband <- ifelse(datw$isProband == "Proband", 1, 0)

## ---------------------------------------------------------------------------------------
# for each subject, if "Maternal.AJ.Ancestry  == 1" or "Paternal.AJ.Ancestry == 1", then "ethnic = AJ"
# no way of knowing if someone is "Italian" or "Other"

# no missing values, so
datw$ethnic <- ifelse(datw$Maternal.AJ.Ancestry + datw$Paternal.AJ.Ancestry >= 1, "AJ", "nonAJ")

## ---------------------------------------------------------------------------------------
# convert Death to correct format
datw$Death <- as.numeric(datw$isDead)

## ---------------------------------------------------------------------------------------
# assuming CurAge is age at which they died
datw$AgeDeath <- datw$CurAge

## ---------------------------------------------------------------------------------------
# recode isHispanic
datw <- 
  datw %>%
  mutate(isHispanic = recode(isHispanic, '2' = '0', '1' = '1', .default = NA_character_),
         isHispanic = as.numeric(isHispanic))

## ---------------------------------------------------------------------------------------
# BayesMendel: Asian, Black, Hispanic, NativeAmerican and White
# PanelPRO: one of All_Races, AIAN, Asian, Black, White, Hispanic, WH, WNH, NA
# recode into 0/1/NA
datw <-
  datw %>%
  mutate(Race..AI.or.AN = recode(Race..AI.or.AN, 'N' = '0', 'Y' = '1', .default = NA_character_),
         #Race..AI.or.AN = as.numeric(Race..AI.or.AN),
         RaceCol_AIAN = as.numeric(Race..AI.or.AN),
        
         Race..Asian = recode(Race..Asian, 'N' = '0', 'Y' = '1', .default = NA_character_),
         #Race..Asian = as.numeric(Race..Asian),
         RaceCol_Asian = as.numeric(Race..Asian),
         
         Race..Bl.or.AA = recode(Race..Bl.or.AA, 'N' = '0', 'Y' = '1', .default = NA_character_),
         #Race..Bl.or.AA = as.numeric(Race..Bl.or.AA),
         RaceCol_Black = as.numeric(Race..Bl.or.AA),
         
         Race..NH.or.other.PI = recode(Race..NH.or.other.PI, 'N' = '0', 'Y' = '1', .default = NA_character_),
         #Race..NH.or.other.PI = as.numeric(Race..NH.or.other.PI),
         RaceCol_NHPI = as.numeric(Race..NH.or.other.PI),
         
         Race..Other = recode(Race..Other, 'N' = '0', 'Y' = '1', .default = NA_character_),
         #Race..Other = as.numeric(Race..Other),
         RaceCol_Other = as.numeric(Race..Other),
         
         Race..White = recode(Race..White, 'N' = '0', 'Y' = '1', .default = NA_character_),
         #Race..White = as.numeric(Race..White),
         RaceCol_White = as.numeric(Race..White),
         
         RaceCol_Hispanic = isHispanic
  )


## ---------------------------------------------------------------------------------------
# NHPI -> Asian
datw$RaceCol_Asian <- ifelse(datw$RaceCol_Asian == 1 | datw$RaceCol_NHPI == 1, 1, 0)
datw <- subset(datw, select = -RaceCol_NHPI)

## ---------------------------------------------------------------------------------------
# if hispanic == 1 and (native american, asian, black, or white) == 1, then make hispanic == 0?
datw$BMRaceCol_Hispanic <- datw$isHispanic
datw$BMRaceCol_Hispanic[which((  datw$RaceCol_Asian == 1  | 
                                 datw$RaceCol_AIAN == 1   | 
                                 datw$RaceCol_Black == 1  | 
                                 datw$RaceCol_White == 1  ) & 
                               datw$RaceCol_Hispanic == 1)] <- 0

## ---------------------------------------------------------------------------------------
# BayesMendel: Asian, Black, Hispanic, NativeAmerican and White
datw <-
  datw %>%
  mutate(BMRaceCol_NativeAmerican = RaceCol_AIAN,
         BMRaceCol_Asian = RaceCol_Asian,
         BMRaceCol_Black = RaceCol_Black,
         BMRaceCol_White = RaceCol_White)

## ---------------------------------------------------------------------------------------
# get race col names
Names <- colnames(datw)
Names_bm <- Names[grepl("^BMRaceCol", Names)]

## ---------------------------------------------------------------------------------------
# BayesMendel: Asian, Black, Hispanic, NativeAmerican and White
# create BM_NA column
datw$sums_bm <- rowSums(datw[, Names_bm], na.rm = TRUE)
datw$BMRaceCol_Unknown <- ifelse(datw$sums_bm != 1 | datw$RaceCol_Other == 1, 1, 0)

## ---------------------------------------------------------------------------------------
for(x in Names_bm){
  datw[which(datw$BMRaceCol_Unknown == 1), x] <- 0
}
datw$sums_bm <- rowSums(datw[, Names_bm], na.rm = TRUE)

## ---------------------------------------------------------------------------------------
# race_bm column
copy_dat <- datw
copy_dat <- copy_dat %>%
 pivot_longer(
   cols = starts_with("BMRaceCol"),
   names_to = "race_bm",
   names_prefix = "BMRaceCol_",
   values_to = "count"
 )

## ---------------------------------------------------------------------------------------
# with race column
dat_new <- copy_dat[which(copy_dat$count == 1), ]


## ---------------------------------------------------------------------------------------
# PanelPRO: one of All_Races, AIAN, Asian, Black, White, Hispanic, WH, WNH, NA
# create WH
dat_new$RaceCol_WH <- ifelse(dat_new$RaceCol_White == 1 & dat_new$RaceCol_Hispanic == 1, 1, 0)

# create WNH
dat_new$RaceCol_WNH <- ifelse(dat_new$RaceCol_White == 1 & dat_new$RaceCol_Hispanic == 0, 1, 0)

# set RaceCol_White = 0 and RaceCol_Hispanic = 0 if RaceCol_WH == 1 or RaceCol_WNH == 1
dat_new[which(dat_new$RaceCol_WH == 1 | dat_new$RaceCol_WNH == 1), c("RaceCol_White", "RaceCol_Hispanic")] <- 0

## ---------------------------------------------------------------------------------------
# PanelPRO: one of All_Races, AIAN, Asian, Black, White, Hispanic, WH, WNH, NA
# get race col names
Names_pp <- colnames(dat_new)
Names_pp <- Names_pp[grepl("^RaceCol", Names_pp)]

## ---------------------------------------------------------------------------------------
# PanelPRO: one of All_Races, AIAN, Asian, Black, White, Hispanic, WH, WNH, NA
# create dummy variable for NA
dat_new$sums <- rowSums(dat_new[, Names_pp], na.rm = TRUE)
dat_new$RaceCol_NA <- ifelse(dat_new$sums == 0, 1, 0)

## ---------------------------------------------------------------------------------------
# PanelPRO: one of All_Races, AIAN, Asian, Black, White, Hispanic, WH, WNH, NA
# create All Races label
# if subject is more than 1 race or 'Other', mark as All Races
# new column for sums to see and check
dat_new$RaceCol_All_Races <- ifelse(dat_new$sums > 1 | dat_new$RaceCol_Other == 1, 1, 0)

## ---------------------------------------------------------------------------------------
for(x in Names_pp){
  dat_new[which(dat_new$RaceCol_All_Races == 1), x] <- 0
}
dat_new$sums <- rowSums(dat_new[, Names_pp], na.rm = TRUE)
dat_new <- subset(dat_new, select = -RaceCol_Other)

## ---------------------------------------------------------------------------------------
# race_pp column
copy_datw <- dat_new
copy_datw <- copy_datw %>%
 pivot_longer(
   cols = starts_with("RaceCol"),
   names_to = "race_pp",
   names_prefix = "RaceCol_",
   values_to = "count_pp")

## ---------------------------------------------------------------------------------------
# with race column
dat_new <- copy_datw[which(copy_datw$count_pp == 1), ]

## ---------------------------------------------------------------------------------------
# recode 'NA' into missing
dat_new <- dat_new %>%
  mutate(race_pp = na_if(race_pp, 'NA'))

## ---------------------------------------------------------------------------------------
col_names <- colnames(dat_new)
col_drop <- c("Maternal.AJ.Ancestry", "Paternal.AJ.Ancestry", "year.of.DOB", "sums", "count", "sums_bm", "count_pp")
col_drop <- c(col_drop, col_names[grepl("^Race..", col_names)])
dat.demographics.clean <- dat_new[- which(col_names %in% col_drop)]
