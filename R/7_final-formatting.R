# final formatting and error checking for BayesMendel and PanelPRO

#### BayesMendel Format ####
# select relevant columns
bm.clean <-
  dat.clean %>%
  select(
    SubjectID,
    PedigreeID,
    Region,
    isHispanic,
    isProband,
    ID,
    Gender,
    FatherID,
    MotherID,
    Twins,
    ethnic,
    race_bm,
    Death,
    AgeDeath,
    CarrierStatus,
    isAffBC,
    isAffOC,
    isAffCBC,
    isAffBC2,
    AgeBC,
    AgeOC,
    AgeCBC,
    AgeBC2,
    BRCA1,
    BRCA2,
    ER,
    PR,
    HER2,
    Oophorectomy,
    AgeOophorectomy,
    Mastectomy,
    AgeMastectomy,
    isAffCOL,
    isAffENDO,
    AgeCOL,
    AgeENDO,
    MLH1,
    MSH2,
    MSH6,
    MSI,
    isAffPANC,
    AgePANC,
    isAffMELA,
    AgeMELA,
    P16
  ) %>%
  add_column(CK14 = rep(0,nrow(dat.clean)), 
             CK5.6 = rep(0,nrow(dat.clean)),
             .after = "ER") %>%
  add_column(location = rep(0,nrow(dat.clean)), .after = "MSI") %>%
  rename('AffectedBreast' = isAffBC,
         'AffectedOvary' = isAffOC,
         'AgeBreast' = AgeBC,
         'AgeOvary' = AgeOC,
         'AffectedBreastContralateral' = isAffCBC,
         'AffectedBreast2' = isAffBC2,
         'AgeBreastContralateral' = AgeCBC,
         'AgeBreast2' = AgeBC2,
         'AffectedColon' = isAffCOL,
         'AffectedEndometrium' = isAffENDO,
         'AgeColon' = AgeCOL,
         'AgeEndometrium' = AgeENDO,
         'AffectedPancreas' = isAffPANC,
         'AgePancreas' = AgePANC,
         'AffectedSkin' = isAffMELA,
         'AgeSkin' = AgeMELA,
         'race' = race_bm)

##### Recode Gene Columns ####
bm.clean <-
  bm.clean %>%
  mutate(across(.cols = c(BRCA1,BRCA2,MLH1,MSH2,MSH6,P16,CarrierStatus),
                ~ recode(.,
                         "1" = "1",
                         "2" = "3",
                         "0" = "2"))) %>%
  mutate(across(.cols = c(BRCA1,BRCA2,MLH1,MSH2,MSH6,P16,CarrierStatus),
                ~ replace_na(., 0))) %>%
  mutate(across(.cols = c(BRCA1,BRCA2,MLH1,MSH2,MSH6,P16,CarrierStatus),
                ~ as.numeric(.)))

##### Race ####
# in BayesMendel each family can have only one race, therefore 
# iterate through each family
peds <- unique(bm.clean$PedigreeID)
for(fam in 1:length(peds)){
  if(fam%%1000 == 0){ print(paste0("Unique family race: ",fam," of ",length(peds))) }
  
  # get the temp pedigree
  tmp.ped.id <- peds[fam]
  tmp.fam <- bm.clean[which(bm.clean$PedigreeID == tmp.ped.id),
                      c("PedigreeID","race")]
  
  # find unique races, but exclude "Unknown"
  all.races <- unique(tmp.fam$race)
  all.races <- all.races[!all.races %in% "Unknown"]
  num.races <- length(all.races)
  
  # if 1 race, assume all family members are the same race
  if(num.races == 1){
    bm.clean$race[which(bm.clean$PedigreeID == tmp.ped.id)] <- all.races
  } else {
    bm.clean$race[which(bm.clean$PedigreeID == tmp.ped.id)] <- "Unknown"
  }
}

##### Affected Age Columns ####

# update affected age column with death age, if person is unaffected
# iterate through each BM cancer type
bm.cancers <- c("Breast","Ovary","Colon","Endometrium","Pancreas","Skin","BreastContralateral","Breast2")
bm.clean <- as.data.frame(bm.clean)
for(cncr in bm.cancers){
  
  # loop vars
  age.col <- paste0("Age",cncr)
  aff.col <- paste0("Affected",cncr)
  
  # find unaffected people and iterate through them
  unaff <- bm.clean$SubjectID[which(bm.clean[,aff.col] == 0)]
  
  # update the affected age with death age
  bm.clean[which(bm.clean$SubjectID %in% unaff),age.col] <- bm.clean$AgeDeath[which(bm.clean$SubjectID %in% unaff)]
}


##### Final Formatting ####
bm.clean <-
  bm.clean %>%
  arrange(PedigreeID,ID)

# change to data frame
bm.clean <- as.data.frame(bm.clean)

#### PanelPRO ####
# select relevant columns
pp.clean <-
  dat.clean %>%
  select(
    SubjectID,
    PedigreeID,
    Region,
    ID,
    Gender,
    MotherID,
    FatherID,
    isProband,
    CurAge,
    isDead,
    race_pp,
    ethnic,
    isHispanic,
    CarrierStatus,
    Twins,
    Oophorectomy,
    AgeOophorectomy,
    Mastectomy,
    AgeMastectomy,
    ER,
    PR,
    HER2,
    MSI,
    starts_with("isAff"),
    all_of(paste0("Age",pp.cancers)),
    all_of(pp.genes)
  ) %>%
  add_column(CK14 = rep(NA,nrow(dat.clean)), 
             CK5.6 = rep(NA,nrow(dat.clean)),
             .after = "ER") %>%
  rename('Sex' = Gender,
         'Ancestry' = ethnic,
         'race' = race_pp)

##### Marker Tests ####
# change from BayesMendel coding to PanelPRO coding
# NA = no test, 0 = neg, 1 = pos
pp.clean <-
  pp.clean %>%
  mutate(across(.cols = c(CK14,CK5.6), ~ as.numeric(.))) %>%
  mutate(across(.cols = c(ER,PR,HER2,CK14,CK5.6,MSI), ~ recode(.,
                                                               "0"="100",
                                                               "2"="0",
                                                               "1"="1"))) %>%
  mutate(across(.cols = c(ER,PR,HER2,CK14,CK5.6,MSI), ~ na_if(.,"100"))) %>%
  mutate(across(.cols = c(ER,PR,HER2,CK14,CK5.6,MSI), ~ as.numeric(.)))
           

##### Prophylactic Surgeries ####
pro.surg <- 
  pp.clean %>%
  select(SubjectID,Mastectomy,Oophorectomy,AgeMastectomy,AgeOophorectomy) %>%
  filter(if_any(.cols = c(Mastectomy,Oophorectomy), ~ . == 1))

# storage
riskmod <- list()
interAge <- list()

# iterate through each subject
ids <- pro.surg$SubjectID
for(id in 1:length(ids)){
  if(id%%1000 == 0){print(paste0("Prophylactic Surgery Patient ",id," of ",nrow(pro.surg)))}
  mast <- pro.surg$Mastectomy[id]
  ooph <- pro.surg$Oophorectomy[id]
  mastAge <- as.numeric(pro.surg$AgeMastectomy[id])
  oophAge <- as.numeric(pro.surg$AgeOophorectomy[id])
  if(mast == 1 & ooph == 1){
    tmp.mod <- list("Mastectomy","Oophorectomy")
    tmp.age <- list(mastAge,oophAge)
  } else if(mast == 1){
    tmp.mod <- list("Mastectomy")
    tmp.age <- list(mastAge)
  } else if(ooph == 1){
    tmp.mod <- list("Oophorectomy")
    tmp.age <- list(oophAge)
  } else {
    tmp.mod <- NA
    tmp.age <- NA
  }
  
  # store results
  riskmod[[id]] <- tmp.mod
  interAge[[id]] <- tmp.age
}

# create empty columns
pp.clean$riskmod <- rep(NA,nrow(pp.clean))
pp.clean$interAge <- rep(NA,nrow(pp.clean))

# populate empty columns
pp.clean$riskmod[which(pp.clean$SubjectID %in% pro.surg$SubjectID)] <- riskmod
pp.clean$interAge[which(pp.clean$SubjectID %in% pro.surg$SubjectID)] <- interAge

# move interAge and riskmod
pp.clean <- 
  pp.clean %>%
  relocate(riskmod,interAge,.after = "AgeMastectomy") %>%
  select(-contains("Mastectomy"),-contains("Oophorectomy"))

##### Final Formatting ####
pp.clean <-
  pp.clean %>%
  arrange(PedigreeID,ID)

pp.clean <- as.data.frame(pp.clean)
