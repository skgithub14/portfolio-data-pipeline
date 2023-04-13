# This script combines the demographics and non-demographics data
# it also cleans and imputes the sexes and creates the test order column

#### Missing Person ####
## -----------------------------------------------------------------------------------------------------------------------------------------------------
# load a missing person
# manually type in the row information
missing.subject <- c("REDACTED VALUE", # SubjectID
                     "REDACTED VALUE", # CurAge
                     "REDACTED VALUE", # PedigreeID
                     "REDACTED VALUE", # Region
                     "REDACTED VALUE", # ID
                     "REDACTED VALUE", # FatherID
                     "REDACTED VALUE", # MotherID
                     "REDACTED VALUE", # Gender
                     "REDACTED VALUE", # isProband
                     "REDACTED VALUE", # isDead
                     "REDACTED VALUE", # isHispanic
                     "REDACTED VALUE", # Twins
                     "REDACTED VALUE", # ethnic
                     "REDACTED VALUE", # Death
                     "REDACTED VALUE", # AgeDeath
                     "REDACTED VALUE", # race_bm
                     "REDACTED VALUE"  # race_pp
)
missing.subject <- as.data.frame(t(missing.subject))
colnames(missing.subject) <- colnames(dat.demographics.clean)

# add the missing person back into the clean demographics data
dat.demographics.clean <- rbind(dat.demographics.clean,missing.subject)

# free up RAM
remove(missing.subject)

#### Join Cancer History and Demographics ####
# combine cancer history, surgeries and demographics
datw <- 
  datw %>%
  select(SubjectID, Mastectomy, Oophorectomy, starts_with("isAff"), starts_with("Age")) %>%
  relocate(AgeMastectomy, AgeOophorectomy, .after = "Oophorectomy")
setDT(datw)
setDT(dat.demographics.clean)
datw <- dat.demographics.clean[datw, on = "SubjectID"]
remove(dat.demographics.clean)

#### Impute/Clean Sex ####
## -----------------------------------------------------------------------------------------------------------------------------------------------------
# if a female cancer or prophylactic surgery is found, assign the sex as female
datw <-
  datw %>%
  mutate(Gender = ifelse(is.na(Gender) & (isAffOC == 1 | isAffENDO == 1), 0, Gender)) %>%
  mutate(Gender = ifelse(is.na(Gender) & (Oophorectomy == 1 | Mastectomy == 1), 0, Gender))

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# iterate through each family with a missing gender
missing.sex <-
  datw %>%
  filter(is.na(Gender))

# check if they had a child with someone from opposite sex
for(person in 1:nrow(missing.sex)){
  if(person%%1000 == 0){ print(paste0("Imputing Gender: ",person," of ",nrow(missing.sex))) }
  tmp.sub <- as.character(missing.sex[person,"SubjectID"])
  tmp.id <- as.numeric(missing.sex[person,"ID"])
  tmp.fam <- as.character(missing.sex[person,"PedigreeID"])
  fam.moms <- as.numeric(datw$MotherID[which(datw$PedigreeID == tmp.fam)])
  fam.dads <- as.numeric(datw$FatherID[which(datw$PedigreeID == tmp.fam)])
  mom.flag <- sum(fam.moms == tmp.id, na.rm = T)
  dad.flag <- sum(fam.dads == tmp.id, na.rm = T)
  
  # if it was found in one column or the other, assign the appropriate sex
  if(dad.flag > 0 & mom.flag == 0){
    tmp.sex <- 1
  } else if(dad.flag == 0 & mom.flag > 0){
    tmp.sex <- 0
  } else {
    tmp.sex <- NA
  }
  
  # store the result
  datw[which(datw$SubjectID == tmp.sub),"Gender"] <- tmp.sex
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------
#drop un-imputable genders 
datw <-
  datw %>%
  filter(!is.na(Gender))

##-------------------------------------------------------------------------------
# check males are not assigned female cancers or prophylactic surgeries
bad.sexes <-
  datw %>%
  select(SubjectID,PedigreeID,ID,Gender,isAffOC,isAffENDO,Mastectomy,Oophorectomy) %>%
  filter(Gender == 1) %>%
  filter(isAffOC == 1 | isAffENDO == 1 | Oophorectomy == 1 | Mastectomy == 1)

# check if this person is a mother or father
for(row in 1:nrow(bad.sexes)){
  tmp.sub <- as.character(bad.sexes[row,"SubjectID"])
  tmp.id <- as.numeric(bad.sexes[row,"ID"])
  tmp.dads <- as.numeric(datw$FatherID[which(datw$PedigreeID == as.character(bad.sexes[row,"PedigreeID"]))])
  tmp.moms <- as.numeric(datw$MotherID[which(datw$PedigreeID == as.character(bad.sexes[row,"PedigreeID"]))])
  is.father <- sum(tmp.dads == tmp.id, na.rm = T)
  is.mother <- sum(tmp.moms == tmp.id, na.rm = T)
  
  # if they are a father (and not a mother as a check), then change their female cancers and surgery statuses to 0
  if(is.father > 0 & is.mother == 0){
    datw[which(datw$SubjectID == tmp.sub),c("isAffOC","isAffENDO","Mastectomy","Oophorectomy")] <- 0
    
    # also update surgery ages to CurAge to match existing data (for compatibility with BayesMendel)
    datw[which(datw$SubjectID == tmp.sub),c("AgeMastectomy","AgeOophorectomy")] <- datw$CurAge[which(datw$SubjectID == tmp.sub)]
    
    # if they are not a father, change the sex to female
  } else if(is.father == 0 & is.mother >= 0){
    datw[which(datw$SubjectID == tmp.sub),"Gender"] <- 0
    
    # if they are labelled as both a mother and a father, remove them
  } else {
    datw <- datw[which(datw$SubjectID != tmp.sub),]
  }
}

# check sexes of mothers and fathers
# change sex to match mother or father assignment

# create unique mother/father ids for the dataset
mf.sexes <-
  datw %>%
  select(SubjectID,PedigreeID,Gender,MotherID,FatherID) %>%
  mutate(u_MotherID = paste0(PedigreeID,"_",MotherID),
         u_FatherID = paste0(PedigreeID,"_",FatherID))
u.mothers <- unique(mf.sexes$u_MotherID)
u.fathers <- unique(mf.sexes$u_FatherID)
mother.sexes <- datw[which(datw$SubjectID %in% u.mothers),c("SubjectID","Gender")]
father.sexes <- datw[which(datw$SubjectID %in% u.fathers),c("SubjectID","Gender")]
wrong.mom.sexes <- mother.sexes %>% filter(Gender != 0)
wrong.dad.sexes <- father.sexes %>% filter(Gender != 1)
datw[which(datw$SubjectID %in% wrong.mom.sexes$SubjectID),"Gender"] <- 0
datw[which(datw$SubjectID %in% wrong.dad.sexes$SubjectID),"Gender"] <- 1

#### Join in Test Data ####
# Note this can only be done in a high RAM environment (at least 6.8GB available)
datw2 <- datw2 %>% select(-(CurAge:AgeOophorectomy))
setDT(datw2)
dat.clean <- merge(datw, datw2, by = c("SubjectID","PedigreeID"))
remove(datw)
remove(datw2)

#### Add Carrier/Non-Carrier Identifier ####
# only the all genes version of the dataset can find the true non-carrier/carriers for all genes
pos.carriers <-
  dat.clean %>%
  filter(if_any(.cols = ABL1:ZRSR2, ~ !is.na(.) & .==1)) %>%
  select(SubjectID) %>%
  as.matrix() %>%
  as.vector()
vus.carriers <-
  dat.clean %>%
  filter(!SubjectID %in% pos.carriers) %>%
  filter(if_any(.cols = ABL1:ZRSR2, ~ !is.na(.) & .==2)) %>%
  select(SubjectID) %>%
  as.matrix() %>%
  as.vector()
untested <-
  dat.clean %>%
  filter(if_all(.cols = ABL1:ZRSR2, ~ is.na(.))) %>%
  select(SubjectID) %>%
  as.matrix() %>%
  as.vector()
dat.clean$CarrierStatus <- rep(0,nrow(dat.clean))
dat.clean <- dat.clean %>% relocate(CarrierStatus, .after = "Twins")
dat.clean$CarrierStatus[which(dat.clean$SubjectID %in% pos.carriers)] <- 1
dat.clean$CarrierStatus[which(dat.clean$SubjectID %in% vus.carriers)] <- 2
dat.clean$CarrierStatus[which(dat.clean$SubjectID %in% untested)] <- NA

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# convert columns to correct class
dat.clean <-
  dat.clean %>%
  mutate(across(.cols = c(isProband,ID,Gender,MotherID,FatherID,CurAge,isDead,
                          Death,AgeDeath,isHispanic,Twins,CarrierStatus,Mastectomy,Oophorectomy,
                          starts_with("Age"),starts_with("isAff"),ABL1:ZRSR2),
                ~ as.numeric(.)))

#### Remove Select Marker Tests ####
# marker tests associated with breast and colorectal cancer have to be removed if
# the individual does not have those cancers
bc.del <-
  dat.clean %>%
  filter(if_any(.cols = c(ER,PR,HER2), ~ .!=0)) %>%
  filter(isAffBC == 0) %>%
  select(SubjectID) %>%
  as.matrix() %>%
  as.vector()
dat.clean[which(dat.clean$SubjectID %in% bc.del),c("ER","PR","HER2")] <- 0
col.del <-
  dat.clean %>%
  filter(MSI != 0 & isAffCOL == 0) %>%
  select(SubjectID) %>%
  as.matrix() %>%
  as.vector()
dat.clean[which(dat.clean$SubjectID %in% col.del),"MSI"] <- 0
