# this script takes the original data set and splits it into demographics
# and non demographics subsets for parallel work

#### Remove Duplicate Records ####
dups <- names(table(dat$DatID)[which(table(dat$DatID) > 1)])
dup.df <- 
  dat[which(dat$DatID %in% dups),] %>%
  arrange(DatID) %>%
  group_by(DatID) %>%
  mutate(IDNum = row_number(), .before = DatID)

# loop through duplicate IDs
for(id in dups){
  tmp.df <- 
    dup.df %>%
    filter(DatID == id)
  conditions <- c("diff")
  values <- c(3)
  names(values) <- "IDNum"
  
  # inner loop to check column data
  for(col in colnames(tmp.df)[2:ncol(tmp.df)]){
    
    # if either value is NA
    if(is.na(tmp.df[1,col]) | is.na(tmp.df[2,col])){
      
      # if both values are NA
      if(is.na(tmp.df[1,col]) & is.na(tmp.df[2,col])){
        cond <- "same"
        val <- tmp.df[1,col]
        
        # if any value is not NA
      } else {
        cond <- "diff"
        
        # if value 1 is NA, use value 2
        if(is.na(tmp.df[1,col])){
          val <- tmp.df[2,col]
          
          # otherwise use value 1
        } else {
          val <- tmp.df[1,col]
        }
      }
      
      # if neither value is NA
    } else {
      
      # if the values are equal
      if(tmp.df[1,col] == tmp.df[2,col]){
        cond <- "same"
        val <- tmp.df[1,col]
        
        # if the values are different, review them manually
      } else {
        cond <- "diff"
        val <- "REVIEW"
      }
    }
    conditions <- c(conditions,cond)
    values <- c(values,val)
  }
  
  # append to data frame of dupicate results
  dup.df <- rbind(dup.df,as.data.frame(values))
}

# remove duplicates from the master
dat <- 
  dat %>%
  filter(!(DatID %in% dups)) %>%
  bind_rows(dup.df %>% filter(IDNum == 3))

#### Rename Columns ####
dat.mod <-
  dat %>%
  select(
    -c(
      PanelPro,
      X_affected.status,
      starts_with("BR.Location"),
      Age.at.menarche:Highest.level.of.education.completed,
      starts_with("FDR"),
      starts_with("SDR"),
      starts_with("TDR"))) %>%
  rename("SubjectID" = DatID,
         "PedigreeID" = Pedigree.name,
         "isHispanic" = Ethnicity_BITS,
         "Region" = region,
         "ID" = UPN,
         "FatherID" = Father.ID,
         "MotherID" = Mother.ID,
         "isProband" = Relationship.to.Proband,
         "CurAge" = Age.at.entrance,
         "isDead" = Deceased.status,
         "AgeMastectomy" = Date.of.Mastectomy,
         "Oophorectomy" = CR_oophorectomy,
         "AgeOophorectomy" = Computed.age.of.oophorectomy) %>%
  relocate(Region, .after = PedigreeID)

#### Fix Ages ####
# Take abs() of all ages, round to nearest integer, cap at 94 for ages from 95 to 125, make NA for ages above 125.  
dat.mod$CurAge <- abs(dat.mod$CurAge)
dat.mod$CurAge <- round(dat.mod$CurAge)
dat.mod$CurAge[which(dat.mod$CurAge > 94 & dat.mod$CurAge <= 125)] <- 94
dat.mod$CurAge[which(dat.mod$CurAge > 94)] <- NA


#### De-identify ####
# remove 2 bad/duplicate pedigrees
dat.mod <- dat.mod[which(dat.mod$PedigreeID != "REDACTED VALUE HERE"),]
dat.mod <- dat.mod[which(dat.mod$PedigreeID != "REDACTED VALUE HERE"),]

# check pedigree IDs for characters
char.peds <- grep(pattern = "[a-zA-Z]", unique(dat.mod$PedigreeID), value = T)

# reduce to only identifiable PedigreeIDs
char.peds <- char.peds[!char.peds %in% c("REDACTED VALUE HERE","REDACTED VALUE HERE")]
char.peds <- char.peds[!char.peds %in% grep(pattern = "^REDACTED VALUE HERE", unique(dat.mod$PedigreeID), value = T)]

# assign new ids
new.ids <- seq(1000, 1000 +length(char.peds))

# replace pedigree and subject ids
for(bad.ped in 1:length(char.peds)){
  dat.mod$SubjectID[which(dat.mod$PedigreeID == char.peds[bad.ped])]  <- paste0(new.ids[bad.ped],"_",dat.mod$ID[which(dat.mod$PedigreeID == char.peds[bad.ped])])
  dat.mod$PedigreeID[which(dat.mod$PedigreeID == char.peds[bad.ped])] <- as.character(new.ids[bad.ped])
}

#### Split the Data ####
# the two data sets can be cleaned separately
dat.demographics <-
  dat.mod %>%
  select(SubjectID:year.of.DOB)
dat.non.demographics <-
  dat.mod %>%
  select(SubjectID,PedigreeID,CurAge,year.of.DOB:Cancer.History.Age.at.Diagnosis11)
