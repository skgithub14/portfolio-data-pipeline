# Script for cleaning the non-demographics and non-gene test data prior

#### Drop Columns ####
# drop all tumor cols
datw <-
  datw %>%
  select(-contains("IHC"))

#### Prophylactic Surgeries ####
## Mastectomy
# recode values
datw <-
  datw %>%
  mutate(Mastectomy = recode(Mastectomy,'N'='0','Y'='1'),
         Mastectomy = replace_na(Mastectomy, '0'),
         Mastectomy = as.numeric(Mastectomy))

# extract the year of the surgery only
datw$AgeMastectomy <- str_sub(datw$AgeMastectomy,-4,-1)

# remove partial year value
datw$AgeMastectomy[which(datw$AgeMastectomy == '/201')] <- NA

# convert to correct class
datw$AgeMastectomy <- as.numeric(datw$AgeMastectomy)

# convert year of DOB to correct format
datw$year.of.DOB <- as.numeric(datw$year.of.DOB)

# calculate age of surgery
datw <- 
  datw %>%
  
  # fix erroneous year of birth values
  mutate(year.of.DOB = ifelse(year.of.DOB < 1800, NA, year.of.DOB),
         year.of.DOB = ifelse(year.of.DOB > 2021, NA, year.of.DOB),
         
         # calculate age at mastectomy
         AgeMastectomy = AgeMastectomy - year.of.DOB,
         
         # for BayesMendel, the age of surgery should be the current age if no surgery occurred
         AgeMastectomy = ifelse(Mastectomy == 0, CurAge, AgeMastectomy),
         
         # round age to the year
         AgeMastectomy = round(AgeMastectomy, 0)
  )

## Oophorectomy
# recode values and change the format
datw <-
  datw %>%
  mutate(Oophorectomy = recode(Oophorectomy,'y'='1'),
         Oophorectomy = replace_na(Oophorectomy, '0'),
         Oophorectomy = as.numeric(Oophorectomy))

# calculate the age of surgery
datw <-
  datw %>%
  mutate(
         # correct age of surgery if out of range
         AgeOophorectomy = ifelse(AgeOophorectomy < 1, NA, AgeOophorectomy),
         AgeOophorectomy = ifelse(AgeOophorectomy > 120, NA, AgeOophorectomy),
         
         # per BayesMendel documentation, those without the surgery should have the age of surgery set to CurAge
         AgeOophorectomy = ifelse(Oophorectomy == 0, CurAge, AgeOophorectomy),
         
         # round the surgery age to the yearest year
         AgeOophorectomy = round(AgeOophorectomy,0))

#### Marker Tests ####
## ER
# recode and change the column class
datw <-
  datw %>%
  mutate(ER = recode(ER,
                     'Indeterminant'='0',
                     'not done'='0',
                     'Unknown'='0',
                     'neg'='2',
                     'negative'='2',
                     'Negative'='2',
                     'pos'='1',
                     'positive'='1',
                     'Positive'='1',
                     'Y'='1'),
         ER = replace_na(ER, '0'),
         ER = as.numeric(ER))

## PR
# recode and adjust the class
datw <-
  datw %>%
  mutate(PR = recode(PR,
                     'Indeterminant'='0',
                     'not done'='0',
                     'Unknown'='0',
                     'neg'='2',
                     'negative'='2',
                     'Negative'='2',
                     'N'='2',
                     'pos'='1',
                     'positive'='1',
                     'Positive'='1',
                     'Y'='1'),
         PR = replace_na(PR, '0'),
         PR = as.numeric(PR))

## HER2
# recode all column values and adjust the classes
datw <-
  datw %>%
  mutate(across(starts_with("Her2"),
                ~recode(.,
                        '<1999'='0',
                        'Indetermine'='0',
                        'na'='0',
                        'Indeterminant'='0',
                        'Indeterminant (2+)'='0',
                        'not done'='0',
                        'Not done'='0',
                        'Unknown'='0',
                        'Amplified'='0',
                        'neg'='2',
                        'negative'='2',
                        'Negative'='2',
                        'Neagative'='2',
                        'Negative (1+)'='2',
                        'N'='2',
                        'pos'='1',
                        'positive'='1',
                        'Positive'='1',
                        'Poitive'='1',
                        'Positive '='1',
                        'Positive (3+)'='1',
                        'Y'='1'))) %>%
  mutate(across(starts_with("Her2"),
                ~replace_na(., '0'))) %>%
  mutate(across(starts_with("Her2"),
                ~as.numeric(.)))

# merge results and resolve discrepancies
her.df <- 
  datw %>%
  select(starts_with("Her2"))
HER2 <- rep(NA,nrow(her.df))

# check each row
for(row in 1:nrow(her.df)){
  tmp.vec <- as.numeric(her.df[row,])
  
  # condition if there were no tests
  if(sum(tmp.vec) == 0) { 
    val <- 0
    
    # condition if there was at least one high test 
  } else if( sum(tmp.vec == 1) ){
    val <- 1
    
    # condition if there was at least one low/stable test and no high tests
  } else {
    val <- 2
  }
  
  # save merged value in the storage vector
  HER2[row] <- val
}

# add combined column into clean df and drop dirty columns
datw <-
  datw %>%
  add_column(HER2,
             .after = "Her2.Neu.by.unknown.method") %>%
  select(-starts_with("Her2.Neu"))

## MSI
# recode
datw <-
  datw %>%
  mutate(across(contains("MSI"),
                ~as.character(.))) %>%
  mutate(across(contains("MSI"),
                ~recode(.,
                        '1'='0',
                        '11/14/2017'='0',
                        '6/6/2018'='0',
                        'Commercial'='0',
                        'Did not amplify'='0',
                        'Inconclusive'='0',
                        'High'='1',
                        'Not done'='0',
                        'Stable'='2',
                        'Low'='2',
                        'Y'='2'))) %>%
  mutate(across(contains("MSI"),
                ~replace_na(., '0'))) %>%
  mutate(across(contains("MSI"),
                ~as.numeric(.)))

# merge results
msi.df <- 
  datw %>%
  select(contains("MSI"))
MSI <- rep(NA,nrow(msi.df))

# check each row
for(row in 1:nrow(msi.df)){
  tmp.vec <- as.numeric(msi.df[row,])
  
  # condition if there were no tests
  if(sum(tmp.vec) == 0) { 
    val <- 0
    
    # condition if there was at least one high test 
  } else if( sum(tmp.vec == 1) ){
    val <- 1
    
    # condition if there was at least one low/stable test and no high tests
  } else {
    val <- 2
  }
  
  # save val in the storage vector
  MSI[row] <- val
}

# add combined column into clean df and drop dirty columns
datw <-
  datw %>%
  select(-contains("MSI")) %>%
  add_column(MSI,
             .after = "HER2")

#### Cancer History ####
# classify each unique cancer in the raw data into a PanelPRO cancer category
raw.cancers <- list(
  Brain                 = c(
    "Astrocytoma",
    "CNS (Central Nervous Sys",
    "CNS (Central Nervous System)",
    "Glioblastoma",
    "Medulloblastoma",
    "Brain",
    "Glioma",
    "Pituitary"
  ),
  Breast                = c(
    "Breast",
    "Breast\t",      
    "DCIS (Breast: Ductal carcinoma in situ )" ,
    "DCIS (Breast: Ductal carcinoma in situ )\t",
    "DCIS (Breast: Ductal carcinoma in situ)",
    "New Primary DCIS",
    "New Primary Breast"
  ),
  Cervical              = c(
    "Cervical"
  ),
  Colorectal            = c(
    "Colorectal",
    "Cecum"
  ),
  Endometrial           = c(
    "Uterine (Endometrial)\t",
    "Uterine",
    "Uterine (Endometrial)",
    "Choriocarcinoma"
  ),
  Gastric               = c(
    "Gastric (stomach)"
  ),
  Kidney                = c(
    "Renal (Kidney)",
    "Wilms' Tumor",
    "Renal-Clear Cell",
    "Renal-Papillary Type I",
    "Renal-Papillary Type II"
  ),
  Leukemia              = c(
    "ALL (Acute Lymphocytic Leukemia)",
    "AML (Acute Myeloid Leukemia)",
    "Leukemia",
    "CLL (Chronic Lymphocytic Leukemia)",
    "CML (Chronic Myelogenous Leukemia)"
  ),
  Melanoma              = c(
    "Melanoma"
  ),
  Ovarian               = c(
    "Ovarian"
  ),  
  Osteosarcoma          = c(
    "Osteosarcoma"
  ),
  Pancreas              = c(
    "Pancreatic",
    "Pancreatic Neuroendocrin",
    "Pancreatic Neuroendocrine Tumor"
  ),
  Prostate              = c(
    "Prostate",
    "Prostate\t"
  ),
  'Small Intestine'     = c(
    "Duodenum",
    "Small Bowel"
  ),
  'Soft Tissue Sarcoma' = c(
    "Soft tissue sarcoma",
    "Lymph Sarcoma",
    "Hodgkin's Lymphoma",
    "Non Hodgkin's Lymphoma",
    "Liposarcoma",
    "Lymphoma",
    "Lymphoma\t",
    "Angiosarcoma",
    "Fibrosarcoma",
    "Leiomyosarcoma",
    "Neuroblastoma",
    "Paraganglioma",
    "Rhabdomyosarcoma",
    "Schwannoma"
  ),
  Thyroid               = c(
    "Thyroid",
    "Thyroid medullary type",
    "Thyroid papillary type",
    "Thyroid\t",
    "Parathyroid",
    "Thyroid follicular type"
  ),
  'Urinary Bladder'     = c(
    "Bladder",
    "Ureter"
  ),
  Hepatobiliary         = c(
    "Gallbladder",
    "Gallbladder\t",
    "Liver",
    "Liver\t",
    "Bile Duct (cholangiocarc",
    "Bile Duct (cholangiocarcinoma)"
  ),
  Contralateral         = c()
)

# store unclassified cancers for posterity
unclass.cancers <- c(
  "1"                                   ,      "191"                                 ,
  "221"                                 ,      "224"                                 ,
  "236"                                 ,      "238"                                 ,
  "243"                                 ,      "252"                                 ,
  "38"                                  ,      "62"                                  ,
  "63"                                  ,      "740"                                 ,
  "75"                                  ,      "76"                                  ,
  "Abdominal-unknown type"              ,      "Abdominal - unknown type"            ,
  "affected"                            ,      "carrier"                             ,
  "Cancer-Unknown type"                 ,      "Cancer-Unknown type\t"               ,
  "Colon Polyps"                        ,      "Lung"                                ,
  "Lung\t"                              ,      "MEN1"                                ,  
  "NF1"                                 ,      "Cancer-Unk type"                     ,
  "184"                                 ,      "192"                                 ,
  "196"                                 ,      "226"                                 ,
  "255"                                 ,      "60"                                  ,
  "Adrenocortical"                      ,      "Carcinoid"                           ,  
  "Chondrosarcoma"                      ,      "Desmoid"                             ,
  "Esophageal"                          ,      "Fallopian tube"                      ,
  "Head & Neck"                         ,
  "Laryngea9"                           ,      "Laryngeal"                           ,
  "LCIS"                                ,      "LCIS (Lobular carcinoma in situ)"    ,
  "Mesothelioma"                        ,
  "Multiple Myeloma"                    ,      "Neurofibroma"                        ,
  "Peritoneal"                          ,      "Pheochromocytoma"                    ,
  "Polycythemia"                        ,
  "Sebaceous adenoma"                   ,      "Testicular"                          ,
  "Eye"                                 ,      "Eczema"                              ,
  "Gynecological-Unknown ty"            ,      "Gynecological-Unknown type"          ,
  "Mouth (Oral)"                        ,      "Nasopharyngeal"                      ,
  "Nose"                                ,      "Optic Nerve"                         ,
  "Penile"                              ,      "Pharyngeal"                          ,  
  "Throat"                              ,      "Tongue"                              ,
  "Acoustic Neuroma"                    ,      "Acoustic Neuroma\t"                  ,
  "Acoustic  Neuroma"                   ,      "ADH"                                 ,
  "Adrenal"                             ,      "Anal"                                ,
  "Anal\t"                              ,            
  "Aplastic Anemia"                     ,      "Appendix"                            ,  
  "Bone"                                ,      "Bone-Unknown type"                   , 
  "Bone Marrow Failure"                 ,      "CF"                                  ,
  "Chondroblastoma"                     ,      "Cleidocranial Dysostosis"            ,
  "Dupuytren's contracture"             ,      "Familial Adenomatous Pol"            ,
  "Familial Adenomatous Polyposis"      ,      "Fanconi Anemia"                      ,
  "Gardner Syndrome"                    ,      "Pinealblastoma"                      ,
  "GIST"                                ,      "GU-NOS"                              ,
  "Heart"                               ,      "Hemangioma"                          ,
  "Hyperparathyroidism"                 ,      "Keratoacanthoma"                     ,
  "Lip"                                 ,      "Meningioma"                          ,
  "Lipoma"                              ,      "Myelodysplastic Syndrome"            ,
  "Phyllodes"                           ,      "Prolactinoma"                        ,
  "Retinoblastoma"                      ,      "Sarcoma"                             ,
  "Seminoma"                            ,      "Skin - Basal Cell"                   ,
  "Skin - Basal Cell\t"                 ,      "Skin - Squamous Cell"                ,
  "Skin Ca-unk type"                    ,      "Skin Ca-unknown type"                , 
  "Spleen"                              ,      "Thrombocytopenia"                    ,
  "Thymoma"                             ,      "Vaginal"                             ,
  "VHL"                                 ,      "Vulvar"                              ,
  "Xanthosarcoma"                       
)

# retrieve relevant raw data
cancer.hx  <- 
  datw %>% 
  select(SubjectID,contains("Cancer.History.Type")) %>%
  pivot_longer(cols = contains("Cancer.History.Type"), values_to = "cancer") %>%
  select(SubjectID,cancer)
cancer.age <- 
  datw %>% 
  select(SubjectID,contains("Cancer.History.Age")) %>%
  pivot_longer(cols = contains("Cancer.History.Age"), values_to = "age") %>%
  select(SubjectID,age)

# Drop entries where no cancer is listed (NA or blank)
isEmptyCancer = is.na(cancer.hx$cancer) | cancer.hx$cancer == ""
cancer.hx = cancer.hx %>% subset(!isEmptyCancer)
cancer.age = cancer.age %>% subset(!isEmptyCancer)

# Order diagnoses and ages by increasing age, so that earlier diagnoses 
# will get iterated on first
order_cancer_ages = order(cancer.age$age)
cancer.hx = cancer.hx[order_cancer_ages,]
cancer.age = cancer.age[order_cancer_ages,]

# create storage template for final data
pp.cans <- c("BRA",  "BC",   "CER",  "COL",  "ENDO", "GAS",  "KID",  
             "LEUK", "MELA", "OC",   "OST",  "PANC", "PROS", "SMA",  
             "STS",  "THY",  "UB",   "HEP",  "CBC")
pp.cans.long <- c("Brain", "Breast", "Cervical", "Colorectal", "Endometrial", "Gastric", "Kidney",
                  "Leukemia", "Melanoma", "Ovarian", "Osteosarcoma", "Pancreas", "Prostate", "Small Intestine",    
                  "Soft Tissue Sarcoma", "Thyroid", "Urinary Bladder", "Hepatobiliary", "Contralateral")

# data frame of cancer diagnoses
pp.cancer.hx <-  as.data.frame(matrix(
  rep(0,nrow(datw)*length(pp.cans)),
  nrow = nrow(datw),
  ncol = length(pp.cans)))
colnames(pp.cancer.hx)  <- paste0("isAff",pp.cans)

# Add column for second BC
pp.cancer.hx$isAffBC2 = 0

# data frame of cancer ages
pp.cancer.age <- as.data.frame(matrix(
  NA,
  nrow = nrow(datw),
  ncol = length(pp.cans)))
colnames(pp.cancer.age) <- paste0("Age",pp.cans)

# Add column for second BC
pp.cancer.age$AgeBC2 = NA

#loop cancer history by row
progress <- 0
cancer.hx.cnt <- 1:nrow(cancer.hx)
for(cnt in cancer.hx.cnt){
  
  # progress
  pcnt <- floor((cnt / nrow(cancer.hx)) * 100)
  if(pcnt > progress){ print(paste0("Cancer Hx Progress: ",pcnt,"%")) }
  progress <- pcnt
  
  # record this iteration's data points
  tmp.cancer <- cancer.hx$cancer[cnt]
  tmp.subject <- cancer.hx$SubjectID[cnt]
  tmp.age <- cancer.age$age[cnt]
  
  # clean the diagnosis age
  tmp.age <- round(abs(tmp.age),0)
  if(!is.na(tmp.age) & tmp.age > 94 & tmp.age <= 125){ tmp.age <- 94 }
  if(!is.na(tmp.age) & tmp.age > 94){ tmp.age <- NA }
  
  # look-up clean cancer label based on raw cancer entry
  cancer.translator <- names(grep(pattern = 1, lapply(raw.cancers, match, tmp.cancer), value = T))
  
  # check that multiple categories are not being found for a cancer type
  if(length(cancer.translator) > 1){
    print(tmp.cancer)
    print(cancer.translator)
  }
  
  # get the short PanelPRO cancer name
  if(length(cancer.translator) > 0){
    pp.cancer <- pp.cans[which(pp.cans.long == cancer.translator)]
  } else {next} # if a code was read that is not a part of PanelPRO, then next
  
  # assign cancer affected status and age to storage data frames
  pp.cancer.col <- paste0("isAff",pp.cancer)
  pp.age.col <- paste0("Age",pp.cancer)
  
  # Assign cancer status if person was not already diagnosed with cancer in the 
  # past (age should be age of first cancer diagnosis)
  if (pp.cancer.hx [which(datw$SubjectID == tmp.subject),pp.cancer.col] == 0) {
    pp.cancer.hx [which(datw$SubjectID == tmp.subject),pp.cancer.col] <- 1
    pp.cancer.age[which(datw$SubjectID == tmp.subject),pp.age.col   ] <- tmp.age
  } else if (pp.cancer == "BC") {
    # Exception is BC, in which case the second BC should be recorded under BC2 
    # (if not already recorded)
    if (pp.cancer.hx [which(datw$SubjectID == tmp.subject),"isAffBC2"] == 0) {
      pp.cancer.hx [which(datw$SubjectID == tmp.subject),"isAffBC2"] <- 1
      pp.cancer.age[which(datw$SubjectID == tmp.subject),"AgeBC2"   ] <- tmp.age
    }
  }
}

# drop old cancer history cols
datw <-
  datw %>%
  select(-starts_with("Cancer.History"))

# add new cols
datw <- cbind(datw, pp.cancer.hx, pp.cancer.age)

#### Remove Cancer and Surgery Ages > CurAge ####
# remove surgery and cancer ages greater than current age
# subset relevant pedigree and columns
age.diffs <-
  datw %>%
  select(SubjectID,CurAge,starts_with("Age")) %>%
  mutate(across(.cols = starts_with("Age"), ~ .-CurAge)) %>%
  select(-CurAge) %>%
  pivot_longer(cols = starts_with("Age"), names_to = "Event", names_prefix = "Age", values_to = "delta") %>%
  filter(!is.na(delta)) %>%
  filter(delta > 0) %>%
  filter(delta != Inf)

# remove future events
for(row in 1:nrow(age.diffs)){
  tmp.sub <- age.diffs$SubjectID[row]
  tmp.event <- age.diffs$Event[row]
  tmp.event.age <- paste0("Age",tmp.event)
  if(tmp.event != "Mastectomy" & tmp.event != "Oophorectomy"){
    tmp.event <- paste0("isAff",tmp.event)
    age.val <- NA
  } else {
    age.val <- datw$CurAge[which(datw$SubjectID == tmp.sub)]
  }
  datw[which(datw$SubjectID == tmp.sub), tmp.event] <- 0
  datw[which(datw$SubjectID == tmp.sub), tmp.event.age] <- age.val
}

# for any BC removed, ensure the BC2 is also removed
datw[which(datw$isAffBC == 0 & datw$AgeBC2 == 1), "AgeBC2"] <- NA
datw[which(datw$isAffBC == 0 & datw$isAffBC2 == 1), "isAffBC2"] <- 0

##### Reconcile Prophylactic Surgeries with Cancer History ####
# the models uses prophylactic surgeries only and the data set contains non prophylactics
# assume surgery was prophylactic if the person is unaffected by the relevant cancers
# or if affected by the cancer, then only prophylactic if the cancer age and surgery age are present 
# AND the cancer age is greater than the surgery age

# note that a woman may recieve a bilateral prophylactic mastectomy after diagnosis of a 1st BC
# to prevent contralateral breast cancer, however given the dataset does not differentiate between 
# bilateral and unilateral mastectomies we cannot assume mastecomies that occur after a 1st BC are 
# prophylactic and therefore these mastecomies are removed from the data set.
surgeries <- c("Mastectomy","Oophorectomy")
surg.cans <- c("BC","OC")
reconcile.surg.and.cancer <- function(df){
  
  # loop surgeries
  for(surg in surgeries){
    can <- surg.cans[which(surgeries == surg)]
    surg.ages <- df[which(df[,surg] == 1),
                      c("SubjectID",paste0("Age",c(surg,can)),paste0("isAff",can))
                      ]
    colnames(surg.ages)[2:length(colnames(surg.ages))] <- c("AgeSurg","AgeCan","isAffCan")
    
    # find all suspected non-prophylactic surgery subjects
    bad.surg.ages.na.ages <-
      surg.ages %>%
      filter(isAffCan == 1) %>%
      filter(is.na(AgeSurg) | is.na(AgeCan))
    bad.surg.ages.present.ages <-
      surg.ages %>%
      filter(isAffCan == 1) %>%
      filter(!is.na(AgeSurg) & !is.na(AgeCan) & AgeSurg >= AgeCan)
    bad.surg.ages <- rbind(bad.surg.ages.na.ages, bad.surg.ages.present.ages)
    
    # remove non-prophylactic surgeries from master (ages are replaced with CurAge for compatibility with BayesMendel)
    df[which(df$SubjectID %in% bad.surg.ages$SubjectID),surg] <- 0
    df[which(df$SubjectID %in% bad.surg.ages$SubjectID),paste0("Age",surg)] <- df$CurAge[which(df$SubjectID %in% bad.surg.ages$SubjectID)]
  }
  df
}

datw <- reconcile.surg.and.cancer(datw)
