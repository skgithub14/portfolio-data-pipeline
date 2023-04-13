# Script for cleaning gene test data prior to panel wrangling

#### Drop Columns ####
# drop extra testing cols
datw <-
  datw %>%
  select(-c(
    contains("specimen"),
    contains("nucleotide"),
    contains("Protein"),
    contains("year.of.DOB")
  ))

#### Test Dates ####
# convert all test dates to the right format
datw <-
  datw %>%
  mutate(across(.cols = contains("Date"), ~as.Date(., "%m/%d/%y")))

# assign easier col names
date.cols.names <-
  datw %>%
  select(contains("Date")) %>%
  colnames()
date.cols.names.new <- sub(pattern = "Test.Result.tr_Date.of.Test.Result", replacement = "Date", date.cols.names)
datw <-
  datw %>%
  rename_with(.cols = contains("Date"), ~ date.cols.names.new[which(date.cols.names == .x)])

#### Test Results ####
# assign easier col names
result.cols.names <-
  datw %>%
  select(matches("Test.Result\\d$|Test.Result\\d\\d$")) %>%
  colnames()
result.cols.names.new <- sub(pattern = "Test.Result.tr_Test.Result", replacement = "Result", result.cols.names)
datw <-
  datw %>%
  rename_with(.cols = matches("Result\\d$|Result\\d\\d$"), ~ result.cols.names.new[which(result.cols.names == .x)])

# recode results
datw <-
  datw %>%
  mutate(across(.cols = matches("Result\\d$|Result\\d\\d$"),
                ~ recode(.,
                         'mutation'='2',
                         'Negative'='0',
                         'no mutation'='0',
                         'not certain'='2',
                         'Other (specify in test comment box)'='2',
                         'Polymorphism'='0',
                         'Positive'='1',
                         'VUS'='2',
                         'VUS favor polymorphism'='2',
                         'VUS suspected deleterious'='2',
                         'VUS, favor polymorphism'='2',
                         'VUS, suspected deleterious'='2'
                )),
         across(.cols = matches("Result\\d$|Result\\d\\d$"),
                ~ na_if(.,
                        "Single Site")),
         across(.cols = matches("Result\\d$|Result\\d\\d$"),
                ~ as.numeric(.)))

# Result Gene
# assign easier col names
gene.cols.names <-
  datw %>%
  select(matches("gene\\d$|gene\\d\\d$")) %>%
  colnames()
gene.cols.names.new <- sub(pattern = "Test.Result.tr_gene", replacement = "Result.Gene", gene.cols.names)
datw <-
  datw %>%
  rename_with(.cols = matches("gene\\d$|gene\\d\\d$"), ~ gene.cols.names.new[which(gene.cols.names == .x)])

# recode incorrect gene names
datw <-
  datw %>%
  
  # correct types and alternative names
  mutate(across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ recode(.,
                         'BRCAZ'='BRCA2',
                         'CHEKZ'='CHEK2',
                         'GALN712'='GALNT12',
                         'GALNTIZ'='GALNT12',
                         'MRE11A'='MRE11',
                         'MYH'='MUTYH',
                         'NFL'='NEFL',
                         'NTHLl'='NTHL1',
                         'p16/CDKN2A'='CDKN2A',
                         'SMARC4'='SMARCA4',
                         'WRV'='WRN')),
         
         # remove unidentifiable/invalid entries
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'other')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'BRCA')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'C.44SA>G')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'C67bdelT')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'CDKN2A/B')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'EX61')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'L147F(560')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'RNFY3')),
         across(matches("Result.Gene\\d$|Result.Gene\\d\\d$"),
                ~ na_if(.,
                        'VH2'))
  )

#### Create Gene columns ####
valid.genes <- panel.genes.template$Gene

# add new gene columns
gene.cols <- as.data.frame(matrix(,
                                  nrow = nrow(datw),
                                  ncol = length(valid.genes)
                                  )
                           )
colnames(gene.cols) <- valid.genes
datw <- cbind(datw,gene.cols)

#### Labs ####
# assign easier col names
lab.cols.names <-
  datw %>%
  select(starts_with("Test.Result.tr_Lab.Name")) %>%
  colnames()
lab.cols.names.new <- sub(pattern = "Test.Result.tr_Lab.Name", replacement = "Lab", lab.cols.names)
datw <-
  datw %>%
  rename_with(.cols = contains("Test.Result.tr_Lab.Name"), ~ lab.cols.names.new[which(lab.cols.names == .x)])

# select all labs and view variations
labs.df <- 
  datw %>%
  select(SubjectID, starts_with("Lab")) %>%
  pivot_longer(cols = starts_with("Lab"), values_to = "Lab", names_to = "Col") %>%
  mutate(across(.cols = starts_with("Lab"),
                ~ tolower(.)))

# previous work identified all spelling variants of the 4 labs that provide genetic panel testing
genedx.variants <- c("gemedx","geme dx","genedxq","genedex","gene-dx","genedy","genecx","gen6dx","gene dx")
ambry.variants <- c("ambry             1","ambrygenetics","ambrygenet1cs","ambryygenetics","ambryr","ampry","ambby","ambrey","amby","amry","ambr","anbry","^mbry","aambry","ambryy"," ambry","1ambry","ambrygenetics")
invitae.variants <- c("inuitae","invitaegenetics","inviate","invitea","invitaae","invitae            -","invitre","invltae","inv1tae","invit","inivtae","invite","invtae","invitae  ","inviatae","^nvitae","lnvitae","1nvitae","invitaeae","invitaee","invitae  ","itae","invinvitae genetics","invinvitae")
myriad.variants <- c("myriad            1","mvriad","mtriad","myraid","myr1ad","myriat","myriadd","myriadmedicalservice","myriadmyriskgenetics","mryaid","myrlad","nyriad","myriad ","myriadmy","myriadgenetic","myriadgenetics","myriadgenetlcs","myriadrisk","myriads")
lab.variants <- list(genedx=genedx.variants,ambry=ambry.variants,invitae=invitae.variants,myriad=myriad.variants)
lab.names <- c("genedx","ambry","invitae","myriad")

# fix lab misspellings
# outer loop to correct all 4 lab misspellings
for(lab in lab.names){
  variant.vec <- lab.variants[[lab]]
  
  # inner loop for each variant of a lab
  for(variant in variant.vec){
    labs.df$Lab <- gsub(pattern = variant, replacement = lab, x = labs.df$Lab)
  }
}

# fix combo name
labs.df$Lab <- gsub(pattern = "myriadambry", replacement = "myriad ambry", x = labs.df$Lab)

# standardize labs
labs.df <-
  labs.df %>%
  mutate(Lab.ambry   = ifelse(grepl(pattern = "ambry"    , Lab), "ambry", NA),
         Lab.invitae = ifelse(grepl(pattern = "invitae"  , Lab), "invitae", NA),
         Lab.genedx  = ifelse(grepl(pattern = "genedx"   , Lab), "genedx", NA),
         Lab.myriad  = ifelse(grepl(pattern = "myriad"   , Lab), "myriad", NA))
labs.df <-
  labs.df %>%
  mutate(Lab = trimws(gsub(pattern = "NA", replacement = "", paste0(Lab.ambry," ",Lab.invitae," ",Lab.genedx," ",Lab.myriad))),
         Lab = ifelse(Lab == "",NA,Lab)) %>%
  select(SubjectID,Col,Lab) %>%
  pivot_wider(id_cols = SubjectID, values_from = Lab, names_from = Col)

# replace dirty columns with clean ones in the master data frame
# order clean data the same way as dirty data
labs.df <- labs.df[order(match(labs.df$SubjectID, datw$SubjectID)),]
for(cnt in 1:12){
  tmp.col <- paste0("Lab",cnt)
  datw[,tmp.col] <- labs.df[,tmp.col]
}

#### Panel Names ####
# assign easier col names
panel.cols.names <-
  datw %>%
  select(contains("Multi")) %>%
  colnames()
panel.cols.names.new <- sub(pattern = "Test.Result.tr_Multi.gene.panel.specify", replacement = "Panel", panel.cols.names)
datw <-
  datw %>%
  rename_with(.cols = contains("Test.Result.tr_Multi.gene.panel.specify"), ~ panel.cols.names.new[which(panel.cols.names == .x)])

# assign easier col names
panel.cols.names <-
  datw %>%
  select(contains("Test.Performed")) %>%
  colnames()
panel.cols.names.new <- sub(pattern = "Test.Result.tr_Test.Performed", replacement = "Panel.Alt1.", panel.cols.names)
datw <-
  datw %>%
  rename_with(.cols = contains("Test.Result.tr_Test.Performed"), ~ panel.cols.names.new[which(panel.cols.names == .x)])

# assign easier col names
panel.cols.names <-
  datw %>%
  select(contains("Single.site")) %>%
  colnames()
panel.cols.names.new <- sub(pattern = "Test.Result.tr_Single.site.gene", replacement = "Panel.Alt2.", panel.cols.names)
datw <-
  datw %>%
  rename_with(.cols = contains("Test.Result.tr_Single.site.gene"), ~ panel.cols.names.new[which(panel.cols.names == .x)])

# combine based on condition
alt.panels.df <- as.data.frame(matrix(,nrow = 0, ncol = 4))
colnames(alt.panels.df) <- c("SubjectID","Panel.Alt1","Panel.Alt2","Col")

# pivot longer, manually
for(group in 1:12){
  tmp.df <- 
    datw %>%
    select(SubjectID, (starts_with("Panel.Alt") & ends_with(paste0(group))))
  if(group == 1 | group == 2){
    tmp.df <-
      tmp.df %>%
      select(-ends_with("11"),-ends_with("12"))
  }
  tmp.df$Col <- rep(group,nrow(tmp.df))
  colnames(tmp.df) <- c("SubjectID","Panel.Alt1","Panel.Alt2","Col")
  alt.panels.df <- rbind(alt.panels.df,tmp.df)
}

# find Panel.Alt1 values that need to be updated and filter
replacement.alt.panels <-
  alt.panels.df %>%
  filter(if_any(.cols = starts_with("Panel.Alt"), ~!is.na(.))) %>%
  mutate(Panel.Alt1 = tolower(Panel.Alt1)) %>%
  filter(Panel.Alt1 %in% c("single site","other","single gene",NA))

# loop through above df
for(row in 1:nrow(replacement.alt.panels)){
  
  # this iteration's data values
  tmp.sub <- replacement.alt.panels[row,"SubjectID"]
  tmp.panel <- replacement.alt.panels[row,"Panel.Alt2"]
  tmp.col <- replacement.alt.panels[row,"Col"]
  
  # replace in the master df
  datw[which(datw$SubjectID == tmp.sub),paste0("Panel.Alt1.",tmp.col)] <- tmp.panel
}

# remove all Panel.Alt2 columns
datw <-
  datw %>%
  select(-starts_with("Panel.Alt2"))

# convert to all lower case and remove invalid values
datw <-
  datw %>%
  mutate(
    across(.cols = starts_with("Panel.Alt1"),
           ~ tolower(.)),
    across(.cols = starts_with("Panel.Alt1"),
           ~ na_if(., "genedx")),
    across(.cols = starts_with("Panel.Alt1"),
           ~ na_if(., "multi-gene panel")),
    across(.cols = starts_with("Panel.Alt1"),
           ~ na_if(., "none"))
  )

# Primary Panel column
# fix typos and standardize abbreviations as full words
breast.ovarian.variants <- c("brlov","br/ov")
ova.variants <- c("0va","dva")
hereditary.variants <- c("herdgpia","here.","hereotary","heriditary","heveditary","inherited")
brca.variants <- c("br ac","brac","bral","brc","broca","rrca")
breast.variants <- c("brea5t","bbeast","bbreas","br ","brea","brease","breas","bc ","br\\-","br\\&")
myrisk.variants <- c("myris","myr1sk","myrlsk","nyrisk")
next.variants <- c("nex ","mext","nex7")
expanded.variants <- c("expandt")
analysis.variants <- c("analyses","analysi")
plus.variants <- c("ptus","plos")
gyn.variants <- c("gny"," and g","byn")
management.variants <- c("mgmt")
gastric.variants <- c("gastr")
cancernext.variants <- c("cancecnext","cancenext","cancernekt","cancernet","canvrtnext")
colaris.variants <- c("co1aris","colarls")
colorectal.variants <- c("coldrectal","crc","colorecta1")
colon.variants <- c("color 30")
endometrial.variants <- c("endomdtrial","endometria1")
high.variants <- c("hiqh")
multi.variants <- c("mulit")
prostate.variants <- c("pprosate")
keywords <- c("ova","hereditary","brca","breast","myrisk","next","expanded","analysis",
              "plus","gyn","management","gastric","breast ovarian","cancernext",
              "colaris","colorectal","colon","endometrial","high","multi","prostate")
keyword.variants <- list(
  'breast ovarian' = breast.ovarian.variants,
  'ova'            = ova.variants,
  'hereditary'     = hereditary.variants,
  'brca'           = brca.variants,
  'breast'         = breast.variants,
  'myrisk'         = myrisk.variants,
  'next'           = next.variants,
  'expanded'       = expanded.variants,
  'analysis'       = analysis.variants,
  'plus'           = plus.variants,
  'gyn'            = gyn.variants,
  'management'     = management.variants,
  'gastric'        = gastric.variants,
  'cancernext'     = cancernext.variants,
  'colaris'        = colaris.variants,
  'colorectal'     = colorectal.variants,
  'colon'          = colon.variants,
  'endometrial'    = endometrial.variants,
  'high'           = high.variants,
  'multi'          = multi.variants,
  'prostate'       = prostate.variants
)
primary.panel.cols <- paste0("Panel",1:12)

# outer loop to go through all 12 panel columns
for(col in primary.panel.cols){
  
  # convert to all lower case
  datw[,col] <- tolower(datw[,col])
  
  # middle loop to go through each keyword
  for(keyword in keywords){
    variant.vec <- keyword.variants[[keyword]]
    
    # inner loop for each variant of a keyword
    for(variant in variant.vec){
      datw[,col] <- gsub(pattern = variant, replacement = keyword, x = datw[,col])
    }
  }
}

# combine based on condition
alt.panels.df <- as.data.frame(matrix(,nrow = 0, ncol = 4))
colnames(alt.panels.df) <- c("SubjectID","Panel.Alt1","Panel","Col")

# pivot longer, manually
for(group in 1:12){
  tmp.df <-
    datw %>%
    select(SubjectID, (starts_with("Panel") & ends_with(paste0(group))))
  if(group == 1 | group == 2){
    tmp.df <-
      tmp.df %>%
      select(-ends_with("11"),-ends_with("12"))
  }
  tmp.df$Col <- rep(group,nrow(tmp.df))
  colnames(tmp.df) <- c("SubjectID","Panel.Alt1","Panel","Col")
  alt.panels.df <- rbind(alt.panels.df,tmp.df)
}

# remove invalid primary panel entries
alt.panels.df <-
  alt.panels.df %>%
  mutate( Panel = na_if(Panel, "single gene"),
          Panel = na_if(Panel, "other"),
          Panel = na_if(Panel, "ambry"),
          Panel = na_if(Panel, "carriers panel "),
          Panel = na_if(Panel, "custom"))

# remove alternate entries if there is a conflict
alt.panels.df <-
  alt.panels.df %>%
  mutate( Panel.Alt1 = ifelse(!is.na(Panel) & !is.na(Panel.Alt1), NA, Panel.Alt1))

# combine primary and alternate columns
alt.panels.df <-
  alt.panels.df %>%
  mutate( Panel = ifelse(is.na(Panel), Panel.Alt1, Panel) )

# remove na values for primary panel and remove alternate panel column
alt.panels.df <- 
  alt.panels.df %>%
  filter(!is.na(Panel) & !is.na(Panel.Alt1)) %>%
  select(-Panel.Alt1)

# remove alternate panel column in master data frame
datw <-
  datw %>%
  select(-starts_with("Panel.Alt"))

for(row in 1:nrow(alt.panels.df)){
  tmp.subject <- alt.panels.df[row,"SubjectID"]
  tmp.col     <- alt.panels.df[row,"Col"]
  tmp.panel   <- alt.panels.df[row,"Panel"]
  datw[which(datw$SubjectID == tmp.subject), paste0("Panel",tmp.col)] <- tmp.panel
}

#### Manual Panel and Lab Validation ####
# #the below code is commented out because it was used to generate a .csv which 
# #was manually adjusted and read back in
# # master list of all unique lab and panel combos
# labs.panels <- as.data.frame(matrix(,nrow = 0,ncol = 4))
# colnames(labs.panels) <- c("SubjectID","Lab","Panel","Col")
# 
# for(set in 1:12){
#   
#   tmp.df <- 
#     datw %>%
#     select(SubjectID,matches(paste0("Lab",set,"$")),matches(paste0("Panel",set,"$")))
#   
#   tmp.df$Col <- rep(set,nrow(tmp.df))
#   
#   colnames(tmp.df) <- colnames(labs.panels)
#   
#   labs.panels <- rbind(labs.panels,tmp.df)
#   
# }
# 
# labs.panels.u <-
#   labs.panels %>%
#   group_by(Lab,Panel) %>%
#   summarise(Count = n())
# 
# labs.panels.u
# 
# # export as a .csv for editing
# write.csv(labs.panels.u, "./R/data/interim-master-data/unique_panels.csv", row.names = F)
# 
# read back-in manually categorized panel/lab combos
valid.lab.panels <- read.csv(paste0("./R/data/reference-data/current-version/unique_panels_categorized.csv"))
valid.lab.panels <- valid.lab.panels %>% select(Lab,Panel,Count,Category)
bad.labs <- 
  valid.lab.panels %>%
  filter(Category == "wrong lab") %>%
  mutate(Combined = paste0(Lab,"_",Panel))

# loop through each of the 12 tests
for(col.num in 1:12){
  
  # subset relevant data and format df
  panel.col <- paste0("Panel",col.num)
  lab.col <- paste0("Lab",col.num)
  tmp.df <- datw[which(!is.na(datw[,panel.col])), c("SubjectID",lab.col,panel.col)]
  tmp.df$Col <- rep(col.num,nrow(tmp.df))
  colnames(tmp.df)[2:3] <- c("Lab","Panel")
  
  # concat lab and panel
  tmp.df <-
    tmp.df %>%
    mutate(Combined = paste0(Lab,"_",Panel))
  
  # progress
  print(paste0("Starting Column: ",col.num))
  print(paste0("Initial Lab Entries: ", sum(!is.na(datw[,lab.col]))))
  
  # check which panels match the bad panels
  tmp.df <-
    tmp.df %>%
    mutate(Flag = ifelse(Combined %in% bad.labs$Combined == T, 1, 0)) %>%
    filter(Flag == 1)
  
  # remove the lab entry if the lab was bad
  datw[which(datw$SubjectID %in% tmp.df$SubjectID),lab.col] <- NA
  
  # checks
  print(paste0("Final Lab Entries: ", sum(!is.na(datw[,lab.col]))))
}
# additional manipulation was manually done the csv and it was read back-in

# read in manually validated .csv
specifics <-
  specifics %>%
  select(Lab,Panel,Specific,Category) %>%
  filter(Category == "specific") %>%
  select(-Category) %>%
  mutate(Specific = na_if(Specific, "")) %>%
  mutate(Combined = paste0(Lab,"_",Panel))

# loop through each of the 12 tests
for(col.num in 1:12){
  
  # subset relevant data and format df
  panel.col <- paste0("Panel",col.num)
  lab.col <- paste0("Lab",col.num)
  tmp.df <- datw[which(!is.na(datw[,panel.col])), c("SubjectID",lab.col,panel.col)]
  tmp.df$Col <- rep(col.num,nrow(tmp.df))
  colnames(tmp.df)[2:3] <- c("Lab","Panel")
  
  # concat lab and panel
  tmp.df <-
    tmp.df %>%
    mutate(Combined = paste0(Lab,"_",Panel))
  
  # progress
  print(paste0("Starting Column: ",col.num))
  
  # check which panels match specific list
  tmp.df <-
    tmp.df %>%
    mutate(Flag = ifelse(Combined %in% specifics$Combined == T, 1, 0)) %>%
    filter(Flag == 1)
  
  # create the specific flag column in the master data frame
  datw[,paste0("Specific",col.num)] <- rep(0,nrow(datw))
  
  # populate the specific flag column based on condition
  datw[which(datw$SubjectID %in% tmp.df$SubjectID),paste0("Specific",col.num)] <- 1
  
  # check
  print(paste0("Specific Count: ",sum(datw[,paste0("Specific",col.num)]))
}

# relocate Specific Columns to between lab and panel
datw <-
  datw %>%
  relocate(Specific1, .after = "Lab1") %>%
  relocate(Specific2, .after = "Lab2") %>%
  relocate(Specific3, .after = "Lab3") %>%
  relocate(Specific4, .after = "Lab4") %>%
  relocate(Specific5, .after = "Lab5") %>%
  relocate(Specific6, .after = "Lab6") %>%
  relocate(Specific7, .after = "Lab7") %>%
  relocate(Specific8, .after = "Lab8") %>%
  relocate(Specific9, .after = "Lab9") %>%
  relocate(Specific10, .after = "Lab10") %>%
  relocate(Specific11, .after = "Lab11") %>%
  relocate(Specific12, .after = "Lab12")

#### Clean Specific Genes Tested ####
# replace dirty specific panel entries with clean ones
# loop through each of the 12 tests
for(col.num in 1:12){
  
  # subset relevant data and format df
  panel.col <- paste0("Panel",col.num)
  lab.col <- paste0("Lab",col.num)
  tmp.df <- datw[which(!is.na(datw[,panel.col])), c("SubjectID",lab.col,panel.col)]
  tmp.df$Col <- rep(col.num,nrow(tmp.df))
  colnames(tmp.df)[2:3] <- c("Lab","Panel")
  
  # concat lab and panel
  tmp.df <-
    tmp.df %>%
    mutate(Combined = paste0(Lab,"_",Panel))
  
  # progress
  print(paste0("Starting Column: ",col.num))
  
  # check which panels match specific list
  tmp.df <-
    tmp.df %>%
    mutate(Flag = ifelse(Combined %in% specifics$Combined == T, 1, 0)) %>%
    filter(Flag == 1)
  
  # replace dirty panel values with clean ones by Subject
  for(row in 1:nrow(specifics)){
    this.match <- specifics[row,"Combined"]
    this.clean <- specifics[row,"Specific"]
    
    # check if match is present in this column
    if(sum(tmp.df$Combined %in% this.match) > 0){
      
      # get relevant subjects
      rel.subs <- tmp.df$SubjectID[which(tmp.df$Combined %in% this.match)]
      
      # replace
      datw[which(datw$SubjectID %in% rel.subs), paste0("Panel",col.num)] <- this.clean
    }
  }
}