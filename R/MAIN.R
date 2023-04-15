#### Set-up ####

# libraries
library(tidyverse)
library(lubridate)
library(stringr)
library(data.table)
library(rmarkdown)

# create dated folders for storing this version's data
source("./R/0_utils.R")
check.and.create.date.dir(dir.name = paste0(interim.master.date.dir))
check.and.create.date.dir(dir.name = paste0(final.date.dir))

#### 1 Split Data ####
# raw data
load(paste0(raw.dir,"US_LATAM_data.RData"))
remove(dat.pro)
remove(dat.rel)

# split the data into demographics and non-demographics and preliminary formatting
source(paste0(root.dir,"1_split-data.R"))

# checkpoint save
save(dat.mod, file = paste0(interim.master.date.dir,"validation-data-for-cleaning.RData"))
save(dat.demographics, file = paste0(interim.master.date.dir,"validation-data-for-cleaning-demographics-cols.RData"))
save(dat.non.demographics, file = paste0(interim.master.date.dir,"validation-data-for-cleaning-non-demographics-cols.RData"))

#### 2 Clean Non-demographics, Non-Test Data ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

# non-demographics data
load(paste0(interim.master.date.dir,"validation-data-for-cleaning-non-demographics-cols.RData"))
datw <- dat.non.demographics
remove(dat.non.demographics)

# demographics data
load(paste0(interim.master.date.dir,"validation-data-for-cleaning-demographics-cols.RData"))

# clean the non-demographics and non-gene testing data
source(paste0(root.dir,"2_non-demographics-non-test-data-cleaning.R"))

# checkpoint save
save(dat.demographics, file = paste0(interim.master.date.dir,"validation-data-for-cleaning-demographics-cols.RData"))
save(datw, file = paste0(interim.master.date.dir,"clean-non-demo-non-test-data.RData"))

#### 3 Clean Test Data ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

# reload data from previous step
load(paste0(interim.master.date.dir,"clean-non-demo-non-test-data.RData"))

# template of possible genes (specific to REDACTED data)
load(paste0(REDACTED.panel.cur.dir,"REDACTED_Panel_Genes_Template.RData"))

# manually categorized panel/lab combos and issues
valid.lab.panels <- read.csv(paste0(ref.cur.dir,"unique_panels_categorized.csv"))
specifics <- read.csv(paste0(ref.cur.dir,"unique_panels_categorized_specific.csv"))

# clean gene testing columns
source(paste0(root.dir,"3_gene-test-data-cleaning.R"))

# checkpoint save
save(datw, file = paste0(interim.master.date.dir,"clean-non-demo-data.RData"))

#### 4 Wrangle Test Data ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

# master data set
load(paste0(interim.master.date.dir,"clean-non-demo-data.RData"))

# template of possible genes (specific to REDACTED data)
load(paste0(REDACTED.panel.cur.dir,"REDACTED_Panel_Genes_Template.RData"))
common.genes <- panel.genes.template$Gene

# REDACTED panel genes and eras tables
load(paste0(REDACTED.panel.cur.dir,"REDACTED_Panel_Eras.RData"))
load(paste0(REDACTED.panel.cur.dir,"REDACTED_Panel_Genes.RData"))

# load panel and gene references from the MGH project
load(paste0(ref.cur.dir,"Panel_Genes.RData"))
load(paste0(ref.cur.dir,"Panel_Eras.RData"))

# functions and keyword object for translating panels to genes
panel.name.keywords <- read.csv(paste0(REDACTED.panel.cur.dir,"panel_name_keywords_with_REDACTED.csv"))
source(paste0(root.dir,"4a_panel-to-gene-mapping-functions.R"))

# wrangle the panels and test results
source(paste0(root.dir,"4_gene-test-wrangling.R"))

# checkpoint save
save(test.dat, file = paste0(interim.master.date.dir,"clean-data-by-test.RData"))
save(datw2, file = paste0(interim.master.date.dir,"clean-wrangled-test-and-non-demo-data.RData"))

#### 5 Clean Demographics Cols ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

# load raw demographics data
load(paste0(interim.master.date.dir,"validation-data-for-cleaning-demographics-cols.RData"))
dat <- dat.demographics
remove(dat.demographics)

# clean demographics
source(paste0(root.dir,"5_demographics-cleaning.R"))

# checkpoint save
save(dat.demographics.clean, file = paste0(interim.master.date.dir,"final-cleaned-demographics-data.RData"))

#### 6 Merge Demographics, Cancer History, Surgeries, and Testing ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

#demographics cols
load(paste0(interim.master.date.dir,"final-cleaned-demographics-data.RData"))

#cancer history
load(paste0(interim.master.date.dir,"clean-non-demo-non-test-data.RData"))

#testing columns
load(paste0(interim.master.date.dir,"clean-wrangled-test-and-non-demo-data.RData"))

# clean all genes version and split into two different formats
source(paste0(root.dir,"6_merge.R"))

# Final Save for All Genes Variant of Data Set
save(dat.clean, file = paste0(final.date.dir,"clean-REDACTED-data-all-genes.RData"))

#### 7 Final Formatting: PanelPRO and BayesMendel ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

# clean master data set
load(paste0(final.date.dir,"clean-REDACTED-data-all-genes.RData"))

# test dates
load(paste0(interim.master.date.dir,"clean-non-demo-data.RData"))

# split into BayesMendel and PanelPRO formats
source(paste0(root.dir,"7_final-formatting.R"))

# Save Final for BayesMendel and PanelPRO variants of the data set
save(bm.clean, file = paste0(final.date.dir,"clean-REDACTED-data-BayesMendel.RData"))
save(pp.clean, file = paste0(final.date.dir,"clean-REDACTED-data-PanelPRO.RData"))

#### Version Documentation ####
# clear the environment and load utilities (includes directories)
remove(list = ls())
source("./R/0_utils.R")

# sessionInfo
si.name <- paste0("sessionInfo-",todays.date,".txt")
si <- capture.output(sessionInfo())
writeLines(si, paste0(final.date.dir, si.name))

# summary stats
input.rmd <- "data-version-summary-stats.Rmd"
output.rmd <- paste0("data-version-summary-stats-",todays.date,".html")
render(input = input.rmd,
       output_dir = final.date.dir,
       output_file = output.rmd,
       output_format = "html_document")

#' reminds user to manually create a README and summary statistics for new data version
reminders <- function(){
  status <- readline(prompt = paste0("Did you manually create a README.md in with version specific changes and an updated summary statistics report using data-version-summary-stats.Rmd (copy the knitted version into ",final.date.dir," (y/n)?"))
  if(status == "y"){
    print("Thank you!")
  } else if(status == "n"){
    print(paste0("Please create these files and put them in ",final.date.dir,"."))
  } else {
    print(paste0("Warning: invalid user response. Please create these files and place them in ",final.date.dir,"."))
  }
}
reminders()
