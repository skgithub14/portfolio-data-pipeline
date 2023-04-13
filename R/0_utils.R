# key directories
root.dir <- "./R/"
data.dir <- paste0(root.dir,"data/")
raw.dir <- paste0(data.dir,"raw-data/")
interim.master.dir <- paste0(data.dir,"interim-master-data/")
final.dir <- paste0(data.dir,"final-data/")
todays.date <- gsub(pattern = "-", replacement = ".", Sys.Date())
interim.master.date.dir <- paste0(interim.master.dir,todays.date,"/")
final.date.dir <- paste0(final.dir,todays.date,"/")

# reference material directories
REDACTED.panel.dir <- paste0(root.dir,"REDACTED-specific-panels/")
REDACTED.panel.cur.dir <- paste0(REDACTED.panel.dir,"current-version/")
REDACTED.panel.date.dir <- paste0(REDACTED.panel.dir,todays.date,"/")

ref.dir <- paste0(data.dir,"reference-data/")
ref.cur.dir <- paste0(ref.dir,"current-version/")
ref.date.dir <- paste0(ref.dir,todays.date,"/")

#' check if dated directory already exists, if not create it
check.and.create.date.dir <- function(dir.name){
  if(!dir.exists(dir.name)){ dir.create(dir.name) }
}

# cancers and genes with BC2 added to cancers
pp.genes <- c("APC","ATM","BARD1","BMPR1A","BRCA1","BRCA2","BRIP1","CDH1","CDK4","CDKN2A","CHEK2","EPCAM",
              "MLH1","MSH2","MSH6","MUTYH","NBN","PALB2","PMS2","PTEN","RAD51C","RAD51D","STK11","TP53")
pp.cancers <- c("BRA","BC","CER","COL","ENDO","GAS","KID","LEUK","MELA","OC","OST","PANC","PROS","SMA",
                "STS","THY","UB","HEP","CBC","BC2")