# This script contains the functions necessary to read a panel column entry from
# the COH data set to identify the set of genes tested.

#-------------------------------------------------------------------------------
# many missing dates exists therefore need to assume a date if missing
# will assume each missing date is the first era of the panel
impute.missing.panel.date <- function(panel){
  #' imputes the missing date of the panel by assuming the panel date
  #' is the date of the first era of the panel
  #' @param panel.  character string of name of the panel (all lower case)
  #' returns a date in YYYY-MM-DD format
  
  # retrieve the matching Eras table
  era.tbl.name <- paste0(panel,".Eras")
  tmp.era.tbl <- get(era.tbl.name)
  
  # look up the 1st Era start date
  date <- tmp.era.tbl[1,"Start.date"]
  
  # ensure the imputed date is at least 1990-01-01
  if(as.Date(date) < as.Date("1990-01-01")){
    
    date <- "1990-01-01"
    
  }
  date
}

#-------------------------------------------------------------------------------
# function to find all genes tested in a panel based on date
# 2 inputs:
# panel: character string of the panel name.  It must match panel names that are 
#        found in the data files Panel_Eras.RData and Panel_Genes.RData without 
#        the .Eras or .Genes suffixes on the data frames in those files
# date: character string in the format YYYY/MM/DD.  Date of the panel.
# Outputs: returns a 1 column data frame of the genes that were tested
# Errors: are produced if the test date is before the panel was released and
#         if the test date was after the panel was retired.
find.genes.in.panel <- function(panel, date) {
  era.df <- get(paste0(panel,".Eras"))
  gene.df.all <- get(paste0(panel,".Genes"))
  
  # impute the date if missing
  if(is.na(date)){
    date <- as.Date(impute.missing.panel.date(panel = panel))
  }
  
  # warn if test date before 1990/01/01
  if(date < "1990/01/01") {
    print("Error: test date is before 1990/01/01 and should be verified")
    print(paste0("Panel: ",panel,", Date: ",date))
    gene.df <- "Before 1990/01/01"
  } else {
    era.match <- sum(date >= era.df[,"Start.date"])
    
    # warn if test date before test release date
    if(era.match < 1) {
      print("Error: test date is before the panel was released")
      print(paste0("Panel: ",panel,", Date: ",date))
      gene.df <- "Pre-release"
    } else {
      
      # subset to relevant era only
      era.name <- paste0("Era",era.match)
      gene.df <- gene.df.all[,c("Gene",era.name)]
      
      # check if test was retired on the given date
      total.genes.tested <- sum(gene.df[,era.name])
      if(total.genes.tested == 0) {
        print("Error: panel retired at date of the test")
        print(paste0("Panel: ",panel,", Date: ",date))
        gene.df <- "Retired"
      } else {
        
        # check if the list of genes tested was in a month that a change of tested genes occurred 
        # change dates were only MM/YYYY so all dates were assumed to be the 1st when era tables were created
        year <- lubridate::year(date)
        era.year <- lubridate::year(era.df[era.match,"Start.date"])
        month <- lubridate::month(date)
        era.month <- lubridate::month(era.df[era.match,"Start.date"])
        if(year == era.year & month == era.month & era.match > 1){
          
          print("Warning: Test date is during the month the list of genes in this panel was changed. Only genes common in both lists will be evaluated.")
          print(paste0("Panel: ",panel,", Date: ",date))
          
          # compare the gene list from the two eras
          past.era <- era.match - 1
          past.era.name <- paste0("Era",past.era)
          two.era.gene.df <- gene.df.all[,c("Gene",past.era.name,era.name)]
          two.era.gene.df$Add <- gene.df.all[,2] + gene.df.all[,3]
          
          # only keep the genes that are common between the two eras
          gene.df <- subset(two.era.gene.df, two.era.gene.df[,"Add"] == 2, select = Gene)
          
        } else { # if there was no special case the do a simple subset to get the list of genes tested
          gene.df <- subset(gene.df, gene.df[,era.name] == 1, select = Gene)
        }
      }
    }
  }
  gene.df
}


#------------------------------------------------------------------------------
# save keywords table as .RData
# panel.name.keywords <- read.csv("C:/Users/steph/Dropbox (Partners HealthCare)/CCGCRN Hispanic Cohort Data/Data Cleaning/COH Specific Panels/current-version/panel_name_keywords_with_coh.csv")
# save(panel.name.keywords, file = "C:/Users/steph/Dropbox (Partners HealthCare)/CCGCRN Hispanic Cohort Data/Data Cleaning/COH Specific Panels/current-version/panel_name_keywords.RData")


#-------------------------------------------------------------------------------
upgrade.panel.based.on.date <- function(panel, date){
  #' Checks for test date and panel name agreement and if the panel was retired
  #' at the date of the test then the current version of the panel is returned.
  #' If the panel was valid at the date of the test then it just returns the 
  #' original panel name.
  #' @param panel character. Panel name (all lower case).
  #' @param date character. Date in the format YYYY-MM-DD.
  #' returns a validated panel name (either the original if it was valid, or
  #' a the upgraded panel name)
  
  # upgrade translator
  upgrades <- list(
    brca.plus.ambry                    = "brca.plus.expanded.ambry",
    breast.cancer.management.genedx    = "breast.cancer.high.moderate.risk.genedx",
    breast.gyn.cancer.genedx           = c("breast.ovarian.cancer.genedx","endometrial.cancer.genedx"),
    comprehensive.common.cancer.genedx = "comprehensive.cancer.genedx",
    common.cancer.management.genedx    = "high.moderate.risk.cancer.genedx"
  )
  
  # check the panel/date combo
  check <- find.genes.in.panel(panel = panel, date = date)
  if(is.data.frame(check)){
    
    # if a data frame was found then it's a valid panel/date combo
    return(panel)
    
    # if "retired" was found in the error message then it needs to be upgraded
  } else if(grepl(pattern = "Retired", check)){
    
    new.panel <- names(grep(pattern = 1, lapply(upgrades, match, panel), value = T))
    return(new.panel)
    
    # unconsidered case occurred, so return a message
  } else { # date is bad
    return(NA)
  }
}

#-------------------------------------------------------------------------------
find.panels.that.match <- function(panel.keywords, lab.name, test.date){
  #' finds all relevant panels that match a vector of keywords
  #' 
  #' @param panel.keywords character.  character vector of keywords that must match 
  #'  the column names of "panel_name_keywords.RData" (excluding the first column)
  #' @param test.date character.  test date string in form of YYYY-MM-DD
  #' @param gene.quant optional integer.  number of genes in the panel, default is NULL
  #'
  #' returns character vector of panel(s) that matched the most keywords
  
  # if lab name provided, add it to the keywords
  if(!is.na(lab.name)){
    panel.keywords <- c(panel.keywords,lab.name)
  }
  
  # special cases
  prelim.flag     <- sum(grepl(pattern = "prelim",     panel.keywords)) > 0
  manage.flag     <- sum(grepl(pattern = "manage",     panel.keywords)) > 0
  breast.flag     <- sum(grepl(pattern = "breast",     panel.keywords)) > 0
  colaris.flag    <- sum(grepl(pattern = "colaris",    panel.keywords)) > 0
  colaris.ap.flag <- sum(grepl(pattern = "colaris.ap", panel.keywords)) > 0
  high.flag       <- sum(grepl(pattern = "high",       panel.keywords)) > 0
  mod.flag        <- sum(grepl(pattern = "mod",        panel.keywords)) > 0
  ova.flag        <- sum(grepl(pattern = "ova",        panel.keywords)) > 0
  colo.flag       <- sum(grepl(pattern = "colo",       panel.keywords)) > 0
  common.flag     <- sum(grepl(pattern = "common",     panel.keywords)) > 0
  lynch.flag      <- sum(grepl(pattern = "lynch",      panel.keywords)) > 0
  next.flag       <- sum(grepl(pattern = "next",       panel.keywords)) > 0
  custom.flag     <- sum(grepl(pattern = "custom",     panel.keywords)) > 0
  brca.analy.flag <- sum(grepl(pattern = "brca.analy", panel.keywords)) > 0
  cancernext.flag <- sum(grepl(pattern = "cancernext", panel.keywords)) > 0
  guide.flag      <- sum(grepl(pattern = "guide"     , panel.keywords)) > 0
  hered.flag      <- sum(grepl(pattern = "hered"     , panel.keywords)) > 0
  my.risk.flag    <- sum(grepl(pattern = "my risk"   , panel.keywords)) > 0
  gyn.flag        <- sum(grepl(pattern = "gyn"       , panel.keywords)) > 0
  canpan.flag     <- sum(grepl(pattern = "canpan",     panel.keywords)) > 0
  n.634.flag      <- sum(grepl(pattern = "634",        panel.keywords)) > 0
  hem.flag        <- sum(grepl(pattern = "hem",        panel.keywords)) > 0
  
  # remove certain keywords for certain cases
  if(breast.flag == F & colo.flag == F)   {panel.keywords <- panel.keywords[!panel.keywords %in% "prelim"]} # only consider prelim if breast or colo
  if(breast.flag == F & common.flag == F) {panel.keywords <- panel.keywords[!panel.keywords %in% "manage"]} # only consider manage if breast or common
  
  # only consider a certain subset of panels if one keyword
  skip.flag <- F
  if(breast.flag == T  & length(panel.keywords) == 1) {
    matching.panels <- c("breast.cancer.invitae","breast.cancer.myriad")
    skip.flag <- T
  } # if breast is only keyword then only consider 2 panels
  if(colo.flag   == T  & length(panel.keywords) == 1) {
    matching.panels <- c("colorectal.cancer.genedx","colorectal.cancer.panel.invitae")
    skip.flag <- T
  } # if colo is only keyword then only consider 2 panels
  if(next.flag   == T  & length(panel.keywords) == 1) {
    matching.panels <- "cancer.next.ambry"
    skip.flag <- T
  } # if next is only keyword then only consider 1
  if(custom.flag == T  & length(panel.keywords) == 1) {
    matching.panels <- c("custom.cancer.genedx","custom.next.cancer.ambry")
    skip.flag <- T
  } # if next is only keyword then only consider 1
  if(brca.analy.flag == T  & length(panel.keywords) == 1) {
    matching.panels <- "brac.analysis.myriad"
    skip.flag <- T
  }
  if(cancernext.flag == T  & next.flag == T & length(panel.keywords) == 2) {
    matching.panels <- "cancer.next.ambry"
    skip.flag <- T
  }
  if(hered.flag  == T  & length(panel.keywords) == 1) {
    matching.panels <- "common.hereditary.cancers.invitae"
    skip.flag <- T
  }
  if(my.risk.flag    == T  & length(panel.keywords) == 1) {
    matching.panels <- "my.risk.myriad"
    skip.flag <- T
  }
  # canpan634 panel gene list is unknown so assume it is canpan700
  if(canpan.flag == T & n.634.flag == T) {
    matching.panels <- "canpan.700.coh"
    skip.flag <- T
  }
  # only consider a certain subset of panels in certain conditions
  if(colaris.flag == T & colaris.ap.flag  == F) {
    matching.panels <- "colaris.myriad"
    skip.flag <- T
  } # if colaris is only keyword then only consider 1
  if(high.flag == T    | mod.flag         == T) {
    if(breast.flag == F & colo.flag == F) {
      matching.panels <- "high.moderate.risk.cancer.genedx"
      skip.flag <- T
    }
  }
  
  # hem panels before 2/2020 are Hem I, otherwise Hem II
  if(hem.flag == T){
    if(is.na(test.date)){
      matching.panels <- "hem.1.coh"
    } else if(test.date < "2020-02-01"){
      matching.panels <- "hem.1.coh"
    } else {
      matching.panels <- "hem.2.coh"
    }
    skip.flag <- T
  }
  
  if(skip.flag == F){
    
    # only look at the columns with relevant keywords
    extra <- rep(FALSE,nrow(panel.name.keywords))
    panel.keywords <- gsub(pattern = "next", replacement = "next.", panel.keywords) # special case where next cannot be a column name
    panel.keywords <- gsub(pattern = "700", replacement = "X700", panel.keywords) # special case where number cannot be a column name
    panel.keywords <- gsub(pattern = "760", replacement = "X760", panel.keywords) # special case where number cannot be a column name
    panel.keywords <- sub(pattern = "cancernext.", replacement = "cancernext", panel.keywords)
    panel.and.keyword.cols <- panel.name.keywords %>% select(c(panel.name,panel.keywords))
    panel.and.keyword.cols$extra <- extra  # this dummy column prevents one keyword cases from turning into a vector instead of data frame
    
    # which panel(s) have the most matches
    cnt.matches <- rowSums(panel.and.keyword.cols[,-1])
    max.indices <- which(cnt.matches==max(cnt.matches))
    matching.panels <- panel.and.keyword.cols$panel.name[max.indices]
    
    # remove panels with diseases not specified in note
    if(ova.flag    == T & lynch.flag == F) {matching.panels <- matching.panels[!matching.panels %in% "hboc.and.lynch.syndrome.myriad"]}
    if(breast.flag == T & gyn.flag   == F) {matching.panels <- matching.panels[!matching.panels %in% grep(pattern = "gyn"  , matching.panels, value = T)]}
    if(breast.flag == T & ova.flag   == F) {matching.panels <- matching.panels[!matching.panels %in% grep(pattern = "ova"  , matching.panels, value = T)]}
    if(guide.flag  == F )                  {matching.panels <- matching.panels[!matching.panels %in% grep(pattern = "guide", matching.panels, value = T)]}
    
    # verify date range is valid and upgrade panels if needed, remove if invalid
    if(!is.na(test.date)){
      
      # loop through each panel to check if the date is valid for it
      for(pan in matching.panels){
        
        # check if the panel is one that is on the possible upgrade list
        upgrade.candidates <- c("brca.plus.expanded.ambry",
                                "breast.cancer.high.moderate.risk.genedx",
                                "breast.ovarian.cancer.genedx",
                                "endometrial.cancer.genedx",
                                "comprehensive.cancer.genedx",
                                "high.moderate.risk.cancer.genedx")
        if(sum(upgrade.candidates == pan) > 0){
          tmp.result <- upgrade.panel.based.on.date(panel = pan, date = test.date)
          
          # if a different panel was returned by the above function, replace it
          if(tmp.result != pan){
            matching.panels <- matching.panels[!matching.panels %in% pan]
            matching.panels <- c(matching.panels, tmp.result)
          } 
          
          # the panel was not on the potential upgrade list
        } else {
          
          # confirm a the date/panel combo is valid, if not remove the panel
          if(!is.null(pan)){
            if(!is.data.frame(find.genes.in.panel(panel = pan, date = test.date))){
              matching.panels <- matching.panels[!matching.panels %in% pan]
            }
          }
        }
      }
    }
  }
  matching.panels
}



#-------------------------------------------------------------------------------
find.common.genes <- function(panel.names.vec, test.date, result.genes, panel.keywords) {
  #' find set of genes common in several panels on a given date
  #' 
  #' @param panel.names.vec character.  character vector of panel names
  #' @param test.date character.  character string that contains a date in format YYYY-MM-DD
  #' @param result.genes string.  text string of the genes in the result
  #' @param panel.keywords character. character vector of panel keywords.  Current purpose
  #' is only to indentify if this vector contains "prelim" to add those panels
  #'
  #' @return character vector of genes
  
  # initialize list of common genes with the template
  common.genes <- panel.genes.template[,1]
  
  # if there is text in result.genes then attempt to extract all valid genes
  result.genes.vec <- NULL
  if(!is.na(result.genes[1])){
    result.genes.split <- str_split(result.genes, pattern = " ")[[1]]
    
    # find all the genes in the result genes entry
    for(gene in result.genes.split){
      if(sum(common.genes %in% gene) == 1){ result.genes.vec <- c(result.genes.vec,gene) }
    }
  }
  
  # if valid genes were found in the result genes string
  if(!is.null(result.genes.vec)){
    
    # loop through each possible panel and narrow down possible genes
    false.flags <- NULL
    for(panel.name in panel.names.vec){
      
      # get the list of genes in this panel
      tmp.panel.genes <- find.genes.in.panel(panel = panel.name, date = test.date)
      # edge cases of find.genes.in.panels produce a character string instead of a data frame, so rule them out
      if(is.data.frame(tmp.panel.genes)) {
        tmp.panel.genes <- tmp.panel.genes[,1]
      } else { next } # if an edge case occurred for this panel then go to the next panel
      
      # append invitae prelim panel genes if needed
      prelim.flag <- sum(grepl(pattern = "prelim", panel.keywords)) > 0
      if(prelim.flag == T){
        prelim.panels <- panel.name.keywords[which(panel.name.keywords$prelim == TRUE),"panel.name"]
        prelim.panel.flag <- sum(grepl(pattern = panel.name, prelim.panels)) > 0
        if(length(prelim.panel.flag) > 0){
          prelim.genes <- find.genes.in.panel(panel = paste0("prelim.",panel.name), date = test.date)
          tmp.panel.genes <- c(tmp.panel.genes,prelim.genes)
        }
      }
      
      # check if all of the result genes are in the panel or not
      for(gene in result.genes.vec){
        
        # stop searching when a gene is found not to be in the panel
        flag <- sum(grepl(pattern = gene, tmp.panel.genes)) > 0
        if(flag == F){ 
          false.flags <- c(false.flags,flag)
          break
        }
      }
      
      # if all the result genes were present then the flag will still be true
      # if that is the case then all the result genes are in the current panel
      # and therefore this list of genes in this panel should be considered
      if(flag == T){ common.genes <- common.genes[common.genes %in% tmp.panel.genes] }
    }
    
    # if none of the possible panels contained the result genes, print warning
    if(!is.null(false.flags)){
      all.false.check <- length(false.flags) == length(panel.names.vec)
      if(all.false.check == T){
        print("Warning: None of the possible panels contained all of the result genes.")
      }
    }
    
  } else { # there was no text in result.genes or there was text but no valid genes
    
    # loop through each possible panel and narrow down possible genes
    for(panel.name in panel.names.vec){
      
      # get the list of genes in this panel
      tmp.panel.genes <- find.genes.in.panel(panel = panel.name, date = test.date)
      # edge cases of find.genes.in.panels produce a character string instead of a data frame, so rule them out
      if(is.data.frame(tmp.panel.genes)) {
        tmp.panel.genes <- tmp.panel.genes[,1]
      } else { next } # if an edge case occurred for this panel then go to the next panel
      common.genes <- common.genes[common.genes %in% tmp.panel.genes]
    }
  }
  as.character(common.genes)
}


#-------------------------------------------------------------------------------
get.specific.genes.tested <- function(panel.entry, result.genes){
  #' checks the panel entry and result genes for individual genes tested
  #' @param panel.entry. string of gene
  #' returns a character vector of genes tested
  
  # check for specific genes listed in the panel entry
  panel.entry.vec <- NULL
  if(!is.na(panel.entry[1])){
    panel.entry.split <- str_split(panel.entry, pattern = " ")[[1]]
    
    # find all the genes in the panel entry
    for(gene in panel.entry.split){
      if(sum(common.genes %in% gene) == 1){ panel.entry.vec <- c(panel.entry.vec,gene) }
    }
  }
  
  # check for result genes
  result.genes.vec <- NULL
  if(!is.na(result.genes[1])){
    result.genes.split <- str_split(result.genes, pattern = " ")[[1]]
    
    # find all the genes in the result genes entry
    for(gene in result.genes.split){
      if(sum(common.genes %in% gene) == 1){ result.genes.vec <- c(result.genes.vec,gene) }
    }
  }
  
  # consolidate
  genes.vec <- unique(panel.entry.vec,result.genes.vec)
  genes.vec
}


#------------------------------------------------------------------------------
compare.results <- function(gene.vec, date.vec, test.num.vec){
  #' function for resolving different results for the same gene across multiple
  #' genetic tests.  Rules:
  #' 1) any non-NA result trumps all NA results
  #' 2) neg, pos, and vus discrepancies are resolved by the most recent test date
  #' 3) if test dates are not available, then it uses the highest test number (1-12)
  #' @param gene.vec. numeric vector of the genes results (NA = no test, 0 = neg, 1 = pos, 2 = vus)
  #' @param date.vec. date vector of test dates
  #' @param test.num.vec. numeric vector of the test number (1-12)
  #' returns NA or numeric value of 0, 1, or 2
  
  # check if there is a single unique type of non-NA value, if so, return that value
  if(length(table(gene.vec)) == 1){ return(as.numeric(names(table(gene.vec)))) }
  
  # RESOLVE CONFLICTS
  # remove NA results
  date.vec.tmp <-     date.vec    [which(!is.na(gene.vec))]
  test.num.vec.tmp <- test.num.vec[which(!is.na(gene.vec))]
  gene.vec.tmp <-     gene.vec    [which(!is.na(gene.vec))]
  
  # remove neg results with NA dates
  if(sum(is.na(date.vec.tmp) & gene.vec.tmp == 0) > 0){
    test.num.vec.tmp <- test.num.vec.tmp[which(!(is.na(date.vec.tmp) & gene.vec.tmp == 0))]
    date.vec.TEMP    <- date.vec.tmp
    date.vec.tmp     <- date.vec.tmp    [which(!(is.na(date.vec.tmp) & gene.vec.tmp == 0))]
    gene.vec.tmp     <- gene.vec.tmp    [which(!(is.na(date.vec.TEMP) & gene.vec.tmp == 0))]
  }
  
  # special case: ignore neg results if they are the only ones with dates
  valid.date.0 <- sum(!is.na(date.vec.tmp[which(gene.vec.tmp == 0)]))
  valid.date.1 <- sum(!is.na(date.vec.tmp[which(gene.vec.tmp == 1)]))
  valid.date.2 <- sum(!is.na(date.vec.tmp[which(gene.vec.tmp == 2)]))
  if(valid.date.0 > 0 & valid.date.1 == 0 & valid.date.2 == 0){
    
    # remove neg results
    test.num.vec.tmp <- test.num.vec.tmp[which(gene.vec.tmp != 0)]
    gene.vec.tmp     <- gene.vec.tmp    [which(gene.vec.tmp != 0)]
    
    # return most recent result by test order
    return(gene.vec.tmp[which(test.num.vec.tmp == max(test.num.vec.tmp))])
  }
  
  # check if any dates remain
  if(sum(!is.na(date.vec.tmp)) > 1){
    
    # sort by date
    gene.vec.tmp <- gene.vec.tmp[order(date.vec.tmp, decreasing = T)]
    
    # if there is a tie for first remove the 0 results as options
    if(sum(date.vec.tmp == max(date.vec.tmp, na.rm = T), na.rm = T) > 1){
      gene.vec.tmp <- gene.vec.tmp[which(date.vec.tmp == max(date.vec.tmp))]
      gene.vec.tmp <- gene.vec.tmp[which(gene.vec.tmp != 0)]
    }
    
    # return most recent and relevant result
    return(gene.vec.tmp[1])
    
    # if all remaining dates are NA, use test number after removing negs
  } else {
    
    # remove negs
    test.num.vec.tmp <- test.num.vec.tmp[which(gene.vec.tmp != 0)]
    gene.vec.tmp     <- gene.vec.tmp    [which(gene.vec.tmp != 0)]
    
    # return most recent result by test order
    return(gene.vec.tmp[which(test.num.vec.tmp == max(test.num.vec.tmp))])
  }
}