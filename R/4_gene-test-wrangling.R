# This script will process all panel entries and will encode them into the individual gene columns

#### Re-stack Data ####
test.cols <- c("Date","Lab","Specific","Panel","Result","Result.Gene")
gene.col.cnt <- ncol(datw %>% select(ABL1:ZRSR2))
test.dat <- as.data.frame(matrix(, nrow = 0, ncol = 1 + length(test.cols) + gene.col.cnt))
for(test.num in 1:12){
  tmp.test.cols <- paste0(test.cols,test.num)
  tmp.test.dat <- 
    datw %>%
    select(SubjectID,one_of(tmp.test.cols),ABL1:ZRSR2) %>%
    rename_with(~ sub(pattern = test.num, replacement = "", .x), .cols = one_of(tmp.test.cols)) %>%
    filter(if_any(.cols = one_of(test.cols[c(2,4,6)]),
                  ~ !is.na(.))) %>%
    add_column(test.num, .after = "SubjectID")
  test.dat <- rbind(test.dat,tmp.test.dat)
}

#### Create Keyword Matrix ####
# convert panel entry to keywords
keywords <- colnames(panel.name.keywords)[-1]
keywords <- sub(pattern = "next.", replacement = "next", keywords)
keywords <- sub(pattern = "X700", replacement = "700", keywords)
keywords <- sub(pattern = "X760", replacement = "760", keywords)
note.keywords <- as.data.frame(matrix(, nrow = nrow(test.dat), ncol = length(keywords)))
note.keywords <- cbind(test.dat$SubjectID, test.dat$test.num, note.keywords)
colnames(note.keywords) <- c("SubjectID","test.num",keywords)

# populate keyword matching data in the note.keywords data frame
for(word in keywords){
  note.keywords[,word] <- grepl(pattern = word, test.dat$Panel)
}

#### Record Test Results by Test  ####
# loop through each test
for(row in 1:nrow(test.dat)){
  
  # retrieve row's critical information
  tmp.sub         <- test.dat[row,"SubjectID"]
  tmp.test.num    <- test.dat[row,"test.num"]
  tmp.date        <- test.dat[row,"Date"]
  tmp.lab         <- test.dat[row,"Lab"]
  tmp.specific    <- test.dat[row,"Specific"]
  tmp.panel       <- test.dat[row,"Panel"]
  tmp.result      <- test.dat[row,"Result"]
  tmp.result.gene <- test.dat[row,"Result.Gene"]
  
  ##### Specific Genes Tested ####
  tmp.tested.genes <- NULL
  
  # check if specific genes were tested or a panel
  if(tmp.specific == 1){
    tmp.tested.genes <- get.specific.genes.tested(panel.entry = tmp.panel, 
                                                  result.genes = tmp.result.gene)
    
  ##### Panel Tested ####  
  # not a specific set of genes ordered, therefore it was a panel
  } else {
    
    # find keyword hits
    tmp.keys.row <- note.keywords[row,keywords]
    tmp.keys <- names(which(colSums(tmp.keys.row) == 1))
    
    # if there are keywords found attempt to determine the possible panels and genes tested
    if(length(tmp.keys) != 0) { 
    
      # get list of possible panels
      possible.panels.vec <- find.panels.that.match(panel.keywords = tmp.keys, lab.name = tmp.lab, test.date = tmp.date)
      
      # get list of genes tested
      tmp.tested.genes <- find.common.genes(panel.names.vec = possible.panels.vec, 
                                            test.date = tmp.date, 
                                            result.genes = tmp.result.gene, 
                                            panel.keywords = tmp.keys)
      
      # if there were no keywords found, at least check if there was a result gene recorded
    } else {
      
      # find all the genes in the result genes entry
      if(!is.na(tmp.result.gene)){
        result.genes.split.1 <- str_split(tmp.result.gene, pattern = " ")[[1]]
        for(gene in result.genes.split.1){
          if(sum(common.genes %in% gene) == 1){ tmp.tested.genes <- unique(c(tmp.tested.genes,gene)) }
        }
      }
    }
  }
  
  ##### Record Tested Genes ####
  # initialize by recording all tested genes as negative (0)
  # genes with pos and vus results will be updated from 0 in the next section
  test.dat[row,tmp.tested.genes] <- 0
  
  ##### Record Result Genes ####
  # if there is an entry in the result and that result is non-negative then record
  if(is.na(tmp.result)){ next }
  if(tmp.result != 0){
    result.genes.vec <- NULL
    result.genes.split <- str_split(tmp.result.gene, pattern = " ")[[1]]

    # find all the genes in the result genes entry
    for(gene in result.genes.split){
      if(sum(common.genes %in% gene) == 1){ result.genes.vec <- unique(c(result.genes.vec,gene)) }
    }

    # if valid result genes were found then update the result
    if(length(result.genes.vec) > 0){
      test.dat[row,result.genes.vec] <- tmp.result
    }
  }
}


#### Combine Subjects' Tests ####

# create new test data frame with single test result for each subject
subs.1.test <- names(which(table(test.dat$SubjectID) == 1))
new.test.dat <- 
  test.dat %>%
  select(SubjectID,ABL1:ZRSR2) %>%
  filter(SubjectID %in% subs.1.test)

# template df for subjects with multiple tests
subs.mult.tests <- names(which(table(test.dat$SubjectID) > 1))
combined.test.dat <- as.data.frame(matrix(,nrow = length(subs.mult.tests), ncol = gene.col.cnt))
combined.test.dat <- as.data.frame(cbind(subs.mult.tests,combined.test.dat))
colnames(combined.test.dat) <- colnames(new.test.dat)

# populate the template by looping through all subjects with multiple tests
for(loop.sub in subs.mult.tests){
  
  # subset test data for this individual
  tmp.df <- test.dat[which(test.dat$SubjectID == loop.sub),]
  
  # remove columns filled with NAs
  na.col.cnt <- sapply(tmp.df %>% select(ABL1:ZRSR2), function(x) sum(is.na(x)), simplify = T)
  non.na.cols <- names(which(na.col.cnt != nrow(tmp.df)))
  if(length(na.col.cnt) == 0){ next } # skip subject if all genes for all tests are NA
  tmp.df.1 <- tmp.df[,c("SubjectID","test.num",test.cols,non.na.cols)]
  
  # resolve conflicts by looping through each gene
  for(tmp.col in non.na.cols){
    best.result <- compare.results(gene.vec = tmp.df.1[,tmp.col], date.vec = tmp.df.1[,"Date"], test.num.vec = tmp.df.1[,"test.num"])
    combined.test.dat[which(combined.test.dat$SubjectID == loop.sub),tmp.col] <- best.result
  }
}

# merge the subjects with one test and multiple tests
new.test.dat <- rbind(new.test.dat,combined.test.dat)

# merge back into master data frame
datw2 <- 
  datw %>%
  select(SubjectID:MSI)

# join data frames
datw2 <- left_join(datw2, new.test.dat, by = "SubjectID")
