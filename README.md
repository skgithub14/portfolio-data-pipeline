# Showcase Data Pipeline for Portfolio

This repository contains a modified copy of a complex data cleaning pipeline I wrote in R. Due to confidentialy, it has modified to remove all identifying data and source information. 

## Portfolio Information

See the included `presentation.pdf` file for an overview with visuals of the project's highlights.

### Description

**Purpose:** clean and curate a messy clinical data set from genetic clinic sites from across the US and Latin America. Format it for analysis by two different hereditary cancer syndrome models, BayesMendel and PanelPRO.

**Content:** The dataset consists of pedigrees (family trees in tabular format) with demographics, genetics, surgical, cancer history data, and other clinical data. 

**Initial Number of Rows:** _____ (each row represents one person in one family)

**Initial Number of Columns:** ______

**Source:** data was collected from genetics clinics. The exact source and dataset are confidential. 

### Complex Pipeline Features:

**Linked Records:** 
 - Records within the same family must be linked for the models to run.
 - Example 1: race of 1 relative is based on a mix of their mother’s and father’s race.
 - Example 2: a child must have both a male and female biological parent who are older than the child.
 - _Solution:_ compared and imputed missing records within families for consistency.

**Gene Tests:**
 - Understanding which genes were tested and which were not is critical data for the models. 
 - Each test consists of a panel of one to 800+ genes. Panel names and the genes they include are non-standard and the genes in each panel change over time. Most genes with “negative” results were not recorded which makes the model results less accurate.
 - _Solution:_ automated the comparison of the unstructured panel names, dates of testing, recorded results, and published panel information from gene testing labs. Selected the most likely panel for each record and use it to impute all missing genes.

**Cancer history:**
 - Cancer history is also critical for the model but…
 - Cancer diagnoses were unstructured, inconsistently coded, and often misspelled.
 - _Solution:_ created a dictionary of word roots, common medical codings, and common misspellings.

This repository showcases my ability to write a complex data pipeline in R however, due to the intentionally omitted data, the pipeline unfortunately cannot be run.



## Pipeline Documentation

This pipeline documentation has four sub-sections:
  - **Overview**
  - **Variable Codings**: description of variable codings for variables not used in BayesMendel or PanelPRO.
  - **Scripts**: description of each step of the data cleaning process by script.
  - **Assumptions**: description of each assumption used in the cleaning process.

### Overview

This repo contains the files required to clean and wrangle the REDACTED genetics pedigree data set.  The pipeline produces three clean variants of the same raw data set every time the pipeline runs.  
  - The "all genes" variant (`clean-REDACTED-data-all-genes.RData`) contain all 800+ genes in the REDACTED data. This file is >6Gb and won't open on most local machines. You can see the data summary for all 800+ genes in the `data-version-summary-stats-[VERSION DATE].html` file though. 
  - The BayesMendel variants (`clean-REDACTED-data-BayesMendel.RData`) have the clean data formatted for the [BayesMendel R package](https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package) with some additional columns for identifying important features of each subject. 
  - The PanelPRO variant (`clean-REDACTED-data-PanelPRO.RData`) have the data set formatted for the [PanelPRO R package](https://github.com/bayesmendel/PanelPRO) with some additional columns for identifying important features of each subject. 

There are also three documentation files generated each time the pipeline runs: 

  - `data-version-summary-stats-[VERSION DATE].html`: basic summaries of each variable.
  - `sessionInfo-[VERSION DATE].txt`: a report of the package versions and system used.
  - `README-[VERSION DATE].md`: contains a list of changes made to the pipeline since the last run.

All versions of the data set can be found on REDACTED Dropbox: 

`REDACTED VALUE` 

The directory above has a sub-directory for each date the pipeline ran.

### Variable Codings

For an explanation of any variables that are utilized by either [PanelPRO R package](https://github.com/bayesmendel/PanelPRO) or [BayesMendel R package](https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package), please see the documentation for those R packages.

The following variables are included in the clean data but are not used by PanelPRO of BayesMendel.  These variables can be utilized by the user to stratify or filter subjects and pedigrees based on the needs of the researchers.

- `PedigreeID`: an anonymized character vector indicating which pedigree (family) a subject belongs to.
- `SubjectID`: a concatenation of the `PedigreeID` and `ID` columns with an underscore between them.  This allows a subject to be uniquely identified in the data set quickly and easily.
- Region: `"US"` or `"LATAM"`. Indicates if a subject was drawn from United States REDACTED data or LATin AMerica REDACTED data.
- isHispanic: `1` or `0`.  Indicates if a subject has Hispanic ethnicity.  Due to the unique race category requirements of BayesMendel and PanelPRO, determining Hispanic ethnicity based on those categories may not be possible for mixed race / mixed ethnicity individuals.  For example, PanelPRO has a race category of `"WH"` (White-hispanic) but no unique categories for a subject that reports as being both Hispanic and any other Non-white race.  An individual who reports being Hispanic and Black is therefore categorized in PanelPRO as `"All_Races"`.  In BayesMendel, an individual reporting both Hispanic ethnicity and any race would be categorized as `"Unknown"`.
- `isProband`: `1` or `0`. This is a required variable in PanelPRO, see the PanelPRO documentation for an explanation. It is a useful variable so it was also included in the BayesMendel versions, although BayesMendel itself does not use it.
- `CarrierStatus`: indicates if a subject is a positive carrier, vus carrier, non-carrier, or was never tested.  It uses the same codings for genes as the BayesMendel and PanelPRO for the respective variants of the data set.  The carrier statuses are determined by the "all-genes" variants of the data set, meaning over 800 genes are considered.  A subject is considered positive if they had at least one positive gene, VUS if they had no positive genes and at least one VUS gene, negative if all genes tested were negative.  If a subject was never tested for any gene, they are coded as such.
- Note on the codings of gene tests for the "all-gene" versions of the data set: these follow the PanelPRO codings.


### Scripts

The pipeline consists of 7 main steps divided into 8 scripts.  The `MAIN.R` script runs all 8 scripts and loads and saves the required data for each step.

#### 1_split-data.R:

Splits the data into demographics data and non-demographics data (which will be rejoined later in the pipeline).

1. Starts with raw data
2. Remove duplicate subjects
3. Rename columns (note `Age.at.entrance` becomes `CurAge` here)
4. Clean ages (absolute value, round to ones place, ages from `95` to `125` censored to `94`, ages above `125` become `NA`)
5. Replace identifiable PedigreeIDs
6. Split data into demographics columns (races, ethnicities, ages, etc.) and non-demographics columns (testing, surgeries, cancer history)

#### 2_non-demographics-non-test-data.cleaning.R

1. Drop unneeded tumor columns (`IHC` dropped but `MSI`, `HER2`, `PR`, and `ER` are kept; note `CK14` and `CK5.6` are not in the raw data).
2. Mastectomy: recode and compute age of surgery using date of birth (note `AgeMastectomy` for subjects without `Mastectomy` set to `CurAge` / `AgeDeath` per BayesMendel format)
3. Oophorectomy: recode and compute surgery age (`NA` if < 1 or > 120)
(note `AgeMastectomy` and `AgeOophorectomy` for subjects without those surgeries set to `CurAge` / `AgeDeath` per BayesMendel format)
4. Marker Tests: `ER` and `PR` recoded. `HER2` and `MSI` wrangle multiple columns of data and recode them.
5. Cancer History: classifies cancer history in raw data into PanelPRO cancer types and wrangles into PanelPRO cancer history format.
6. Drop all surgeries and cancers occurring after current age.
7. Eliminate non-prophylactic surgeries.  Assumes any surgery conducted before onset of related cancer is prophylactic.  Assumes any surgery for unaffected person for relevant cancer is prophylactic.

#### 3_gene-test-data-cleaning.R

1. Drop variants and DOB columns
2. Assign easy and intuitive column names
3. Format test dates
4. Recode genetic test results
5. Clean result gene names
6. Create storage columns for results by gene (empty)
7. Clean lab names
8. Clean panel names
9. Combine 3 different testing/panel columns based on condition
10. Manual panel name cleaning integration
11. Clean specific genes tested

#### 4a_panel-to-gene-mapping-functions.R

Contains the following functions (see function doc strings for details):
  - `impute.missing.panel.date()`
  - `find.genes.in.panel()`
  - `upgrade.panel.based.on.date()`
  - `find.panels.that.match()`
  - `find.common.genes()`
  - `get.specific.genes.tests()`
  - `compare.results()`

#### 4_gene-test-wrangling.R

Produces two data frames with test results, one by test and one by subject.

1. Reshapes data to index by unique test instead of subject
2. Populate panel keyword matrix for each panel
3. Record specific genes tested
4. For panels, search keyword matrix for clues, find all possible panels, and find genes those panels share
5. Record tests for entries without discernible specific test or panel test but that have result genes
6. Record test results by gene
7. Resolve conflicts and consolidate test results for individuals with multiple genetic tests

#### 5_demographics-cleaning.R

1. Recode gender
2. Drop subjects with missing `MotherID` and `FatherID`
3. Identify `Twins`
4. Recode `isProband`
5. `AJ` versus `non-AJ` ethnicity
6. Recode `isHispanic`
7. Race codings (for both BayesMendel and PanelPRO)

#### 6_merge.R

Creates the clean master data frame with all tested genes with ages, cancers, and surgeries.

1. Join cancer history, surgeries, and demographics
2. Impute and clean sex assignments, female surgeries, and female cancers
3. Join in test data
4. Create carrier identifier column (`1` means at least 1 positive gene, `2` means at least 1 VUS gene and no positive genes, `0` means testing but no positive or VUS genes, `NA` means no tests conducted)
5. Convert columns to correct classes
6. Remove marker tests for subjects without associated cancers

#### 7_final-formatting.R

This script takes the master clean data frame with all genes tested and formats it into variants for BayesMendel and PanelPRO.  

BayesMendel Steps:
1. Select and rename columns
2. Create `CK14` and `CK5.6` dummy columns
3. Recode gene columns for BayesMendel
4. Assign one race per family
5. Change affected ages for unaffected people to `CurAge`/`AgeDeath` (BayesMendel requirement)
6. Order by `PedigreeID` and `ID` and convert to data frame

PanelPRO Steps:
1. Select and rename columns
2. Create `CK14` and `CK5.6` dummy columns
3. Recode marker tests for PanelPRO format
4. Wrangle surgery columns into PanelPRO format
5. Order by `PedigreeID` and `ID` and convert to data frame

Note the differences in gene test coding between BayesMendel and PanelPRO:
  - BayesMendel: `0` = no test, `1` = positive, `2` = negative, `3` = VUS
  - PanelPRO: `NA` = no test, `0` = negative, `1` = positive, `2` = VUS

### Assumptions

#### 1. Age cleaning and adjustments:

  - Ages above `125` and less than `1` were assumed to be misreported and turned to `NA`.
  - Ages of surgeries and cancers that were older than the `CurAge` / `AgeDeath` were assumed to be reported after the family entered the study so those surgeries and cancers were removed from the data set.

#### 2. Imputation of sex: 
  - if a female cancer or prophylactic surgery was found the sex was assumed to be female. 
  - if subject had a child with a person with a known sex, then the opposite sex was assumed. 
  - if males were assigned female cancers or surgeries: if they were fathers, remove the female cancers and surgeries, otherwise change the subjects sex to female
  - if mother assigned male sex or father assigned female sex then the sexes were changed to match parental title

#### 3. Race, ethnic and Ancestry assumptions:
  - Native Hawaiian / Pacific Islander was grouped into `"Asian"` because NHPI is not a race group in either BayesMendel or PanelPRO. 
  - If a subject's mother or father was `"AJ"`, the subject was assumed to be `"AJ"`, otherwise `"non-AJ"` was assumed.  There was no data on `"Italian"` or `"other"` in the raw data. 
  - For multi-race subjects in BayesMendel (which includes Hispanic as its own race category) those subjects had their race listed as `"Unknown"`.
  - For multi-race subjects in PanelPRO (which has `"WH"` (White-Hispanic) `"WNH"` (White Non-Hispanic), and `"White"` categories), Hispanic only subjects were assumed to be `"WH"` if they did not report anything for White, White subjects without reporting on Hispanic ethnicity were just listed as `"White"`, subjects that reported being both White and Hispanic were listed as such, subjects reporting White and Non-Hispanic were listed as such.  Any other combinations of races were assigned `"All_Races"`.
  - to account for the fact BayesMendel can only handle one race label per family, all races in a family were checked for a single unique value.  If a single unique values was found all members of the family were changed to that race.  If multiple races labels were found in the family, all family members' races were changed to `"Unknown"`.

#### 4. Surgery cleaning and adjustments:

  - `Year.of.DOB`, which is used to help calculate age of oophorectomy, was assumed to be mis-entered if it was less than `1800` or greater than `2021`.  These entries were set to `NA`.
  - Non-prophylactic surgeries were reported in this raw data set.  To eliminate those the following rationale was used: assume surgery was prophylactic if the person is unaffected by the relevant cancers or if affected by the relevant cancers, then assume the surgery was prophylactic if the cancer age and surgery age are present AND the cancer age is greater than the surgery age.

#### 5. Cancers cleaning and adjustments:

  - See the lists named `raw.cancers` and `unclass.cancers` in script `2_non-demographics-non-test-data-cleaning.R` for how cancer names from the medical records were classified into the 19 PanelPRO cancer types (and which could not be classified).

#### 6. Biomarker codings and adjustments:

  - The following code excerpts were used to translate the corresponding marker test result labels to the BayesMendel codings (`0`=no test, `1`=pos, `2`=neg). VUS as coded as `3`.
	```
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
	```
  - Note: CK14 and CK5/6 were not reported in the raw data and therefore the columns are filled with the equivalent of no test, depending on the format.
  - Marker tests reported for people without the corresponding cancers (Breast Cancer for `HER2`, `ER`, `PR` and Colorectal Cancer for `MSI`) were assumed to be invalid and dropped.  For `MSI` and `HER2`, if multiple test results were reported: if at least one was positive then final result assumed positive, if no positives and at least one negative then assumed to be negative.
        
#### 7. Gene codings and adjustments:

  - The following original gene test result labels were translated using the PanelPRO codings (`NA`=no test, `0`=neg, `1`=pos).  VUS's were coded as `2`.
	```
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
	```
  - Positive, Negative, and VUS status: due to limited reporting on variant informations, the status of the gene at the time of the test was assumed to have not changed since the original test despite the highly likely chance that the variant classifications changed over time.
  - The following code excerpt was used to assume the standardized gene names and remove gene names that could not be assumed.
	```
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
	```
	                   
#### 8. Panel assignments and testing:

  - The following panel name keywords were assumed to relate to the following misspellings and short hand.
	```
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
	```
  - The following keywords were used to help classify panel labels into standard panel names:
	```
	keywords <- c("ova","hereditary","brca","breast","myrisk","next","expanded","analysis",
	              "plus","gyn","management","gastric","breast ovarian","cancernext",
		      "colaris","colorectal","colon","endometrial","high","multi","prostate")
	```     
  - These two files contain lists of manually assigned panel names and genes based based on the provided panel / gene label from the raw data:
    - `./R/data/reference-data/current-version/unique_panels_categorized.csv`
    - `./R/data/reference-data/current-version/unique_panels_categorized_specific.csv`
  - Determining which subset of specific genes were tested for individuals based on panel name notes in the raw data was complex and relied on an algorithm to look at all possible standard panel names, the dates those panels were run, and test result genes.  Because the algorithm relied on messy notes, it is imperfect and therefore the resulting clean data may not reflect reality. The results of the algorithm were manually evaluated and adjusted accordingly.  For a complete understanding of how the list of genes tested for each subject was determined, review script `4a_panel-to-gene-mapping-functions.R`.  The following is a list of assumptions made in that script:
    - if a panel's date was missing, it was assumed to be the same date as the first date that panel was offered.  If that first offered date was also unavailable then the test date was set to an arbitrary date of `1990-01-01`.
    - If a test date was before 1990-01-01 or the date fell before or after the panel was created or retired then it was assumed the proposed panel was wrong and the function will not return a list of genes tested.  There was an exception for panels that were upgraded to a new panel name but had the old panel name in the raw data.  If this occurred it was assumed the upgraded panel was tested.  The newer panels are the list names and the older ones are the values:
		```
		upgrades <- list(brca.plus.ambry                    = "brca.plus.expanded.ambry",
				 breast.cancer.management.genedx    = "breast.cancer.high.moderate.risk.genedx",
				 breast.gyn.cancer.genedx           = c("breast.ovarian.cancer.genedx","endometrial.cancer.genedx"),
				 comprehensive.common.cancer.genedx = "comprehensive.cancer.genedx",
				 common.cancer.management.genedx    = "high.moderate.risk.cancer.genedx")
		```
    - The menu of genes in a panel usually changes over time.  For most of the change dates we only had access to the month and year of the change.  Any panel test date that fell during a month of change in menu assumed the genes tested were only the ones that were common between the previous version of the menu and the next version of the menu.
    - If the all possible panels had 0 common genes, then only positive or vus genes listed in the result were considered to be tested.
    - If results for the same gene in the same individual conflict, then the most recent result was used.
