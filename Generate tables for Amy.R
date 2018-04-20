#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.04.19.
#' 
#' Prepare data tables for Amy Herring.  
#' 
#' She is a Bayesian statistician and we are interested in working with her 
#' to identify recurrent patterns of multiple congenital anomalies.
#' 
#' As a first step she requested some data tables, containing information 
#' like:
#' - The number and proportion of children with 1,2,...,n defects 
#' stratified by cancer status.
#' - For each defect, a table listing the median number of co-occurring 
#' birth defects in children with that defect, stratified by cancer status
#' (she didn't ask for this, but I think it would be helpful). 
#' - The prevalence of each individual defect, plus the prevalence of 
#' combinations of 2,...,n defects conditional on having the first defect, 
#' also stratified by cancer status.
#' - The number of children with multiple defects in each organ system.
#' 
#' For now, I'll do this in the non-chromosomal dataset.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180419.rdata')



# N and % of children by number of defects --------------------------------

require(dplyr); require(tictoc)

#' Overall in non-chromosomal set.
tic()
tab <- as.numeric(table(goback.nochrom$defect.total))
n <- as.numeric(nrow(goback.nochrom))
num <- sort(as.numeric(unique(goback.nochrom$defect.total)))
toc()

defect.total.overall <- data.frame(n.defects = num, n.obs.overall = tab, prop.obs.overall = tab/n)

#' Stratified by cancer status in non-chromosomal set.
tab <- table(goback.nochrom$cancer, goback.nochrom$defect.total, useNA = 'ifany')
tab.controls <- as.numeric(tab[1,])
tab.cases <- as.numeric(tab[2,])

n.controls <- 10144203
n.cases <- 14773

defect.total.controls <- data.frame(n.defects = num, n.obs.controls = tab.controls, prop.obs.controls = tab.controls/n.controls)
defect.total.cases <- data.frame(n.defects = num, n.obs.cases = tab.cases, prop.obs.cases = tab.cases/n.cases)

defect.totals <- full_join(
                            full_join(defect.total.overall, defect.total.cases, by = 'n.defects'),
                            defect.total.controls, by = 'n.defects')

write.csv(defect.totals, file = 'Z:/Jeremy/Multiple anomalies project/number.of.defects.by.cancer.status.csv', row.names = FALSE)

rm(list = ls()); gc()



# The median number of co-occurring defects stratified by cancer --------

require(dplyr)

#' Initialize an empty data frame to hold new data.
comorbid.counts <- data.frame(defect = as.character(NA), 
                              comorbid.noncancer = as.numeric(NA),
                              comorbid.cancer = as.numeric(NA))

#' A vector of column indices for all specific non-chromosomal defects to iterate over.
cols <- c(23:29,31:34,36,37,40:44,46:49,51,52,54:60,62:64,66,67,69:75,77:82,84:94)

for (i in cols){

  tmp <- filter(goback.nochrom, goback.nochrom[ ,i] == 1)

  flag <- as.numeric(length(unique(tmp$cancer)))
  
  if (flag > 1){ 
    
    counts <- aggregate(defect.total ~ cancer, data = tmp, median)
    counts.noncan <- counts[1,2]
    counts.can <- counts[2,2]
    defect.name <- names(goback.nochrom[i])
    
    new.data <- data.frame(defect = defect.name, comorbid.noncancer = counts.noncan, comorbid.cancer = counts.can)
    
    comorbid.counts <- rbind(comorbid.counts, new.data)

  }
  
  else {
    
    rm(tmp, flag)
    
  }
}

comorbid.counts <- comorbid.counts[2:58, ]

write.csv(comorbid.counts, file = 'Z:/Jeremy/Multiple anomalies project/median.defect.count.by.cancer.status.csv', row.names = FALSE)

rm(list = ls()); gc()



# Frequency of each individual defect -------------------------------------

require(tictoc)

#' A vector of column indices for all specific non-chromosomal defects to iterate over.
cols <- c(23:29,31:34,36,37,40:44,46:49,51,52,54:60,62:64,66,67,69:75,77:82,84:94)

#' Initialize an empty data frame to house results.
defect.freq <- data.frame(defect = as.character(NA), 
                          n.defect.total = as.numeric(NA),
                          prop.defect.total = as.numeric(NA),
                          n.defect.nocancer = as.numeric(NA),
                          prop.defect.nocancer = as.numeric(NA),
                          n.defect.cancer = as.numeric(NA),
                          prop.defect.cancer = as.numeric(NA))

tic()

for (i in cols){
  
  print(paste('Starting iteration', i))
  
  defect <- names(goback.nochrom[i])
  
  tab <- table(goback.nochrom$cancer, goback.nochrom[, i])
  
  overall.defect <- as.numeric(sum(tab[, 2]))
  overall.total <- as.numeric(sum(tab[, 1] + tab[, 2]))
  nocan.defect <- as.numeric(tab[1,2])
  nocan.total <- as.numeric(tab[1,1] + tab[2,2])
  can.defect <- as.numeric(tab[2,2])
  can.total <- as.numeric(tab[2,1] + tab[2,2])
  
  new.data <- data.frame(defect = defect, 
                            n.defect.total = overall.defect,
                            prop.defect.total = overall.defect/overall.total,
                            n.defect.nocancer = nocan.defect,
                            prop.defect.nocancer = nocan.defect/nocan.total,
                            n.defect.cancer = can.defect,
                            prop.defect.cancer = can.defect/can.total)
  
  defect.freq <- rbind(defect.freq, new.data)
  
  rm(can.defect, can.total, defect, i, nocan.defect, nocan.total, overall.defect, overall.total, tab, new.data)
  
}

toc()

defect.freq <- defect.freq[2:59, ]

write.csv(defect.freq, file = 'Z:/Jeremy/Multiple anomalies project/proportion.with.each.defect.csv', row.names = FALSE)

rm(list = ls()); gc()



# What are the most common comorbid defects by defect ---------------------

require(dplyr); require(tictoc)

#' A vector of column indices for all specific non-chromosomal defects to iterate over.
cols <- c(23:29,31:34,36,37,40:44,46:49,51,52,54:60,62:64,66,67,69:75,77:82,84:94)

count.data <- data.frame(index.defect = as.character(NA),
                         comorbid.defect = as.character(NA),
                         comorbid.count = as.numeric(NA))

tic()
for (i in cols){
  
  tmp <- filter(goback.nochrom, goback.nochrom[, i] == 1)
  
  for (j in cols){
    
    tab <- as.numeric(sum(table(tmp[, i], tmp[, j])[1, ]))
    
    comorbid.defect <- names(goback.nochrom[j])
    
    new.count <- data.frame(index.defect = names(goback.nochrom[i]),
                            comorbid.defect = names(goback.nochrom[j]),
                            comorbid.count = tab)
    
    count.data <- rbind(count.data, new.count)
    
  }
  
}
toc()

count.data <- count.data[2:3601, ]
count.data <- count.data[count.data$index.defect != count.data$comorbid.defect, ]

write.csv(count.data, file = 'Z:/Jeremy/Multiple anomalies project/conditional.defect.frequencies.csv', row.names = FALSE)
