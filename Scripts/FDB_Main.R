##### Set up #####
setwd("~/Documents/GitHub/FDB_Freeland/")

library(doParallel)
library(doRNG)
library(missForest)
library(dplyr)
library(mixOmics)

#### Imputation: CRISPR (NA's in Data) ####

## pull in data
file.crispr <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/CRISPRGeneEffect.csv"

CRISPR <- read.delim(
  file = file.crispr,
  row.names = 1,
  stringsAsFactors = F,
  sep = ","
)

table(colSums(is.na(CRISPR)))

## miss forest  
doParallel::registerDoParallel(cores = detectCores() - 2)
doRNG::registerDoRNG(seed = 999)
set.seed(999)

CRISPR_mf <- missForest(
  xmis = CRISPR,
  parallelize = "variables",
  verbose = T
)

CRISPR_mf.imp <- CRISPR_mf$ximp
CRISPR_mf.imp_t <- as.data.frame(
  t(CRISPR_mf.imp),
  stringsAsFactors = F
)

## save 
write.table(
  x = CRISPR_mf.imp, 
  file = gsub(".csv$","_MFImputed.txt",file.crispr), 
  quote = F, 
  sep = "\t",
  col.names = NA
)

write.table(
  x = CRISPR_mf.imp_t, 
  file = gsub(".csv$","_MFImputed_sg.txt",file.crispr), 
  quote = F, 
  sep = "\t",
  col.names = NA
)

##### Imputation: CTRP (NA's in Data)#####

## provide initial paths
path.dm   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/"
path.ctrp <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/CTRPv2/"

## load in cell line info from depamp
models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F)

## load in CTRP Data
ctrp.expt   <- read.delim(paste0(path.ctrp,"v20.meta.per_experiment.txt"), sep = "\t", stringsAsFactors = F)
ctrp.cell   <- read.delim(paste0(path.ctrp,"v20.meta.per_cell_line.txt"), sep = "\t", stringsAsFactors = F)
ctrp.inform <- read.delim(paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F)
ctrp.curves <- read.delim(paste0(path.ctrp,"v20.data.curves_post_qc.txt"), sep = "\t", stringsAsFactors = F)

# temp <- as.data.frame(table(ctrp.expt$master_ccl_id))

## add ID's and name to curves
ctrp.curves$master_ccl_id <- ctrp.expt$master_ccl_id[match(ctrp.curves$experiment_id, ctrp.expt$experiment_id)]
ctrp.curves$ccl_name      <- ctrp.cell$ccl_name[match(ctrp.curves$master_ccl_id, ctrp.cell$master_ccl_id)]
ctrp.curves$DepMap_ID     <- models$ModelID[match(ctrp.curves$ccl_name, models$StrippedCellLineName)]
ctrp.curves$cpd_name      <- ctrp.inform$cpd_name[match(ctrp.curves$master_cpd_id, ctrp.inform$master_cpd_id)]

## move important columns to front
ctrp.curves <- ctrp.curves %>% 
  dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, area_under_curve, experiment_id, everything())

not.mapped.celllines <- ctrp.curves[is.na(ctrp.curves$DepMap_ID),]
not.mapped.celllines <- not.mapped.celllines[,1:3] %>% dplyr::distinct()

write.table(
  x = ctrp.curves,
  file = paste0(path.ctrp,"ctrp.curves.txt"),
  quote = F, col.names=T, row.names = F, sep = "\t")

# ctrp.rsl3 <-  ctrp.curves[ctrp.curves$cpd_name == "1S,3R-RSL-3",]
# table(table(ctrp.rsl3$DepMap_ID))
# write.table(ctrp.rsl3, paste0(path.ctrp,"ctrp.rsl3.txt"), quote = F, col.names=T, row.names = F, sep = "\t")

## trim and average data frame
ctrp.curves.abr <- ctrp.curves %>% 
  dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, master_cpd_id, area_under_curve)

ctrpv2.ave <- ctrp.curves.abr %>%
  dplyr::group_by(DepMap_ID, ccl_name,master_ccl_id, cpd_name,master_cpd_id) %>%
  dplyr::summarise(avg = mean(area_under_curve)) %>%
  dplyr::ungroup()

# ctrp.ave.rsl3 <- ctrpv2.ave[ctrpv2.ave$cpd_name == "1S,3R-RSL-3",]
# names(ctrpv2.ave)

ctrpv2.ave.wide <- ctrpv2.ave %>% 
  dplyr::select(DepMap_ID, ccl_name, cpd_name, avg)
ctrpv2.ave.wide.pre <- as.data.frame(ctrpv2.ave %>% dplyr::select(DepMap_ID, cpd_name, avg))

ctrpv2.ave.wide <- reshape(
  ctrpv2.ave.wide.pre,
  idvar = "DepMap_ID",
  v.names= c("avg"),
  timevar = "cpd_name",
  direction = "wide")

names(ctrpv2.ave.wide) <- gsub("^avg\\.","",names(ctrpv2.ave.wide))

file.ctrpv2.wide = paste0(path.ctrp,"ctrpv2.wide.txt")

write.table(
  x = ctrpv2.ave.wide, 
  file = file.ctrpv2.wide,
  quote = F, col.names=T, row.names = F, sep = "\t")

ctrpv2 <- ctrpv2.ave.wide

## create a culled versio - rename one entry with no DepMap_ID - keep only drugs with data for > 20% of the cell lines
ctrpv3 <-  ctrpv2

## rename one entry with no DepMap_ID
ctrpv3$DepMap_ID <-  ifelse(is.na(ctrpv3$DepMap_ID),"no.depmap.match",ctrpv3$DepMap_ID)
row.names(ctrpv3) <- ctrpv3$DepMap_ID
ctrpv3 <-  ctrpv3[,-1]

## remove drugs with > 80% NAs
percent.nas <- as.data.frame(colMeans(is.na(ctrpv3)) * 100)
names(percent.nas) <- "percent.nas"

## keep only drugs with data for > 20% of the cell lines
percent.nas$eighty.percent.keep.flag <- ifelse(percent.nas$percent.nas > 80, 0, 1)
print("cut drugs")
print(percent.nas[percent.nas$eighty.percent.keep.flag == 0, ])

ctrpv3.culled <- ctrpv3[,names(ctrpv3) %in% row.names(percent.nas[percent.nas$eighty.percent.keep.flag==1,])]

write.table(
  x = ctrpv3.culled,
  file = gsub(".(csv|txt)$","_culled80.\\1",file.ctrpv2.wide),
  quote = F, col.names=T, row.names = T, sep = "\t")

## run imputation
doParallel::registerDoParallel(cores = detectCores() - 2)
doRNG::registerDoRNG(seed = 999)
set.seed(999)

ctrpv3.culled_mf = missForest(
  xmis = ctrpv3.culled,
  parallelize = "variables",
  verbose = T
  )

ctrpv3.culled_mf.imp = ctrpv3.culled_mf$ximp
ctrpv3.culled_mf.imp_t = as.data.frame(
  t(ctrpv3.culled_mf.imp),
  stringsAsFactors = F
  )
    
write.table(
  x = ctrpv3.culled_mf.imp,
  file = gsub(".(csv|txt)$","_culled80_MFImputed.txt",file.ctrpv2.wide),
  col.names = NA
  )
    
write.table(
  x = ctrpv3.culled_mf.imp_t,
  file = gsub(".(csv|txt)$","_culled80_MFImputed_sg.txt",file.ctrpv2.wide),
  quote = F,sep = "\t",
  col.names = NA
  )

##### PLSR: CRISPR & CTRP #####

## set paths
path.dm   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/"
path.ctrp <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/CTRPv2/"

## read in data and set row names
CRISPR <- read.delim(file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"), sep = "\t", stringsAsFactors = F) %>%
  tibble::column_to_rownames(var = "X")
CTRP <- read.delim(file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"), sep = "\t", stringsAsFactors = F) %>%
  tibble::column_to_rownames(var = "X")

## filter for shared cell lines, make matrix (for mixomics), ensure numeric
ids <- intersect(rownames(CRISPR), rownames(CTRP))

X <- CRISPR[ids, , drop = FALSE]
Y <- CTRP[ids, , drop = FALSE]

X[] <- lapply(X, base::as.numeric)
Y[] <- lapply(Y, base::as.numeric)

X <- as.matrix(X)
Y <- as.matrix(Y)

## run PLS
ncomp <- 15

pls_fit <- mixOmics::pls(
  X = X,
  Y = Y,
  ncomp = ncomp,
  scale = TRUE,
  mode  = "regression"    # default
)

## extract from pls_fit object
print(pls_fit$prop_expl_var$X)

x.variates <- data.frame(pls_fit$variates$X)
y.variates <- data.frame(pls_fit$variates$Y)

x.loadings <- data.frame(pls_fit$loadings$X)
y.loadings <- data.frame(pls_fit$loadings$Y)

dim(x.variates);dim(x.loadings)
dim(y.variates);dim(y.loadings)

x.exp_variance <- data.frame(pls_fit$prop_expl_var$X)
y.exp_variance <- data.frame(pls_fit$prop_expl_var$Y)

variates.X <- cbind(Score = rownames(pls_fit$variates$X), x.variates)
variates.Y <- cbind(Score = rownames(pls_fit$variates$Y), y.variates)

loadings.X <- cbind(Loading = rownames(pls_fit$loadings$X), x.loadings)
loadings.Y <- cbind(Loading = rownames(pls_fit$loadings$Y), y.loadings)

rownames(x.exp_variance) = paste0("comp.",seq(1,nrow(x.exp_variance)))

loadings.X.Y = merge(loadings.X,loadings.Y,by="Loading",suffixes = c(".geneexp",".crispr"))
variates.X.Y = merge(variates.X,variates.Y,by="Score",suffixes = c(".geneexp",".crispr"))

# save files
file.gs.gs = "PLSR geneexp vs crispr.txt"

write.table(as.data.frame(variates.X.Y), paste0(gsub(".txt", "", file.gs.gs),"_ncomp",ncomp,"_weight.type-",weight.type, "_PLSR_X.Y.variates.txt"), sep = "\t", row.names = F, quote = F)
write.table(as.data.frame(loadings.X.Y), paste0(gsub(".txt", "", file.gs.gs),"_ncomp",ncomp,"_weight.type-",weight.type, "_PLSR_X.Y.loadings.txt"), sep = "\t", row.names = F, quote = F)

write.table(as.data.frame(x.exp_variance), paste0(gsub(".txt", "", file.gs.gs),"_weight.type-",weight.type, "_PLSR_Xpve.txt"), sep = "\t", row.names = T, quote = F)