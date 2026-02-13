##### Set up #####
if (1) {
  library(BiocManager)
  library(doParallel)
  library(doRNG)
  library(missForest)
  library(dplyr)
  library(mixOmics)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(GO.db)
  library(impute)
  library(preprocessCore)
  library(WGCNA)
  library(ComplexHeatmap)
  library(circlize)
  library(enrichplot)
  library(openxlsx)
  library(doParallel)
  library(GSVA)
}

#### Imputation: CRISPR (NA's in Data) ####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd     <- paste0(path.OS, "WD_FDB_Freeland/")
file.crispr <- paste0(path.wd, "DataSets/DepMap_25Q3/CRISPRGeneEffect.csv")

## Pull in data
CRISPR <- read.delim(
  file = file.crispr,
  row.names = 1,
  stringsAsFactors = F,
  sep = ",",
  check.names = F
)

table(colSums(is.na(CRISPR)))

# 0     1    71    72   228   309   441   524   541   716   769   970   974  
# 9     2    93     9     3     1     2     4     7   155     4     1     6

# 1059  1583  2053  3929  3930 3931  3932  3947  3949  3950  4560  4587  4601
# 1     6     40    60    12   3     7     3     1     1     2     1     1
  
# 4610  4639  5402  5404  5434  6173  6218  8922  8923  8924 
# 2     1     1     1     1     96    11    2     18    57

# 8925  8926  8927  8928  8943  8967 11130 11132 11220 11630 13697 14311 
# 54    13    12     1     1     1     1     1     1     1     1     1 

## Random forrest  
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

## Save 
write.table(
  x = CRISPR_mf.imp, 
  file = gsub(".csv$","_MFImputed.txt", file.crispr), 
  quote = F, 
  sep = "\t",
  col.names = NA
)

write.table(
  x = CRISPR_mf.imp_t, 
  file = gsub(".csv$","_MFImputed_sg.txt", file.crispr), 
  quote = F, 
  sep = "\t",
  col.names = NA
)

#### Imputation: RNAi (NA's in Data) ####

## Set OS (for swapping between personal and workstation)
OS <- "Linux" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd   <- paste0(path.OS, "WD_FDB_Freeland/")
file.rnai <- paste0(path.wd, "DataSets/DepMap_25Q3/D2_combined_gene_dep_scores.csv")

## Pull in data
RNAi <- read.delim(
  file = file.rnai,
  row.names = 1,
  stringsAsFactors = F,
  sep = ",",
  check.names = F
)

table(colSums(is.na(RNAi)))

# 0     1    71    72   228   309   441   524   541   716   769   970   974
# 9     2    93     9     3     1     2     4     7   155     4     1     6

# 1059  1583  2053  3929  3930  3931  3932  3947  3949  3950  4560  4587  4601  
# 1     6     40    60    12    3     7     3     1     1     2     1     1     

# 4610  4639  5402  5404  5434  6173  6218  8922  8923  8924  8925  8926  8927
# 2     1     1     1     1     96    11     2    18    57    54    13    12

# 8928  8943  8967 11130 11132 11220 11630 13697 14311 
# 1     1     1    1     1     1     1     1     1

## Random forrest
doParallel::registerDoParallel(cores = detectCores() - 2)
doRNG::registerDoRNG(seed = 999)
set.seed(999)

RNAi_mf <- missForest(
  xmis = RNAi,
  parallelize = "variables",
  verbose = T
)

RNAi_mf.imp <- RNAi_mf$ximp
RNAi_mf.imp_t <- as.data.frame(
  t(RNAi_mf.imp),
  stringsAsFactors = F
)

## Save 
write.table(
  x = RNAi_mf.imp, 
  file = gsub(".csv$","_MFImputed.txt", file.rnai), 
  quote = F, 
  sep = "\t",
  col.names = NA
)

write.table(
  x = RNAi_mf.imp_t, 
  file = gsub(".csv$","_MFImputed_sg.txt",file.rnai), 
  quote = F, 
  sep = "\t",
  col.names = NA
)

##### Imputation: CTRP (NA's in Data) #####

## Set OS (for swapping between personal and workstation)
OS <- "Linux" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd   <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm   <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp <- paste0(path.wd, "DataSets/CTRPv2/")

## Load in cell line info from depamp
models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F)

## Load in CTRP Data
ctrp.expt   <- read.delim(paste0(path.ctrp,"v20.meta.per_experiment.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
ctrp.cell   <- read.delim(paste0(path.ctrp,"v20.meta.per_cell_line.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
ctrp.inform <- read.delim(paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
ctrp.curves <- read.delim(paste0(path.ctrp,"v20.data.curves_post_qc.txt"), sep = "\t", stringsAsFactors = F, check.names = F)

## Add ID's and name to curves
ctrp.curves$master_ccl_id <- ctrp.expt$master_ccl_id[match(ctrp.curves$experiment_id, ctrp.expt$experiment_id)]
ctrp.curves$ccl_name      <- ctrp.cell$ccl_name[match(ctrp.curves$master_ccl_id, ctrp.cell$master_ccl_id)]
ctrp.curves$DepMap_ID     <- models$ModelID[match(ctrp.curves$ccl_name, models$StrippedCellLineName)]
ctrp.curves$cpd_name      <- ctrp.inform$cpd_name[match(ctrp.curves$master_cpd_id, ctrp.inform$master_cpd_id)]

## Move important columns to front
ctrp.curves <- ctrp.curves %>% 
  dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, area_under_curve, experiment_id, everything())

not.mapped.celllines <- ctrp.curves[is.na(ctrp.curves$DepMap_ID),]
not.mapped.celllines <- not.mapped.celllines[,1:3] %>% dplyr::distinct()

write.table(
  x = ctrp.curves,
  file = paste0(path.ctrp, "ctrp.curves.txt"),
  quote = F, col.names=T, row.names = F, sep = "\t"
  )

## Trim and average data frame
ctrp.curves.abr <- ctrp.curves %>% 
  dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, master_cpd_id, area_under_curve)

ctrpv2.ave <- ctrp.curves.abr %>%
  dplyr::group_by(DepMap_ID, ccl_name,master_ccl_id, cpd_name,master_cpd_id) %>%
  dplyr::summarise(avg = mean(area_under_curve)) %>%
  dplyr::ungroup()

ctrpv2.ave.wide <- ctrpv2.ave %>% 
  dplyr::select(DepMap_ID, ccl_name, cpd_name, avg)
ctrpv2.ave.wide.pre <- as.data.frame(ctrpv2.ave %>% dplyr::select(DepMap_ID, cpd_name, avg))

ctrpv2.ave.wide <- reshape(
  ctrpv2.ave.wide.pre,
  idvar = "DepMap_ID",
  v.names= c("avg"),
  timevar = "cpd_name",
  direction = "wide"
  )

names(ctrpv2.ave.wide) <- gsub("^avg\\.","",names(ctrpv2.ave.wide))

file.ctrpv2.wide = paste0(path.ctrp,"ctrpv2.wide.txt")

write.table(
  x = ctrpv2.ave.wide, 
  file = file.ctrpv2.wide,
  quote = F, col.names=T, row.names = F, sep = "\t"
  )

ctrpv2 <- ctrpv2.ave.wide

## Create a culled version
ctrpv3 <-  ctrpv2

## Rename one entry with no DepMap_ID
ctrpv3$DepMap_ID  <- ifelse(is.na(ctrpv3$DepMap_ID),"no.depmap.match",ctrpv3$DepMap_ID)
row.names(ctrpv3) <- ctrpv3$DepMap_ID
ctrpv3 <- ctrpv3[,-1]

## Remove drugs with > 80% NAs
percent.nas <- as.data.frame(colMeans(is.na(ctrpv3)) * 100)
names(percent.nas) <- "percent.nas"

## Keep only drugs with data for > 20% of the cell lines
percent.nas$eighty.percent.keep.flag <- ifelse(percent.nas$percent.nas > 80, 0, 1)
print("cut drugs")
print(percent.nas[percent.nas$eighty.percent.keep.flag == 0, ])

ctrpv3.culled <- ctrpv3[,names(ctrpv3) %in% row.names(percent.nas[percent.nas$eighty.percent.keep.flag==1,])]

write.table(
  x = ctrpv3.culled,
  file = gsub(".(csv|txt)$","_culled80.\\1",file.ctrpv2.wide),
  quote = F, col.names=T, row.names = T, sep = "\t")

## Run imputation
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

##### PLS: CRISPR & CTRP #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.pls     <- paste0(path.wd, "DataSets/PLS/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Set PLS parameters
X_source <- "CRISPR" # CRISPR or CTRP
Y_source <- "CTRP"   # CRISPR or CTRP

ncomp <- 15
mode  <- "canonical" # default = regression, symmetric = canonical

#### 1. Execute to prep for PLS
if(1) {

  ## For saving files later
  file_tag <- paste0("PLS_Mode.", mode, "_X.", X_source, "_Y.", Y_source)
  
  ## Read in data and set row names
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
    ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
    )
  
  ## Filter for shared cell lines, make matrix (for mixomics), ensure numeric
  if (X_source == "CRISPR") X_data <- CRISPR
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "CRISPR") Y_data <- CRISPR
  if (Y_source == "CTRP")   Y_data <- CTRP
  
  ids <- intersect(rownames(X_data), rownames(Y_data))
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
}

#### 2. Execute to run PLS and save output files (requires Step 1)
if(1){
  
  ## Run PLS
  pls_fit <- mixOmics::pls(
    X = X,
    Y = Y,
    ncomp = ncomp,
    scale = TRUE,
    mode  = mode
  )
  
  print(pls_fit$prop_expl_var$X)
  print(pls_fit$prop_expl_var$Y)
  
  ## Extract from pls_fit object
  x.variates <- data.frame(pls_fit$variates$X) %>%
    tibble::rownames_to_column(var = "Score")
  y.variates <- data.frame(pls_fit$variates$Y) %>%
    tibble::rownames_to_column(var = "Score")
  
  x.loadings <- data.frame(pls_fit$loadings$X) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(comp1)
  y.loadings <- data.frame(pls_fit$loadings$Y) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(comp1)
  
  dim(x.variates);dim(x.loadings)
  dim(y.variates);dim(y.loadings)
  
  x.exp_variance <- data.frame(pls_fit$prop_expl_var$X)
  y.exp_variance <- data.frame(pls_fit$prop_expl_var$Y)
  
  variates.X.Y <- merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = (c(paste0(".", X_source), paste0(".", Y_source)))
  )
  
  ## Save files
  write.table(
    x = x.variates,
    file = paste0(path.pls, file_tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.variates,
    file = paste0(path.pls, file_tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = variates.X.Y,
    file = paste0(path.pls, file_tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.loadings,
    file = paste0(path.pls, file_tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.loadings,
    file = paste0(path.pls, file_tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.exp_variance,
    file = paste0(path.pls, file_tag, "_X.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.exp_variance,
    file = paste0(path.pls, file_tag, "_Y.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

}

#### 3. Execute to plot PLS (requires Step 1)
if(1) {

  ## Load saved loading files
  X_loadings <- read.delim(
    file = paste0(path.pls, file_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- read.delim(
    file = paste0(path.pls, file_tag, "_Y.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )

  ## Bring in raw matrices to compute %NA later
  if (!exists("CRISPR_mat") || !exists("CTRP_mat")) {
    
    CRISPR_mat <- read.delim(
      file = paste0(path.dm, "CRISPRGeneEffect.csv"),
      sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub(" .*", "", .))
    
    CTRP_mat <- read.delim(
      file = paste0(path.ctrp, "ctrpv2.wide.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
  }
  
  ## Helper function for NA-safe pattern detection (useful when labeling for plotting)
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for CTRP loadings file (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    ## Map compound name -> target info
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 MDM2i",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        ),
        group.atp5 = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## Target-category bucketing
    df <- df %>%
      dplyr::mutate(target.category = NA_character_) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA damage"), "DNA.damage", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "(micro|mi)rotubule"), "microtubule", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 1|\\bPLK1\\b"), "PLK1", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 2|\\bPLK2\\b"), "PLK2", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "aurora kinase"), "aurora", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA methyltransferase"), "DNA meth", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA replication"), "DNA rep", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "nicotinamide phosphoribosyltransferase|\\bNAMPT\\b"), "NAMPT", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "dihydrofolate reductase|\\bDHFR\\b"), "DHFR", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "BCL2"), "BCL2.", target.category))
    
    ## %NA per compound using non-imputed CTRP matrix
    percent.nas <- as.data.frame(colMeans(is.na(CTRP_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ### Annotation for CRISPR loadings file 
  annotate_crispr <- function(df, side_label) {
    
    ## Adjusting gene nomenclature
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^MDM[24]$")                                  ~ "04 MDM2.MDM4",
          stringr::str_detect(Loading, "^ATP5")                                      ~ "05 ATP5",
          stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                         ~ "06 dasa targets",
          stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                ~ "07 BCL2+",
          stringr::str_detect(Loading, "^MYC(|N|L)")                                 ~ "08 MYC.",
          stringr::str_detect(Loading, "^(GRB2|CRKL)$")                              ~ "09 SRC-related",
          stringr::str_detect(Loading, "^TP53$")                                     ~ "10 TP53",
          stringr::str_detect(Loading, "^MED12$")                                    ~ "11 MED12",
          TRUE ~ NA_character_
        ),
        group.atp5        = dplyr::if_else(stringr::str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_),
        group.na          = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na     = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na      = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_)
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## %NA using CRISPR matrix
    percent.nas <- as.data.frame(colMeans(is.na(CRISPR_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  
  ## annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_crispr(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_crispr(Y_loadings, "Y")
  
  ## Plotting colors (always plot both sides)
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    # Create alpha flag based on the same thing you use for labeling
    df <- df %>%
      dplyr::mutate(
        label_flag = dplyr::if_else(
          is.na(.data[[color_col]]),
          "Unlabeled", 
          "Labeled"
        )
      )
    
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(
          x = comp1,
          y = comp2,
          color = color_col,
          alpha = "label_flag"
        )
      ) +
        geom_point(size = 2.5) +
        
        geom_text_repel(
          data = df %>% dplyr::filter(label_flag == "Labeled"),
          aes_string(label = label_col),
          size = 2
        ) +
        
        scale_alpha_manual(
          values = c(Labeled = 1, Unlabeled = 0.2),
          guide  = "none"
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  
  ## If side is CTRP: color by target.category, label by Loading
  ## If side is CRISPR: color by group, label by Loading
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  }
}

##### PLS: RNAi & CTRP #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.pls     <- paste0(path.wd, "DataSets/PLS/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Set PLS parameters
X_source <- "CTRP"  # RNAi or CTRP
Y_source <- "RNAi"  # RNAi or CTRP

ncomp <- 15
mode  <- "regression" # default = regression, symmetric = canonical

#### 1. Execute to prep for PLS
if(1) {
  
  ## For saving files later
  file_tag <- paste0("PLS_Mode.", mode, "_X.", X_source, "_Y.", Y_source)
  
  ## Read in data
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  ## Convert RNAi sample nomenclature (CCLEName -> ModelID via Model.csv)
  models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F) %>%
    dplyr::select(ModelID, CCLEName)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## filter for shared cell lines, make matrix (for mixomics), ensure numeric
  if (X_source == "RNAi") X_data <- RNAi_t_ModelID
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP")   Y_data <- CTRP
  
  ids <- intersect(rownames(X_data), rownames(Y_data))
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  }

#### 2. Execute to run PLS and save output files (requires Step 1)
if(1){
  
  ## Run PLS
  pls_fit <- mixOmics::pls(
    X = X,
    Y = Y,
    ncomp = ncomp,
    scale = TRUE,
    mode  = mode
  )
  
  print(pls_fit$prop_expl_var$X)
  print(pls_fit$prop_expl_var$Y)
  
  ## Extract from pls_fit object
  x.variates <- data.frame(pls_fit$variates$X) %>%
    tibble::rownames_to_column(var = "Score")
  y.variates <- data.frame(pls_fit$variates$Y) %>%
    tibble::rownames_to_column(var = "Score")
  
  x.loadings <- data.frame(pls_fit$loadings$X) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(comp1)
  y.loadings <- data.frame(pls_fit$loadings$Y) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(comp1)
  
  dim(x.variates);dim(x.loadings)
  dim(y.variates);dim(y.loadings)
  
  x.exp_variance <- data.frame(pls_fit$prop_expl_var$X)
  y.exp_variance <- data.frame(pls_fit$prop_expl_var$Y)
  
  variates.X.Y <- merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = (c(paste0(".", X_source), paste0(".", Y_source)))
  )
  
  ## Save files
  write.table(
    x = x.variates,
    file = paste0(path.pls, file_tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.variates,
    file = paste0(path.pls, file_tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = variates.X.Y,
    file = paste0(path.pls, file_tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.loadings,
    file = paste0(path.pls, file_tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.loadings,
    file = paste0(path.pls, file_tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.exp_variance,
    file = paste0(path.pls, file_tag, "_X.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.exp_variance,
    file = paste0(path.pls, file_tag, "_Y.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
}

#### 3. Execute to plot PLS (requires Step 1)
if(1) {
  
  ## Load saved loading files
  X_loadings <- read.delim(
    file = paste0(path.pls, file_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- read.delim(
    file = paste0(path.pls, file_tag, "_Y.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Bring in raw matrices to compute %NA later
  if (!exists("RNAi_mat") || !exists("CTRP_mat")) {
    
    RNAi_mat <- read.delim(
      file = paste0(path.dm, "D2_combined_gene_dep_scores.csv"),
      sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub(" .*", "", .))
    
    CTRP_mat <- read.delim(
      file = paste0(path.ctrp, "ctrpv2.wide.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE #, row.names = 1
    )
    
  }
  
  ## Helper function for NA-safe pattern detection (useful when labeling for plotting)
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for CTRP loadings file (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    ## Map compound name -> target info
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 MDM2i",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        ),
        group.atp5 = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## Target-category bucketing
    df <- df %>%
      dplyr::mutate(target.category = NA_character_) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA damage"), "DNA.damage", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "(micro|mi)rotubule"), "microtubule", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 1|\\bPLK1\\b"), "PLK1", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 2|\\bPLK2\\b"), "PLK2", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "aurora kinase"), "aurora", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA methyltransferase"), "DNA meth", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA replication"), "DNA rep", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "nicotinamide phosphoribosyltransferase|\\bNAMPT\\b"), "NAMPT", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "dihydrofolate reductase|\\bDHFR\\b"), "DHFR", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "BCL2"), "BCL2.", target.category))
    
    ## %NA per compound using non-imputed CTRP matrix
    percent.nas <- as.data.frame(colMeans(is.na(CTRP_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ### Annotation for a RNAi loadings frame (gene metadata buckets)
  annotate_rnai <- function(df, side_label) {
    
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^MDM[24]$")                                  ~ "04 MDM2.MDM4",
          stringr::str_detect(Loading, "^ATP5")                                      ~ "05 ATP5",
          stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                         ~ "06 dasa targets",
          stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                ~ "07 BCL2+",
          stringr::str_detect(Loading, "^MYC(|N|L)")                                 ~ "08 MYC.",
          stringr::str_detect(Loading, "^(GRB2|CRKL)$")                              ~ "09 SRC-related",
          stringr::str_detect(Loading, "^TP53$")                                     ~ "10 TP53",
          stringr::str_detect(Loading, "^MED12$")                                    ~ "11 MED12",
          TRUE ~ NA_character_
        ),
        group.atp5     = dplyr::if_else(stringr::str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_),
        group.na       = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na  = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na   = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_)
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## %NA per compound using non-imputed RNAi matrix
    percent.nas <- as.data.frame(colMeans(is.na(RNAi_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ## annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_rnai(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_rnai(Y_loadings, "Y")
  
  ## Plotting colors (always plot both sides)
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    comp_cols <- grep("comp", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(x = comp1, y = comp2, color = color_col)  # <- no label here
      ) +
        geom_point(size = 2.5) +
        
        # Only label rows where color_col is NOT NA
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col),  # inherits x, y, color from main ggplot
          size = 2
        ) +
        
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## If side is CTRP: color by target.category, label by Loading
  ## If side is CRISPR: color by group, label by Loading
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  }
}

##### rCCA: CRISPR & CTRP #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.rcca    <- paste0(path.wd, "DataSets/rCCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Set RCCA parameters
X_source <- "CTRP"    # "CRISPR" or "CTRP"
Y_source <- "CRISPR"  # "CRISPR" or "CTRP"

ncomp <- 15

## Regularization controls
mode_rcca      <- "shrinkage" # ridge (default) requires parameters or tuning, shrinkage

tune_lambda    <- FALSE   # set TRUE to run automatic tuning
lambda1_manual <- 0.20    # default penalty on X (CRISPR side if X_source == "CRISPR")
lambda2_manual <- 0.10    # default penalty on Y

#### 1. Execute to prep for RCCA
if (1) {
  
  ## Read in data and set row names
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  ## Filter for shared cell lines, make matrix (for mixOmics), ensure numeric
  if (X_source == "CRISPR") X_data <- CRISPR
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "CRISPR") Y_data <- CRISPR
  if (Y_source == "CTRP")   Y_data <- CTRP
  
  ids <- intersect(rownames(X_data), rownames(Y_data))
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  ## Choose lambdas: either tuned or manual
  if (mode_rcca == "ridge") {
    
    if (tune_lambda) {
      
      grid1      <- c(0.10, 0.20, 0.30)  # candidate lambdas for X
      grid2      <- c(0.05, 0.10, 0.20)  # candidate lambdas for Y
      ncomp_tune <- min(5L, ncomp)       # tune only first few components
      
      set.seed(999)
      tune_time <- system.time({
        tune.out <- mixOmics::tune.rcc(
          X          = X,
          Y          = Y,
          grid1      = grid1,
          grid2      = grid2,
          ncomp      = ncomp_tune,
          validation = "loo"
        )
      })
      
      print(tune_time)
      print(tune.out$opt.lambda1)
      print(tune.out$opt.lambda2)
      
      lambda1 <- tune.out$opt.lambda1
      lambda2 <- tune.out$opt.lambda2
      
    } else {
      
      ## Use predefined defaults
      lambda1 <- lambda1_manual
      lambda2 <- lambda2_manual
    }
    
    ## Common file tag for ridge mode
    file_tag <- paste0(
      "RCCA_ridge",
      "_lambda1.", format(lambda1, digits = 3),
      "_lambda2.", format(lambda2, digits = 3),
      "_X.", X_source, "_Y.", Y_source
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    ## Shrinkage mode: no lambda parameters
    file_tag <- paste0(
      "RCCA_shrinkage",
      "_X.", X_source, "_Y.", Y_source
    )
    
  }
  
}

#### 2. Execute to run RCCA and save output files (requires Step 1)
if (1) {

  ## Run RCCA with mode-specific call
  if (mode_rcca == "ridge") {
    
    message("Running rCCA in ridge mode with lambda1 = ", lambda1,
            ", lambda2 = ", lambda2)
    
    rcca_fit <- mixOmics::rcc(
      X       = X,
      Y       = Y,
      ncomp   = ncomp,
      lambda1 = lambda1,
      lambda2 = lambda2,
      method  = "ridge"
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    message("Running rCCA in shrinkage mode (automatic lambda estimation).")
    
    rcca_fit <- mixOmics::rcc(
      X      = X,
      Y      = Y,
      ncomp  = ncomp,
      method = "shrinkage"
    )
    
  } else {
    
    stop("mode_rcca must be 'ridge' or 'shrinkage', not: ", mode_rcca)
  }
  
  ## Canonical correlations per component
  print(rcca_fit$cor)
  
  ## Extract from rcca_fit object
  x.variates <- data.frame(rcca_fit$variates$X) %>%
    tibble::rownames_to_column(var = "Score")
  y.variates <- data.frame(rcca_fit$variates$Y) %>%
    tibble::rownames_to_column(var = "Score")
  
  x.loadings <- data.frame(rcca_fit$loadings$X) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(X1)
  y.loadings <- data.frame(rcca_fit$loadings$Y) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(X1)
  
  variates.X.Y <- merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = c(paste0(".", X_source), paste0(".", Y_source))
  )
  
  ## Canonical correlations data.frame
  cancor.df <- data.frame(
    comp                   = seq_along(rcca_fit$cor),
    canonical_correlation  = rcca_fit$cor
  )
  
  ## Save files
  if (!dir.exists(path.rcca)) dir.create(path.rcca, recursive = TRUE)
  
  write.table(
    x = x.variates,
    file = paste0(path.rcca, file_tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.variates,
    file = paste0(path.rcca, file_tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = variates.X.Y,
    file = paste0(path.rcca, file_tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.loadings,
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.loadings,
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = cancor.df,
    file = paste0(path.rcca, file_tag, "_canonical_correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

#### 3. Execute to plot RCCA (requires Step 1)
if(1) {
  
  ## Load saved loading files
  X_loadings <- read.delim(
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- read.delim(
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Bring in raw matrices to compute %NA later
  if (!exists("CRISPR_mat") || !exists("CTRP_mat")) {
    
    CRISPR_mat <- read.delim(
      file = paste0(path.dm, "CRISPRGeneEffect.csv"),
      sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub(" .*", "", .))
    
    CTRP_mat <- read.delim(
      file = paste0(path.ctrp, "ctrpv2.wide.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
  }
  
  ## Helper function for NA-safe pattern detection (useful when labeling for plotting)
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for CTRP loadings file (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    ## Map compound name -> target info
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 MDM2i",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        ),
        group.atp5 = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## Target-category bucketing
    df <- df %>%
      dplyr::mutate(target.category = NA_character_) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA damage"), "DNA.damage", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "(micro|mi)rotubule"), "microtubule", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 1|\\bPLK1\\b"), "PLK1", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 2|\\bPLK2\\b"), "PLK2", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "aurora kinase"), "aurora", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA methyltransferase"), "DNA meth", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA replication"), "DNA rep", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "nicotinamide phosphoribosyltransferase|\\bNAMPT\\b"), "NAMPT", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "dihydrofolate reductase|\\bDHFR\\b"), "DHFR", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "BCL2"), "BCL2.", target.category))
    
    ## %NA per compound using non-imputed CTRP matrix
    percent.nas <- as.data.frame(colMeans(is.na(CTRP_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ### Annotation for CRISPR loadings file 
  annotate_crispr <- function(df, side_label) {
    
    ## Adjusting gene nomenclature
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^MDM[24]$")                                  ~ "04 MDM2.MDM4",
          stringr::str_detect(Loading, "^ATP5")                                      ~ "05 ATP5",
          stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                         ~ "06 dasa targets",
          stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                ~ "07 BCL2+",
          stringr::str_detect(Loading, "^MYC(|N|L)")                                 ~ "08 MYC.",
          stringr::str_detect(Loading, "^(GRB2|CRKL)$")                              ~ "09 SRC-related",
          stringr::str_detect(Loading, "^TP53$")                                     ~ "10 TP53",
          stringr::str_detect(Loading, "^MED12$")                                    ~ "11 MED12",
          TRUE ~ NA_character_
        ),
        group.atp5     = dplyr::if_else(stringr::str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_),
        group.na       = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na  = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na   = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_)
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## %NA using CRISPR matrix
    percent.nas <- as.data.frame(colMeans(is.na(CRISPR_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ## annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_crispr(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_crispr(Y_loadings, "Y")
  
  ## Plotting colors (always plot both sides)
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    comp_cols <- grep("X", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot(
        df,
        aes_string(x = comp1, y = comp2, color = color_col)  # <- no label here
      ) +
        geom_point(size = 2.5) +
        
        # Only label rows where color_col is NOT NA
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col),  # inherits x, y, color from main ggplot
          size = 2
        ) +
        
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## Color by group, label by Loading
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  }
}

##### rCCA: RNAi & CTRP #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.rcca    <- paste0(path.wd, "DataSets/rCCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Set RCCA parameters
X_source <- "RNAi"   # "RNAi" or "CTRP"
Y_source <- "CTRP"   # "RNAi" or "CTRP"

ncomp <- 15

## Regularization controls
mode_rcca      <- "ridge" # ridge (default) requires parameters or tuning, shrinkage

tune_lambda    <- FALSE   # set TRUE to run automatic tuning
lambda1_manual <- 0.20    # penalty on X (RNAi side if X_source == "RNAi")
lambda2_manual <- 0.10    # penalty on Y

#### 1. Execute to prep for RCCA
if (1) {
  
  ## Read in data
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  ## Convert RNAi sample nomenclature (CCLEName -> ModelID via Model.csv)
  models <- read.delim(
    file = paste0(path.dm,"Model.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::select(ModelID, CCLEName)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName")
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## Filter for shared cell lines, make matrix (for mixOmics), ensure numeric
  if (X_source == "RNAi") X_data <- RNAi_t_ModelID
  if (X_source == "CTRP") X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP") Y_data <- CTRP
  
  ids <- intersect(rownames(X_data), rownames(Y_data))
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (mode_rcca == "ridge") {
    
    if (tune_lambda) {
      
      grid1      <- c(0.10, 0.20, 0.30)  # candidate lambdas for X
      grid2      <- c(0.05, 0.10, 0.20)  # candidate lambdas for Y
      ncomp_tune <- min(5L, ncomp)       # tune only first few components
      
      set.seed(999)
      tune_time <- system.time({
        tune.out <- mixOmics::tune.rcc(
          X          = X,
          Y          = Y,
          grid1      = grid1,
          grid2      = grid2,
          ncomp      = ncomp_tune,
          validation = "loo"
        )
      })
      
      print(tune_time)
      print(tune.out$opt.lambda1)
      print(tune.out$opt.lambda2)
      
      lambda1 <- tune.out$opt.lambda1
      lambda2 <- tune.out$opt.lambda2
      
    } else {
      
      ## Use predefined defaults
      lambda1 <- lambda1_manual
      lambda2 <- lambda2_manual
    }
    
    ## Common file tag for ridge mode
    file_tag <- paste0(
      "RCCA_ridge",
      "_lambda1.", format(lambda1, digits = 3),
      "_lambda2.", format(lambda2, digits = 3),
      "_X.", X_source, "_Y.", Y_source
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    ## Shrinkage mode: no lambda parameters
    file_tag <- paste0(
      "RCCA_shrinkage",
      "_X.", X_source, "_Y.", Y_source
    )
    
  }
}

#### 2. Execute to run RCCA and save output files (requires Step 1)
if (1) {
  
  if (mode_rcca == "ridge") {
    
    message("Running rCCA in ridge mode with lambda1 = ", lambda1,
            ", lambda2 = ", lambda2)
    
    rcca_fit <- mixOmics::rcc(
      X       = X,
      Y       = Y,
      ncomp   = ncomp,
      lambda1 = lambda1,
      lambda2 = lambda2,
      method  = "ridge"
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    message("Running rCCA in shrinkage mode (automatic lambda estimation).")
    
    rcca_fit <- mixOmics::rcc(
      X      = X,
      Y      = Y,
      ncomp  = ncomp,
      method = "shrinkage"
    )
    
  } else {
    
    stop("mode_rcca must be 'ridge' or 'shrinkage', not: ", mode_rcca)
  }
  
  ## Canonical correlations per component (full spectrum; first ncomp used)
  print(rcca_fit$cor[1:ncomp])

  ## Extract from rcca_fit object
  x.variates <- data.frame(rcca_fit$variates$X) %>%
    tibble::rownames_to_column(var = "Score")
  y.variates <- data.frame(rcca_fit$variates$Y) %>%
    tibble::rownames_to_column(var = "Score")
  
  x.loadings <- data.frame(rcca_fit$loadings$X) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(X1)
  y.loadings <- data.frame(rcca_fit$loadings$Y) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(X1)
  
  variates.X.Y <- merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = c(paste0(".", X_source), paste0(".", Y_source))
  )
  
  ## Canonical correlations data.frame
  cancor.df <- data.frame(
    comp                   = seq_along(rcca_fit$cor),
    canonical_correlation  = rcca_fit$cor
  )
  
  ## Save files
  if (!dir.exists(path.rcca)) dir.create(path.rcca, recursive = TRUE)
  
  write.table(
    x = x.variates,
    file = paste0(path.rcca, file_tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.variates,
    file = paste0(path.rcca, file_tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = variates.X.Y,
    file = paste0(path.rcca, file_tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.loadings,
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.loadings,
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = cancor.df,
    file = paste0(path.rcca, file_tag, "_canonical_correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

#### 3. Execute to plot RCCA (requires Step 2)
if (1) {
  
  ## Load saved loading files
  X_loadings <- read.delim(
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- read.delim(
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
 
  ## Bring in raw matrices to compute %NA later
  if (!exists("RNAi_mat") || !exists("CTRP_mat")) {
    
    RNAi_mat <- read.delim(
      file = paste0(path.dm, "D2_combined_gene_dep_scores.csv"),
      sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub(" .*", "", .))
    
    CTRP_mat <- read.delim(
      file = paste0(path.ctrp, "ctrpv2.wide.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE #, row.names = 1
    )
    
  }
  
  ## Helper function for NA-safe pattern detection (useful when labeling for plotting)
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for CTRP loadings file (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    ## Map compound name -> target info
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 MDM2i",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        ),
        group.atp5 = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## Target-category bucketing
    df <- df %>%
      dplyr::mutate(target.category = NA_character_) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA damage"), "DNA.damage", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "(micro|mi)rotubule"), "microtubule", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 1|\\bPLK1\\b"), "PLK1", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 2|\\bPLK2\\b"), "PLK2", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "aurora kinase"), "aurora", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA methyltransferase"), "DNA meth", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA replication"), "DNA rep", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "nicotinamide phosphoribosyltransferase|\\bNAMPT\\b"), "NAMPT", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "dihydrofolate reductase|\\bDHFR\\b"), "DHFR", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "BCL2"), "BCL2.", target.category))
    
    ## %NA per compound using non-imputed CTRP matrix
    percent.nas <- as.data.frame(colMeans(is.na(CTRP_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ### Annotation for a RNAi loadings frame (gene metadata buckets)
  annotate_rnai <- function(df, side_label) {
    
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^MDM[24]$")                                  ~ "04 MDM2.MDM4",
          stringr::str_detect(Loading, "^ATP5")                                      ~ "05 ATP5",
          stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                         ~ "06 dasa targets",
          stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                ~ "07 BCL2+",
          stringr::str_detect(Loading, "^MYC(|N|L)")                                 ~ "08 MYC.",
          stringr::str_detect(Loading, "^(GRB2|CRKL)$")                              ~ "09 SRC-related",
          stringr::str_detect(Loading, "^TP53$")                                     ~ "10 TP53",
          stringr::str_detect(Loading, "^MED12$")                                    ~ "11 MED12",
          TRUE ~ NA_character_
        ),
        group.atp5       = dplyr::if_else(stringr::str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_),
        group.na         = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na    = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na     = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_)
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
    ## %NA using RNAi matrix
    percent.nas <- as.data.frame(colMeans(is.na(RNAi_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ## annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_rnai(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_rnai(Y_loadings, "Y")
  
  ## Plotting colors (always plot both sides)
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    comp_cols <- grep("X", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot(
        df,
        aes_string(x = comp1, y = comp2, color = color_col)  # <- no label here
      ) +
        geom_point(size = 2.5) +
        
        # Only label rows where color_col is NOT NA
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col),  # inherits x, y, color from main ggplot
          size = 2
        ) +
        
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## Color by group, label by Loading
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  }
}

##### Max loading: Scatter & GSEA #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.pls     <- paste0(path.wd, "DataSets/PLS/")
path.rcca    <- paste0(path.wd, "DataSets/rCCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.max     <- paste0(path.wd, "DataSets/MaxLoading/")
path.scripts <- paste0(path.OS, "FDB_Freeland/Scripts/")

source("/Users/jack/Documents/GitHub/FDB_Freeland/Scripts/FGSEA_functions.R")

## Set dim red technique. PLS, rCCA
DimRedTec <- "PLS"

## Set parameters. CRISPR, RNAi, CTRP
X1_source <- "CRISPR"
Y1_source <- "CTRP"

X2_source <- "RNAi"
Y2_source <- "CTRP"

mode  <- "canonical" # default = regression, symmetric = canonical

### Create scatter plot and generate table of distances and theta
if(1){
  
  if(DimRedTec == "PLS"){
    ## file tag
    file1_tag <- paste0("PLS_Mode.", mode, "_X.", X1_source, "_Y.", Y1_source)
    file2_tag <- paste0("PLS_Mode.", mode, "_X.", X2_source, "_Y.", Y2_source)
  }
  
  if(DimRedTec == "rCCA"){
    ## file tag
    file1_tag <- paste0("rCCA_X.", X1_source, "_Y.", Y1_source)
    file2_tag <- paste0("rCCA_X.", X2_source, "_Y.", Y2_source)
  }
  
  ## read in data
  X1_loadings <- read.delim(
    file = paste0(path.pls, file1_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading)) %>%
    dplyr::select(Loading, paste0("comp", 1:10))
  
  X2_loadings <- read.delim(
    file = paste0(path.pls, file2_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading)) %>%
    dplyr::select(Loading, paste0("comp", 1:10))
  
  ## find max abs()
  X1_loadings_max <- X1_loadings %>%
    tidyr::pivot_longer(
      cols = paste0("comp", 1:10),
      names_to  = "component",
      values_to = "loading"
    ) %>%
    dplyr::mutate(abs_loading = abs(loading)) %>%
    dplyr::group_by(Loading) %>%
    dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      component_CRISPR    = component,
      loading_CRISPR      = loading,
      abs_loading_CRISPR  = abs_loading
    )
  
  X2_loadings_max <- X2_loadings %>%
    tidyr::pivot_longer(
      cols = paste0("comp", 1:10),
      names_to  = "component",
      values_to = "loading"
    ) %>%
    dplyr::mutate(abs_loading = abs(loading)) %>%
    dplyr::group_by(Loading) %>%
    dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      component_RNAi    = component,
      loading_RNAi      = loading,
      abs_loading_RNAi  = abs_loading
    )
  
  ## merge
  Max <- merge(X1_loadings_max, X2_loadings_max, by = "Loading")
  
  xcol <- names(Max)[7]   # 7th column name
  ycol <- names(Max)[4]   # 4th column name
  
  Max <- Max %>%
    dplyr::mutate(
      theta_rad = atan2(.data[[ycol]], .data[[xcol]]),
      theta_deg = theta_rad * 180 / pi,
      r         = sqrt(.data[[xcol]]^2 + .data[[ycol]]^2)
    )
  
  write.table(
    x = Max,
    file = paste0(path.max, "MaxLoadingsDF_", mode, "_X1_", X1_source, "_vs_", Y1_source, "_X2_", X2_source, "_vs_", Y2_source, ".txt"),
    quote = F, sep = "\t", col.names = T, row.names = F)
  
  ## plot 
  plot_df <- data.frame(
    x      = Max[[7]],   # 4th column = abs_loading_CRISPR
    y      = Max[[4]],   # 7th column = abs_loading_RNAi
    gene   = Max$Loading
  )
  
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(size = 0.075, alpha = 0.3) +
    geom_abline(
      slope = tan(pi/6), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    geom_abline(
      slope = tan(pi/3), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    # theme_bw(base_size = 10) +
    labs(
      x = names(Max)[7],
      y = names(Max)[4],
      title = "Max Absolute Loadings per Gene comp 1-10"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    theme_classic(base_size = 10)
  
  ggsave(
    filename = paste0(path.plots, "MaxLoadingsDF_", mode, "_X1_", X1_source, "_vs_", Y1_source, "_X2_", X2_source, "_vs_", Y2_source, "_Scatter.pdf"),
    plot = p, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
}

#### Create GSEA Plot

## Load msigdb pathways
msig_df <- load.MSigDB(species = 'Homo sapiens')

gsea_list <- get.MSigDB.genesets(
  msig_df = rbind(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "C5")
  ),
  genesets = c("BP")
)

keyword_groups <- list(
  # OX =
  #   c("OXIDATIVE_PHOSPHORYLATION", "ELECTRON_TRANSPORT_CHAIN", "MITOCHONDRIAL_COMPLEX",
  #     "NADH_DEHYDROGENASE", "MITOCHONDRIAL_LARGE_RIBOSOMAL"),
  # DNA =
  #   c("DNA_REPAIR", "DNA_DAMAGE_STIMULUS", "FANCONI", "END_JOINING"),
  # NEURO = 
  #   c("NEURO", "NEUROTRANSMITTER", "SYNAPTIC", "VOLTAGE", "AXON",
  #     "CEREBRAL", "CORTEX", "DENDRITE", "GLUTAMATE"),
  IMMUNE =
    c("INFLAME", "IMMUNE", "INTERLEUKIN", "LEUKOCYTE", "CD4",
      "MACROPHAGE", "NEUTROPHILE"),
  # TISSUE_DEVELOPMENT =
  #   c("MORPHOGENESIS", "VESSEL", "TISSUE_DEVELOPMENT", "TISSUE"),
  # SENSORY_PERCEPTION =
  #   c("SENSORY", "AUDITORY", "SMELL"),
  # KINASE_ACTIVITY =
  #   c("MAPK", "KINASE", "GTP", "TYROSINE"),
  # CELL_DIFFERENTIATION =
  #   c("KERATINOCYTE", "DIFFERENTIATION"),
  # CELL_CELL_INTERACTION =
  #   c("ADHESION", "ADHERENS", "CELL-CELL", "COMMUNICATION"),
  # PEROXISOME =
  #   c("PEROXIDE", "PEROXISOME"),
  # CELL_PROLIFERATION =
  #   c("PROLIFERATION"),
  PROTEIN_PROCESSING =
    c("PEPTIDE", "AMINO_ACID", "UBIQUITIN", "UBIQUITINATION"),
  VIRAL_PROCESSES =
    c("VIRAL", "SYMBIOTIC", "DSRNA"),
  # MICRO_RNA =
  #   c("MIRNA"),
  STRESS_RESPONSE =
    c("DNA_DAMAGE", "APOPTOTIC", "REPAIR", "HYPOXIA", "STRESS"),
  METABOLIC_PATHWAY =
    c("CATABOLIC", "ATP", "POLYSACCHARIDE", "FRUCTOSE",
      "GLYCOSYLATION", "GLYCOGEN", "BIOSYNTHESIS", "LIPID"),
  MITOCHONDRIA =
    c("MITOCHONDRIAL", "MITOCHONDRION"),
  # ORGANELLE_TRANSPORT =
  #   c("ENDOPLASMIC_RETICULUM", "GOLGI", "VACUOLE"),
  # ORGANELLE_ORGANIZATION =
  #   c("ORGANELLE"),
  # EPIGENETIC =
  #   c("HISTONE", "NUCLEOSIDE", "DEMETHYLATION", "METHYLATION", "EPIGENETIC"),
  # SPLICEOSOME =
  #   c("SNRNA", "SPLICING", "SPLICEOSOME"),
  TRANSLATION =
    c("RIBOSOME", "RRNA", "TRNA", "TRANSLATION", "RIBONUCLEOPROTEIN"),
  DNA_TRANSCRIPTION =
    c("MRNA", "TRANSCRIPTION", "POLYMERASE", "TRANSCRIBED", "GENE_EXPRESSION"),
  CELL_CYCLE =
    c("CELL_CYCLE", "MITOTIC", "DNA_REPLICATION", "CHROMOSOME_SEGREGATION",
      "CHROMATID_SEGREGATION", "SPINDLE", "CELL_DIVISION",
      "KINETOCHORE", "CENTRIOLE", "ANAPHASE")
)

pathway_names <- names(gsea_list)

## For each keyword group, find all pathways whose name contains ANY of the strings
keyword_to_genes <- purrr::imap(
  keyword_groups,
  function(pattern_vec, kw_name) {
    
    # single regex that ORs all patterns in the keyword group
    pattern_regex <- paste(pattern_vec, collapse = "|")
    
    # which pathways match any of these patterns (case-insensitive)
    hit_idx <- grepl(pattern_regex, pathway_names, ignore.case = TRUE)
    
    # union of all genes in those matched pathways
    genes <- unique(unlist(gsea_list[hit_idx], use.names = FALSE))
    
    genes
  }
)

sapply(keyword_to_genes, length)

## Long data frame: one row per (gene, keyword_group)
keyword_gene_df <- purrr::imap_dfr(
  keyword_to_genes,
  ~ dplyr::tibble(
    Loading       = .x,   # gene symbol
    keyword_group = .y    # name of the list element (e.g. "OX", "DNA", ...)
  )
)

## Make sure groups have a fixed order
keyword_gene_df$keyword_group <- factor(
  keyword_gene_df$keyword_group,
  levels = names(keyword_groups)
)

## Join to Max on the gene column "Loading"
Max_kw <- Max %>%
  dplyr::inner_join(keyword_gene_df, by = "Loading") %>%
  dplyr::mutate(
    keyword_group = factor(keyword_group, levels = names(keyword_groups))
  )

group_order <- Max_kw %>%
  dplyr::group_by(keyword_group) %>%
  dplyr::summarise(mean_theta = mean(theta_deg, na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(mean_theta)) %>%
  dplyr::pull(keyword_group)

Max_kw$keyword_group <- factor(Max_kw$keyword_group, levels = group_order)

## Set custom colors
my_group_colors <- c(
  IMMUNE              = "#fea605",
  PROTEIN_PROCESSING  = "#0f34fe",
  VIRAL_PROCESSES     = "#fefb39",
  METABOLIC_PATHWAY   = "#006500",
  MITOCHONDRIA        = "#5fe2d1",
  TRANSLATION         = "#aa3337",
  DNA_TRANSCRIPTION   = "#fd2600",
  CELL_CYCLE          = "#fec0cc"
)

## Plot and save
p <- ggplot2::ggplot(
  Max_kw,
  aes(x = r, y = theta_deg, color = keyword_group)
) +
  ## density wave contours — colored per group
  geom_density_2d(
    aes(group = keyword_group),   # ensures independent estimation per facet
    linewidth = 0.5,
    alpha = 0.8
  ) +
  geom_point(size = 0.8, alpha = 0.6) +
  ## reference lines
  geom_hline(
    yintercept = c(30, 60),
    linetype   = "dotted",
    color      = "grey40",
    linewidth  = 0.4
  ) +
  ## facets
  facet_grid(
    . ~ keyword_group,
    scales = "free_x",
    space  = "free_x"
  ) +
  ## apply your palette
  scale_color_manual(values = my_group_colors, guide = "none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 0.07))

print(p)

ggsave(
  filename = paste0(path.plots, "MaxLoadingsDF_", mode, "_X1_", X1_source, "_vs_", Y1_source, "_X2_", X2_source, "_vs_", Y2_source, "_GSEA.pdf"),
  plot = p, width = 13, height = 9, units = "in", device = cairo_pdf
)

##### GSEA #####

## Set OS (for swapping between personal and workstation)
OS <- "Linux" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.max     <- paste0(path.wd, "DataSets/MaxLoading/")
path.scripts <- "/Users/jack/Documents/GitHub/FDB_Freeland/Scripts/"

source(paste0(path.scripts, "FGSEA_functions.R"))

### Load msigdb pathways
msig_df <- load.MSigDB(species = 'Homo sapiens')

gsea_list <- get.MSigDB.genesets(
  msig_df = rbind(
    msigdbr::msigdbr(species = "Homo sapiens", category = "C2"),
    msigdbr::msigdbr(species = "Homo sapiens", category = "C5"),
    msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  ), # restart R session if error (C2, C5, H)
  genesets = c("CP", "GO", "H$")
)

# gsea_list <- get.MSigDB.genesets(
#   msig_df = rbind(
#     msigdbr(species = "Homo sapiens", category = "H"),
#     msigdbr(species = "Homo sapiens", category = "C1"),
#     msigdbr(species = "Homo sapiens", category = "C2"),
#     msigdbr(species = "Homo sapiens", category = "C3"),
#     msigdbr(species = "Homo sapiens", category = "C4"),
#     msigdbr(species = "Homo sapiens", category = "C5"),
#     msigdbr(species = "Homo sapiens", category = "C6"),
#     msigdbr(species = "Homo sapiens", category = "C7"),
#     msigdbr(species = "Homo sapiens", category = "C8"),
#     msigdbr(species = "Homo sapiens", category = "C8")
#   ),
#   genesets = c()
# )

### Read input data
input.path <- paste0(path.max, "MaxLoadingsDF_canonical_X1_CRISPR_vs_CTRP_X2_RNAi_vs_CTRP.txt")

rankings <- read.delim(
  file = input.path,
  stringsAsFactors = F,
  sep = "\t",
  check.names = F
) %>%
  tibble::column_to_rownames(var = "Loading")

rnk <- get.rnk.vector(DE_results = rankings, column_name = "theta_deg")

### Run GSEA
FGSEA_results <- run.FGSEA(rnk, gsea_list, nproc = 2, minGenes = 3, maxGenes = 5000, reformat = T, filename = F, minP = 1e-20)

FGSEA_results_rnk <- FGSEA_results %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(rank = 1:n())

### Save file
write.table(FGSEA_results_rnk,
            file = paste0(path.max, "FGSEA_", basename(input.path)),
            sep = "\t", quote = F,
            row.names = F)

##### GSEA^2 #####

## Set OS (for swapping between personal and workstation)
OS <- "Linux" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd    <- paste0(path.OS, "WD_FDB_Freeland/")
path.max   <- paste0(path.wd, "DataSets/MaxLoading/")
path.plots <- paste0(path.wd, "Plots/")

## Load input data
path.input <- paste0(
  path.max,
  "FGSEA_MaxLoadingsDF_canonical_X1_CRISPR_vs_CTRP_X2_RNAi_vs_CTRP.txt"
)

FGSEA_results_rnk <- read.table(
  file   = path.input,
  sep    = "\t",
  header = TRUE
)

### Prep for GSEA^2 

# Define keyword groups
keyword_groups <- list(
  OX =
    c("OXIDATIVE_PHOSPHORYLATION", "ELECTRON_TRANSPORT_CHAIN", "MITOCHONDRIAL_COMPLEX",
      "NADH_DEHYDROGENASE", "MITOCHONDRIAL_LARGE_RIBOSOMAL"),
  DNA =
    c("DNA_REPAIR", "DNA_DAMAGE_STIMULUS", "FANCONI", "END_JOINING"),
  NEURO = 
    c("NEURO", "NEUROTRANSMITTER", "SYNAPTIC", "VOLTAGE", "AXON",
      "CEREBRAL", "CORTEX", "DENDRITE", "GLUTAMATE"),
  IMMUNE =
    c("INFLAME", "IMMUNE", "INTERLEUKIN", "LEUKOCYTE", "CD4",
      "MACROPHAGE", "NEUTROPHILE"),
  TISSUE_DEVELOPMENT =
    c("MORPHOGENESIS", "VESSEL", "TISSUE_DEVELOPMENT", "TISSUE"),
  SENSORY_PERCEPTION =
    c("SENSORY", "AUDITORY", "SMELL"),
  KINASE_ACTIVITY =
    c("MAPK", "KINASE", "GTP", "TYROSINE"),
  CELL_DIFFERENTIATION =
    c("KERATINOCYTE", "DIFFERENTIATION"),
  CELL_CELL_INTERACTION =
    c("ADHESION", "ADHERENS", "CELL-CELL", "COMMUNICATION"),
  PEROXISOME =
    c("PEROXIDE", "PEROXISOME"),
  CELL_PROLIFERATION =
    c("PROLIFERATION"),
  PROTEIN_PROCESSING =
    c("PEPTIDE", "AMINO_ACID", "UBIQUITIN", "UBIQUITINATION"),
  VIRAL_PROCESSES =
    c("VIRAL", "SYMBIOTIC", "DSRNA"),
  MICRO_RNA =
    c("MIRNA"),
  STRESS_RESPONSE =
    c("DNA_DAMAGE", "APOPTOTIC", "REPAIR", "HYPOXIA", "STRESS"),
  METABOLIC_PATHWAY =
    c("CATABOLIC", "ATP", "POLYSACCHARIDE", "FRUCTOSE",
      "GLYCOSYLATION", "GLYCOGEN", "BIOSYNTHESIS", "LIPID"),
  MITOCHONDRIA =
    c("MITOCHONDRIAL", "MITOCHONDRION"),
  ORGANELLE_TRANSPORT =
    c("ENDOPLASMIC_RETICULUM", "GOLGI", "VACUOLE"),
  ORGANELLE_ORGANIZATION =
    c("ORGANELLE"),
  EPIGENETIC =
    c("HISTONE", "NUCLEOSIDE", "DEMETHYLATION", "METHYLATION", "EPIGENETIC"),
  SPLICEOSOME =
    c("SNRNA", "SPLICING", "SPLICEOSOME"),
  TRANSLATION =
    c("RIBOSOME", "RRNA", "TRNA", "TRANSLATION", "RIBONUCLEOPROTEIN"),
  DNA_TRANSCRIPTION =
    c("MRNA", "TRANSCRIPTION", "POLYMERASE", "TRANSCRIBED", "GENE_EXPRESSION"),
  CELL_CYCLE =
    c("CELL_CYCLE", "MITOTIC", "DNA_REPLICATION", "CHROMOSOME_SEGREGATION",
      "CHROMATID_SEGREGATION", "SPINDLE", "CELL_DIVISION",
      "KINETOCHORE", "CENTRIOLE", "ANAPHASE")
)

# Get ranks per keyword group 
get_ranks_for_keywords <- function(results_df, keywords) {
  pattern <- paste(keywords, collapse = "|")  # OR across all keywords
  results_df %>%
    dplyr::filter(stringr::str_detect(pathway, pattern)) %>%
    dplyr::pull(rank)
}

# Build GSEA^2 data frame
data_df <- purrr::imap_dfr(
  keyword_groups,
  ~ tibble::tibble(
    Category = .y,
    Value    = get_ranks_for_keywords(FGSEA_results_rnk, .x)
  )
)

print(table(data_df$Category))

# KS Test (enrichment + signed ordering)
max_rank <- max(FGSEA_results_rnk$rank, na.rm = TRUE)

ks_results <- data_df %>%
  dplyr::group_by(Category) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    mean_rank = mean(Value, na.rm = TRUE),
    p_ks      = if (n > 1) {
      # Scale ranks to [0,1] and test vs Uniform(0,1)
      scaled_vals <- Value / max_rank
      stats::ks.test(scaled_vals, "punif")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # '+' = enriched toward low ranks (top of GSEA ranking)
    # '-' = enriched toward high ranks (bottom of GSEA ranking)
    direction = dplyr::if_else(mean_rank <= max_rank / 2, "+", "-"),
    # signed significance: large negative = strong '-' enrichment,
    # large positive      = strong '+' enrichment
    logp        = dplyr::if_else(is.na(p_ks), NA_real_, -log10(p_ks)),
    signed_logp = dplyr::if_else(direction == "+", logp, -logp),
    # significance stars
    sig_flag = dplyr::case_when(
      is.na(p_ks)        ~ "",
      p_ks <= 0.0001     ~ "****",
      p_ks <= 0.001      ~ "***",
      p_ks <= 0.01       ~ "**",
      p_ks <= 0.05       ~ "*",
      TRUE               ~ ""
    ),
    label_suffix = dplyr::if_else(
      is.na(p_ks),
      "",
      paste0(
        " p = ",
        format(p_ks, scientific = TRUE, digits = 3),
        " (", direction, ") ",
        sig_flag
      )
    )
  ) %>%
  dplyr::arrange(signed_logp)

print(ks_results)

# Use this ordering for the factor, so along the axis you go:
# most significant '-'  --> weaker '-' --> weaker '+' --> most significant '+'
ordered_levels <- ks_results$Category

data_df$Category <- factor(
  data_df$Category,
  levels = ordered_levels
)

# Left-side category labels (edit right-hand side later if you like)
custom_labels <- c(
  OX                     = "OX",
  DNA                    = "DNA",
  NEURO                  = "Neuro",
  IMMUNE                 = "Immune",
  TISSUE_DEVELOPMENT     = "Tissue Development",
  SENSORY_PERCEPTION     = "Sensory Perception",
  KINASE_ACTIVITY        = "Kinase Activity",
  CELL_DIFFERENTIATION   = "Cell Differentiation",
  CELL_CELL_INTERACTION  = "Cell–Cell Interaction",
  PEROXISOME             = "Peroxisome",
  CELL_PROLIFERATION     = "Cell Proliferation",
  PROTEIN_PROCESSING     = "Protein Processing",
  VIRAL_PROCESSES        = "Viral Processes",
  MICRO_RNA              = "Micro RNA",
  STRESS_RESPONSE        = "Stress Response",
  METABOLIC_PATHWAY      = "Metabolic Pathway",
  MITOCHONDRIA           = "Mitochondria",
  ORGANELLE_TRANSPORT    = "Organelle Transport",
  ORGANELLE_ORGANIZATION = "Organelle Organization",
  EPIGENETIC             = "Epigenetic",
  SPLICEOSOME            = "Spliceosome",
  TRANSLATION            = "Translation",
  DNA_TRANSCRIPTION      = "DNA Transcription",
  CELL_CYCLE             = "Cell Cycle"
)

# Right-side labels: p-value + direction + stars
right_labels <- ks_results %>%
  dplyr::mutate(
    right_lab = 
      # paste0(
      # "p = ", format(p_ks, scientific = TRUE, digits = 3),
      # " (", direction, ") ", sig_flag)
    paste0(
      " p = ",
      signif(p_ks, 3),
      " (", direction, ") ", sig_flag
    )
  ) %>%
  dplyr::select(Category, right_lab) %>%
  tibble::deframe()

# Plot (x = rank, y = category; labels on left and right)
plt <- ggplot(data_df, ggplot2::aes(x = Value, y = Category)) +
  geom_jitter(
    height = 0.2,
    width  = 0,
    aes(color = Category),
    size   = 1,
    shape  = 16
  ) +
  scale_y_discrete(
    labels   = custom_labels,
    sec.axis = dup_axis(
      labels = right_labels[levels(data_df$Category)],
      name   = ""
    )
  ) +
  scale_x_continuous(
    breaks = c(max_rank * 0.15, max_rank * 0.85),
    labels = c("Enriched in\nRNAi", "Enriched in\nCRISPR")
  ) +
  labs(x = "Rank", y = "") +
  theme_minimal() +
  theme(
    axis.text.y.left  = element_text(size = 7),
    axis.text.y.right = element_text(size = 7, hjust = 0),
    axis.text.x       = element_text(size = 7),
    legend.position   = "none",
    panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank()
  )

ggsave(
  filename = paste0(
    path.plots,
    "GSEA_sq_",
    gsub(".txt", ".pdf", basename(path.input))
  ),
  plot   = plt,
  width  = 6,
  height = 4,
  units  = "in",
  device = cairo_pdf
)

##### Investigating drug data for better classification #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.general <- paste0(path.wd, "DataSets/General/")

ctrp.inform  <- read.delim(paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F, check.names = F)

View(table(ctrp.inform$target_or_activity_of_compound))
##### Group distance diagnostics from loadings #####
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.pls     <- paste0(path.wd, "DataSets/PLS/")
path.rcca    <- paste0(path.wd, "DataSets/rCCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Which side of the PLS to analyse
side_to_use <- "X"        # "X" or "Y"

## What type of data this side corresponds to
source_type <- "CTRP"   # "CRISPR" or "CTRP"

## Path to the loadings file you want to analyse
loadings_path <- paste0(path.pls, "PLS_Mode.canonical_X.CRISPR_Y.CTRP_Y.loadings.txt")

if (1) {

  ## Which column encodes the groups you used for color. CRISPR = group, CTRP = target.category
  color_col <- if (source_type == "CTRP") "target.category" else "group"
  
  ## Read and annotate loadings
  loadings_raw <- read.delim(
    file = loadings_path,
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    ## Map compound name -> target info
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 MDM2i",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        )
      )
    
    ## Target-category bucketing
    df <- df %>%
      dplyr::mutate(target.category = NA_character_) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA damage"), "DNA.damage", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "(micro|mi)rotubule"), "microtubule", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 1|\\bPLK1\\b"), "PLK1", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "polo\\-like kinase 2|\\bPLK2\\b"), "PLK2", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "aurora kinase"), "aurora", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA methyltransferase"), "DNA meth", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "DNA replication"), "DNA rep", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "nicotinamide phosphoribosyltransferase|\\bNAMPT\\b"), "NAMPT", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "dihydrofolate reductase|\\bDHFR\\b"), "DHFR", target.category)) %>%
      dplyr::mutate(target.category = dplyr::if_else(detect(drug.target, "BCL2"), "BCL2.", target.category))
    
    df
  }
  
  ### Annotation for CRISPR loadings file 
  annotate_crispr <- function(df, side_label) {
    
    ## Adjusting gene nomenclature
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    ## Groupings
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^MDM[24]$")                                  ~ "04 MDM2.MDM4",
          stringr::str_detect(Loading, "^ATP5")                                      ~ "05 ATP5",
          stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                         ~ "06 dasa targets",
          stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                ~ "07 BCL2+",
          stringr::str_detect(Loading, "^MYC(|N|L)")                                 ~ "08 MYC.",
          stringr::str_detect(Loading, "^(GRB2|CRKL)$")                              ~ "09 SRC-related",
          stringr::str_detect(Loading, "^TP53$")                                     ~ "10 TP53",
          stringr::str_detect(Loading, "^MED12$")                                    ~ "11 MED12",
          TRUE ~ NA_character_
        )
      )
    
    df
  }
  
  ## Re-use your existing annotation helpers to bring in group info
  loadings_annot <- if (source_type == "CTRP") {
    annotate_ctrp(loadings_raw, side_to_use)
  } else {
    annotate_crispr(loadings_raw, side_to_use)
  }
  
  ## Keep only rows with a defined group (for distance summaries)
  df_use <- loadings_annot %>%
    dplyr::filter(!is.na(.data[["group"]]))
  
  ## 1D group distance metrics along each component separately
  compute_group_distances <- function(df, group) {
    
    ## Find component columns (comp1, comp2, ...)
    comp_cols <- grep("comp", names(df), value = TRUE)

    results <- list()
    
    for (comp in comp_cols) {
      
      tmp <- df %>%
        dplyr::filter(!is.na(.data[[group]])) %>%
        dplyr::mutate(val = .data[[comp]]) %>%  # 1D values on this component
        dplyr::group_by(.data[[group]]) %>%
        dplyr::summarise(
          center          = mean(val, na.rm = TRUE),
          dist_origin     = abs(center),
          mean_dist_center = mean(abs(val - mean(val, na.rm = TRUE)), na.rm = TRUE),
          n               = dplyr::n(),
          .groups         = "drop"
        ) %>%
        dplyr::mutate(component = comp)
      
      ## Standardize group column name
      names(tmp)[names(tmp) == group] <- "group_var"
      results[[length(results) + 1]] <- tmp
    }
    
    dplyr::bind_rows(results)
  }
  
  ## Compute distances
  group_distances <- compute_group_distances(df_use, "group")
  
  ## Ensure component facets appear in numeric order (comp1, comp2, comp3, ...)
  group_distances$component <- factor(
    group_distances$component,
    levels = paste0("comp", sort(unique(as.numeric(gsub("comp", "", group_distances$component)))))
  )
  
  ## save
  write.table(
    x = group_distances, 
    file = paste0(path.pls, file_tag, "_", side_to_use, ".GroupDistanceMetrics.txt"), 
    quote = F, 
    sep = "\t"
  )
  
  ## Plot: distance from origin vs cluster tightness
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  ## Scatter: distance from origin vs 1D cluster size, per component
  p_dist <- ggplot(
    group_distances,
    aes(
      x     = dist_origin,
      y     = mean_dist_center,
      color = group_var,
      label = group_var,
      size  = n
    )
  ) +
    geom_point(alpha = 0.8) +
    geom_text_repel(size = 2, show.legend = FALSE) +
    facet_wrap(~ component) +
    labs(
      x     = "Distance from origin",
      y     = "Mean Cluster Spread",
      color = "Group",
      size  = "n per group",
      title = basename(loadings_path)
    ) +
    scale_color_manual(values = my_colors, na.value = "grey80") +
    theme_bw(base_size = 10)
  
  print(p_dist)
  
  ggsave(
    filename = paste0(path.plots, "Plot_", gsub(".txt", "", basename(loadings_path)), "clustering.pdf"),
    plot   = p_dist,
    width  = 8,
    height = 7,
    units  = "in",
    device = cairo_pdf
  )
  
}

##### WGCNA: CRISPR #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.plots   <- paste0(path.wd, "Plots/")

## WGCNA parameters (tune as needed)
soft_power    <- 6L
min_module_sz <- 5L

#### Prep for WGCNA by creating shared RNAi and CRISPR files
if (1) {
  
  ## Read in CRISPR data
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  ## Read in and format RNAi data
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F) %>%
    dplyr::select(ModelID, CCLEName)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## Filter data frames to common genes and cell lines
  common_genes    <- intersect(colnames(CRISPR), colnames(RNAi_t_ModelID))
  common_cells    <- intersect(rownames(CRISPR), rownames(RNAi_t_ModelID))
  
  CRISPR_common   <- CRISPR[common_cells, common_genes, drop = FALSE]
  RNAi_common     <- RNAi_t_ModelID[common_cells, common_genes, drop = FALSE]
  
  CRISPR_common[] <- lapply(CRISPR_common, as.numeric)
  RNAi_common[]   <- lapply(RNAi_common, as.numeric)
  
  CRISPR_common   <- as.data.frame(CRISPR_common)
  RNAi_common     <- as.data.frame(RNAi_common)

}

#### Run WGCNA (1) or Read in WGCNA object (0)
if (1) {
  
  ## Allow multi-threading for WGCNA
  WGCNA::enableWGCNAThreads()
  
  ## Run WGCNA on CRISPR dependencies
  net_CRISPR <- WGCNA::blockwiseModules(
    CRISPR_common,
    power              = soft_power,
    minModuleSize      = min_module_sz,
    networkType        = "signed", # anti correlate genes are not emphasized
    TOMType            = "signed",
    reassignThreshold  = 0,
    mergeCutHeight     = 0.25, # increase to merge similar modules
    numericLabels      = FALSE,
    pamRespectsDendro  = TRUE,
    verbose            = 3,
    deepSplit = 2
  )
  
  ## Save WGCNA object
  saveRDS(
    net_CRISPR,
    file = paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".rds"))

} else {
  
  ## Load in WGCNA object
  net_CRISPR <- readRDS(paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
}

#### Perform correlation on all non grey modules and prep for plot
if (1) {
  
  ## Extract module colors and gene tree
  moduleColors_CRISPR <- net_CRISPR$colors
  # table(moduleColors_CRISPR)
  
  ## Get names of genes assigned to clusters (remove "grey" genes)
  non_grey_genes <- colnames(CRISPR_common)[moduleColors_CRISPR != "grey"]
  length(non_grey_genes)
  # 4277, soft_power - 6L, min_module_sz - 5L
  # 1465, soft_power - 10L, min_module_sz - 5L
  
  CRISPR_ng <- CRISPR_common[, non_grey_genes, drop = FALSE]
  
  ## Perform correlation without grey genes
  cor_CRISPR <- stats::cor(
    CRISPR_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Reorder cor matrix so genes are grouped by their WGCNA module color
  mod_ng <- moduleColors_CRISPR[non_grey_genes]
  gene_order <- order(mod_ng)
  
  cor_CRISPR_ord <- cor_CRISPR[gene_order, gene_order]
  mod_ng_ord     <- mod_ng[gene_order]
  
  col_fun <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  ) # fix the legend scale
  
  ## Make module color mapping that matches WGCNA names exactly
  module_levels <- unique(mod_ng_ord)
  module_col <- stats::setNames(module_levels, module_levels)  # names == values == R colors

}

#### Plot 
if (1) {
  
  p_crispr <- ComplexHeatmap::Heatmap(
    cor_CRISPR_ord,
    name = "Pearson r",
    col  = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_ng_ord,
    column_split = mod_ng_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    
    use_raster = FALSE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(
      path.plots,
      "HEATMAP_WGCNA_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_crispr)
  grDevices::dev.off()
  
}

#### Repeat CRISPR plot but only on top # of modules
top_k <- 5L # Number of clusters

if (1) {

  ## Filter for top # of modules
  mod_sizes <- sort(table(mod_ng), decreasing = TRUE)
  top_modules <- names(mod_sizes)[seq_len(min(top_k, length(mod_sizes)))]
  
  ## Genes in those top modules
  top_genes <- non_grey_genes[mod_ng %in% top_modules]
  
  ## Subset CRISPR to the top modules and cor()
  CRISPR_top <- CRISPR_ng[, top_genes, drop = FALSE]
  
  cor_top <- stats::cor(
    CRISPR_top,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Module labels for these genes
  mod_top <- mod_ng[mod_ng %in% top_modules]
  
  ## Order genes by module
  gene_order_top <- order(mod_top)
  
  cor_top_ord <- cor_top[gene_order_top, gene_order_top, drop = FALSE]
  mod_top_ord <- mod_top[gene_order_top]

}

#### Plot
if (1) {
  
  p_crispr_top <- ComplexHeatmap::Heatmap(
    cor_top_ord,
    name = "Pearson r",
    col  = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_top_ord,
    column_split = mod_top_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(path.plots, "HEATMAP_WGCNA_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_Top", top_k, "Clusters.png"),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_crispr_top)
  grDevices::dev.off()

}

#### Look at all CRISPR modules now in RNAi
if (1) {
  ## Filter for all non grey modules
  RNAi_ng <- RNAi_common[, non_grey_genes, drop = FALSE]
  
  ## Gene–gene correlation for RNAi on the same genes
  cor_RNAi <- stats::cor(
    RNAi_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Use the SAME module labels/order you used for CRISPR
  cor_RNAi_ord <- cor_RNAi[gene_order, gene_order, drop = FALSE]

  ## Plot
  p_RNAi <- ComplexHeatmap::Heatmap(
    cor_RNAi_ord,
    name = "Pearson r",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_ng_ord,
    column_split = mod_ng_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 6
  )

  Cairo::CairoPNG(
    filename = paste0(
      path.plots,"HEATMAP_WGCNA_RNAi_OrderedByCRISPR_SoftPower_", soft_power,"_MinModuleSize_", min_module_sz, "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_RNAi)
  grDevices::dev.off()
  
}

#### Look at top CRISPR modules now in RNAi
if (1) {
  RNAi_top <- RNAi_common[, top_genes, drop = FALSE]
  
  cor_RNAi_top <- stats::cor(
    RNAi_top,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  cor_RNAi_top_ord <- cor_RNAi_top[gene_order_top, gene_order_top, drop = FALSE]
  
  p_RNAi_top <- ComplexHeatmap::Heatmap(
    cor_RNAi_top_ord,
    name = "Pearson r",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_top_ord,
    column_split = mod_top_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(path.plots,"HEATMAP_WGCNA_RNAi_OrderedByCRISPR_SoftPower_", soft_power,"_MinModuleSize_", min_module_sz, "_Top", top_k, "Clusters.png"),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_RNAi_top)
  grDevices::dev.off()
  
}

#### GO:BP Enrichment per CRISPR module 
if (1) {
  
  ## Get all unique modules (excluding grey)
  moduleColors_CRISPR <- net_CRISPR$colors
  unique_modules <- unique(moduleColors_CRISPR[moduleColors_CRISPR != "grey"])
  
  ## Create a list to store enrichment results for each module
  enrich_results <- list()
  
  ## Loop through each module and perform ORA
  for (module in unique_modules) {
    
    # Get genes in this module
    module_genes <- names(moduleColors_CRISPR)[moduleColors_CRISPR == module]
    
    # Gene Ontology enrichment
    ego <- enrichGO(
      gene          = module_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    # Store results
    enrich_results[[module]] <- list(
      GO      = ego,
      n_genes = length(module_genes)
    )
    
    cat("Module:", module, "- Genes:", length(module_genes), 
        "- GO terms:", nrow(ego@result), "\n")
  }
  
  ## Save RDS Object
  saveRDS(enrich_results, 
          file = paste0(path.wd, "DataSets/WGCNA/Enrichment_Results_CRISPR_SoftPower_", 
                        soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  ## Write file and sort by module size
  wb <- createWorkbook()
  
  # Order modules by size (largest to smallest)
  module_sizes <- sapply(names(enrich_results), function(m) enrich_results[[m]]$n_genes)
  modules_ordered <- names(sort(module_sizes, decreasing = TRUE))
  
  for (module in modules_ordered) {
    if (!is.null(enrich_results[[module]]$GO) && 
        nrow(enrich_results[[module]]$GO@result) > 0) {
      
      go_df <- enrich_results[[module]]$GO@result
      
      # Add sheet for this module
      addWorksheet(wb, sheetName = module)
      writeData(wb, sheet = module, x = go_df)
    }
  }
  
  saveWorkbook(wb, 
               file = paste0(path.wd, "DataSets/WGCNA/GO_Enrichment_CRISPR_AllModules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".xlsx"),
               overwrite = TRUE)
}

#### Visualize enrichment for top modules with significant results
n_modules_to_plot <- 5 # Number of modules

if (1) {
  
  ## Get all modules sorted by size
  all_modules <- names(sort(table(moduleColors_CRISPR[moduleColors_CRISPR != "grey"]), 
                            decreasing = TRUE))
  
  ## Filter to only modules with significant GO terms
  modules_with_sig_results <- c()
  for (module in all_modules) {
    if (!is.null(enrich_results[[module]]$GO) && 
        nrow(enrich_results[[module]]$GO@result) > 0 &&
        sum(enrich_results[[module]]$GO@result$p.adjust < 0.05) > 0) {
      modules_with_sig_results <- c(modules_with_sig_results, module)
    }
  }
  
  ## Take top n modules that have significant results
  top_modules <- head(modules_with_sig_results, n_modules_to_plot)
  
  cat("Plotting", length(top_modules), "modules with significant GO enrichment:\n")
  cat(paste(top_modules, collapse = ", "), "\n\n")
  
  ## Loop through and create plots for each
  for (target_module in top_modules) {
    
    ## Dotplot for GO terms
    p_go_dot <- dotplot(enrich_results[[target_module]]$GO, 
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(
      paste0(path.plots, "WGCGO_Dotplot_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
      p_go_dot,
      width = 10,
      height = 8)
    
    ## Barplot for GO terms
    p_go_bar <- barplot(enrich_results[[target_module]]$GO,
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(
      paste0(path.plots, "WGCNA_GO_Barplot_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
      p_go_bar,
      width = 10,
      height = 8)
    
    ## Enrichment map to show GO term relationships (with error handling)
    tryCatch({
      p_emap <- emapplot(pairwise_termsim(enrich_results[[target_module]]$GO),
                         showCategory = 30)
      
      ggsave(
        paste0(path.plots, "WGCNA_EnrichmentMap_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
        p_emap,
        width = 12,
        height = 10)
      
    }, error = function(e) {
      cat("Could not create enrichment map for module:", target_module, 
          "(not enough similar terms)\n")
    })
    
    cat("Plots saved for module:", target_module, "\n")
  }
  
  # Report if fewer than requested modules had significant results
  if (length(top_modules) < n_modules_to_plot) {
    cat("\nNote: Only", length(top_modules), "modules had significant GO enrichment (requested", n_modules_to_plot, ")\n")
  }
  
}

#### Visualize enrichment for n = # of modules (OG)
n_modules_to_plot <- 5 # Number of modules

if (1) {
  
  top_modules <- names(sort(table(moduleColors_CRISPR[moduleColors_CRISPR != "grey"]), 
                            decreasing = TRUE))[1:n_modules_to_plot]
  
  ## Loop through and create plots for each
  for (target_module in top_modules) {
    
    # Check if GO results exist
    if (!is.null(enrich_results[[target_module]]$GO) && 
        nrow(enrich_results[[target_module]]$GO@result) > 0) {
      
      ## Dotplot for GO terms
      p_go_dot <- dotplot(enrich_results[[target_module]]$GO, 
                          showCategory = 15,
                          title = paste0(target_module, " module - GO:BP enrichment"))
      
      ggsave(
        paste0(path.plots, "WGCGO_Dotplot_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
        p_go_dot,
        width = 10,
        height = 8)
      
      ## Barplot for GO terms
      p_go_bar <- barplot(enrich_results[[target_module]]$GO,
                          showCategory = 15,
                          title = paste0(target_module, " module - GO:BP enrichment"))
      
      ggsave(
        paste0(path.plots, "WGCNA_GO_Barplot_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
        p_go_bar,
        width = 10,
        height = 8)
      
      ## Enrichment map to show GO term relationships (with error handling)
      tryCatch({
        p_emap <- emapplot(pairwise_termsim(enrich_results[[target_module]]$GO),
                           showCategory = 30)
        
        ggsave(
          paste0(path.plots, "WGCNA_EnrichmentMap_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
          p_emap,
          width = 12,
          height = 10)
        
      }, error = function(e) {
        cat("Could not create enrichment map for module:", target_module, 
            "(not enough similar terms)\n")
      })
      
      cat("Plots saved for module:", target_module, "\n")
      
    } else {
      cat("No significant GO terms for module:", target_module, "\n")
    }
  }

}

#### Checking for conservation between CRISPR and RNAi
if (1) {
  
  moduleColors_CRISPR <- net_CRISPR$colors
  
  multiExpr <- list(
    CRISPR = list(data = CRISPR_common),
    RNAi   = list(data = RNAi_common)
  )
  
  multiColor <- list(
    CRISPR = moduleColors_CRISPR
  )
  
  ## Set up to run in parallel
  n_cores <- parallel::detectCores() - 1
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  WGCNA::enableWGCNAThreads(nThreads = n_cores)
  
  ## Run module preservation with parallelization
  set.seed(999)
  
  mp <- WGCNA::modulePreservation(
    multiExpr,
    multiColor,
    referenceNetworks = 1,      # CRISPR is reference (index 1)
    nPermutations = 200,        # increase to 200+ for publication
    randomSeed = 999,
    quickCor = 0,               # 0 = use WGCNA cor, 1 = use cor()
    verbose = 3,
    maxGoldModuleSize = 1000,   # modules larger than this use approximations
    maxModuleSize = 1000
  )
  
  ## Stop the cluster when done
  stopCluster(cl)
  
  ## Save the results
  saveRDS(mp, 
          file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", 
                        soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  ## Extract preservation statistics
  ref <- 1  # CRISPR
  test <- 2 # RNAi
  
  stats <- mp$preservation$Z$ref.CRISPR$inColumnsAlsoPresentIn.RNAi
  
  ## Interpretation thresholds (Langfelder & Horvath)
  # Zsummary < 2: no preservation
  # 2 < Zsummary < 10: weak to moderate preservation  
  # Zsummary > 10: strong preservation
  # Note: gold module (all genes) and grey (unassigned) are not informative
  
  write.table(
    x = stats,
    file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".txt"),
    row.names = T,
    sep = "\t",
    quote = F
  )

}

#### Visualize preservation statistics
if (1) {
  
  mp <- readRDS(file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", 
                              soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  stats <- mp$preservation$Z$ref.CRISPR$inColumnsAlsoPresentIn.RNAi
  obsStats <- mp$preservation$observed$ref.CRISPR$inColumnsAlsoPresentIn.RNAi
  
  plotData <- data.frame(
    module     = rownames(stats),
    size       = stats$moduleSize,
    Zsummary   = stats$Zsummary.pres,
    medianRank = obsStats$medianRank.pres
  )
  
  ## Remove gold and grey
  plotData <- plotData[!plotData$module %in% c("gold", "grey"), ]
  
  ## Add preservation category
  plotData$preservation <- cut(
    plotData$Zsummary,
    breaks = c(-Inf, 2, 10, Inf),
    labels = c("No preservation", "Weak-Moderate", "Strong preservation")
  )
  
  ## Plot 1: Zsummary vs module size
  p_preservation <- ggplot(plotData, aes(x = size, y = Zsummary, color = module, label = module)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x = "Module Size (number of genes)",
      y = "Preservation Z-summary",
      title = paste0("Module Preservation: CRISPR modules in RNAi data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    
    annotate("text", x = max(plotData$size) * 0.7, y = 2, 
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$size) * 0.7, y = 10, 
             label = "Z = 10 (strong)", vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_Zsummary_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".png"), p_preservation, width = 8, height = 7)
  
  ## Plot 2: Median rank vs Zsummary
  p_rank <- ggplot(plotData, aes(x = medianRank, y = Zsummary, color = module, label = module)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x = "Median Rank",
      y = "Preservation Z-summary",
      title = paste0("Module Preservation: CRISPR modules in RNAi data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    annotate("text", x = max(plotData$size) * 0.2, y = 2, 
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$size) * 0.2, y = 10, 
             label = "Z = 10 (strong)", vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_MedianRank_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".png"), p_rank, width = 8, height = 7)

}

##### WGCNA: RNAi #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.plots   <- paste0(path.wd, "Plots/")

## WGCNA parameters (tune as needed)
soft_power    <- 6L
min_module_sz <- 5L

#### Prep for WGCNA by creating shared RNAi and CRISPR files
if (1) {
  
  ## Read in CRISPR data
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  ## Read in and format RNAi data
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F) %>%
    dplyr::select(ModelID, CCLEName)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## Filter data frames to common genes and cell lines
  common_genes    <- intersect(colnames(CRISPR), colnames(RNAi_t_ModelID))
  common_cells    <- intersect(rownames(CRISPR), rownames(RNAi_t_ModelID))
  
  CRISPR_common   <- CRISPR[common_cells, common_genes, drop = FALSE]
  RNAi_common     <- RNAi_t_ModelID[common_cells, common_genes, drop = FALSE]
  
  CRISPR_common[] <- lapply(CRISPR_common, as.numeric)
  RNAi_common[]   <- lapply(RNAi_common, as.numeric)
  
  CRISPR_common   <- as.data.frame(CRISPR_common)
  RNAi_common     <- as.data.frame(RNAi_common)
  
}

#### Run WGCNA (1) or Read in WGCNA object (0)
if (1) {
  
  ## Allow multi-threading for WGCNA
  WGCNA::enableWGCNAThreads()
  
  ## Run WGCNA on RNAi dependencies
  net_RNAi <- WGCNA::blockwiseModules(
    RNAi_common,
    power              = soft_power,
    minModuleSize      = min_module_sz,
    networkType        = "signed", # anti correlate genes are not emphasized
    TOMType            = "signed",
    reassignThreshold  = 0,
    mergeCutHeight     = 0.25, # increase to merge similar modules
    numericLabels      = FALSE,
    pamRespectsDendro  = TRUE,
    verbose            = 3,
    deepSplit = 2
  )
  
  ## Save WGCNA object
  saveRDS(
    net_RNAi,
    file = paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
} else {
  
  ## Load in WGCNA object
  net_RNAi <- readRDS(paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
}

#### Perform correlation on all non grey modules and prep for plot
if (1) {
  
  ## Extract module colors and gene tree
  moduleColors_RNAi <- net_RNAi$colors
  # table(moduleColors_RNAi)
  
  ## Get names of genes assigned to clusters (remove "grey" genes)
  non_grey_genes <- colnames(RNAi_common)[moduleColors_RNAi != "grey"]
  length(non_grey_genes)
  # , soft_power - 6L, min_module_sz - 5L
  # , soft_power - 10L, min_module_sz - 5L
  
  RNAi_ng <- RNAi_common[, non_grey_genes, drop = FALSE]
  
  ## Perform correlation without grey genes
  cor_RNAi <- stats::cor(
    RNAi_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Reorder cor matrix so genes are grouped by their WGCNA module color
  mod_ng <- moduleColors_RNAi[non_grey_genes]
  gene_order <- order(mod_ng)
  
  cor_RNAi_ord <- cor_RNAi[gene_order, gene_order]
  mod_ng_ord     <- mod_ng[gene_order]
  
  col_fun <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  ) # fix the legend scale
  
  ## Make module color mapping that matches WGCNA names exactly
  module_levels <- unique(mod_ng_ord)
  module_col <- stats::setNames(module_levels, module_levels)  # names == values == R colors
  
}

#### Plot 
if (1) {
  
  p_RNAi <- ComplexHeatmap::Heatmap(
    cor_RNAi_ord,
    name = "Pearson r",
    col  = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_ng_ord,
    column_split = mod_ng_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    
    use_raster = FALSE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(
      path.plots,
      "HEATMAP_WGCNA_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_RNAi)
  grDevices::dev.off()
  
}

#### Repeat RNAi plot but only on top # of modules
top_k <- 5L # Number of clusters

if (1) {
  
  ## Filter for top # of modules
  mod_sizes <- sort(table(mod_ng), decreasing = TRUE)
  top_modules <- names(mod_sizes)[seq_len(min(top_k, length(mod_sizes)))]
  
  ## Genes in those top modules
  top_genes <- non_grey_genes[mod_ng %in% top_modules]
  
  ## Subset RNAi to the top modules and cor()
  RNAi_top <- RNAi_ng[, top_genes, drop = FALSE]
  
  cor_top <- stats::cor(
    RNAi_top,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Module labels for these genes
  mod_top <- mod_ng[mod_ng %in% top_modules]
  
  ## Order genes by module
  gene_order_top <- order(mod_top)
  
  cor_top_ord <- cor_top[gene_order_top, gene_order_top, drop = FALSE]
  mod_top_ord <- mod_top[gene_order_top]
  
}

#### Plot
if (1) {
  
  p_RNAi_top <- ComplexHeatmap::Heatmap(
    cor_top_ord,
    name = "Pearson r",
    col  = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_top_ord,
    column_split = mod_top_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(path.plots, "HEATMAP_WGCNA_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_Top", top_k, "Clusters.png"),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_RNAi_top)
  grDevices::dev.off()
  
}

#### Look at all RNAi modules now in CRISPR
if (1) {
  ## Filter for all non grey modules
  CRISPR_ng <- CRISPR_common[, non_grey_genes, drop = FALSE]
  
  ## Gene–gene correlation for RNAi on the same genes
  cor_CRISPR <- stats::cor(
    CRISPR_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Use the SAME module labels/order you used for RNAi
  cor_CRISPR_ord <- cor_CRISPR[gene_order, gene_order, drop = FALSE]
  
  ## Plot
  p_CRISPR <- ComplexHeatmap::Heatmap(
    cor_CRISPR_ord,
    name = "Pearson r",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_ng_ord,
    column_split = mod_ng_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(
      path.plots,"HEATMAP_WGCNA_CRISPR_OrderedByRNAi_SoftPower_", soft_power,"_MinModuleSize_", min_module_sz, "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_CRISPR)
  grDevices::dev.off()
  
}

#### Look at top RNAi modules now in CRISPR
if (1) {
  
  CRISPR_top <- CRISPR_common[, top_genes, drop = FALSE]
  
  cor_CRISPR_top <- stats::cor(
    CRISPR_top,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  cor_CRISPR_top_ord <- cor_CRISPR_top[gene_order_top, gene_order_top, drop = FALSE]
  
  p_CRISPR_top <- ComplexHeatmap::Heatmap(
    cor_CRISPR_top_ord,
    name = "Pearson r",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = mod_top_ord,
    column_split = mod_top_ord,
    row_title = NULL,
    column_title = NULL,
    row_gap = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_top_ord,
      col = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(path.plots,"HEATMAP_WGCNA_CRISPR_OrderedByRNAi_SoftPower_", soft_power,"_MinModuleSize_", min_module_sz, "_Top", top_k, "Clusters.png"),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_CRISPR_top)
  grDevices::dev.off()
  
}

#### GO:BP Enrichment per RNAi module 
if (1) {
  
  ## Get all unique modules (excluding grey)
  moduleColors_RNAi <- net_RNAi$colors
  unique_modules <- unique(moduleColors_RNAi[moduleColors_RNAi != "grey"])
  
  ## Create a list to store enrichment results for each module
  enrich_results <- list()
  
  ## Loop through each module and perform ORA
  for (module in unique_modules) {
    
    # Get genes in this module
    module_genes <- names(moduleColors_RNAi)[moduleColors_RNAi == module]
    
    # Gene Ontology enrichment
    ego <- enrichGO(
      gene          = module_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    # Store results
    enrich_results[[module]] <- list(
      GO      = ego,
      n_genes = length(module_genes)
    )
    
    cat("Module:", module, "- Genes:", length(module_genes), 
        "- GO terms:", nrow(ego@result), "\n")
  }
  
  ## Save RDS Object
  saveRDS(enrich_results, 
          file = paste0(path.wd, "DataSets/WGCNA/Enrichment_Results_RNAi_SoftPower_", 
                        soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  ## Write file and sort by module size
  wb <- createWorkbook()
  
  # Order modules by size (largest to smallest)
  module_sizes <- sapply(names(enrich_results), function(m) enrich_results[[m]]$n_genes)
  modules_ordered <- names(sort(module_sizes, decreasing = TRUE))
  
  for (module in modules_ordered) {
    if (!is.null(enrich_results[[module]]$GO) && 
        nrow(enrich_results[[module]]$GO@result) > 0) {
      
      go_df <- enrich_results[[module]]$GO@result
      
      # Add sheet for this module
      addWorksheet(wb, sheetName = module)
      writeData(wb, sheet = module, x = go_df)
    }
  }
  
  saveWorkbook(wb, 
               file = paste0(path.wd, "DataSets/WGCNA/GO_Enrichment_RNAi_AllModules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".xlsx"),
               overwrite = TRUE)
}

#### Visualize enrichment for top modules with significant results
n_modules_to_plot <- 5 # Number of modules

if (1) {
  
  ## Get all modules sorted by size
  all_modules <- names(sort(table(moduleColors_RNAi[moduleColors_RNAi != "grey"]), 
                            decreasing = TRUE))
  
  ## Filter to only modules with significant GO terms
  modules_with_sig_results <- c()
  for (module in all_modules) {
    if (!is.null(enrich_results[[module]]$GO) && 
        nrow(enrich_results[[module]]$GO@result) > 0 &&
        sum(enrich_results[[module]]$GO@result$p.adjust < 0.05) > 0) {
      modules_with_sig_results <- c(modules_with_sig_results, module)
    }
  }
  
  ## Take top n modules that have significant results
  top_modules <- head(modules_with_sig_results, n_modules_to_plot)
  
  cat("Plotting", length(top_modules), "modules with significant GO enrichment:\n")
  cat(paste(top_modules, collapse = ", "), "\n\n")
  
  ## Loop through and create plots for each
  for (target_module in top_modules) {
    
    ## Dotplot for GO terms
    p_go_dot <- dotplot(enrich_results[[target_module]]$GO, 
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(paste0(path.plots, "WGCGO_Dotplot_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
           p_go_dot, width = 10, height = 8)
    
    ## Barplot for GO terms
    p_go_bar <- barplot(enrich_results[[target_module]]$GO,
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(paste0(path.plots, "WGCNA_GO_Barplot_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
           p_go_bar, width = 10, height = 8)
    
    ## Enrichment map to show GO term relationships (with error handling)
    tryCatch({
      p_emap <- emapplot(pairwise_termsim(enrich_results[[target_module]]$GO),
                         showCategory = 30)
      
      ggsave(paste0(path.plots, "WGCNA_EnrichmentMap_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, target_module, "_Module.png"),
             p_emap, width = 12, height = 10)
    }, error = function(e) {
      cat("Could not create enrichment map for module:", target_module, 
          "(not enough similar terms)\n")
    })
    
    cat("Plots saved for module:", target_module, "\n")
  }
  
  ## Report if fewer than requested modules had significant results
  if (length(top_modules) < n_modules_to_plot) {
    cat("\nNote: Only", length(top_modules), "modules had significant GO enrichment (requested", n_modules_to_plot, ")\n")
  }
  
}

#### Checking for conservation between RNAi and CRISPR
if (1) {
  
  moduleColors_RNAi <- net_RNAi$colors
  
  multiExpr <- list(
    RNAi   = list(data = RNAi_common),
    CRISPR = list(data = CRISPR_common)
  )
  
  multiColor <- list(
    RNAi = moduleColors_RNAi
  )
  
  ## Set up to run in parallel
  n_cores <- parallel::detectCores() - 1
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  WGCNA::enableWGCNAThreads(nThreads = n_cores)
  
  ## Run module preservation with parallelization
  set.seed(999)
  
  mp <- WGCNA::modulePreservation(
    multiExpr,
    multiColor,
    referenceNetworks = 1,      # RNAi is reference (index 1)
    nPermutations = 200,        # increase to 200+ for publication
    randomSeed = 999,
    quickCor = 0,               # 0 = use WGCNA cor, 1 = use cor()
    verbose = 3,
    maxGoldModuleSize = 1000,   # modules larger than this use approximations
    maxModuleSize = 1000
  )
  
  ## Stop the cluster when done
  stopCluster(cl)
  
  ## Save the results
  saveRDS(mp, 
          file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_RNAi_in_CRISPR_SoftPower_", 
                        soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  ## Extract preservation statistics
  ref <- 1  # RNAi
  test <- 2 # CRISPR
  
  stats <- mp$preservation$Z$ref.RNAi$inColumnsAlsoPresentIn.CRISPR
  
  ## Interpretation thresholds (Langfelder & Horvath)
  # Zsummary < 2: no preservation
  # 2 < Zsummary < 10: weak to moderate preservation  
  # Zsummary > 10: strong preservation
  # Note: gold module (all genes) and grey (unassigned) are not informative
  
  write.table(
    x = stats,
    file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_RNAi_in_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".txt"),
    row.names = T,
    sep = "\t",
    quote = F
  )
  
}

#### Visualize preservation statistics
if (1) {
  
  mp <- readRDS(file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_RNAi_in_CRISPR_SoftPower_", 
                              soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  stats <- mp$preservation$Z$ref.RNAi$inColumnsAlsoPresentIn.CRISPR
  obsStats <- mp$preservation$observed$ref.RNAi$inColumnsAlsoPresentIn.CRISPR
  
  plotData <- data.frame(
    module     = rownames(stats),
    size       = stats$moduleSize,
    Zsummary   = stats$Zsummary.pres,
    medianRank = obsStats$medianRank.pres
  )
  
  ## Remove gold and grey
  plotData <- plotData[!plotData$module %in% c("gold", "grey"), ]
  
  ## Add preservation category
  plotData$preservation <- cut(
    plotData$Zsummary,
    breaks = c(-Inf, 2, 10, Inf),
    labels = c("No preservation", "Weak-Moderate", "Strong preservation")
  )
  
  ## Plot 1: Zsummary vs module size
  p_preservation <- ggplot(plotData, aes(x = size, y = Zsummary, color = module, label = module)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x = "Module Size (number of genes)",
      y = "Preservation Z-summary",
      title = paste0("Module Preservation: RNAi modules in CRISPR data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    
    annotate("text", x = max(plotData$size) * 0.7, y = 2, 
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$size) * 0.7, y = 10, 
             label = "Z = 10 (strong)", vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_Zsummary_RNAi_in_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".png"), p_preservation, width = 8, height = 7)
  
  ## Plot 2: Median rank vs Zsummary
  p_rank <- ggplot(plotData, aes(x = medianRank, y = Zsummary, color = module, label = module)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x = "Median Rank",
      y = "Preservation Z-summary",
      title = paste0("Module Preservation: RNAi modules in CRISPR data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    annotate("text", x = max(plotData$size) * 0.2, y = 2, 
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$size) * 0.2, y = 10, 
             label = "Z = 10 (strong)", vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_MedianRank_RNAi_in_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".png"), p_rank, width = 8, height = 7)
  
}

##### Gene-Only PCA + WGCNA Cluster Investigation #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")
path.pca     <- paste0(path.wd, "DataSets/PCA/")
path.max     <- paste0(path.wd, "DataSets/MaxLoading/")

## WGCNA Parameters
soft_power    <- 10L
min_module_sz <- 5L

#### Run to created shared CRISPR/RNAi files & save (only needs to be ran once)
if (1) {
  
  ## Read in CRISPR data
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  ## Read in and format RNAi data
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F) %>%
    dplyr::select(ModelID, CCLEName)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## Filter data frames to common genes and cell lines
  common_genes    <- intersect(colnames(CRISPR), colnames(RNAi_t_ModelID))
  common_cells    <- intersect(rownames(CRISPR), rownames(RNAi_t_ModelID))
  
  CRISPR_common   <- CRISPR[common_cells, common_genes, drop = FALSE]
  RNAi_common     <- RNAi_t_ModelID[common_cells, common_genes, drop = FALSE]
  
  CRISPR_common[] <- lapply(CRISPR_common, as.numeric)
  RNAi_common[]   <- lapply(RNAi_common, as.numeric)
  
  CRISPR_common_t   <- as.data.frame(t(CRISPR_common))
  RNAi_common_t     <- as.data.frame(t(RNAi_common))
  
  write.table(
    CRISPR_common_t,
    file = paste0(path.pca, "CRISPR_common_PCA.txt"),
    row.names = T,
    sep = "\t",
    quote = F
  )
  
  write.table(
    RNAi_common_t,
    file = paste0(path.pca, "RNAi_common_PCA.txt"),
    row.names = T,
    sep = "\t",
    quote = F
  )
  
}

#### Run PCA 
if (1) {
  
  source("/Users/jack/Documents/GitHub/FDB_Freeland/Scripts/PCA_from_file.R")
  
  CRISPR_path <- paste0(path.pca, "CRISPR_common_PCA.txt")
  RNAi_path   <- paste0(path.pca, "RNAi_common_PCA.txt")
  
  PCA_from_file(
    file = CRISPR_path,
    center = TRUE,
    scale = FALSE,
    fread = FALSE
  )
  
  PCA_from_file(
    file = RNAi_path,
    center = TRUE,
    scale = FALSE,
    fread = FALSE
  )
}

### Generate table of distances and theta and create plots
if (1) {
  
  ## Read in PCA loadings
  CRISPR_loadings <- read.delim(
    file = paste0(path.pca, "CRISPR_common_PCA_prcomp_loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading)) %>%
    dplyr::select(Loading, paste0("PC", 1:10))
  
  RNAi_loadings <- read.delim(
    file = paste0(path.pca, "RNAi_common_PCA_prcomp_loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading)) %>%
    dplyr::select(Loading, paste0("PC", 1:10))
  
  ## Read in WGCNA objects for clusters
  WGCNA_RNAi <- readRDS(paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  WGCNA_CRISPR <- readRDS(paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".rds"))
  
  ## Extract module colors
  moduleColors_RNAi <- WGCNA_RNAi$colors
  moduleColors_CRISPR <- WGCNA_CRISPR$colors
  
  ## Find max abs()
  CRISPR_loadings_max <- CRISPR_loadings %>%
    tidyr::pivot_longer(
      cols = paste0("PC", 1:10),
      names_to  = "component",
      values_to = "loading"
    ) %>%
    dplyr::mutate(abs_loading = abs(loading)) %>%
    dplyr::group_by(Loading) %>%
    dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      component_CRISPR    = component,
      loading_CRISPR      = loading,
      abs_loading_CRISPR  = abs_loading
    )
  
  RNAi_loadings_max <- RNAi_loadings %>%
    tidyr::pivot_longer(
      cols = paste0("PC", 1:10),
      names_to  = "component",
      values_to = "loading"
    ) %>%
    dplyr::mutate(abs_loading = abs(loading)) %>%
    dplyr::group_by(Loading) %>%
    dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      component_RNAi    = component,
      loading_RNAi      = loading,
      abs_loading_RNAi  = abs_loading
    )
  
  ## Merge
  Max_PCA <- merge(CRISPR_loadings_max, RNAi_loadings_max, by = "Loading")
  
  ## Add WGCNA module assignments
  Max_PCA$CRISPR_module <- moduleColors_CRISPR[Max_PCA$Loading]
  Max_PCA$RNAi_module <- moduleColors_RNAi[Max_PCA$Loading]
  
  ## Create color columns: module color if non-grey, otherwise grey
  Max_PCA <- Max_PCA %>%
    dplyr::mutate(
      CRISPR_module_color = ifelse(!is.na(CRISPR_module) & CRISPR_module != "grey", 
                                   CRISPR_module, 
                                   "grey80"),
      RNAi_module_color = ifelse(!is.na(RNAi_module) & RNAi_module != "grey", 
                                 RNAi_module, 
                                 "grey80")
    )
  
  xcol <- names(Max_PCA)[7]   # 7th column name abs_loading_RNAi
  ycol <- names(Max_PCA)[4]   # 4th column name abs_loading_CRISPR
  
  Max_PCA <- Max_PCA %>%
    dplyr::mutate(
      theta_rad = atan2(.data[[ycol]], .data[[xcol]]),
      theta_deg = theta_rad * 180 / pi,
      r         = sqrt(.data[[xcol]]^2 + .data[[ycol]]^2)
    )
  
  write.table(
    x = Max_PCA,
    file = paste0(path.max, "MaxPCALoadings_CRISPR_vs_RNAi.txt"),
    quote = F, sep = "\t", col.names = T, row.names = F
  )
  
  ## Plot scatter
  plot_df <- data.frame(
    x    = Max_PCA[[xcol]],
    y    = Max_PCA[[ycol]],
    gene = Max_PCA$Loading
  )
  
  p_scatter <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(size = 0.075, alpha = 0.3) +
    geom_abline(
      slope = tan(pi/6), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    geom_abline(
      slope = tan(pi/3), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    labs(
      x = xcol,
      y = ycol,
      title = "Max Absolute PCA Loadings per Gene PC 1-10"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    theme_classic(base_size = 10)
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_CRISPR_vs_RNAi_Scatter.pdf"),
    plot = p_scatter, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
  ## Plot polar - CRISPR modules (all genes, colored by CRISPR module)
  Max_PCA_CRISPR_ordered <- Max_PCA %>%
    dplyr::arrange(desc(CRISPR_module_color == "grey80"))
  
  p_polar_CRISPR <- ggplot(Max_PCA_CRISPR_ordered, aes(x = r, y = theta_deg, color = CRISPR_module_color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(
      yintercept = c(30, 60),
      linetype   = "dotted",
      color      = "grey40",
      linewidth  = 0.4
    ) +
    scale_color_identity() +
    labs(
      x = "r (magnitude)",
      y = "θ (degrees)",
      title = paste0("Max PCA Loadings: CRISPR WGCNA Modules: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, 90)) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_Polar_CRISPR_Modules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".pdf"),
    plot = p_polar_CRISPR, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
  ## Plot polar - RNAi modules (all genes, colored by RNAi module)
  Max_PCA_RNAi_ordered <- Max_PCA %>%
    dplyr::arrange(desc(RNAi_module_color == "grey80"))
  
  p_polar_RNAi <- ggplot(Max_PCA_RNAi_ordered, aes(x = r, y = theta_deg, color = RNAi_module_color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(
      yintercept = c(30, 60),
      linetype   = "dotted",
      color      = "grey40",
      linewidth  = 0.4
    ) +
    scale_color_identity() +
    labs(
      x = "r (magnitude)",
      y = "θ (degrees)",
      title = paste0("Max PCA Loadings: RNAi WGCNA Modules: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, 90)) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_Polar_RNAi_Modules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, ".pdf"),
    plot = p_polar_RNAi, width = 5, height = 4, units = "in", device = cairo_pdf
  )
}

##### Melanoma / Differentiation Signature #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.plots   <- paste0(path.wd, "Plots/")
path.mel     <- paste0(path.wd, "DataSets/Melanoma/")

#### Get Hallmark EMT gene set from MSigDB
if (1) {

  ## Get Hallmark EMT gene set
  hallmark_emt <- msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "H"
  ) %>%
    dplyr::filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
    dplyr::pull(gene_symbol) %>%
    unique()
  
  cat("Number of genes in Hallmark EMT gene set:", length(hallmark_emt), "\n")
  
  ## Save gene list
  write.table(
    x = data.frame(Gene = hallmark_emt),
    file = paste0(path.mel, "Hallmark_EMT_GeneSet.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
}

#### ssGSEA
if (1) {
  
  ## Load in data and format
  counts <- read.table(
    file = paste0(path.dm, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
    sep = ",",
    header = TRUE)
  
  counts_trim <- counts %>%
    dplyr::filter(IsDefaultEntryForModel == "Yes") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .)) %>%
    dplyr::select(-X, -SequencingID, -IsDefaultEntryForModel, -ModelConditionID, -IsDefaultEntryForMC) %>%
    tibble::column_to_rownames(var = "ModelID") %>%
    t() %>%
    data.frame()
  
  ## Calculate Signature via ssGSEA
  run.ssGSEA2 <- function(exp_mat, gene_list, method = "ssgsea", norm = F){
    gsvaPar <- GSVA::ssgseaParam(exp_mat, gene_list, normalize = norm)
    as.data.frame(t(GSVA::gsva(gsvaPar)))
  }
  
  EMT_scores <- run.ssGSEA2(
    exp_mat = as.matrix(counts_trim),
    gene_list = list(EMT = hallmark_emt)
  )

}








### MOVE FROM HERE


#### Step 5: Add cell line metadata
if (1) {
  
  cat("\n=== Step 5: Adding Cell Line Metadata ===\n")
  
  ## Load model information
  models <- read.delim(
    file = paste0(path.dm, "Model.csv"),
    sep = ",",
    stringsAsFactors = FALSE
  )
  
  cat("Model metadata loaded:", nrow(models), "cell lines\n")
  
  ## Merge EMT signature with metadata
  emt_with_metadata <- emt_signature %>%
    dplyr::left_join(
      models %>% dplyr::select(
        ModelID, 
        StrippedCellLineName,
        OncotreeLineage, 
        OncotreePrimaryDisease,
        OncotreeSubtype
      ),
      by = c("DepMap_ID" = "ModelID")
    )
  
  cat("Merged dataset:", nrow(emt_with_metadata), "cell lines with metadata\n")
  
  ## Save merged data
  write.table(
    emt_with_metadata,
    file = paste0(path.fig3, "EMT_Signature_with_Metadata.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  ## Show top and bottom cell lines by EMT signature
  cat("\nTop 10 cell lines with highest EMT (most mesenchymal-like):\n")
  print(
    emt_with_metadata %>%
      dplyr::arrange(desc(EMT_Signature)) %>%
      dplyr::select(StrippedCellLineName, OncotreeLineage, EMT_Signature) %>%
      head(10)
  )
  
  cat("\nTop 10 cell lines with lowest EMT (most epithelial-like):\n")
  print(
    emt_with_metadata %>%
      dplyr::arrange(EMT_Signature) %>%
      dplyr::select(StrippedCellLineName, OncotreeLineage, EMT_Signature) %>%
      head(10)
  )
  
  ## Plot by lineage
  cat("\n=== Creating lineage-specific EMT plots ===\n")
  
  ## Get top lineages by count
  top_lineages <- emt_with_metadata %>%
    dplyr::filter(!is.na(OncotreeLineage)) %>%
    dplyr::count(OncotreeLineage, sort = TRUE) %>%
    head(10) %>%
    dplyr::pull(OncotreeLineage)
  
  emt_top_lineages <- emt_with_metadata %>%
    dplyr::filter(OncotreeLineage %in% top_lineages)
  
  p_lineage <- ggplot(
    emt_top_lineages,
    aes(x = reorder(OncotreeLineage, EMT_Signature, FUN = median), 
        y = EMT_Signature)
  ) +
    geom_boxplot(aes(fill = OncotreeLineage), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3") +
    labs(
      x = "Cancer Lineage",
      y = "EMT Signature Score",
      title = "EMT Signature by Cancer Lineage",
      subtitle = "Top 10 lineages by cell line count"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(
    filename = paste0(path.plots, "EMT_Signature_by_Lineage.pdf"),
    plot = p_lineage,
    width = 10,
    height = 8,
    device = cairo_pdf
  )
  
  print(p_lineage)
  
}
