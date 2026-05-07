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
  library(MASS)
  library(fitdistrplus)
  library(ggrepel)
  library(WGCNA)
  library(patchwork)
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
OS <- "Mac" # Linux or Mac

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

##### Imputation: GDSC (NA's in Data) #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd   <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm   <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.gdsc <- paste0(path.wd, "DataSets/GDSC/")

## Load in cell line info from depamp
gdsc <- read.delim(paste0(path.gdsc,"sanger-dose-response.csv"), sep = ",", stringsAsFactors = F, check.names = F)











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
path.stat    <- paste0(path.wd, "DataSets/Stats/")

## Set PLS parameters
X_source <- "CRISPR" # CRISPR or CTRP
Y_source <- "CTRP"   # CRISPR or CTRP

ncomp <- 15
mode  <- "canonical" # default = regression, symmetric = canonical

## Derived label for plot titles
mode_label <- if (mode == "canonical") "PLS-C" else "PLS-R"

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- character(0)  # e.g. c("Myeloid", "Lymphoid") or character(0)

## For plot iterations 
Special_string <- "_MED12" # character(0) or "_VALUE"

## Filtered for all three data sets shared lines?
FilteredAll3 <- TRUE # TRUE or FALSE

#### 1. Execute to prep for PLS
if (1) {
  
  ## For saving files later
  excl_tag <- if (length(exclude_lineages) > 0) {
    paste0("_excl.", paste(exclude_lineages, collapse = "."))
  } else {
    ""
  }
  file_tag <- paste0("PLS_Mode.", mode, "_X.", X_source, "_Y.", Y_source, excl_tag)
  
  if (FilteredAll3 == TRUE) {
    Filtered_Tag <- "_Filtered3"
  } else {
    Filtered_Tag <- character(0)
  }
  
  ## Read in CRISPR
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  ## Read in CTRP
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  if (FilteredAll3 == TRUE) {
    
    ## Read in RNAi
    RNAi <- read.delim(
      file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
      sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
    )
    
    models <- read.delim(paste0(path.dm, "Model.csv"), sep = ",", stringsAsFactors = F, check.names = F) %>%
      dplyr::select(ModelID, CCLEName, OncotreeLineage)
    
    RNAi_t <- RNAi %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "CCLEName") %>%
      dplyr::rename_with(~ sub("\\.\\..*", "", .))
    
    RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
      dplyr::select(-CCLEName, -OncotreeLineage) %>%
      tibble::column_to_rownames(var = "ModelID")
    
    ## Filter CRISPR columns to genes shared with RNAi
    shared_genes <- intersect(colnames(CRISPR), colnames(RNAi_t_ModelID))
    message("Shared genes between CRISPR and RNAi: ", length(shared_genes))
    CRISPR <- CRISPR[, shared_genes, drop = FALSE]
    
  }
  
  ## Filter for shared cell lines, make matrix (for mixomics), ensure numeric
  if (X_source == "CRISPR") X_data <- CRISPR
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "CRISPR") Y_data <- CRISPR
  if (Y_source == "CTRP")   Y_data <- CTRP
  
  if (FilteredAll3 == TRUE) {
    # Three-way intersection: CRISPR, CTRP, and RNAi samples
    ids <- Reduce(intersect, list(
      rownames(X_data),
      rownames(Y_data), 
      rownames(RNAi_t_ModelID) 
    ))
  } else {
    # Two-way intersection: CRISPR and CTRP only
    ids <- intersect(rownames(X_data), rownames(Y_data))
  }
  
  ## Filter out cell lines belonging to excluded lineages
  if (length(exclude_lineages) > 0) {
    models_filt <- read.csv(paste0(path.dm, "Model.csv"))
    keep_ids <- models_filt$ModelID[!(models_filt$OncotreeLineage %in% exclude_lineages)]
    ids <- intersect(ids, keep_ids)
  }
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
}

#### 2. Execute to run PLS and save output files (requires Step 1)
if (0) {
  
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
  
  ## Canonical correlations between X and Y variates for comps 1-10
  n_cancor <- min(10L, ncomp)
  cancor.df <- data.frame(
    comp                  = seq_len(n_cancor),
    canonical_correlation = sapply(
      seq_len(n_cancor),
      function(i) cor(pls_fit$variates$X[, i], pls_fit$variates$Y[, i])
    )
  )
  print(cancor.df)
  
  ## Save files
  write.table(
    x = x.variates,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.variates,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = variates.X.Y,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.loadings,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.loadings,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.exp_variance,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.exp_variance,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = cancor.df,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_canonical_correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
}

#### 3. Execute to plot PLS loadings (requires Step 1)
if (1) {
  
  ## Load saved loading files
  X_loadings <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.loadings.txt"),
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
        Group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|dabrafenib|PLX\\-4720|PLX\\-4032|dabrafenib|GDC\\-0879)$") ~ "BRAF & MEK\nInhibitors", # sorafenib, regorafenib, RAF265
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "EGFR & HER2\nInhibitors",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "Ferroptosis\nInducers",
          # stringr::str_detect(Loading, "^(clofarabine|procarbazine|carboplatin|cytarabine hydrochloride|bleomycin A2)$") ~ "DNA Damage\nInducers",
          # stringr::str_detect(Loading, "^(AZD8055|XL765|OSI-027|sirolimus|temsirolimus|KU-0063794)$") ~ "mTOR Inhibitors",
          # stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "DDR Pathway\nInhibitors",
          # stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          # stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          # detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ "Other"
        ),
        # group.atp5 = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na = dplyr::if_else(is.na(Group), 1L, 0L),
        # group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(Group), Loading, NA_character_),
        # label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
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
    
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    df <- df %>%
      dplyr::mutate(
        Group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)\\b") ~ "BRAF Signalling",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|ERBB2)$")     ~ "EGFR Signalling",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "Ferroptosis",
          # stringr::str_detect(Loading, "^(RRM1|RRM2|POLA1|POLE|POLD1|DCK|SLC29A1|ERCC1|ERCC2|XPA|XPC|BRCA1|BRCA2|PALB2|RAD51|BLMH|TOP2A|CYP1A2|CYP2B6)$") ~ "DNA Damage Repair Genes",
          # stringr::str_detect(Loading, "^(MTOR|PIK3CA|PIK3CB|PTEN|TSC1|TSC2|AKT1|AKT2)$") ~ "mTOR Signalling",
          # stringr::str_detect(Loading, "^(MDM2|PPM1D|USP7|MDM4|CDKN1A|ATM|ERBB3|TP53|CHEK2|TP53BP1|USP28)$") ~ "DNA Damage\nResponse",
          # stringr::str_detect(Loading, "^MDM[24]$")                                  ~ "04 MDM2.MDM4",
          # stringr::str_detect(Loading, "^ATP5")                                      ~ "05 ATP5",
          # stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                         ~ "06 dasa targets",
          # stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                ~ "07 BCL2+",
          # stringr::str_detect(Loading, "^MYC(|N|L)")                                 ~ "08 MYC.",
          # stringr::str_detect(Loading, "^(GRB2|CRKL)$")                              ~ "09 SRC-related",
          # stringr::str_detect(Loading, "^TP53$")                                     ~ "10 TP53",
          stringr::str_detect(Loading, "^(MED12)$")
          ~ "MED12", # MED12
          # ~ "CDK8 Kinase Module", # MED12|MED13|CDK8|CCNC
          # ~ "Wnt/B-catenin Sig", # MED12|MED13|CDK8|CCNC|CTNNB1|APC|AXIN1|TCF7L2
          # ~ "Super-Enhancer Axis", # MED12|MED13|CDK8|CCNC|MED1|BRD4|CDK9|CCNT1
          TRUE ~ "Other"
        ),
        group.atp5        = dplyr::if_else(stringr::str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_),
        group.na          = dplyr::if_else(is.na(Group), 1L, 0L),
        group.atp5.na     = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na      = dplyr::if_else(!is.na(Group), Loading, NA_character_),
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
  # my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
  #                "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  # 
  
  group_colors <- c(
    "BRAF Signalling"  = "#F8766D",
    "BRAF & MEK\nInhibitors"    = "#F8766D",
    "EGFR Signalling"  = "#DE8C00",
    "EGFR & HER2\nInhibitors"   = "#DE8C00",
    "Ferroptosis"      = "#B79F00",
    "Ferroptosis\nInducers"       = "#B79F00",
    "MED12"            = "#00BA38",
    "DNA Damage Repair Genes" = "purple",
    "DNA Damage\nInducers" = "purple",
    "mTOR Inhibitors" = "black",
    "mTOR Signalling" = "black",
    "Other"            = "grey80"
    # "MDM2\nInhibitors" = "#00BA38",
    # "DNA Damage\nResponse" = "#00BA38"
    
  )
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    df <- df %>%
      dplyr::mutate(
        label_flag = dplyr::if_else(
          is.na(.data[[color_col]]) | .data[[color_col]] == "Other",
          "Unlabeled",
          "Labeled"
        )
      ) %>%
      dplyr::arrange(desc(.data[[color_col]] == "Other" | is.na(.data[[color_col]])))
    
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1,
          y     = comp2,
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
        scale_color_manual(
          values = group_colors,
          na.value = "grey80",
          breaks = c(names(group_colors)[names(group_colors) != "Other"], "Other")
        ) +
        labs(
          # title = paste0(mode_label, " | ", source_label, " loadings: ", comp1, " vs ", comp2)
          title = paste0(mode_label, ": ", source_label," Loadings")
          ) +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".loadings_", comp1, "vs", comp2, Special_string, ".pdf"
        ),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, X_source, "Group", "Loading") # paste0("X.", X_source)
  } else {
    plot_loadings_side(X_plot, X_source, "Group", "Loading") # paste0("X.", X_source)
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, Y_source, "Group", "Loading") # paste0("Y.", Y_source)
  } else {
    plot_loadings_side(Y_plot, Y_source, "Group", "Loading") # paste0("Y.", Y_source)
  }
}

#### 4. Execute to plot PLS scores colored by cancer type (requires Step 1)
if (1) {
  
  ## Load model metadata
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  ## Load saved variates files
  x.variates.plot <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.variates.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  y.variates.plot <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.variates.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Annotate with cancer type
  x.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(x.variates.plot$Score, model$ModelID)]
  y.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(y.variates.plot$Score, model$ModelID)]
  
  ## Get top N lineages by cell line count for coloring
  top_lineages_n <- 15
  top_lineages <- names(sort(table(x.variates.plot$OncotreeLineage), decreasing = TRUE))[1:top_lineages_n]
  
  x.variates.plot <- x.variates.plot %>%
    dplyr::mutate(lineage_label = dplyr::if_else(OncotreeLineage %in% top_lineages, OncotreeLineage, "Other"))
  
  y.variates.plot <- y.variates.plot %>%
    dplyr::mutate(lineage_label = dplyr::if_else(OncotreeLineage %in% top_lineages, OncotreeLineage, "Other"))
  
  ## Color palette
  lineage_colors <- c(
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(7, "Set2"),
    "grey70"  # for "Other"
  )
  names(lineage_colors) <- c(top_lineages, "Other")
  
  ## Helper: scatter plots of scores
  plot_scores_side <- function(df, source_label) {
    
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    n_lineages <- length(unique(na.omit(df$OncotreeLineage)))
    lineage_pal <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_lineages),
      sort(unique(na.omit(df$OncotreeLineage)))
    )
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "comp1"
      comp2_col <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1_col,
          y     = comp2_col,
          color = "OncotreeLineage"
        )
      ) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_pal, name = "Lineage") +
        labs(
          # title = paste0(mode_label, " | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
          title = paste0(mode_label, ": ", source_label," Variate Scores"),
          x     = comp1_col,
          y     = comp2_col
        ) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_",
          comp1_col, "vs", comp2_col, ".pdf"
        ),
        plot = p, width = 7, height = 6, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, X_source) # paste0("X.", X_source)
  plot_scores_side(y.variates.plot, Y_source) # paste0("Y.", Y_source)
  
  ## Helper: boxplots of scores per lineage
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^comp([1-9]|10)$", names(df), value = TRUE)
    
    ## For the multi-facet overview, order lineages by comp1 median
    lineage_order_comp1 <- df %>%
      dplyr::filter(!is.na(OncotreeLineage)) %>%
      dplyr::group_by(OncotreeLineage) %>%
      dplyr::summarise(med = median(comp1, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(med) %>%
      dplyr::pull(OncotreeLineage)
    
    df_long <- df %>%
      dplyr::select(Score, OncotreeLineage, dplyr::all_of(comp_cols)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(comp_cols),
        names_to  = "Component",
        values_to = "Score_value"
      ) %>%
      dplyr::filter(!is.na(OncotreeLineage)) %>%
      dplyr::mutate(
        Component       = factor(Component, levels = comp_cols),
        OncotreeLineage = factor(OncotreeLineage, levels = lineage_order_comp1)
      )
    
    ## Multi-facet overview PDF (ordered by comp1 median)
    p <- ggplot(
      df_long,
      aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)
    ) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, size = 0.3) +
      facet_wrap(~ Component, scales = "free_y", ncol = 2) +
      scale_fill_manual(
        values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(
          length(levels(df_long$OncotreeLineage))
        ),
        guide = "none"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
      labs(
        title = paste0(mode_label, ": ", source_label, " Scores by Cancer Lineage (comps 1–10)"),
        # title = NULL,
        x     = NULL,
        y     = "Score"
      ) +
      theme_bw(base_size = 9) +
      theme(
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
        strip.text    = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines")
      )
    
    ggsave(
      filename = paste0(
        path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_comps1to10.pdf"
      ),
      plot = p, width = 12, height = 18, units = "in", device = cairo_pdf
    )
    
    ## Individual per-comp PDFs, each ordered by that comp's own median
    for (comp in comp_cols) {
      
      lineage_order_comp <- df %>%
        dplyr::filter(!is.na(OncotreeLineage)) %>%
        dplyr::group_by(OncotreeLineage) %>%
        dplyr::summarise(med = median(.data[[comp]], na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(med) %>%
        dplyr::pull(OncotreeLineage)
      
      df_comp <- df_long %>%
        dplyr::filter(Component == comp) %>%
        dplyr::mutate(OncotreeLineage = factor(OncotreeLineage, levels = lineage_order_comp))
      
      p_ind <- ggplot(
        df_comp,
        aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)
      ) +
        geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, size = 0.3) +
        scale_fill_manual(
          values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(
            length(levels(df_comp$OncotreeLineage))
          ),
          guide = "none"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
        labs(
          # title = paste0(mode_label, " | ", source_label, " scores — ", comp, " by cancer lineage"),
          title = paste0(mode_label, ": ", source_label, " Scores by Cancer Lineage (", comp, ")"),
          x     = NULL,
          y     = "Score"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_", comp, ".pdf"
        ),
        plot = p_ind, width = 10, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_boxplot(x.variates.plot, X_source) # paste0("X.", X_source)
  plot_scores_boxplot(y.variates.plot, Y_source) #  paste0("Y.", Y_source)
  
  ## Wilcoxon rank-sum: each lineage vs. all others, per component, FDR-corrected
  wilcox_lineage_tests <- function(df, source_label) {
    
    comp_cols <- grep("^comp([1-9]|10)$", names(df), value = TRUE)
    lineages  <- unique(na.omit(df$OncotreeLineage))
    
    results <- purrr::map_dfr(comp_cols, function(comp) {
      purrr::map_dfr(lineages, function(lin) {
        
        in_group  <- df[[comp]][!is.na(df$OncotreeLineage) & df$OncotreeLineage == lin]
        out_group <- df[[comp]][!is.na(df$OncotreeLineage) & df$OncotreeLineage != lin]
        
        if (length(in_group) < 3 || length(out_group) < 3) return(NULL)
        
        wt <- wilcox.test(in_group, out_group, exact = FALSE)
        
        data.frame(
          Component  = comp,
          Lineage    = lin,
          n_lineage  = length(in_group),
          median_in  = median(in_group),
          median_out = median(out_group),
          W          = wt$statistic,
          p_value    = wt$p.value
        )
      })
    })
    
    results$p_adj_BH <- p.adjust(results$p_value, method = "BH")
    results <- results %>% dplyr::arrange(p_adj_BH)
    
    write.table(
      x         = results,
      file      = paste0(path.pls, file_tag, Filtered_Tag, "_", source_label, "_wilcox_lineage.txt"),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
    
    results
  }
  
  wilcox_x <- wilcox_lineage_tests(x.variates.plot, paste0("X.", X_source)) %>%
    dplyr::arrange(Component, desc(median_in))
  
  wilcox_y <- wilcox_lineage_tests(y.variates.plot, paste0("Y.", Y_source)) %>%
    dplyr::arrange(Component, desc(median_in))
  
  ## Quick view of top hits
  head(wilcox_x, 20)
  head(wilcox_y, 20)
  
  ## Save Wilcoxon results tables
  write.table(
    x         = wilcox_x,
    file      = paste0(path.stat, file_tag, Filtered_Tag, "_X.", X_source, "_wilcox_lineage.txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  write.table(
    x         = wilcox_y,
    file      = paste0(path.stat, file_tag, Filtered_Tag, "_Y.", Y_source, "_wilcox_lineage.txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
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
path.stat    <- paste0(path.wd, "DataSets/Stats/")

## Set PLS parameters
X_source <- "RNAi"  # RNAi or CTRP
Y_source <- "CTRP"  # RNAi or CTRP

ncomp <- 15
mode  <- "canonical" # default = regression, symmetric = canonical

## Derived label for plot titles
mode_label <- if (mode == "canonical") "PLS-C" else "PLS-R"

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- character(0)  # e.g. c("Myeloid", "Lymphoid") or character(0)

## For plot iterations
Special_string <- "_MED12" # character(0) or "_VALUE"

## Filtered for all three data sets shared lines?
FilteredAll3 <- TRUE # TRUE or FALSE

#### 1. Execute to prep for PLS
if (1) {
  
  ## For saving files later
  excl_tag <- if (length(exclude_lineages) > 0) {
    paste0("_excl.", paste(exclude_lineages, collapse = "."))
  } else {
    ""
  }
  file_tag <- paste0("PLS_Mode.", mode, "_X.", X_source, "_Y.", Y_source, excl_tag)
  
  if (FilteredAll3 == TRUE) {
    Filtered_Tag <- "_Filtered3"
  } else {
    Filtered_Tag <- character(0)
  }
  
  ## Read in RNAi
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  ## Read in CTRP
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  )
  
  ## Convert RNAi sample nomenclature (CCLEName -> ModelID via Model.csv)
  models <- read.delim(paste0(path.dm, "Model.csv"), sep = ",", stringsAsFactors = F, check.names = F) %>%
    dplyr::select(ModelID, CCLEName, OncotreeLineage)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName, -OncotreeLineage) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  if (FilteredAll3 == TRUE) {
    
    ## Read in CRISPR (for shared sample and gene filtering)
    CRISPR <- read.delim(
      file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
      sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub("\\.\\..*", "", .))
    
    ## Filter RNAi columns to genes shared with CRISPR
    shared_genes <- intersect(colnames(RNAi_t_ModelID), colnames(CRISPR))
    message("Shared genes between RNAi and CRISPR: ", length(shared_genes))
    RNAi_t_ModelID <- RNAi_t_ModelID[, shared_genes, drop = FALSE]
    
  }
  
  ## Assign X and Y data
  if (X_source == "RNAi") X_data <- RNAi_t_ModelID
  if (X_source == "CTRP") X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP") Y_data <- CTRP
  
  if (FilteredAll3 == TRUE) {
    # Three-way intersection: RNAi, CTRP, and CRISPR samples
    ids <- Reduce(intersect, list(
      rownames(X_data),
      rownames(Y_data),
      rownames(CRISPR)
    ))
  } else {
    # Two-way intersection: RNAi and CTRP only
    ids <- intersect(rownames(X_data), rownames(Y_data))
  }
  
  ## Filter out cell lines belonging to excluded lineages
  if (length(exclude_lineages) > 0) {
    models_filt <- read.csv(paste0(path.dm, "Model.csv"))
    keep_ids <- models_filt$ModelID[!(models_filt$OncotreeLineage %in% exclude_lineages)]
    ids <- intersect(ids, keep_ids)
  }
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
}

#### 2. Execute to run PLS and save output files (requires Step 1)
if (0) {
  
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
  
  ## Canonical correlations between X and Y variates for comps 1-10
  n_cancor <- min(10L, ncomp)
  cancor.df <- data.frame(
    comp                  = seq_len(n_cancor),
    canonical_correlation = sapply(
      seq_len(n_cancor),
      function(i) cor(pls_fit$variates$X[, i], pls_fit$variates$Y[, i])
    )
  )
  print(cancor.df)
  
  ## Save files
  write.table(
    x = x.variates,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.variates,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = variates.X.Y,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.loadings,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.loadings,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = x.exp_variance,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = y.exp_variance,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.expvar.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    x = cancor.df,
    file = paste0(path.pls, file_tag, Filtered_Tag, "_canonical_correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
}

#### 3. Execute to plot PLS loadings (requires Step 1)
if (1) {
  
  ## Load saved loading files
  X_loadings <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.loadings.txt"),
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
        Group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|dabrafenib|PLX\\-4720|PLX\\-4032|dabrafenib|GDC\\-0879)$") ~ "BRAF & MEK\nInhibitors",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "EGFR & HER2\nInhibitors",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "Ferroptosis\nInducers",
          # stringr::str_detect(Loading, "^(clofarabine|procarbazine|carboplatin|cytarabine hydrochloride|bleomycin A2)$") ~ "DNA Damage\nInducers",
          # stringr::str_detect(Loading, "^(AZD8055|XL765|OSI-027|sirolimus|temsirolimus|KU-0063794)$") ~ "mTOR Inhibitors",
          # stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "DDR Pathway\nInhibitors",
          # stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          # stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          # detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ "Other"
        ),
        # group.atp5 = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na = dplyr::if_else(is.na(Group), 1L, 0L),
        # group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(Group), Loading, NA_character_),
        # label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
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
  
  ### Annotation for RNAi loadings file
  annotate_rnai <- function(df, side_label) {
    
    gene.info.all <- read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    df <- df %>%
      dplyr::mutate(
        Group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)\\b") ~ "BRAF Signalling",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|ERBB2)$")               ~ "EGFR Signalling",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "Ferroptosis",
          stringr::str_detect(Loading, "^(MED12)$")                                     ~ "MED12",
          # stringr::str_detect(Loading, "^(RRM1|RRM2|POLA1|POLE|POLD1|DCK|SLC29A1|ERCC1|ERCC2|XPA|XPC|BRCA1|BRCA2|PALB2|RAD51|BLMH|TOP2A|CYP1A2|CYP2B6)$") ~ "DNA Damage Repair Genes",
          # stringr::str_detect(Loading, "^(MTOR|PIK3CA|PIK3CB|PTEN|TSC1|TSC2|AKT1|AKT2)$") ~ "mTOR Signalling",
          # stringr::str_detect(Loading, "^(MDM2|PPM1D|USP7|MDM4|CDKN1A|ATM|TP53|CHEK2|TP53BP1|USP28)$") ~ "DNA Damage\nResponse",
          TRUE ~ "Other"
        ),
        group.atp5        = dplyr::if_else(stringr::str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_),
        group.na          = dplyr::if_else(is.na(Group), 1L, 0L),
        group.atp5.na     = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na      = dplyr::if_else(!is.na(Group), Loading, NA_character_),
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
  
  ## Annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_rnai(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_rnai(Y_loadings, "Y")
  
  ## Plotting colors (always plot both sides)
  group_colors <- c(
    "BRAF Signalling"        = "#F8766D",
    "BRAF & MEK\nInhibitors" = "#F8766D",
    "EGFR Signalling"        = "#DE8C00",
    "EGFR & HER2\nInhibitors"= "#DE8C00",
    "Ferroptosis"            = "#B79F00",
    "Ferroptosis\nInducers"  = "#B79F00",
    "MED12"                  = "#00BA38",
    "DNA Damage Repair Genes" = "purple",
    "DNA Damage\nInducers" = "purple",
    "mTOR Inhibitors" = "black",
    "mTOR Signalling" = "black",
    "Other"                  = "grey80"
    # "MDM2\nInhibitors" = "#00BA38",
    # "DNA Damage\nResponse" = "#619CFF"
  )
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    df <- df %>%
      dplyr::mutate(
        label_flag = dplyr::if_else(
          is.na(.data[[color_col]]) | .data[[color_col]] == "Other",
          "Unlabeled",
          "Labeled"
        )
      ) %>%
      dplyr::arrange(desc(.data[[color_col]] == "Other" | is.na(.data[[color_col]])))
    
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1,
          y     = comp2,
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
        scale_color_manual(
          values = group_colors,
          na.value = "grey80",
          breaks = c(names(group_colors)[names(group_colors) != "Other"], "Other")
        ) +
        labs(
          # title = paste0(mode_label, " | ", source_label, " loadings: ", comp1, " vs ", comp2)
          title = paste0(mode_label, ": ", source_label, " Loadings")
          ) +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".loadings_", comp1, "vs", comp2, Special_string, ".pdf"
        ),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, X_source, "Group", "Loading") # paste0("X.", X_source)
  } else {
    plot_loadings_side(X_plot, X_source, "Group", "Loading") # paste0("X.", X_source)
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot,Y_source, "Group", "Loading") # paste0("Y.", Y_source)
  } else {
    plot_loadings_side(Y_plot, Y_source, "Group", "Loading") # paste0("Y.", Y_source)
  }
  
}

#### 4. Execute to plot PLS scores colored by cancer type (requires Step 1)
if (1) {
  
  ## Load model metadata
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  ## Load saved variates files
  x.variates.plot <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_X.variates.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  y.variates.plot <- read.delim(
    file = paste0(path.pls, file_tag, Filtered_Tag, "_Y.variates.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Annotate with cancer type
  x.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(x.variates.plot$Score, model$ModelID)]
  y.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(y.variates.plot$Score, model$ModelID)]
  
  ## Helper: scatter plots of scores
  plot_scores_side <- function(df, source_label) {
    
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    n_lineages <- length(unique(na.omit(df$OncotreeLineage)))
    lineage_pal <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_lineages),
      sort(unique(na.omit(df$OncotreeLineage)))
    )
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "comp1"
      comp2_col <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1_col,
          y     = comp2_col,
          color = "OncotreeLineage"
        )
      ) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_pal, name = "Lineage") +
        labs(
          title = paste0(mode_label, ": ", source_label, " Variate Scores"),
          x     = comp1_col,
          y     = comp2_col
        ) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_",
          comp1_col, "vs", comp2_col, ".pdf"
        ),
        plot = p, width = 7, height = 6, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, X_source)
  plot_scores_side(y.variates.plot, Y_source)
  
  ## Helper: boxplots of scores per lineage
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^comp([1-9]|10)$", names(df), value = TRUE)
    
    ## For the multi-facet overview, order lineages by comp1 median
    lineage_order_comp1 <- df %>%
      dplyr::filter(!is.na(OncotreeLineage)) %>%
      dplyr::group_by(OncotreeLineage) %>%
      dplyr::summarise(med = median(comp1, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(med) %>%
      dplyr::pull(OncotreeLineage)
    
    df_long <- df %>%
      dplyr::select(Score, OncotreeLineage, dplyr::all_of(comp_cols)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(comp_cols),
        names_to  = "Component",
        values_to = "Score_value"
      ) %>%
      dplyr::filter(!is.na(OncotreeLineage)) %>%
      dplyr::mutate(
        Component       = factor(Component, levels = comp_cols),
        OncotreeLineage = factor(OncotreeLineage, levels = lineage_order_comp1)
      )
    
    ## Multi-facet overview PDF (ordered by comp1 median)
    p <- ggplot(
      df_long,
      aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)
    ) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, size = 0.3) +
      facet_wrap(~ Component, scales = "free_y", ncol = 2) +
      scale_fill_manual(
        values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(
          length(levels(df_long$OncotreeLineage))
        ),
        guide = "none"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
      labs(
        title = paste0(mode_label, ": ", source_label, " Scores by Cancer Lineage (comps 1–10)"),
        x     = NULL,
        y     = "Score"
      ) +
      theme_bw(base_size = 9) +
      theme(
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
        strip.text    = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines")
      )
    
    ggsave(
      filename = paste0(
        path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_comps1to10.pdf"
      ),
      plot = p, width = 12, height = 18, units = "in", device = cairo_pdf
    )
    
    ## Individual per-comp PDFs, each ordered by that comp's own median
    for (comp in comp_cols) {
      
      lineage_order_comp <- df %>%
        dplyr::filter(!is.na(OncotreeLineage)) %>%
        dplyr::group_by(OncotreeLineage) %>%
        dplyr::summarise(med = median(.data[[comp]], na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(med) %>%
        dplyr::pull(OncotreeLineage)
      
      df_comp <- df_long %>%
        dplyr::filter(Component == comp) %>%
        dplyr::mutate(OncotreeLineage = factor(OncotreeLineage, levels = lineage_order_comp))
      
      p_ind <- ggplot(
        df_comp,
        aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)
      ) +
        geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, size = 0.3) +
        scale_fill_manual(
          values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(
            length(levels(df_comp$OncotreeLineage))
          ),
          guide = "none"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
        labs(
          title = paste0(mode_label, ": ", source_label, " Scores by Cancer Lineage (", comp, ")"),
          x     = NULL,
          y     = "Score"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_", comp, ".pdf"
        ),
        plot = p_ind, width = 10, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_boxplot(x.variates.plot, X_source)
  plot_scores_boxplot(y.variates.plot, Y_source)
  
  ## Wilcoxon rank-sum: each lineage vs. all others, per component, FDR-corrected
  wilcox_lineage_tests <- function(df, source_label) {
    
    comp_cols <- grep("^comp([1-9]|10)$", names(df), value = TRUE)
    lineages  <- unique(na.omit(df$OncotreeLineage))
    
    results <- purrr::map_dfr(comp_cols, function(comp) {
      purrr::map_dfr(lineages, function(lin) {
        
        in_group  <- df[[comp]][!is.na(df$OncotreeLineage) & df$OncotreeLineage == lin]
        out_group <- df[[comp]][!is.na(df$OncotreeLineage) & df$OncotreeLineage != lin]
        
        if (length(in_group) < 3 || length(out_group) < 3) return(NULL)
        
        wt <- wilcox.test(in_group, out_group, exact = FALSE)
        
        data.frame(
          Component  = comp,
          Lineage    = lin,
          n_lineage  = length(in_group),
          median_in  = median(in_group),
          median_out = median(out_group),
          W          = wt$statistic,
          p_value    = wt$p.value
        )
      })
    })
    
    results$p_adj_BH <- p.adjust(results$p_value, method = "BH")
    results <- results %>% dplyr::arrange(p_adj_BH)
    
    write.table(
      x         = results,
      file      = paste0(path.pls, file_tag, Filtered_Tag, "_", source_label, "_wilcox_lineage.txt"),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
    
    results
  }
  
  wilcox_x <- wilcox_lineage_tests(x.variates.plot, X_source) %>%
    dplyr::arrange(Component, desc(median_in))
  
  wilcox_y <- wilcox_lineage_tests(y.variates.plot, Y_source) %>%
    dplyr::arrange(Component, desc(median_in))
  
  ## Quick view of top hits
  head(wilcox_x, 20)
  head(wilcox_y, 20)
  
  ## Save Wilcoxon results tables
  write.table(
    x         = wilcox_x,
    file      = paste0(path.stat, file_tag, Filtered_Tag, "_X.", X_source, "_wilcox_lineage.txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  write.table(
    x         = wilcox_y,
    file      = paste0(path.stat, file_tag, Filtered_Tag, "_Y.", Y_source, "_wilcox_lineage.txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
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
mode_rcca      <- "shrinkage" # ridge or shrinkage
tune_lambda    <- FALSE
lambda1_manual <- 0.20
lambda2_manual <- 0.10

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- c("Myeloid", "Lymphoid")  # e.g. c("Myeloid", "Lymphoid") or character(0)

## Filtered for all three data sets shared lines?
FilteredAll3 <- TRUE # TRUE or FALSE

#### 1. Execute to prep for RCCA
if (1) {
  
  if (FilteredAll3 == TRUE) {
    Filtered_Tag <- "_Filtered3"
  } else {
    Filtered_Tag <- ""
  }
  
  ## Read in CRISPR
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  ## Read in CTRP
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  if (FilteredAll3 == TRUE) {
    
    ## Read in RNAi
    RNAi <- read.delim(
      file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    )
    
    models <- read.delim(
      paste0(path.dm, "Model.csv"), sep = ",", stringsAsFactors = FALSE, check.names = FALSE
    ) %>%
      dplyr::select(ModelID, CCLEName, OncotreeLineage)
    
    RNAi_t <- RNAi %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "CCLEName") %>%
      dplyr::rename_with(~ sub("\\.\\..*", "", .))
    
    RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
      dplyr::select(-CCLEName, -OncotreeLineage) %>%
      tibble::column_to_rownames(var = "ModelID")
    
    ## Filter CRISPR columns to genes shared with RNAi
    shared_genes <- intersect(colnames(CRISPR), colnames(RNAi_t_ModelID))
    message("Shared genes between CRISPR and RNAi: ", length(shared_genes))
    CRISPR <- CRISPR[, shared_genes, drop = FALSE]
    
  }
  
  ## Assign X and Y data
  if (X_source == "CRISPR") X_data <- CRISPR
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "CRISPR") Y_data <- CRISPR
  if (Y_source == "CTRP")   Y_data <- CTRP
  
  if (FilteredAll3 == TRUE) {
    # Three-way intersection: CRISPR, CTRP, and RNAi samples
    ids <- Reduce(intersect, list(
      rownames(X_data),
      rownames(Y_data),
      rownames(RNAi_t_ModelID)
    ))
  } else {
    # Two-way intersection: CRISPR and CTRP only
    ids <- intersect(rownames(X_data), rownames(Y_data))
  }
  
  ## Filter out cell lines belonging to excluded lineages
  if (length(exclude_lineages) > 0) {
    models_filt <- read.csv(paste0(path.dm, "Model.csv"))
    keep_ids <- models_filt$ModelID[!(models_filt$OncotreeLineage %in% exclude_lineages)]
    ids <- intersect(ids, keep_ids)
  }
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  ## Build file_tag
  excl_tag <- if (length(exclude_lineages) > 0) {
    paste0("_excl.", paste(exclude_lineages, collapse = "."))
  } else {
    ""
  }
  
  if (mode_rcca == "ridge") {
    
    if (tune_lambda) {
      
      grid1      <- c(0.10, 0.20, 0.30)
      grid2      <- c(0.05, 0.10, 0.20)
      ncomp_tune <- min(5L, ncomp)
      
      set.seed(999)
      tune_time <- system.time({
        tune.out <- mixOmics::tune.rcc(
          X = X, Y = Y, grid1 = grid1, grid2 = grid2,
          ncomp = ncomp_tune, validation = "loo"
        )
      })
      
      print(tune_time)
      print(tune.out$opt.lambda1)
      print(tune.out$opt.lambda2)
      
      lambda1 <- tune.out$opt.lambda1
      lambda2 <- tune.out$opt.lambda2
      
    } else {
      
      lambda1 <- lambda1_manual
      lambda2 <- lambda2_manual
    }
    
    file_tag <- paste0(
      "RCCA_ridge",
      "_lambda1.", format(lambda1, digits = 3),
      "_lambda2.", format(lambda2, digits = 3),
      "_X.", X_source, "_Y.", Y_source, excl_tag
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    file_tag <- paste0(
      "RCCA_shrinkage",
      "_X.", X_source, "_Y.", Y_source, excl_tag
    )
    
  }
  
}

#### 2. Execute to run RCCA and save output files (requires Step 1)
if (1) {
  
  if (mode_rcca == "ridge") {
    
    message("Running rCCA in ridge mode with lambda1 = ", lambda1, ", lambda2 = ", lambda2)
    
    rcca_fit <- mixOmics::rcc(
      X = X, Y = Y, ncomp = ncomp,
      lambda1 = lambda1, lambda2 = lambda2, method = "ridge"
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    message("Running rCCA in shrinkage mode (automatic lambda estimation).")
    
    rcca_fit <- mixOmics::rcc(
      X = X, Y = Y, ncomp = ncomp, method = "shrinkage"
    )
    
  } else {
    
    stop("mode_rcca must be 'ridge' or 'shrinkage', not: ", mode_rcca)
  }
  
  print(rcca_fit$cor)
  
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
  
  cancor.df <- data.frame(
    comp                  = seq_along(rcca_fit$cor),
    canonical_correlation = rcca_fit$cor
  )
  
  if (!dir.exists(path.rcca)) dir.create(path.rcca, recursive = TRUE)
  
  write.table(x = x.variates,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.variates.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = y.variates,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.variates.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = variates.X.Y, file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.Y.variates.txt"),           sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = x.loadings,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.loadings.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = y.loadings,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.loadings.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = cancor.df,    file = paste0(path.rcca, file_tag, Filtered_Tag, "_canonical_correlations.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
}

#### 3. Execute to plot rCCA loadings (requires Step 1)
if (1) {
  
  X_loadings <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.loadings.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  Y_loadings <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.loadings.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
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
  
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 p53.pathway",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        ),
        group.atp5        = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na          = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na     = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na      = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag          = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
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
    
    percent.nas <- as.data.frame(colMeans(is.na(CTRP_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    df
  }
  
  annotate_crispr <- function(df, side_label) {
    
    gene.info.all <- read.delim(file = paste0(path.general, "Homo_sapiens.gene_info.20251028"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(MDM2|PPM1D|USP7|MDM4|CDKN1A|ATM|TP53|CHEK2|TP53BP1|USP28)$") ~ "04 DDR",
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
    
    percent.nas <- as.data.frame(colMeans(is.na(CRISPR_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    df
  }
  
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_crispr(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_crispr(Y_loadings, "Y")
  
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot(df, aes_string(x = comp1, y = comp2, color = color_col)) +
        geom_point(size = 2.5) +
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col), size = 2
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        labs(title = paste0("rCCA | ", source_label, " loadings: ", comp1, " vs ", comp2)) +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  
}

#### 4. Execute to plot rCCA scores colored by cancer type (requires Step 1 + saved variates)
if (1) {
  
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  x.variates.plot <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.variates.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  y.variates.plot <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.variates.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  x.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(x.variates.plot$Score, model$ModelID)]
  y.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(y.variates.plot$Score, model$ModelID)]
  
  top_lineages_n <- 15
  top_lineages <- names(sort(table(x.variates.plot$OncotreeLineage), decreasing = TRUE))[1:top_lineages_n]
  
  x.variates.plot <- x.variates.plot %>%
    dplyr::mutate(lineage_label = dplyr::if_else(OncotreeLineage %in% top_lineages, OncotreeLineage, "Other"))
  y.variates.plot <- y.variates.plot %>%
    dplyr::mutate(lineage_label = dplyr::if_else(OncotreeLineage %in% top_lineages, OncotreeLineage, "Other"))
  
  lineage_colors <- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(7, "Set2"), "grey70")
  names(lineage_colors) <- c(top_lineages, "Other")
  
  plot_scores_side <- function(df, source_label) {
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "X1"
      comp2_col <- paste0("X", i)
      
      p <- ggplot(df, aes_string(x = comp1_col, y = comp2_col, color = "lineage_label")) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_colors, name = "Lineage") +
        labs(
          title = paste0("rCCA | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
          x = comp1_col, y = comp2_col
        ) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
      
      ggsave(
        filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_", comp1_col, "vs", comp2_col, ".pdf"),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, paste0("X.", X_source))
  plot_scores_side(y.variates.plot, paste0("Y.", Y_source))
  
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^X([1-9]|10)$", names(df), value = TRUE)
    
    df_long <- df %>%
      dplyr::select(Score, OncotreeLineage, dplyr::all_of(comp_cols)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(comp_cols), names_to = "Component", values_to = "Score_value") %>%
      dplyr::filter(!is.na(OncotreeLineage)) %>%
      dplyr::mutate(
        Component = factor(Component, levels = comp_cols),
        OncotreeLineage = factor(
          OncotreeLineage,
          levels = df %>%
            dplyr::filter(!is.na(OncotreeLineage)) %>%
            dplyr::group_by(OncotreeLineage) %>%
            dplyr::summarise(med = median(X1, na.rm = TRUE), .groups = "drop") %>%
            dplyr::arrange(med) %>%
            dplyr::pull(OncotreeLineage)
        )
      )
    
    p <- ggplot(df_long, aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, size = 0.3) +
      facet_wrap(~ Component, scales = "free_y", ncol = 2) +
      scale_fill_manual(
        values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(df_long$OncotreeLineage))),
        guide = "none"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
      labs(
        title = paste0("rCCA | ", source_label, " scores by cancer lineage (comps 1–10)"),
        x = NULL, y = "Score"
      ) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), strip.text = element_text(size = 9, face = "bold"), panel.spacing = unit(0.4, "lines"))
    
    ggsave(
      filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_comps1to10.pdf"),
      plot = p, width = 12, height = 18, units = "in", device = cairo_pdf
    )
    
    for (comp in comp_cols) {
      
      df_comp <- df_long %>% dplyr::filter(Component == comp)
      
      p_ind <- ggplot(df_comp, aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)) +
        geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, size = 0.3) +
        scale_fill_manual(
          values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(df_comp$OncotreeLineage))),
          guide = "none"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
        labs(
          title = paste0("rCCA | ", source_label, " scores — ", comp, " by cancer lineage"),
          x = NULL, y = "Score"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      
      ggsave(
        filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_", comp, ".pdf"),
        plot = p_ind, width = 10, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_boxplot(x.variates.plot, paste0("X.", X_source))
  plot_scores_boxplot(y.variates.plot, paste0("Y.", Y_source))
  
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
mode_rcca      <- "ridge" # ridge or shrinkage
tune_lambda    <- FALSE
lambda1_manual <- 0.20
lambda2_manual <- 0.10

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- character(0)  # e.g. c("Myeloid", "Lymphoid") or character(0)

## Filtered for all three data sets shared lines?
FilteredAll3 <- TRUE # TRUE or FALSE

#### 1. Execute to prep for RCCA
if (1) {
  
  if (FilteredAll3 == TRUE) {
    Filtered_Tag <- "_Filtered3"
  } else {
    Filtered_Tag <- ""
  }
  
  ## Read in RNAi
  RNAi <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  ## Read in CTRP
  CTRP <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  ## Convert RNAi sample nomenclature (CCLEName -> ModelID via Model.csv)
  models <- read.delim(
    paste0(path.dm, "Model.csv"), sep = ",", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::select(ModelID, CCLEName, OncotreeLineage)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName") %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName, -OncotreeLineage) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  if (FilteredAll3 == TRUE) {
    
    ## Read in CRISPR (for shared sample and gene filtering)
    CRISPR <- read.delim(
      file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub("\\.\\..*", "", .))
    
    ## Filter RNAi columns to genes shared with CRISPR
    shared_genes <- intersect(colnames(RNAi_t_ModelID), colnames(CRISPR))
    message("Shared genes between RNAi and CRISPR: ", length(shared_genes))
    RNAi_t_ModelID <- RNAi_t_ModelID[, shared_genes, drop = FALSE]
    
  }
  
  ## Assign X and Y data
  if (X_source == "RNAi") X_data <- RNAi_t_ModelID
  if (X_source == "CTRP") X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP") Y_data <- CTRP
  
  if (FilteredAll3 == TRUE) {
    # Three-way intersection: RNAi, CTRP, and CRISPR samples
    ids <- Reduce(intersect, list(
      rownames(X_data),
      rownames(Y_data),
      rownames(CRISPR)
    ))
  } else {
    # Two-way intersection: RNAi and CTRP only
    ids <- intersect(rownames(X_data), rownames(Y_data))
  }
  
  ## Filter out cell lines belonging to excluded lineages
  if (length(exclude_lineages) > 0) {
    keep_ids <- models$ModelID[!(models$OncotreeLineage %in% exclude_lineages)]
    ids <- intersect(ids, keep_ids)
  }
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- lapply(X, as.numeric)
  Y[] <- lapply(Y, as.numeric)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  ## Build file_tag
  excl_tag <- if (length(exclude_lineages) > 0) {
    paste0("_excl.", paste(exclude_lineages, collapse = "."))
  } else {
    ""
  }
  
  if (mode_rcca == "ridge") {
    
    if (tune_lambda) {
      
      grid1      <- c(0.10, 0.20, 0.30)
      grid2      <- c(0.05, 0.10, 0.20)
      ncomp_tune <- min(5L, ncomp)
      
      set.seed(999)
      tune_time <- system.time({
        tune.out <- mixOmics::tune.rcc(
          X = X, Y = Y, grid1 = grid1, grid2 = grid2,
          ncomp = ncomp_tune, validation = "loo"
        )
      })
      
      print(tune_time)
      print(tune.out$opt.lambda1)
      print(tune.out$opt.lambda2)
      
      lambda1 <- tune.out$opt.lambda1
      lambda2 <- tune.out$opt.lambda2
      
    } else {
      
      lambda1 <- lambda1_manual
      lambda2 <- lambda2_manual
    }
    
    file_tag <- paste0(
      "RCCA_ridge",
      "_lambda1.", format(lambda1, digits = 3),
      "_lambda2.", format(lambda2, digits = 3),
      "_X.", X_source, "_Y.", Y_source, excl_tag
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    file_tag <- paste0(
      "RCCA_shrinkage",
      "_X.", X_source, "_Y.", Y_source, excl_tag
    )
    
  }
  
}

#### 2. Execute to run RCCA and save output files (requires Step 1)
if (1) {
  
  if (mode_rcca == "ridge") {
    
    message("Running rCCA in ridge mode with lambda1 = ", lambda1, ", lambda2 = ", lambda2)
    
    rcca_fit <- mixOmics::rcc(
      X = X, Y = Y, ncomp = ncomp,
      lambda1 = lambda1, lambda2 = lambda2, method = "ridge"
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    message("Running rCCA in shrinkage mode (automatic lambda estimation).")
    
    rcca_fit <- mixOmics::rcc(
      X = X, Y = Y, ncomp = ncomp, method = "shrinkage"
    )
    
  } else {
    
    stop("mode_rcca must be 'ridge' or 'shrinkage', not: ", mode_rcca)
  }
  
  print(rcca_fit$cor[1:ncomp])
  
  x.variates <- data.frame(rcca_fit$variates$X) %>% tibble::rownames_to_column(var = "Score")
  y.variates <- data.frame(rcca_fit$variates$Y) %>% tibble::rownames_to_column(var = "Score")
  
  x.loadings <- data.frame(rcca_fit$loadings$X) %>% tibble::rownames_to_column(var = "Loading") %>% dplyr::arrange(X1)
  y.loadings <- data.frame(rcca_fit$loadings$Y) %>% tibble::rownames_to_column(var = "Loading") %>% dplyr::arrange(X1)
  
  variates.X.Y <- merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = c(paste0(".", X_source), paste0(".", Y_source))
  )
  
  cancor.df <- data.frame(
    comp                  = seq_along(rcca_fit$cor),
    canonical_correlation = rcca_fit$cor
  )
  
  if (!dir.exists(path.rcca)) dir.create(path.rcca, recursive = TRUE)
  
  write.table(x = x.variates,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.variates.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = y.variates,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.variates.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = variates.X.Y, file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.Y.variates.txt"),           sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = x.loadings,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.loadings.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = y.loadings,   file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.loadings.txt"),             sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = cancor.df,    file = paste0(path.rcca, file_tag, Filtered_Tag, "_canonical_correlations.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
}

#### 3. Execute to plot rCCA loadings (requires Step 1)
if (1) {
  
  X_loadings <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.loadings.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  Y_loadings <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.loadings.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!exists("RNAi_mat") || !exists("CTRP_mat")) {
    
    RNAi_mat <- read.delim(
      file = paste0(path.dm, "D2_combined_gene_dep_scores.csv"),
      sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
    ) %>%
      dplyr::rename_with(~ sub(" .*", "", .))
    
    CTRP_mat <- read.delim(
      file = paste0(path.ctrp, "ctrpv2.wide.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
  }
  
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(selumetinib|PD318088|trametinib|RAF265|dabrafenib|regorafenib|PLX\\-4720|PLX\\-4032|sorafenib|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
          stringr::str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib|PD 153035)$") ~ "02 EGFRi.HER2i",
          stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 p53.pathway",
          stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
          stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",
          detect(drug.target, "BCL2") & !stringr::str_detect(Loading, ":") ~ "07 BCL2+i",
          TRUE ~ NA_character_
        ),
        group.atp5        = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na          = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na     = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na      = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag          = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
      ) %>%
      dplyr::arrange(dplyr::desc(group.na))
    
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
    
    percent.nas <- as.data.frame(colMeans(is.na(CTRP_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    df
  }
  
  annotate_rnai <- function(df, side_label) {
    
    gene.info.all <- read.delim(file = paste0(path.general, "Homo_sapiens.gene_info.20251028"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- sub("\\.\\..*$", "", df$Loading)
    df <- merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
    df <- df %>%
      dplyr::mutate(
        group = dplyr::case_when(
          stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)") ~ "01 BRAF sig",
          stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")     ~ "02 EGFR sig",
          stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$") ~ "03 ferropt",
          stringr::str_detect(Loading, "^(MDM2|PPM1D|USP7|MDM4|CDKN1A|ATM|TP53|CHEK2|TP53BP1|USP28)$") ~ "04 DDR",
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
    
    percent.nas <- as.data.frame(colMeans(is.na(RNAi_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    df
  }
  
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_rnai(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_rnai(Y_loadings, "Y")
  
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot(df, aes_string(x = comp1, y = comp2, color = color_col)) +
        geom_point(size = 2.5) +
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col), size = 2
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        labs(title = paste0("rCCA | ", source_label, " loadings: ", comp1, " vs ", comp2)) +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_loadings_side(X_plot, paste0("X.", X_source), "group", "Loading")
  plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "Loading")
  
}

#### 4. Execute to plot rCCA scores colored by cancer type (requires Step 1 + saved variates)
if (1) {
  
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  x.variates.plot <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_X.variates.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  y.variates.plot <- read.delim(file = paste0(path.rcca, file_tag, Filtered_Tag, "_Y.variates.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  x.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(x.variates.plot$Score, model$ModelID)]
  y.variates.plot$OncotreeLineage <- model$OncotreeLineage[match(y.variates.plot$Score, model$ModelID)]
  
  top_lineages_n <- 15
  top_lineages <- names(sort(table(x.variates.plot$OncotreeLineage), decreasing = TRUE))[1:top_lineages_n]
  
  x.variates.plot <- x.variates.plot %>%
    dplyr::mutate(lineage_label = dplyr::if_else(OncotreeLineage %in% top_lineages, OncotreeLineage, "Other"))
  y.variates.plot <- y.variates.plot %>%
    dplyr::mutate(lineage_label = dplyr::if_else(OncotreeLineage %in% top_lineages, OncotreeLineage, "Other"))
  
  lineage_colors <- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(7, "Set2"), "grey70")
  names(lineage_colors) <- c(top_lineages, "Other")
  
  plot_scores_side <- function(df, source_label) {
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "X1"
      comp2_col <- paste0("X", i)
      
      p <- ggplot(df, aes_string(x = comp1_col, y = comp2_col, color = "lineage_label")) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_colors, name = "Lineage") +
        labs(
          title = paste0("rCCA | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
          x = comp1_col, y = comp2_col
        ) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
      
      ggsave(
        filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_", comp1_col, "vs", comp2_col, ".pdf"),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, paste0("X.", X_source))
  plot_scores_side(y.variates.plot, paste0("Y.", Y_source))
  
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^X([1-9]|10)$", names(df), value = TRUE)
    
    df_long <- df %>%
      dplyr::select(Score, OncotreeLineage, dplyr::all_of(comp_cols)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(comp_cols), names_to = "Component", values_to = "Score_value") %>%
      dplyr::filter(!is.na(OncotreeLineage)) %>%
      dplyr::mutate(
        Component = factor(Component, levels = comp_cols),
        OncotreeLineage = factor(
          OncotreeLineage,
          levels = df %>%
            dplyr::filter(!is.na(OncotreeLineage)) %>%
            dplyr::group_by(OncotreeLineage) %>%
            dplyr::summarise(med = median(X1, na.rm = TRUE), .groups = "drop") %>%
            dplyr::arrange(med) %>%
            dplyr::pull(OncotreeLineage)
        )
      )
    
    p <- ggplot(df_long, aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, size = 0.3) +
      facet_wrap(~ Component, scales = "free_y", ncol = 2) +
      scale_fill_manual(
        values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(df_long$OncotreeLineage))),
        guide = "none"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
      labs(
        title = paste0("rCCA | ", source_label, " scores by cancer lineage (comps 1–10)"),
        x = NULL, y = "Score"
      ) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), strip.text = element_text(size = 9, face = "bold"), panel.spacing = unit(0.4, "lines"))
    
    ggsave(
      filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_comps1to10.pdf"),
      plot = p, width = 12, height = 18, units = "in", device = cairo_pdf
    )
    
    for (comp in comp_cols) {
      
      df_comp <- df_long %>% dplyr::filter(Component == comp)
      
      p_ind <- ggplot(df_comp, aes(x = OncotreeLineage, y = Score_value, fill = OncotreeLineage)) +
        geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, size = 0.3) +
        scale_fill_manual(
          values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(df_comp$OncotreeLineage))),
          guide = "none"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.3) +
        labs(
          title = paste0("rCCA | ", source_label, " scores — ", comp, " by cancer lineage"),
          x = NULL, y = "Score"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      
      ggsave(
        filename = paste0(path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".scores_boxplot_", comp, ".pdf"),
        plot = p_ind, width = 10, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_boxplot(x.variates.plot, paste0("X.", X_source))
  plot_scores_boxplot(y.variates.plot, paste0("Y.", Y_source))
  
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

## Cell lines excluded in the upstream PLS/rCCA runs (must match what was used)
## Set to character(0) if no filtering was applied
exclude_lineages_1 <- character(0) # for file1 (X1/Y1) e.g. c("Myeloid", "Lymphoid")
exclude_lineages_2 <- character(0) # for file2 (X2/Y2) e.g. c("Myeloid", "Lymphoid")

## Filtered for all three data sets shared lines? (must match upstream PLS runs)
FilteredAll3_1 <- TRUE  # for file1 (X1/Y1)
FilteredAll3_2 <- TRUE  # for file2 (X2/Y2)

### Create scatter plot and generate table of distances and theta
if(1) {
  
  ## Build excl tags to match upstream file_tag construction
  excl_tag_1 <- if (length(exclude_lineages_1) > 0) {
    paste0("_excl.", paste(exclude_lineages_1, collapse = "."))
  } else {
    ""
  }
  excl_tag_2 <- if (length(exclude_lineages_2) > 0) {
    paste0("_excl.", paste(exclude_lineages_2, collapse = "."))
  } else {
    ""
  }
  
  ## Build Filtered_Tags to match upstream
  filtered_tag_1 <- if (FilteredAll3_1) "_Filtered3" else ""
  filtered_tag_2 <- if (FilteredAll3_2) "_Filtered3" else ""
  
  if (DimRedTec == "PLS") {
    file1_tag <- paste0("PLS_Mode.", mode, "_X.", X1_source, "_Y.", Y1_source, excl_tag_1)
    file2_tag <- paste0("PLS_Mode.", mode, "_X.", X2_source, "_Y.", Y2_source, excl_tag_2)
    path_in <- path.pls
  }
  
  if (DimRedTec == "rCCA") {
    file1_tag <- paste0("RCCA_shrinkage_X.", X1_source, "_Y.", Y1_source, excl_tag_1)
    file2_tag <- paste0("RCCA_shrinkage_X.", X2_source, "_Y.", Y2_source, excl_tag_2)
    path_in <- path.rcca
  }
  
  ## Read in data
  X1_loadings <- read.delim(
    file = paste0(path_in, file1_tag, filtered_tag_1, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading)) %>%
    dplyr::select(Loading, paste0("comp", 1:10))
  
  X2_loadings <- read.delim(
    file = paste0(path_in, file2_tag, filtered_tag_2, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading)) %>%
    dplyr::select(Loading, paste0("comp", 1:10))
  
  ## Find max abs()
  X1_loadings_max <- X1_loadings %>%
    tidyr::pivot_longer(
      cols      = paste0("comp", 1:10),
      names_to  = "component",
      values_to = "loading"
    ) %>%
    dplyr::mutate(abs_loading = abs(loading)) %>%
    dplyr::group_by(Loading) %>%
    dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      component_CRISPR   = component,
      loading_CRISPR     = loading,
      abs_loading_CRISPR = abs_loading
    )
  
  X2_loadings_max <- X2_loadings %>%
    tidyr::pivot_longer(
      cols      = paste0("comp", 1:10),
      names_to  = "component",
      values_to = "loading"
    ) %>%
    dplyr::mutate(abs_loading = abs(loading)) %>%
    dplyr::group_by(Loading) %>%
    dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      component_RNAi   = component,
      loading_RNAi     = loading,
      abs_loading_RNAi = abs_loading
    )
  
  ## Merge
  Max <- merge(X1_loadings_max, X2_loadings_max, by = "Loading")
  
  xcol <- names(Max)[7]   # 7th column name
  ycol <- names(Max)[4]   # 4th column name
  
  Max <- Max %>%
    dplyr::mutate(
      theta_rad = atan2(.data[[ycol]], .data[[xcol]]),
      theta_deg = theta_rad * 180 / pi,
      r         = sqrt(.data[[xcol]]^2 + .data[[ycol]]^2)
    )
  
  ## Build output tag for file/plot names
  out_tag <- paste0(
    mode,
    "_X1_", X1_source, "_vs_", Y1_source, excl_tag_1, filtered_tag_1,
    "_X2_", X2_source, "_vs_", Y2_source, excl_tag_2, filtered_tag_2
  )
  
  write.table(
    x    = Max,
    file = paste0(path.max, "MaxLoadingsDF_", out_tag, ".txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )
  
  ## Plot
  plot_df <- data.frame(
    x     = Max[[7]],
    y     = Max[[4]],
    gene  = Max$Loading,
    theta = Max$theta_deg
  ) %>%
    dplyr::mutate(
      angle_group = dplyr::case_when(
        theta < 30  ~ "RNAi",
        theta <= 60 ~ "Neutral",
        TRUE        ~ "CRISPR"
      ),
      angle_group = factor(angle_group, levels = c("CRISPR", "Neutral", "RNAi"))
    )
  
  top5 <- plot_df %>%
    dplyr::mutate(r = sqrt(x^2 + y^2)) %>%
    dplyr::group_by(angle_group) %>%
    dplyr::slice_max(r, n = 5, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  p <- ggplot(plot_df, aes(x = x, y = y, color = angle_group)) +
    geom_point(size = 0.075, alpha = 0.3) +
    geom_abline(
      slope = tan(pi/6), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    geom_abline(
      slope = tan(pi/3), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    geom_text_repel(
      data          = top5,
      aes(label     = gene),
      size          = 2.5,
      fontface      = "italic",
      show.legend   = FALSE,
      max.overlaps  = 20,
      segment.size  = 0.3,
      segment.color = "grey50"
    ) +
    scale_color_manual(
      values = c("RNAi" = "#5E2F80", "Neutral" = "#BDBDBD", "CRISPR" = "#D47D37"),
      name   = "Bias"
    ) +
    labs(
      x = "Max Absolute PLS-C Loading RNAi",
      y = "Max Absolute PLS-C Loading CRISPR"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic(base_size = 10) +
    theme(
      legend.background = element_blank(),
      legend.key        = element_blank()
    )
  
  ggsave(
    filename = paste0(path.plots, "MaxLoadingsDF_", out_tag, "_Scatter.pdf"),
    plot = p, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
}

#### Create Estelles GSEA Plot

if (1) {
  
  ## Load msigdb pathways
  msig_df <- load.MSigDB(species = 'Homo sapiens')
  
  gsea_list <- get.MSigDB.genesets(
    msig_df = rbind(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "C5")
    ),
    genesets = c("BP")
  )
  
  keyword_groups <- list(
    IMMUNE =
      c("INFLAME", "IMMUNE", "INTERLEUKIN", "LEUKOCYTE", "CD4",
        "MACROPHAGE", "NEUTROPHILE"),
    PROTEIN_PROCESSING =
      c("PEPTIDE", "AMINO_ACID", "UBIQUITIN", "UBIQUITINATION"),
    VIRAL_PROCESSES =
      c("VIRAL", "SYMBIOTIC", "DSRNA"),
    STRESS_RESPONSE =
      c("DNA_DAMAGE", "APOPTOTIC", "REPAIR", "HYPOXIA", "STRESS"),
    METABOLIC_PATHWAY =
      c("CATABOLIC", "ATP", "POLYSACCHARIDE", "FRUCTOSE",
        "GLYCOSYLATION", "GLYCOGEN", "BIOSYNTHESIS", "LIPID"),
    MITOCHONDRIA =
      c("MITOCHONDRIAL", "MITOCHONDRION"),
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
      pattern_regex <- paste(pattern_vec, collapse = "|")
      hit_idx <- grepl(pattern_regex, pathway_names, ignore.case = TRUE)
      genes <- unique(unlist(gsea_list[hit_idx], use.names = FALSE))
      genes
    }
  )
  
  sapply(keyword_to_genes, length)
  
  ## Long data frame: one row per (gene, keyword_group)
  keyword_gene_df <- purrr::imap_dfr(
    keyword_to_genes,
    ~ dplyr::tibble(
      Loading       = .x,
      keyword_group = .y
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
    geom_density_2d(
      aes(group = keyword_group),
      linewidth = 0.5,
      alpha     = 0.8
    ) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(
      yintercept = c(30, 60),
      linetype   = "dotted",
      color      = "grey40",
      linewidth  = 0.4
    ) +
    facet_grid(
      . ~ keyword_group,
      scales = "free_x",
      space  = "free_x"
    ) +
    scale_color_manual(values = my_group_colors, guide = "none") +
    theme_classic() +
    scale_x_continuous(limits = c(0, 0.07))
  
  print(p)
  
  ggsave(
    filename = paste0(path.plots, "MaxLoadingsDF_", out_tag, "_GSEA.pdf"),
    plot = p, width = 13, height = 9, units = "in", device = cairo_pdf
  )
  
}

##### GSEA #####

## Set OS (for swapping between personal and workstation)
OS <- "Mac" # Linux or Mac

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
input.path <- paste0(path.max, "MaxLoadingsDF_canonical_X1_CRISPR_vs_CTRP_Filtered3_X2_RNAi_vs_CTRP_Filtered3.txt")

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
OS <- "Mac" # Linux or Mac

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
  "FGSEA_MaxLoadingsDF_canonical_X1_CRISPR_vs_CTRP_Filtered3_X2_RNAi_vs_CTRP_Filtered3.txt"
)

FGSEA_results_rnk <- read.table(
  file   = path.input,
  sep    = "\t",
  header = TRUE
)

### Prep for GSEA^2 

## Define keyword groups
keyword_groups <- list(
  NEURO = 
    c("NEURO", "NEUROTRANSMITTER", "SYNAPTIC", "VOLTAGE", "AXON",
      "CEREBRAL", "CORTEX", "DENDRITE", "GLUTAMATE"),
  # IMMUNE =
  # c("INFLAME", "IMMUNE", "INTERLEUKIN", "LEUKOCYTE", "CD4",
  # "MACROPHAGE", "NEUTROPHILE"),
  # KINASE_ACTIVITY =
  # c("MAPK", "KINASE", "GTP", "TYROSINE"),
  # CELL_DIFFERENTIATION =
  # c("KERATINOCYTE", "DIFFERENTIATION"),
  CELL_CELL_INTERACTION =
    c("ADHESION", "ADHERENS", "CELL-CELL", "COMMUNICATION"),
  # CELL_PROLIFERATION =
    # c("PROLIFERATION"),
  PROTEIN_PROCESSING =
    c("PEPTIDE", "AMINO_ACID", "UBIQUITIN", "UBIQUITINATION"),
  STRESS_RESPONSE =
    c("DNA_DAMAGE", "APOPTOTIC", "REPAIR", "HYPOXIA", "STRESS"),
  # METABOLIC_PATHWAY =
  #   c("CATABOLIC", "ATP", "POLYSACCHARIDE", "FRUCTOSE",
  #     "GLYCOSYLATION", "GLYCOGEN", "BIOSYNTHESIS", "LIPID"),
  # MITOCHONDRIA =
    # c("MITOCHONDRIAL", "MITOCHONDRION"),
  # TRANSLATION =
    # c("RIBOSOME", "RRNA", "TRNA", "TRANSLATION", "RIBONUCLEOPROTEIN"),
  DNA_TRANSCRIPTION =
    c("MRNA", "TRANSCRIPTION", "POLYMERASE", "TRANSCRIBED", "GENE_EXPRESSION"),
  CELL_CYCLE =
    c("CELL_CYCLE", "MITOTIC", "DNA_REPLICATION", "CHROMOSOME_SEGREGATION",
      "CHROMATID_SEGREGATION", "SPINDLE", "CELL_DIVISION",
      "KINETOCHORE", "CENTRIOLE", "ANAPHASE"),
  VIRAL_PROCESSES =
    c("VIRAL", "SYMBIOTIC", "DSRNA"),
  # ORGANELLE_TRANSPORT =
    # c("ENDOPLASMIC_RETICULUM", "GOLGI", "VACUOLE"),
  EPIGENETIC =
    c("HISTONE", "NUCLEOSIDE", "DEMETHYLATION", "METHYLATION", "EPIGENETIC")
)

## Get ranks per keyword group 
get_ranks_for_keywords <- function(results_df, keywords) {
  pattern <- paste(keywords, collapse = "|")
  results_df %>%
    dplyr::filter(stringr::str_detect(pathway, pattern)) %>%
    dplyr::pull(rank)
}

## Build GSEA^2 data frame
data_df <- purrr::imap_dfr(
  keyword_groups,
  ~ tibble::tibble(
    Category = .y,
    Value    = get_ranks_for_keywords(FGSEA_results_rnk, .x)
  )
)

## print(table(data_df$Category))

## KS Test (enrichment + signed ordering)
max_rank <- max(FGSEA_results_rnk$rank, na.rm = TRUE)

ks_results <- data_df %>%
  dplyr::group_by(Category) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    mean_rank = mean(Value, na.rm = TRUE),
    p_ks      = if (n > 1) {
      scaled_vals <- Value / max_rank
      stats::ks.test(scaled_vals, "punif")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    direction = dplyr::if_else(mean_rank <= max_rank / 2, "+", "-"),
    logp        = dplyr::if_else(is.na(p_ks), NA_real_, -log10(p_ks)),
    signed_logp = dplyr::if_else(direction == "+", logp, -logp),
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

## print(ks_results)

ordered_levels <- ks_results$Category

data_df$Category <- factor(
  data_df$Category,
  levels = ordered_levels
)

# Color scale
# Stepwise gradient: endpoints are darkest, middle is lightest.
# Steps are evenly spaced (rank-based), ignoring actual p-value distances.
# non-significant (p > 0.05) : gray #888780
#
# Within each arm, categories are ranked by |signed_logp| (most extreme = rank 1)
# and t_val = (rank - 1) / (n - 1) so rank 1 → t=0 (darkest), rank n → t=1 (lightest).
#
# ── OPTION A: Dark Red / Blue ──────────────────────────────────────────────────
#   RNAi arm   (signed_logp < 0): dark red  #7a1515 → light red  #f0a0a0
#   CRISPR arm (signed_logp > 0): light blue #93b8de → dark blue #1a3a6b
#
# ── OPTION B: Dark Purple / Dark Orange ────────────────────────────────────────
#   RNAi arm   (signed_logp < 0): dark purple #5E2F80 → light purple #c4a0e0
#   CRISPR arm (signed_logp > 0): light orange #f0c080 → dark orange #D47D37

hex_interp <- function(t, r1, g1, b1, r2, g2, b2) {
  sprintf("#%02X%02X%02X",
          round(r1 + t * (r2 - r1)),
          round(g1 + t * (g2 - g1)),
          round(b1 + t * (b2 - b1))
  )
}

color_map <- ks_results %>%
  dplyr::select(Category, signed_logp, p_ks) %>%
  dplyr::mutate(
    sig = !is.na(p_ks) & p_ks <= 0.05,
    # rank within each arm by distance from zero (most extreme = rank 1)
    rank_neg = dplyr::if_else(sig & signed_logp < 0,
                              dplyr::dense_rank(signed_logp),   # most negative gets rank 1
                              NA_integer_),
    rank_pos = dplyr::if_else(sig & signed_logp > 0,
                              dplyr::dense_rank(-signed_logp),  # most positive gets rank 1
                              NA_integer_),
    n_neg = sum(!is.na(rank_neg)),
    n_pos = sum(!is.na(rank_pos)),
    dot_color = dplyr::case_when(
      # Non-significant → gray
      !sig ~ "#888780",
      
      # ── OPTION A: Dark Red / Blue (comment out to disable) ──────────────────
      # # RNAi arm: rank 1 = darkest red #7a1515, rank n = lightest red #f0a0a0
      # signed_logp < 0 ~ hex_interp(
      #   (rank_neg - 1) / pmax(n_neg - 1, 1),
      #   0x7a, 0x15, 0x15,   # dark red
      #   0xf0, 0xa0, 0xa0    # light red
      # ),
      # # CRISPR arm: rank 1 = darkest blue #1a3a6b, rank n = lightest blue #93b8de
      # signed_logp > 0 ~ hex_interp(
      #   (rank_pos - 1) / pmax(n_pos - 1, 1),
      #   0x1a, 0x3a, 0x6b,   # dark blue
      #   0x93, 0xb8, 0xde    # light blue
      # ),
      
      # ── OPTION B: Dark Purple / Dark Orange (comment out to disable) ─────────
      # RNAi arm: rank 1 = darkest purple #5E2F80, rank n = lightest purple #c4a0e0
      signed_logp < 0 ~ hex_interp(
        (rank_neg - 1) / pmax(n_neg - 1, 1),
        0x5E, 0x2F, 0x80,   # dark purple
        0xc4, 0xa0, 0xe0    # light purple
      ),
      # CRISPR arm: rank 1 = darkest orange #D47D37, rank n = lightest orange #f0c080
      signed_logp > 0 ~ hex_interp(
        (rank_pos - 1) / pmax(n_pos - 1, 1),
        0xD4, 0x7D, 0x37,   # dark orange
        0xf0, 0xc0, 0x80    # light orange
      ),
      
      TRUE ~ "#888780"
    )
  )

# Named vector: Category → hex color
category_colors <- tibble::deframe(
  dplyr::select(color_map, Category, dot_color)
)

custom_labels <- c(
  NEURO                  = "Neuro",
  # IMMUNE                 = "Immune",
  # KINASE_ACTIVITY        = "Kinase Activity",
  CELL_DIFFERENTIATION   = "Cell Differentiation",
  CELL_CELL_INTERACTION  = "Cell\u2013Cell Interaction",
  # CELL_PROLIFERATION     = "Cell Proliferation",
  PROTEIN_PROCESSING     = "Protein Processing",
  STRESS_RESPONSE        = "Stress Response",
  # METABOLIC_PATHWAY      = "Metabolic Pathway",
  # MITOCHONDRIA           = "Mitochondria",
  # TRANSLATION            = "Translation",
  DNA_TRANSCRIPTION      = "DNA Transcription",
  CELL_CYCLE             = "Cell Cycle",
  EPIGENETIC             = "Epigenetic",
  # ORGANELLE_TRANSPORT    = "Organelle Transport",
  VIRAL_PROCESSES        = "Viral Processes"
)

right_labels <- ks_results %>%
  dplyr::mutate(
    right_lab = paste0(
      " p = ",
      signif(p_ks, 3),
      " (", direction, ")" #sig_flag
    )
  ) %>%
  dplyr::select(Category, right_lab) %>%
  tibble::deframe()

# Compute y-positions for the two dividing lines.
# ordered_levels runs from most-negative signed_logp (bottom of plot = level 1)
# to most-positive (top). In ggplot discrete y, level 1 = y=1, level 2 = y=2, etc.
# We want lines between:
#   (a) last sig RNAi category and first non-sig category
#   (b) last non-sig category and first sig CRISPR category
sig_rnai  <- ks_results %>% dplyr::filter(!is.na(p_ks) & p_ks <= 0.05 & signed_logp < 0)
sig_crispr <- ks_results %>% dplyr::filter(!is.na(p_ks) & p_ks <= 0.05 & signed_logp > 0)

# Position of each category on the y-axis (1 = bottom level in ordered_levels)
cat_pos <- setNames(seq_along(ordered_levels), ordered_levels)

# Line between last RNAi-sig and first non-sig (or CRISPR-sig)
rnai_top_pos   <- max(cat_pos[sig_rnai$Category])
crispr_bot_pos <- min(cat_pos[sig_crispr$Category])

hline_y <- c(rnai_top_pos + 0.5, crispr_bot_pos - 0.5)

# Plot
plt <- ggplot2::ggplot(data_df, ggplot2::aes(x = Value, y = Category)) +
  ggplot2::geom_hline(
    yintercept = hline_y,
    linetype   = "dotted",
    linewidth  = 0.4,
    color      = "grey40"
  ) +
  ggplot2::geom_jitter(
    height = 0.2,
    width  = 0,
    ggplot2::aes(color = Category),
    size   = 1,
    shape  = 16
  ) +
  ggplot2::scale_color_manual(values = category_colors) +
  ggplot2::scale_y_discrete(
    labels   = custom_labels,
    sec.axis = ggplot2::dup_axis(
      labels = right_labels[levels(data_df$Category)],
      name   = ""
    )
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(max_rank * 0.15, max_rank * 0.85),
    labels = c("Enriched in\nCRISPR (+)", "Enriched in\nRNAi (-)")
  ) +
  ggplot2::labs(x = "Rank", y = "") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y.left   = ggplot2::element_text(size = 7),
    axis.text.y.right  = ggplot2::element_text(size = 7, hjust = 0),
    axis.text.x        = ggplot2::element_text(size = 8),
    axis.ticks.x       = ggplot2::element_blank(),
    legend.position    = "none",
    panel.border       = ggplot2::element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid.major   = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = paste0(
    path.plots,
    "GSEA_sq_",
    gsub(".txt", "V2.pdf", basename(path.input))
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
soft_power      <- 4L     # transforms correlation matrix & determines overall connectivity. the higher the value, the more strong correlations are emphasized and weaker are suppressed. This is determined by the scale-free topology fit (below)
deep_Split      <- 4      # [0:4], determines how aggressively the dendogram is cut into initial clusters. higher = more aggressive splitting = more modules detected, less in grey
min_module_sz   <- 30L    # modules bellow this size get assigned to grey
merge_CutHeight <- 0.25   # after modules are built, any two modules who correlate above this threshold get merged. Higher = more aggressive merging = fewer final modules. 0.25 = modules >75% similar get collapsed.

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

#### Run to investigate soft power option (CRISPR = 4, RNAi = 3)
if (0) {
  
  ## Set range of powers and run
  powers <- c(1:20)
  
  sft_CRISPR <- pickSoftThreshold(CRISPR_common,
                                  powerVector = powers,
                                  verbose = 5)
  
  ## Plot
  pdf(paste0(path.plots, "soft_power_selection_CRISPR.pdf"), width = 10, height = 5)
  
  par(mfrow = c(1,2))
  
  plot(sft_CRISPR$fitIndices[,1], -sign(sft_CRISPR$fitIndices[,3])*sft_CRISPR$fitIndices[,2],
       xlab="Soft Power", ylab="Scale Free Topology R²",
       main="Scale independence")
  abline(h=0.8, col="red")
  
  plot(sft_CRISPR$fitIndices[,1], sft_CRISPR$fitIndices[,5],
       xlab="Soft Power", ylab="Mean Connectivity",
       main="Mean connectivity")
  
  dev.off()
}

#### Run WGCNA (1) or Read in WGCNA object (0)
if (0) {
  
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
    mergeCutHeight     = merge_CutHeight, # increase to merge similar modules
    numericLabels      = FALSE,
    pamRespectsDendro  = TRUE,
    verbose            = 3,
    deepSplit          = deep_Split
  )
  
  ## Save WGCNA object
  saveRDS(
    net_CRISPR,
    file = paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  table(net_CRISPR$colors)
  
} else {
  
  ## Load in WGCNA object
  net_CRISPR <- readRDS(paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  table(net_CRISPR$colors)
}

#### Perform correlation on all non grey modules and plot
if (1) {
  
  ## Extract module colors and gene tree
  moduleColors_CRISPR <- net_CRISPR$colors
  table(moduleColors_CRISPR)
  
  ## Manual color/name remapping
  all_modules <- unique(moduleColors_CRISPR[moduleColors_CRISPR != "grey"])
  
  name_remap <- c(
    # "green" = "Cell Cycle",
    # "black" = "DNA Repair",
    # "red" = ""
  )
  
  color_remap <- c(
    "green" = "#228B22",
    "red" = "#bb0a1e",
    "yellow" = "#FFDA03",
    "brown" = "#481F01"
  )
  
  ## Fill unmapped modules with original name
  for (m in all_modules) {
    if (!m %in% names(name_remap)) name_remap[m] <- m
  }
  
  ## Get all new names and fill unmapped colors with original R color
  all_new_names <- name_remap[all_modules]
  for (n in all_new_names) {
    if (!n %in% names(color_remap)) color_remap[n] <- n
  }
  
  ## Get names of genes assigned to clusters (remove "grey" genes)
  non_grey_genes <- colnames(CRISPR_common)[moduleColors_CRISPR != "grey"]
  length(non_grey_genes)
  
  CRISPR_ng <- CRISPR_common[, non_grey_genes, drop = FALSE]
  
  ## Perform correlation without grey genes
  cor_CRISPR <- stats::cor(
    CRISPR_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Apply remaps
  mod_ng_orig  <- moduleColors_CRISPR[non_grey_genes]  # original WGCNA color strings
  mod_ng_name  <- name_remap[mod_ng_orig]              # remapped display names
  mod_ng_color <- color_remap[mod_ng_name]             # remapped hex colors (keyed off new names)
  names(mod_ng_name)  <- non_grey_genes
  names(mod_ng_color) <- non_grey_genes
  
  ## Reorder by display name
  gene_order     <- order(mod_ng_name)
  cor_CRISPR_ord <- cor_CRISPR[gene_order, gene_order]
  
  mod_ng_name_ord  <- mod_ng_name[gene_order]
  mod_ng_color_ord <- mod_ng_color[gene_order]
  
  ## Build color mapping: display name -> hex color
  module_col <- setNames(mod_ng_color_ord, mod_ng_name_ord)
  module_col <- module_col[!duplicated(names(module_col))]
  
  col_fun <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  )
  
  p_crispr <- ComplexHeatmap::Heatmap(
    cor_CRISPR_ord,
    name = "\nPearson r\n",
    col  = col_fun,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_split    = mod_ng_name_ord,
    column_split = mod_ng_name_ord,
    row_title    = NULL,
    column_title = NULL,
    row_gap    = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    heatmap_legend_param = list(
      title_gp      = grid::gpar(fontsize = 20, fontface = "bold"),
      labels_gp     = grid::gpar(fontsize = 20),
      title_gap     = grid::unit(20, "mm"),
      legend_height = grid::unit(40, "mm")
    ),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_name_ord,
      col    = list(Module = module_col),
      show_legend          = TRUE,
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        Module = list(
          title_gp    = grid::gpar(fontsize = 20, fontface = "bold"),
          labels_gp   = grid::gpar(fontsize = 16),
          grid_height = grid::unit(10, "mm"),
          grid_width  = grid::unit(8, "mm")
        )
      )
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_name_ord,
      col    = list(Module = module_col),
      show_legend          = FALSE,
      show_annotation_name = FALSE
    ),
    use_raster    = FALSE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(
      path.plots,
      "HEATMAP_WGCNA_CRISPR_SoftPower_", soft_power,
      "_MinModuleSize_", min_module_sz,
      "_mergeCutHeight_", merge_CutHeight,
      "_deepSplit_", deep_Split,
      "_AllClusters.png"
    ),
    width  = 2500,
    height = 2500,
    res    = 200
  )
  ComplexHeatmap::draw(
    p_crispr,
    merge_legends   = TRUE,
    column_title    = "CRISPR WGCNA Co-dependency Modules",
    column_title_gp = grid::gpar(fontsize = 24, fontface = "bold")
  )
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
      path.plots,"HEATMAP_WGCNA_RNAi_OrderedByCRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_RNAi)
  grDevices::dev.off()
  
}

#### GO:BP Enrichment per CRISPR module (1 = run, 0 = read)
if (1) {
  
  ## Get all unique modules (excluding grey)
  moduleColors_CRISPR <- net_CRISPR$colors
  unique_modules <- unique(moduleColors_CRISPR[moduleColors_CRISPR != "grey"])
  
  ## Create a list to store enrichment results for each module
  enrich_results <- list()
  
  ## Loop through each module and perform ORA
  for (module in unique_modules) {
    
    module_genes <- names(moduleColors_CRISPR)[moduleColors_CRISPR == module]
    
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
    
    enrich_results[[module]] <- list(
      GO      = ego,
      n_genes = length(module_genes)
    )
    
    cat("Module:", module, "- Genes:", length(module_genes),
        "- GO terms:", nrow(ego@result), "\n")
  }
  
  ## Save RDS object with original names (needed for downstream lookups)
  saveRDS(enrich_results,
          file = paste0(path.wd, "DataSets/WGCNA/Enrichment_Results_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  ## Rename list using display names for Excel output
  enrich_results_renamed <- enrich_results
  names(enrich_results_renamed) <- name_remap[names(enrich_results_renamed)]
  
  ## Write file sorted by module size
  wb <- createWorkbook()
  
  module_sizes    <- sapply(names(enrich_results_renamed), function(m) enrich_results_renamed[[m]]$n_genes)
  modules_ordered <- names(sort(module_sizes, decreasing = TRUE))
  
  for (module in modules_ordered) {
    if (!is.null(enrich_results_renamed[[module]]$GO) &&
        nrow(enrich_results_renamed[[module]]$GO@result) > 0) {
      go_df <- enrich_results_renamed[[module]]$GO@result
      addWorksheet(wb, sheetName = module)
      writeData(wb, sheet = module, x = go_df)
    }
  }
  
  saveWorkbook(wb,
               file = paste0(path.wd, "DataSets/WGCNA/GO_Enrichment_CRISPR_AllModules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".xlsx"),
               overwrite = TRUE)
  
} else {
  
  enrich_results <- readRDS(paste0(path.wd, "DataSets/WGCNA/Enrichment_Results_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
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
      paste0(path.plots, "WGCGO_Dotplot_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", target_module, "_Module.png"),
      p_go_dot,
      width = 10,
      height = 8)
    
    ## Barplot for GO terms
    p_go_bar <- barplot(enrich_results[[target_module]]$GO,
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(
      paste0(path.plots, "WGCNA_GO_Barplot_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", target_module, "_Module.png"),
      p_go_bar,
      width = 10,
      height = 8)
    
    ## Enrichment map to show GO term relationships (with error handling)
    tryCatch({
      p_emap <- emapplot(pairwise_termsim(enrich_results[[target_module]]$GO),
                         showCategory = 30)
      
      ggsave(
        paste0(path.plots, "WGCNA_EnrichmentMap_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", "_Module.png"),
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
    referenceNetworks = 1,
    nPermutations     = 200,
    randomSeed        = 999,
    quickCor          = 0,
    verbose           = 3,
    maxGoldModuleSize = 1000,
    maxModuleSize     = 1000
  )
  
  ## Stop the cluster when done
  stopCluster(cl)
  
  ## Save RDS with original names (needed for visualization block)
  saveRDS(mp,
          file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  ## Extract preservation statistics
  stats <- mp$preservation$Z$ref.CRISPR$inColumnsAlsoPresentIn.RNAi
  
  ## Interpretation thresholds (Langfelder & Horvath)
  # Zsummary < 2: no preservation
  # 2 < Zsummary < 10: weak to moderate preservation  
  # Zsummary > 10: strong preservation
  # Note: gold module (all genes) and grey (unassigned) are not informative
  
  ## Remap row names to display names before writing
  stats_renamed <- stats
  rownames(stats_renamed) <- ifelse(
    rownames(stats_renamed) %in% names(name_remap),
    name_remap[rownames(stats_renamed)],
    rownames(stats_renamed)
  )
  
  write.table(
    x         = stats_renamed,
    file      = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".txt"),
    row.names = TRUE,
    sep       = "\t",
    quote     = FALSE
  )
  
}

#### Visualize preservation statistics
if (1) {
  
  mp <- readRDS(file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", 
                              soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  stats    <- mp$preservation$Z$ref.CRISPR$inColumnsAlsoPresentIn.RNAi
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
  
  ## Map original WGCNA module names -> display names -> hex colors
  ## Note: color_remap is keyed off NEW names, so go via name_remap first
  plotData$display_name <- name_remap[plotData$module]
  plotData$plot_color   <- color_remap[plotData$display_name]
  
  ## Plot 1: Zsummary vs module size
  p_preservation <- ggplot(plotData, aes(x = size, y = Zsummary, color = plot_color, label = display_name)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x     = "Module Size (number of genes)",
      y     = "Preservation Z-summary",
      title = paste0("Module Preservation: CRISPR modules in RNAi data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    annotate("text", x = max(plotData$size) * 0.7, y = 2,
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$size) * 0.7, y = 10,
             label = "Z = 10 (strong)",   vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_Zsummary_CRISPR_in_RNAi_SoftPower_", soft_power,
                "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight,
                "_deepSplit_", deep_Split, "_.png"),
         p_preservation, width = 8, height = 7)
  
  ## Plot 2: Median rank vs Zsummary
  p_rank <- ggplot(plotData, aes(x = medianRank, y = Zsummary, color = plot_color, label = display_name)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x     = "Median Rank",
      y     = "Preservation Z-summary",
      title = paste0("Module Preservation: CRISPR modules in RNAi data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    annotate("text", x = max(plotData$medianRank) * 0.2, y = 2,
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$medianRank) * 0.2, y = 10,
             label = "Z = 10 (strong)",   vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_MedianRank_CRISPR_in_RNAi_SoftPower_", soft_power,
                "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight,
                "_deepSplit_", deep_Split, "_.png"),
         p_rank, width = 8, height = 7)
  
  ## Plot 3: Z-score bar plot
  p_bar <- ggplot(plotData, aes(x = reorder(display_name, -Zsummary), y = Zsummary, fill = plot_color)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept =  0, linetype = "solid",  color = "black",     linewidth = 0.5) +
    geom_hline(yintercept =  2, linetype = "dashed", color = "blue",      linewidth = 0.6) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen", linewidth = 0.6) +
    annotate("text", x = nrow(plotData) * 0.7, y = 2,
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue",      size = 3) +
    annotate("text", x = nrow(plotData) * 0.7, y = 10,
             label = "Z = 10 (strong)",   vjust = -0.5, color = "darkgreen", size = 3) +
    scale_fill_identity() +
    labs(
      x     = "Module",
      y     = "Preservation Z-summary",
      title = "CRISPR Module Preservation in RNAi"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, color = "black"),
      legend.position = "none"
    )
  
  ggsave(paste0(path.plots, "ModulePreservation_Zbar_CRISPR_in_RNAi_SoftPower_", soft_power,
                "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight,
                "_deepSplit_", deep_Split, "_.png"),
         p_bar, width = 8, height = 3.5)
  
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
soft_power      <- 3L
deep_Split      <- 4
min_module_sz   <- 30L
merge_CutHeight <- 0.25

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

#### Run to investigate soft power option (CRISPR = 4, RNAi = 3)
if (0) {
  
  powers <- c(1:20)
  
  sft_RNAi <- pickSoftThreshold(RNAi_common,
                                powerVector = powers,
                                verbose = 5)
  
  pdf(paste0(path.plots, "soft_power_selection_RNAi.pdf"), width = 10, height = 5)
  
  par(mfrow = c(1,2))
  
  plot(sft_RNAi$fitIndices[,1], -sign(sft_RNAi$fitIndices[,3])*sft_RNAi$fitIndices[,2],
       xlab="Soft Power", ylab="Scale Free Topology R²",
       main="Scale independence")
  abline(h=0.8, col="red")
  
  plot(sft_RNAi$fitIndices[,1], sft_RNAi$fitIndices[,5],
       xlab="Soft Power", ylab="Mean Connectivity",
       main="Mean connectivity")
  
  dev.off()
}

#### Run WGCNA (1) or Read in WGCNA object (0)
if (0) {
  
  WGCNA::enableWGCNAThreads()
  
  net_RNAi <- WGCNA::blockwiseModules(
    RNAi_common,
    power              = soft_power,
    minModuleSize      = min_module_sz,
    networkType        = "signed",
    TOMType            = "signed",
    reassignThreshold  = 0,
    mergeCutHeight     = merge_CutHeight,
    numericLabels      = FALSE,
    pamRespectsDendro  = TRUE,
    verbose            = 3,
    deepSplit          = deep_Split
  )
  
  saveRDS(
    net_RNAi,
    file = paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  table(net_RNAi$colors)
  
} else {
  
  net_RNAi <- readRDS(paste0(path.wd, "/DataSets/WGCNA/WGCNA_Object_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  table(net_RNAi$colors)
}

#### Perform correlation on all non grey modules and plot
if (1) {
  
  ## Extract module colors
  moduleColors_RNAi <- net_RNAi$colors
  table(moduleColors_RNAi)
  
  all_modules <- unique(moduleColors_RNAi[moduleColors_RNAi != "grey"])
  
  name_remap <- c(
    "red" = "green",
    "yellow" = "red",
    "green" = "black",
    
    "turquoise" = "orange",
    "brown" = "violet",
    "blue" = "sapphire"
  )
  
  color_remap <- c(
    "green" = "#228B22",
    "red" = "#bb0a1e",
    "black" = "black",
    
    "orange" = "#FF6E00",
    "violet" = "#7F00FF",
    "sapphire" = "#1B4FA8"
  )
  
  ## Fill unmapped modules with original name
  for (m in all_modules) {
    if (!m %in% names(name_remap)) name_remap[m] <- m
  }
  
  ## Get all new names and fill unmapped colors with original R color
  all_new_names <- name_remap[all_modules]
  for (n in all_new_names) {
    if (!n %in% names(color_remap)) color_remap[n] <- n
  }
  
  ## Get names of genes assigned to clusters (remove "grey" genes)
  non_grey_genes <- colnames(RNAi_common)[moduleColors_RNAi != "grey"]
  length(non_grey_genes)
  
  RNAi_ng <- RNAi_common[, non_grey_genes, drop = FALSE]
  
  ## Perform correlation without grey genes
  cor_RNAi <- stats::cor(
    RNAi_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  ## Apply remaps
  mod_ng_orig  <- moduleColors_RNAi[non_grey_genes]
  mod_ng_name  <- name_remap[mod_ng_orig]
  mod_ng_color <- color_remap[mod_ng_name]
  names(mod_ng_name)  <- non_grey_genes
  names(mod_ng_color) <- non_grey_genes
  
  ## Reorder by display name
  gene_order   <- order(mod_ng_name)
  cor_RNAi_ord <- cor_RNAi[gene_order, gene_order]
  
  mod_ng_name_ord  <- mod_ng_name[gene_order]
  mod_ng_color_ord <- mod_ng_color[gene_order]
  
  ## Build color mapping: display name -> hex color
  module_col <- setNames(mod_ng_color_ord, mod_ng_name_ord)
  module_col <- module_col[!duplicated(names(module_col))]
  
  col_fun <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  )
  
  p_RNAi <- ComplexHeatmap::Heatmap(
    cor_RNAi_ord,
    name = "\nPearson r\n",
    col  = col_fun,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_split    = mod_ng_name_ord,
    column_split = mod_ng_name_ord,
    row_title    = NULL,
    column_title = NULL,
    row_gap    = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    heatmap_legend_param = list(
      title_gp      = grid::gpar(fontsize = 20, fontface = "bold"),
      labels_gp     = grid::gpar(fontsize = 20),
      title_gap     = grid::unit(20, "mm"),
      legend_height = grid::unit(40, "mm")
    ),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_name_ord,
      col    = list(Module = module_col),
      show_legend          = TRUE,
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        Module = list(
          title_gp    = grid::gpar(fontsize = 20, fontface = "bold"),
          labels_gp   = grid::gpar(fontsize = 16),
          grid_height = grid::unit(10, "mm"),
          grid_width  = grid::unit(8, "mm")
        )
      )
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_name_ord,
      col    = list(Module = module_col),
      show_legend          = FALSE,
      show_annotation_name = FALSE
    ),
    use_raster    = FALSE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(
      path.plots,
      "HEATMAP_WGCNA_RNAi_SoftPower_", soft_power,
      "_MinModuleSize_", min_module_sz,
      "_mergeCutHeight_", merge_CutHeight,
      "_deepSplit_", deep_Split,
      "_AllClusters.png"
    ),
    width  = 2500,
    height = 2500,
    res    = 200
  )
  ComplexHeatmap::draw(
    p_RNAi,
    merge_legends   = TRUE,
    column_title    = "RNAi WGCNA Co-dependency Modules",
    column_title_gp = grid::gpar(fontsize = 24, fontface = "bold")
  )
  grDevices::dev.off()
  
}

#### Look at all RNAi modules now in CRISPR
if (0) {
  
  CRISPR_ng <- CRISPR_common[, non_grey_genes, drop = FALSE]
  
  cor_CRISPR <- stats::cor(
    CRISPR_ng,
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  cor_CRISPR_ord <- cor_CRISPR[gene_order, gene_order, drop = FALSE]
  
  p_CRISPR <- ComplexHeatmap::Heatmap(
    cor_CRISPR_ord,
    name = "Pearson r",
    col  = col_fun,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_split    = mod_ng_name_ord,
    column_split = mod_ng_name_ord,
    row_title    = NULL,
    column_title = NULL,
    row_gap    = grid::unit(0, "pt"),
    column_gap = grid::unit(0, "pt"),
    rect_gp = grid::gpar(col = NA),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Module = mod_ng_name_ord,
      col    = list(Module = module_col),
      show_legend = TRUE
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Module = mod_ng_name_ord,
      col    = list(Module = module_col),
      show_legend = FALSE
    ),
    use_raster    = TRUE,
    raster_quality = 6
  )
  
  Cairo::CairoPNG(
    filename = paste0(
      path.plots,
      "HEATMAP_WGCNA_CRISPR_OrderedByRNAi_SoftPower_", soft_power,
      "_MinModuleSize_", min_module_sz,
      "_mergeCutHeight_", merge_CutHeight,
      "_deepSplit_", deep_Split,
      "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_CRISPR)
  grDevices::dev.off()
  
}

#### GO:BP Enrichment per RNAi module (1 = run, 0 = read)
if (1) {
  
  moduleColors_RNAi <- net_RNAi$colors
  unique_modules    <- unique(moduleColors_RNAi[moduleColors_RNAi != "grey"])
  
  enrich_results <- list()
  
  for (module in unique_modules) {
    
    module_genes <- names(moduleColors_RNAi)[moduleColors_RNAi == module]
    
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
    
    enrich_results[[module]] <- list(
      GO      = ego,
      n_genes = length(module_genes)
    )
    
    cat("Module:", module, "- Genes:", length(module_genes),
        "- GO terms:", nrow(ego@result), "\n")
  }
  
  ## Save RDS with original names (needed for downstream lookups)
  saveRDS(enrich_results,
          file = paste0(path.wd, "DataSets/WGCNA/Enrichment_Results_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  ## Rename list using display names for Excel output
  enrich_results_renamed <- enrich_results
  names(enrich_results_renamed) <- name_remap[names(enrich_results_renamed)]
  
  ## Write file sorted by module size
  wb <- createWorkbook()
  
  module_sizes    <- sapply(names(enrich_results_renamed), function(m) enrich_results_renamed[[m]]$n_genes)
  modules_ordered <- names(sort(module_sizes, decreasing = TRUE))
  
  for (module in modules_ordered) {
    if (!is.null(enrich_results_renamed[[module]]$GO) &&
        nrow(enrich_results_renamed[[module]]$GO@result) > 0) {
      go_df <- enrich_results_renamed[[module]]$GO@result
      addWorksheet(wb, sheetName = module)
      writeData(wb, sheet = module, x = go_df)
    }
  }
  
  saveWorkbook(wb,
               file = paste0(path.wd, "DataSets/WGCNA/GO_Enrichment_RNAi_AllModules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".xlsx"),
               overwrite = TRUE)
  
} else {
  
  enrich_results <- readRDS(paste0(path.wd, "DataSets/WGCNA/Enrichment_Results_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
}

#### Visualize enrichment for top modules with significant results
n_modules_to_plot <- 5

if (1) {
  
  all_modules <- names(sort(table(moduleColors_RNAi[moduleColors_RNAi != "grey"]),
                            decreasing = TRUE))
  
  modules_with_sig_results <- c()
  for (module in all_modules) {
    if (!is.null(enrich_results[[module]]$GO) &&
        nrow(enrich_results[[module]]$GO@result) > 0 &&
        sum(enrich_results[[module]]$GO@result$p.adjust < 0.05) > 0) {
      modules_with_sig_results <- c(modules_with_sig_results, module)
    }
  }
  
  top_modules <- head(modules_with_sig_results, n_modules_to_plot)
  
  cat("Plotting", length(top_modules), "modules with significant GO enrichment:\n")
  cat(paste(top_modules, collapse = ", "), "\n\n")
  
  for (target_module in top_modules) {
    
    p_go_dot <- dotplot(enrich_results[[target_module]]$GO,
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(
      paste0(path.plots, "WGCGO_Dotplot_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", target_module, "_Module.png"),
      p_go_dot, width = 10, height = 8)
    
    p_go_bar <- barplot(enrich_results[[target_module]]$GO,
                        showCategory = 15,
                        title = paste0(target_module, " module - GO:BP enrichment"))
    
    ggsave(
      paste0(path.plots, "WGCNA_GO_Barplot_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", target_module, "_Module.png"),
      p_go_bar, width = 10, height = 8)
    
    tryCatch({
      p_emap <- emapplot(pairwise_termsim(enrich_results[[target_module]]$GO),
                         showCategory = 30)
      
      ggsave(
        paste0(path.plots, "WGCNA_EnrichmentMap_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", target_module, "_Module.png"),
        p_emap, width = 12, height = 10)
      
    }, error = function(e) {
      cat("Could not create enrichment map for module:", target_module,
          "(not enough similar terms)\n")
    })
    
    cat("Plots saved for module:", target_module, "\n")
  }
  
  if (length(top_modules) < n_modules_to_plot) {
    cat("\nNote: Only", length(top_modules), "modules had significant GO enrichment (requested", n_modules_to_plot, ")\n")
  }
  
}

#### Checking for conservation between RNAi and CRISPR
if (0) {
  
  moduleColors_RNAi <- net_RNAi$colors
  
  multiExpr <- list(
    RNAi   = list(data = RNAi_common),
    CRISPR = list(data = CRISPR_common)
  )
  
  multiColor <- list(
    RNAi = moduleColors_RNAi
  )
  
  n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  WGCNA::enableWGCNAThreads(nThreads = n_cores)
  
  set.seed(999)
  
  mp <- WGCNA::modulePreservation(
    multiExpr,
    multiColor,
    referenceNetworks = 1,
    nPermutations     = 200,
    randomSeed        = 999,
    quickCor          = 0,
    verbose           = 3,
    maxGoldModuleSize = 1000,
    maxModuleSize     = 1000
  )
  
  stopCluster(cl)
  
  ## Save RDS with original names (needed for visualization block)
  saveRDS(mp,
          file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_RNAi_in_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  ## Extract preservation statistics
  stats <- mp$preservation$Z$ref.RNAi$inColumnsAlsoPresentIn.CRISPR
  
  ## Remap row names to display names before writing
  stats_renamed <- stats
  rownames(stats_renamed) <- ifelse(
    rownames(stats_renamed) %in% names(name_remap),
    name_remap[rownames(stats_renamed)],
    rownames(stats_renamed)
  )
  
  write.table(
    x         = stats_renamed,
    file      = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_RNAi_in_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".txt"),
    row.names = TRUE,
    sep       = "\t",
    quote     = FALSE
  )
  
}

#### Visualize preservation statistics
if (1) {
  
  mp <- readRDS(file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_RNAi_in_CRISPR_SoftPower_",
                              soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  
  stats    <- mp$preservation$Z$ref.RNAi$inColumnsAlsoPresentIn.CRISPR
  obsStats <- mp$preservation$observed$ref.RNAi$inColumnsAlsoPresentIn.CRISPR
  
  plotData <- data.frame(
    module     = rownames(stats),
    size       = stats$moduleSize,
    Zsummary   = stats$Zsummary.pres,
    medianRank = obsStats$medianRank.pres
  )
  
  plotData <- plotData[!plotData$module %in% c("gold", "grey"), ]
  
  plotData$preservation <- cut(
    plotData$Zsummary,
    breaks = c(-Inf, 2, 10, Inf),
    labels = c("No preservation", "Weak-Moderate", "Strong preservation")
  )
  
  ## Map original WGCNA module names -> display names -> hex colors
  plotData$display_name <- name_remap[plotData$module]
  plotData$plot_color   <- color_remap[plotData$display_name]
  
  ## Plot 1: Zsummary vs module size
  p_preservation <- ggplot(plotData, aes(x = size, y = Zsummary, color = plot_color, label = display_name)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x     = "Module Size (number of genes)",
      y     = "Preservation Z-summary",
      title = paste0("Module Preservation: RNAi modules in CRISPR data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    annotate("text", x = max(plotData$size) * 0.7, y = 2,
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$size) * 0.7, y = 10,
             label = "Z = 10 (strong)",   vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_Zsummary_RNAi_in_CRISPR_SoftPower_", soft_power,
                "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight,
                "_deepSplit_", deep_Split, "_.png"),
         p_preservation, width = 8, height = 7)
  
  ## Plot 2: Median rank vs Zsummary
  p_rank <- ggplot(plotData, aes(x = medianRank, y = Zsummary, color = plot_color, label = display_name)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 2,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
    geom_text(hjust = -0.2, vjust = -0.2, size = 3, show.legend = FALSE) +
    scale_color_identity() +
    labs(
      x     = "Median Rank",
      y     = "Preservation Z-summary",
      title = paste0("Module Preservation: RNAi modules in CRISPR data: Soft Power ", soft_power, ", Min Mod Size ", min_module_sz)
    ) +
    annotate("text", x = max(plotData$medianRank) * 0.2, y = 2,
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue") +
    annotate("text", x = max(plotData$medianRank) * 0.2, y = 10,
             label = "Z = 10 (strong)",   vjust = -0.5, color = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(path.plots, "ModulePreservation_MedianRank_RNAi_in_CRISPR_SoftPower_", soft_power,
                "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight,
                "_deepSplit_", deep_Split, "_.png"),
         p_rank, width = 8, height = 7)
  
  ## Plot 3: Z-score bar plot
  p_bar <- ggplot(plotData, aes(x = reorder(display_name, -Zsummary), y = Zsummary, fill = plot_color)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept =  0, linetype = "solid",  color = "black",     linewidth = 0.5) +
    geom_hline(yintercept =  2, linetype = "dashed", color = "blue",      linewidth = 0.6) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen", linewidth = 0.6) +
    annotate("text", x = nrow(plotData) * 0.7, y = 2,
             label = "Z = 2 (threshold)", vjust = -0.5, color = "blue",      size = 3) +
    annotate("text", x = nrow(plotData) * 0.7, y = 10,
             label = "Z = 10 (strong)",   vjust = -0.5, color = "darkgreen", size = 3) +
    scale_fill_identity() +
    labs(
      x     = "Module",
      y     = "Preservation Z-summary",
      title = "RNAi Module Preservation in CRISPR"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, color = "black"),
      legend.position = "none"
    )
  
  ggsave(paste0(path.plots, "ModulePreservation_Zbar_RNAi_in_CRISPR_SoftPower_", soft_power,
                "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight,
                "_deepSplit_", deep_Split, "_.png"),
         p_bar, width = 8, height = 3.5)
  
}

##### Manual GO Plot V1 #####
plot_go_enrichment <- function(pathways, padj, modules, colors, FoldEnrichment,
                               gene_ratio = NULL,
                               title = "RNAi GO:BP Results (FDR < 0.05)") {
  
  df <- data.frame(
    pathway        = pathways,
    padj           = padj,
    module         = modules,
    color          = colors,
    log_padj       = -log10(padj),
    FoldEnrichment = FoldEnrichment,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(gene_ratio)) {
    df$gene_ratio <- sapply(gene_ratio, function(x) {
      if (is.character(x)) eval(parse(text = x)) else as.numeric(x)
    })
  } else {
    df$gene_ratio <- 0.2
  }
  
  # rank modules by max FoldEnrichment, then sort within module
  module_rank <- tapply(df$FoldEnrichment, df$module, max)
  df$module_rank <- module_rank[df$module]
  df <- df[order(df$module_rank, df$FoldEnrichment), ]
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))
  
  # set module factor order so facets stack by module_rank
  module_order <- unique(df$module[order(-df$module_rank)])
  df$module <- factor(df$module, levels = rev(module_order))
  
  # named color scale
  color_map <- setNames(df$color, df$module)
  color_map <- color_map[!duplicated(names(color_map))]
  
  p <- ggplot(df, aes(x = FoldEnrichment, y = pathway, color = module, size = gene_ratio)) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, yend = pathway),
                 linewidth = 0.8, alpha = 0.4) +
    geom_point() +
    facet_grid(module ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    scale_color_manual(values = color_map) +
    scale_size_continuous(range = c(3, 8), name = "gene ratio") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x     = "Fold Enrichment",
      y     = NULL,
      title = title,
      color = "module"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major.x  = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor    = element_blank(),
      axis.text.y         = element_text(size = 10),
      legend.position     = "right",
      plot.title          = element_text(size = 13, face = "plain"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.text          = element_blank(),
      strip.background    = element_blank(),
      panel.spacing       = unit(0.3, "lines"),
      strip.placement     = "outside",
      axis.line = element_blank(),
    )
  
  p
}

pdf(file = paste0(path.plots, "GO_AllModules_CRISPR_V1.pdf"), height = 6, width = 10)
plot_go_enrichment(
  pathways = c(
    "adenylate cyclase-activating\nGPCR signaling",
    "adenylate cyclase-modulating\nGPCR signaling",
    "mitochondrial gene expression",
    "respiratory electron transport chain",
    "RNA splicing",
    "protein-RNA complex organization",
    "cell-substrate adhesion",
    "actin filament organization",
    "polysaccharide catabolic process",
    "glucan catabolic process"
  ),
  padj = c(0.01971604, 0.01971604, 5.08303050835317E-147, 1.82718E-45, 2.31869E-30, 5.1e-7, 9.04064E-11, 1.26283E-09, 0.033358824, 0.033358824),
  FoldEnrichment = c(3.579880213, 2.969655892, 28.08579882, 16.36724138, 9.652784388, 9.345594487, 12.37477595, 9.607084124, 36.37164751, 39.06584362),
  modules = c("yellow","yellow","green","green","red","red","black","black",
              "pink","pink"),
  colors  = c("#FFDA03","#FFDA03","#228B22","#228B22","#bb0a1e","#bb0a1e","black","black",
              "pink","pink"),
  gene_ratio = c(0.035, 0.043, 0.25, 0.108, 0.246, 0.115, 0.235, 0.235, 0.056, 0.056)
)
dev.off()

pdf(file = paste0(path.plots, "GO_AllModules_RNAi_V1.pdf"), height = 6, width = 10)
plot_go_enrichment(
  pathways = c(
    "rRNA metabolic process",
    "ribosome biogenesis",
    "substrate adhesion-dependent\ncell spreading",
    "positive regulation of\nlamellipodium assembly",
    "respiratory electron transport chain",
    "mitochondrial gene expression",
    "cytoplasmic translation"
  ),
  padj = c(1.94651E-42, 4.01514E-52, 0.020327341, 0.010704066, 0.003020666, 1.64197E-50, 1.64475E-06),
  FoldEnrichment = c(16.6690079, 17.02832931, 7.786162048, 20.06128487, 17.35862857, 89.87455621, 4.439300412),
  modules = c("yellow","yellow", "green", "green", "red", "red", "brown"),
  colors  = c("#FFDA03","#FFDA03", "green", "green", "red", "red", "brown"),
  gene_ratio = c(0.235, 0.284, 0.0414, 0.029585, 0.114, 0.8, 0.0386)
)
dev.off()


##### Manual GO Plot V2 (more colors) #####
library(ggh4x)

plot_go_enrichment <- function(pathways, padj, modules, colors, FoldEnrichment,
                               gene_ratio = NULL,
                               title = "CRISPR GO:BP Results (FDR > 0.05)") {
  
  df <- data.frame(
    pathway        = pathways,
    padj           = padj,
    module         = modules,
    color          = colors,
    log_padj       = -log10(padj),
    FoldEnrichment = FoldEnrichment,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(gene_ratio)) {
    df$gene_ratio <- sapply(gene_ratio, function(x) {
      if (is.character(x)) eval(parse(text = x)) else as.numeric(x)
    })
  } else {
    df$gene_ratio <- 0.2
  }
  
  # rank modules by max FoldEnrichment, then sort within module
  module_rank <- tapply(df$FoldEnrichment, df$module, max)
  df$module_rank <- module_rank[df$module]
  df <- df[order(df$module_rank, df$FoldEnrichment), ]
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))
  
  # set module factor order so facets stack by module_rank
  module_order <- unique(df$module[order(-df$module_rank)])
  df$module <- factor(df$module, levels = rev(module_order))
  
  # named color scale
  color_map <- setNames(df$color, df$module)
  color_map <- color_map[!duplicated(names(color_map))]
  
  # strip colors in facet level order
  strip_colors <- color_map[levels(df$module)]
  
  p <- ggplot(df, aes(x = FoldEnrichment, y = pathway, color = module, size = gene_ratio)) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, yend = pathway),
                 linewidth = 0.8, alpha = 0.4) +
    geom_point() +
    facet_grid2(module ~ ., scales = "free_y", space = "free_y",
                switch = "y",
                strip = strip_themed(
                  background_y = elem_list_rect(fill = strip_colors),
                  text_y = elem_list_text(size = 1, color = "transparent")
                )) +
    scale_color_manual(values = color_map) +
    scale_size_continuous(range = c(3, 8), name = "gene ratio") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x     = "Fold Enrichment",
      y     = NULL,
      title = title,
      color = "module"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major.x   = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor     = element_blank(),
      axis.text.y          = element_text(size = 10),
      axis.line            = element_blank(),
      legend.position      = "right",
      plot.title           = element_text(size = 13, face = "plain"),
      panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.background     = element_blank(),
      panel.spacing        = unit(0.3, "lines"),
      strip.placement      = "outside",
      strip.switch.pad.grid = unit(0.05, "cm")
    )
  
  p
}

pdf(file = paste0(path.plots, "GO_AllModules_CRISPR_V2.pdf"), height = 6.25, width = 8)
plot_go_enrichment(
  pathways = c(
    "adenylate cyclase-activating\nGPCR signaling",
    "adenylate cyclase-modulating\nGPCR signaling",
    "mitochondrial gene expression",
    "respiratory electron transport chain",
    "RNA splicing",
    "  protein-RNA complex organization",
    "cell-substrate adhesion",
    "actin filament organization",
    "polysaccharide catabolic process",
    "glucan catabolic process"
  ),
  padj = c(0.01971604, 0.01971604, 5.08303050835317E-147, 1.82718E-45, 2.31869E-30, 5.1e-7, 9.04064E-11, 1.26283E-09, 0.033358824, 0.033358824),
  FoldEnrichment = c(3.579880213, 2.969655892, 28.08579882, 16.36724138, 9.652784388, 9.345594487, 12.37477595, 9.607084124, 36.37164751, 39.06584362),
  modules = c("yellow","yellow","green","green","red","red","black","black",
              "pink","pink"),
  colors  = c("#FFDA03","#FFDA03","green","green","red","red","black","black",
              "pink","pink"),
  gene_ratio = c(0.035, 0.043, 0.25, 0.108, 0.246, 0.115, 0.235, 0.235, 0.056, 0.056)
)
dev.off()


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
soft_power_crispr <- 4L
soft_power_rnai   <- 3L
deep_Split        <- 4
min_module_sz     <- 30L
merge_CutHeight   <- 0.25

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

### Plot the PCA
source("/Users/jack/Documents/GitHub/FDB_Freeland/Scripts/PCA_plot.R")

if (1) {
  
  CRISPR_path <- paste0(path.pca, "CRISPR_common_PCA.txt")
  RNAi_path   <- paste0(path.pca, "RNAi_common_PCA.txt")
  
  CRISPR_PCA_input <- read.table(file = CRISPR_path, sep = "\t", header = T)
  
  model <- read.csv(paste0(path.dm, "Model.csv"))
  model$ModelID <- gsub("-", "\\.", model$ModelID)
  
  info_name <- colnames(CRISPR_PCA_input)
  info_type <- model$OncotreeLineage[match(info_name, model$ModelID)]
  
  pdf(file = paste0(path.plots, "PCA_CRISPR_common.pdf"), height = 6, width = 10)
  PCA_plot(
    file = paste0(path.pca, "CRISPR_common_PCA_prcomp_scores.txt"),
    info.name = info_name,
    info.type = as.factor(info_type),
    title = "CRISPR PCA: Shared Genes",
    ellipse = F,
    labels = F,
    colors = NULL,
    info.shape = NULL,
    shape = NULL,
    PCx="PC1", PCy="PC2", conf = 0.95, density = F, fliph = F, flipv = F, show.legend = TRUE
  )
  dev.off()
  
  pdf(file = paste0(path.plots, "PCA_RNAi_common.pdf"), height = 6, width = 10)
  PCA_plot(
    file = paste0(path.pca, "RNAi_common_PCA_prcomp_scores.txt"),
    info.name = info_name,
    info.type = as.factor(info_type),
    title = "RNAi PCA: Shared Genes",
    ellipse = F,
    labels = F,
    colors = NULL,
    info.shape = NULL,
    shape = NULL,
    PCx="PC1", PCy="PC2", conf = 0.95, density = F, fliph = F, flipv = F, show.legend = TRUE
  )
  dev.off()
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
  WGCNA_RNAi <- readRDS(paste0(path.wd, "DataSets/WGCNA/WGCNA_Object_RNAi_SoftPower_", soft_power_rnai, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))
  WGCNA_CRISPR <- readRDS(paste0(path.wd, "DataSets/WGCNA/WGCNA_Object_CRISPR_SoftPower_", soft_power_crispr, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, ".rds"))

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
    x     = Max_PCA[[xcol]],
    y     = Max_PCA[[ycol]],
    gene  = Max_PCA$Loading,
    theta = Max_PCA$theta_deg
  ) %>%
    dplyr::mutate(
      angle_group = dplyr::case_when(
        theta < 30  ~ "RNAi",
        theta <= 60 ~ "Neutral",
        TRUE        ~ "CRISPR"
      ),
      angle_group = factor(angle_group, levels = c("CRISPR", "Neutral", "RNAi"))
    )
  
  top5 <- plot_df %>%
    dplyr::mutate(r = sqrt(x^2 + y^2)) %>%
    dplyr::group_by(angle_group) %>%
    dplyr::slice_max(r, n = 5, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  p_scatter <- ggplot(plot_df, aes(x = x, y = y, color = angle_group)) +
    geom_point(size = 0.075, alpha = 0.3) +
    geom_abline(
      slope = tan(pi/6), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    geom_abline(
      slope = tan(pi/3), intercept = 0,
      linetype = "dotted", linewidth = 0.4, color = "grey40"
    ) +
    geom_text_repel(
      data          = top5,
      aes(label     = gene),
      size          = 2.5,
      fontface      = "italic",
      show.legend   = FALSE,
      max.overlaps  = 20,
      segment.size  = 0.3,
      segment.color = "grey50"
    ) +
    scale_color_manual(
      values = c("RNAi" = "#5E2F80", "Neutral" = "#BDBDBD", "CRISPR" = "#D47D37"),
      name   = "Bias"
    ) +
    labs(
      x = "Max Absolute PCA Loading (RNAi)",
      y = "Max Absolute PCA Loading (CRISPR)",
      title = "Max Absolute PCA Loadings per Gene PC 1-10"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic(base_size = 10) +
    theme(
      legend.background = element_blank(),
      legend.key        = element_blank()
    )
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_CRISPR_vs_RNAi_Scatter.pdf"),
    plot = p_scatter, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
  ## Plot polar - CRISPR modules (all genes, colored by CRISPR module)
  crispr_color_remap <- c(
    "green"  = "#228B22",
    "red"    = "#bb0a1e",
    "yellow" = "#FFDA03",
    "brown"  = "#481F01"
  )
  
  Max_PCA_CRISPR_ordered <- Max_PCA %>%
    dplyr::arrange(desc(CRISPR_module_color == "grey80")) %>%
    dplyr::mutate(
      CRISPR_display_color = dplyr::recode(CRISPR_module_color, !!!crispr_color_remap, .default = CRISPR_module_color)
    )
  
  crispr_legend <- Max_PCA_CRISPR_ordered %>%
    dplyr::filter(CRISPR_module_color != "grey80") %>%
    dplyr::distinct(CRISPR_display_color, CRISPR_module_color)
  
  p_polar_CRISPR <- ggplot(Max_PCA_CRISPR_ordered, aes(x = r, y = theta_deg, color = CRISPR_display_color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(
      yintercept = c(30, 60),
      linetype   = "dotted",
      color      = "grey40",
      linewidth  = 0.4
    ) +
    scale_color_identity(
      guide  = "legend",
      name   = "CRISPR modules",
      breaks = c(crispr_legend$CRISPR_display_color, "grey80"),
      labels = c(crispr_legend$CRISPR_module_color,  "unassigned")
    ) +
    labs(
      x = "r (magnitude)",
      y = "θ (degrees)",
      title = "Gene-Only PCA Max Loadings: CRISPR WGCNA Modules"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, 90)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic(base_size = 8) +
    theme(
      legend.key   = element_blank(),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold")
      # plot.margin  = margin(r = 15)
    ) +
    annotate("text", x = Inf, y = 15, label = "RNAi Biased",   color = "#5E2F80", size = 2.5, hjust = 1.1, fontface = "bold") +
    annotate("text", x = Inf, y = 45, label = "Neutral",        color = "#BDBDBD", size = 2.5, hjust = 1.1, fontface = "bold") +
    annotate("text", x = Inf, y = 75, label = "CRISPR Biased",  color = "#D47D37", size = 2.5, hjust = 1.1, fontface = "bold")
  
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_Polar_CRISPR_Modules_SoftPower_", soft_power_crispr, "_MinModuleSize_", min_module_sz, ".pdf"),
    plot = p_polar_CRISPR, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
  # Apply remapping to RNAi module colors
  name_remap <- c(
    "red"       = "green",
    "yellow"    = "red",
    "green"     = "black",
    "turquoise" = "orange",
    "brown"     = "violet",
    "blue"      = "sapphire"
  )
  
  color_remap <- c(
    "green"   = "#228B22",
    "red"     = "#bb0a1e",
    "black"   = "black",
    "orange"  = "#FF6E00",
    "violet"  = "#7F00FF",
    "sapphire" = "#1B4FA8"
  )
  
  Max_PCA_RNAi_ordered <- Max_PCA %>%
    dplyr::arrange(desc(RNAi_module_color == "grey80")) %>%
    dplyr::mutate(
      RNAi_display_name  = dplyr::recode(RNAi_module_color, !!!name_remap, .default = RNAi_module_color),
      RNAi_display_color = dplyr::recode(RNAi_display_name, !!!color_remap, .default = RNAi_module_color)
    )
  
  # Build legend breaks/labels from the remapped non-grey modules
  rnai_legend <- Max_PCA_RNAi_ordered %>%
    dplyr::filter(RNAi_module_color != "grey80") %>%
    dplyr::distinct(RNAi_display_color, RNAi_display_name)
  
  p_polar_RNAi <- ggplot(Max_PCA_RNAi_ordered, aes(x = r, y = theta_deg, color = RNAi_display_color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(
      yintercept = c(30, 60),
      linetype   = "dotted",
      color      = "grey40",
      linewidth  = 0.4
    ) +
    scale_color_identity(
      guide  = "legend",
      name   = "RNAi modules",
      breaks = c(rnai_legend$RNAi_display_color, "grey80"),
      labels = c(rnai_legend$RNAi_display_name,  "unassigned")
    ) +
    labs(
      x = "r (magnitude)",
      y = "θ (degrees)",
      title = 
        # paste0("Max PCA Loadings: RNAi WGCNA Modules: Soft Power ", soft_power_rnai, ", Min Mod Size ", min_module_sz)
        "Gene-Only PCA Max Loadings: RNAi WGCNA Modules"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, 90)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic(base_size = 8) +
    theme(
      legend.key   = element_blank(),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold")
    ) +
    annotate("text", x = Inf, y = 15, label = "RNAi Biased",    color = "#5E2F80", size = 2.5, hjust = 1.1, fontface = "bold") +
    annotate("text", x = Inf, y = 45, label = "Neutral",  color = "#BDBDBD", size = 2.5, hjust = 1.1, fontface = "bold") +
    annotate("text", x = Inf, y = 75, label = "CRISPR Biased",   color = "#D47D37", size = 2.5, hjust = 1.1, fontface = "bold")
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_Polar_RNAi_Modules_SoftPower_", soft_power_rnai, "_MinModuleSize_", min_module_sz, ".pdf"),
    plot = p_polar_RNAi, width = 5, height = 4, units = "in", device = cairo_pdf
  )
  
  
  ### Add density plots
  
  ## RNAi 

  p_density_RNAi <- ggplot(
    Max_PCA_RNAi_ordered %>% dplyr::filter(RNAi_module_color != "grey80"),
    aes(x = theta_deg, fill = RNAi_display_color, color = RNAi_display_color)
  ) +
    geom_density(alpha = 0.4, linewidth = 0.4) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_flip() +
    scale_x_continuous(limits = c(0, 90), expand = expansion(mult = 0)) +
    scale_y_reverse(expand = expansion(mult = c(0, 0.05))) +  # reverse so it reads left-to-right naturally when flipped
    labs(x = NULL, y = "Density") +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  p_polar_RNAi_combined <- p_density_RNAi + p_polar_RNAi +
    plot_layout(widths = c(1, 4))
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_Polar_RNAi_Modules_SoftPower_", soft_power_rnai, "_MinModuleSize_", min_module_sz, "_density.pdf"),
    plot = p_polar_RNAi_combined, width = 6, height = 4, units = "in", device = cairo_pdf
  )
  
  ## CRISPR
  
  p_density_CRISPR <- ggplot(
    Max_PCA_CRISPR_ordered %>% dplyr::filter(CRISPR_module_color != "grey80"),
    aes(x = theta_deg, fill = CRISPR_module_color, color = CRISPR_module_color)
  ) +
    geom_density(alpha = 0.4, linewidth = 0.4) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_flip() +
    scale_x_continuous(limits = c(0, 90), expand = expansion(mult = 0)) +
    scale_y_reverse(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Density") +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.title.y    = element_blank(),
      legend.position = "none"
    )
  
  p_polar_CRISPR_combined <- p_density_CRISPR + p_polar_CRISPR +
    plot_layout(widths = c(1, 4))
  
  ggsave(
    filename = paste0(path.plots, "MaxPCALoadings_Polar_CRISPR_Modules_SoftPower_", soft_power_crispr, "_MinModuleSize_", min_module_sz, "_desnity.pdf"),
    plot = p_polar_CRISPR_combined, width = 6, height = 4, units = "in", device = cairo_pdf
  )
  
}

### KS Tests: Do modules preferentially fall on one side of theta? ###
if (1) {
  
  theta_cutoff_low  <- 30   # below = RNAi biased
  theta_cutoff_high <- 60   # above = CRISPR biased
  
  run_module_ks_tests <- function(data, module_col, display_name_col, theta_col = "theta_deg", dataset_label) {
    
    modules <- unique(data[[module_col]])
    modules <- modules[!is.na(modules) & modules != "grey80"]
    
    # Background = all non-grey genes
    background_theta <- data %>%
      dplyr::filter(.data[[module_col]] != "grey80") %>%
      dplyr::pull(.data[[theta_col]])
    
    results <- purrr::map_dfr(modules, function(mod) {
      
      module_theta <- data %>%
        dplyr::filter(.data[[module_col]] == mod) %>%
        dplyr::pull(.data[[theta_col]])
      
      display_name <- data %>%
        dplyr::filter(.data[[module_col]] == mod) %>%
        dplyr::pull(.data[[display_name_col]]) %>%
        unique() %>%
        .[1]
      
      n_mod <- length(module_theta)
      if (n_mod < 5) return(NULL)
      
      ks_res   <- ks.test(module_theta, background_theta)
      
      # Directional: is the module shifted toward RNAi (<30) or CRISPR (>60)?
      frac_rnai_biased   <- mean(module_theta < theta_cutoff_low)
      frac_crispr_biased <- mean(module_theta > theta_cutoff_high)
      frac_neutral       <- mean(module_theta >= theta_cutoff_low & module_theta <= theta_cutoff_high)
      median_theta       <- median(module_theta)
      
      direction <- dplyr::case_when(
        median_theta < theta_cutoff_low  ~ "RNAi biased",
        median_theta > theta_cutoff_high ~ "CRISPR biased",
        TRUE                             ~ "Neutral"
      )
      
      data.frame(
        dataset           = dataset_label,
        module_color      = mod,
        module_name       = display_name,
        n_genes           = n_mod,
        median_theta      = round(median_theta, 2),
        direction         = direction,
        frac_rnai_biased  = round(frac_rnai_biased,   3),
        frac_neutral      = round(frac_neutral,        3),
        frac_crispr_biased = round(frac_crispr_biased, 3),
        ks_statistic      = round(ks_res$statistic, 4),
        ks_pval           = ks_res$p.value,
        stringsAsFactors  = FALSE
      )
    })
    
    results %>%
      dplyr::mutate(ks_padj = p.adjust(ks_pval, method = "BH")) %>%
      dplyr::arrange(ks_padj)
  }
  
  ## Run for CRISPR modules
  ks_CRISPR <- run_module_ks_tests(
    data              = Max_PCA_CRISPR_ordered,
    module_col        = "CRISPR_module_color",
    display_name_col  = "CRISPR_module_color",   # CRISPR uses raw color as display name (post-remap via CRISPR_display_color separately)
    dataset_label     = "CRISPR"
  )
  
  ## Run for RNAi modules
  ks_RNAi <- run_module_ks_tests(
    data              = Max_PCA_RNAi_ordered,
    module_col        = "RNAi_module_color",
    display_name_col  = "RNAi_display_name",
    dataset_label     = "RNAi"
  )
  
  ks_combined <- dplyr::bind_rows(ks_CRISPR, ks_RNAi)
  
  print(ks_combined)
  
  write.table(
    x         = ks_combined,
    file      = paste0(path.max, "KS_Test_Module_Theta_Bias.txt"),
    quote     = F,
    sep       = "\t",
    col.names = T,
    row.names = F
  )
  
  ## Visualization: forest-style dot plot of median theta per module
  ks_sig <- ks_combined %>%
    dplyr::mutate(
      label      = paste0(dataset, ": ", module_name),
      sig_marker = dplyr::case_when(
        ks_padj < 0.001 ~ "***",
        ks_padj < 0.01  ~ "**",
        ks_padj < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  p_ks_dot <- ggplot(ks_sig, aes(x = median_theta, y = reorder(label, median_theta), color = direction)) +
    geom_vline(xintercept = c(30, 60), linetype = "dotted", color = "grey40", linewidth = 0.4) +
    geom_point(aes(size = n_genes), alpha = 0.85) +
    geom_text(
      aes(label = sig_marker),
      hjust  = -0.4,
      size   = 3,
      color  = "black"
    ) +
    scale_color_manual(
      values = c("RNAi biased" = "#5E2F80", "Neutral" = "#BDBDBD", "CRISPR biased" = "#D47D37"),
      name   = "Direction"
    ) +
    scale_size_continuous(name = "# genes", range = c(2, 7)) +
    scale_x_continuous(limits = c(0, 90), breaks = c(0, 30, 60, 90)) +
    labs(
      x     = "Median θ (degrees)",
      y     = NULL,
      title = "KS Test: Module Theta Bias vs Background",
      subtitle = "* p<0.05, ** p<0.01, *** p<0.001 (BH adjusted)"
    ) +
    theme_classic(base_size = 9) +
    theme(
      legend.key      = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    filename = paste0(path.plots, "KS_Test_Module_Theta_Bias_DotPlot.pdf"),
    plot     = p_ks_dot,
    width    = 6, height = max(3, nrow(ks_sig) * 0.4),
    units    = "in",
    device   = cairo_pdf
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
path.pls     <- paste0(path.wd, "DataSets/PLS/")

#### Get Hallmark EMT gene set from MSigDB
if (1) {

  ## Get Hallmark EMT gene set
  hallmark_emt <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ) %>%
    dplyr::filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  cat("Number of genes in Hallmark EMT gene set:", length(hallmark_emt))

  ## Save gene list
  write.table(
    x = data.frame(Gene = hallmark_emt),
    file = paste0(path.mel, "Hallmark_EMT_GeneSet.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
}

#### Score EMT signature
score_method <- "ssgsea" # zscore or ssgsea

if (0) {
  
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
  
  colnames(counts_trim) <- gsub("\\.", "-", colnames(counts_trim))
  
  ## Score gene set
  if (score_method == "ssgsea") {
    
    ## ssGSEA: single-sample enrichment score per cell line
    run.ssGSEA2 <- function(exp_mat, gene_list, norm = FALSE) {
      gsvaPar <- GSVA::ssgseaParam(exp_mat, gene_list, normalize = norm)
      as.data.frame(t(GSVA::gsva(gsvaPar)))
    }
    
    EMT_scores <- run.ssGSEA2(
      exp_mat   = as.matrix(counts_trim),
      gene_list = list(EMT = hallmark_emt)
    )

    cat("EMT scores computed via ssGSEA")
    
  } else {
    
    ## Z-score
    
    ## Genes present in both the expression matrix and the gene set
    emt_genes_present <- hallmark_emt[hallmark_emt %in% rownames(counts_trim)]
    cat("EMT genes present in expression matrix:", length(emt_genes_present),
        "of", length(hallmark_emt), "\n")
    
    ## Z-score each gene across cell lines (row-wise z-score)
    counts_mat <- as.matrix(counts_trim)
    counts_z   <- t(scale(t(counts_mat)))   # scale() works on columns, so transpose twice
    
    ## Subset to EMT genes and average per cell line
    emt_z <- counts_z[emt_genes_present, , drop = FALSE]
    emt_sum_z <- colSums(emt_z, na.rm = TRUE)
    
    ## Format to match ssGSEA output structure (cell lines x 1 column named "EMT")
    EMT_scores <- data.frame(
      EMT = emt_sum_z,
      row.names = names(emt_sum_z)
    )
    
    cat("EMT scores computed via sum z-score across", length(emt_genes_present), "genes")
    
  }
  
  cat("EMT score range:", round(range(EMT_scores$EMT), 3), "\n")
  cat("Number of cell lines scored:", nrow(EMT_scores), "\n")
  
  ### Flip for graphic purposes
  EMT_scores <- EMT_scores * -1
  
  write.table(
    x = EMT_scores,
    file = paste0(path.mel, "CCLE_Exhaustion_HALLMARK_EMT_", score_method, ".txt"),
    sep = "\t",
    quote = F,
    row.names =
  )
  
}

#### Co-Rank to Compare Z-score vs ssGSEA score
if(1) {
  
  ssGSEA_rank <- read.table(
    file = paste0(path.mel, "CCLE_Exhaustion_HALLMARK_EMT_ssgsea.txt"),
    header = T, 
    sep = "\t"
  )
  
  ZScore_rank <- read.table(
    file = paste0(path.mel, "CCLE_Exhaustion_HALLMARK_EMT_zscore.txt"),
    header = T, 
    sep = "\t"
  )
  
  shared_samples <- intersect(rownames(ssGSEA_rank), rownames(ZScore_rank))
  
  ssGSEA_rank_shared <- ssGSEA_rank %>%
    dplyr::filter(rownames(.) %in% shared_samples) %>%
    dplyr::mutate(rank_ssGSEA = c(1:nrow(.))) %>%
    dplyr::select(rank_ssGSEA)
  
  ZScore_rank_shared <- ZScore_rank %>%
    dplyr::filter(rownames(.) %in% shared_samples) %>%
    dplyr::mutate(rank_ZScore = c(1:nrow(.))) %>%
    dplyr::select(rank_ZScore)
  
  merged_df_all <- merge(ssGSEA_rank_shared, ZScore_rank_shared, by = "row.names")
  
  length <- dim(merged_df_all)[1]
  
  # point <- merged_df_all %>% dplyr::filter(Row.names == "")
  
  pdf(file = paste0(path.mel, "Plot_CoRank_EMT_Scoring.pdf"), width = 7.25, height = 6.5)
  ggplot(merged_df_all, aes(x = rank_ZScore, y = rank_ssGSEA)) +
    theme(axis.ticks = element_blank()) +
    stat_density_2d(
      aes(fill = ..density..),
      geom = "raster", 
      contour = FALSE
    ) +
    scale_fill_distiller(
      palette= "Spectral",
      name = "Density"
    ) +
    scale_x_continuous(
      breaks = c(length*.5),
      labels = c("Z-Score Rank"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(length*.5),
      labels = c("ssGSEA Rank"),
      expand = c(0, 0)
    ) +
    theme(
      legend.position = "right", 
      panel.border = element_rect(colour = "black", fill=NA, size=.75), 
      text = element_text(size = 20),
      axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 18),
      axis.text.x = element_text(size = 18),
      plot.title = element_text(size = 14)
    ) +
    xlab("") + ylab("") +
    geom_point(size = .25, alpha = 0.2) +
    # geom_point(
    #   data = point,
    #   color = "black",
    #   fill = "#E23F44",
    #   shape = 21,
    #   size = 4,
    #   stroke = 0.6
    # )
  # geom_text_repel(
  #   data = point,
  #   aes(label = Row.names),
  #   color = "#E23F44",
  #   size = 7.5,
  #   fontface = "bold",
  #   nudge_y = 50,
  #   nudge_x = 50,
  #   box.padding = 0.25,
  #   point.padding = 0.3,
  #   max.overlaps = Inf
  # )
  labs(title = "Hallmart EMT Signature Score")
  dev.off()

}

#### Correlate Differentiation Score with PLS Variate Scores (Top 10 Components)
score_method <- "ssgsea"

if (1) {
  
  ## Read in signature score
  EMT_scores_cor <- read.table(
    file   = paste0(path.mel, "CCLE_Exhaustion_HALLMARK_EMT_", score_method, ".txt"),
    sep    = "\t",
    header = TRUE
  )
  
  ## Load PLS-C variates (all components)
  PLS_variates_all <- read.delim(
    file             = paste0(path.pls, "PLS_Mode.canonical_X.CRISPR_Y.CTRP_X.variates.txt"),
    sep              = "\t",
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  
  ## Identify component columns (assumes naming like comp1, comp2, ...)
  comp_cols <- grep("^comp", colnames(PLS_variates_all), value = TRUE)
  comp_cols <- comp_cols[seq_len(min(10, length(comp_cols)))]
  
  ## Build named variate matrix: cell lines x components
  variate_mat <- PLS_variates_all %>%
    tibble::column_to_rownames(var = "Score") %>%
    dplyr::select(dplyr::all_of(comp_cols))
  
  ## Find common cell lines
  common_ids_cor <- intersect(rownames(EMT_scores_cor), rownames(variate_mat))
  cat("Cell lines used for variate-EMT correlation:", length(common_ids_cor), "\n")
  
  EMT_vec        <- EMT_scores_cor[common_ids_cor, "EMT"]
  names(EMT_vec) <- common_ids_cor
  
  ## Correlate each component with the differentiation score
  cor_results <- data.frame(
    component = comp_cols,
    r = sapply(comp_cols, function(comp) {
      cor(variate_mat[common_ids_cor, comp], EMT_vec, use = "pairwise.complete.obs")
    }),
    stringsAsFactors = FALSE
  )
  
  cor_results$r_abs     <- abs(cor_results$r)
  cor_results$direction <- ifelse(cor_results$r >= 0, "Positive", "Negative")
  cor_results           <- cor_results %>% dplyr::arrange(dplyr::desc(r_abs))
  cor_results$rank      <- seq_len(nrow(cor_results))
  cor_results$component <- factor(cor_results$component, levels = cor_results$component)
  
  cat("Component correlations with", score_method, "EMT score (ranked by |r|):\n")
  print(cor_results[, c("component", "r", "r_abs", "rank")])
  
  write.table(
    x         = cor_results,
    file      = paste0(path.mel, "PLS_Variate_EMT_Correlation_Top10_", score_method, ".txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  p_variate_cor <- ggplot2::ggplot(
    cor_results,
    ggplot2::aes(x = component, y = r, fill = direction)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0("r = ", round(r, 3)),
        vjust = ifelse(r >= 0, -0.4, 1.3)
      ),
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c("Positive" = "#1a3f6f", "Negative" = "#7b2d8b"),
      name   = "Direction"
    ) +
    ggplot2::labs(
      x     = "PLS Component (ranked by |r|)",
      y     = paste0("Pearson r (Variate Score vs EMT Score [", score_method, "])"),
      title = paste0("PLS Variate vs Differentiation Status | ", score_method)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.title      = ggplot2::element_text(size = 12)
    )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "PLS_Variate_EMT_Correlation_Top10_", score_method, ".pdf"),
    plot     = p_variate_cor,
    width    = 7,
    height   = 5
  )
  
  cat("Saved variate-EMT correlation plot | score_method:", score_method, "\n")
  
  
  ## Read in Model file for cancer type annotation
  model_cor <- read.csv(paste0(path.dm, "Model.csv"))
  
  df_cor <- data.frame(
    ModelID = common_ids_cor,
    EMT     = EMT_vec
  )
  df_cor <- merge(df_cor, model_cor[, c("ModelID", "OncotreeLineage")], by = "ModelID")
  
  ## Order lineages by median EMT score for clean plotting
  lineage_order <- df_cor %>%
    dplyr::group_by(OncotreeLineage) %>%
    dplyr::summarise(median_EMT = median(EMT, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(median_EMT)) %>%
    dplyr::pull(OncotreeLineage)
  
  df_cor$OncotreeLineage <- factor(df_cor$OncotreeLineage, levels = lineage_order)
  
  ## Boxplot: differentiation score distribution per cancer type
  p_box <- ggplot2::ggplot(df_cor, ggplot2::aes(x = OncotreeLineage, y = EMT, fill = OncotreeLineage)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::labs(
      x     = "",
      y     = paste0("EMT Score [", score_method, "]"),
      title = paste0("Differentiation Score Distribution by Cancer Type | ", score_method)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "none",
      plot.title      = ggplot2::element_text(size = 12)
    )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "EMT_Score_Boxplot_ByCancerType_", score_method, ".pdf"),
    plot     = p_box,
    width    = 14,
    height   = 5
  )
  
  cat("Saved EMT score boxplot by cancer type\n")
  
  #### Partial correlation: variate vs EMT score controlling for cancer type
  ## Residualize both EMT and each variate against OncotreeLineage
  
  ## Restrict variate_mat to cell lines in df_cor (post-annotation)
  common_ids_partial <- df_cor$ModelID
  EMT_vec_partial    <- df_cor$EMT
  names(EMT_vec_partial) <- df_cor$ModelID
  
  lineage_vec <- df_cor$OncotreeLineage
  names(lineage_vec) <- df_cor$ModelID
  
  resid_EMT <- residuals(lm(EMT_vec_partial ~ lineage_vec))
  
  partial_cor_results <- data.frame(
    component = comp_cols,
    r_partial = sapply(comp_cols, function(comp) {
      variate_scores        <- variate_mat[common_ids_partial, comp]
      names(variate_scores) <- common_ids_partial
      resid_variate         <- residuals(lm(variate_scores ~ lineage_vec))
      cor(resid_variate, resid_EMT, use = "pairwise.complete.obs")
    }),
    stringsAsFactors = FALSE
  )
  
  partial_cor_results$r_abs     <- abs(partial_cor_results$r_partial)
  partial_cor_results$direction <- ifelse(partial_cor_results$r_partial >= 0, "Positive", "Negative")
  partial_cor_results           <- partial_cor_results %>% dplyr::arrange(dplyr::desc(r_abs))
  partial_cor_results$rank      <- seq_len(nrow(partial_cor_results))
  partial_cor_results$component <- factor(partial_cor_results$component, levels = partial_cor_results$component)
  
  cat("Partial correlations (controlling for cancer type) with", score_method, "EMT score (ranked by |r|):\n")
  print(partial_cor_results[, c("component", "r_partial", "r_abs", "rank")])
  
  write.table(
    x         = partial_cor_results,
    file      = paste0(path.mel, "PLS_Variate_EMT_PartialCorrelation_Top10_", score_method, ".txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  p_partial_cor <- ggplot2::ggplot(
    partial_cor_results,
    ggplot2::aes(x = component, y = r_partial, fill = direction)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0("r = ", round(r_partial, 3)),
        vjust = ifelse(r_partial >= 0, -0.4, 1.3)
      ),
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c("Positive" = "#1a3f6f", "Negative" = "#7b2d8b"),
      name   = "Direction"
    ) +
    ggplot2::labs(
      x     = "PLS Component (ranked by |r|)",
      y     = paste0("Partial r (Variate vs EMT Score [", score_method, "] | Cancer Type)"),
      title = paste0("PLS Variate vs Differentiation Status\nControlling for Cancer Type | ", score_method)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.title      = ggplot2::element_text(size = 12)
    )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "PLS_Variate_EMT_PartialCorrelation_Top10_", score_method, ".pdf"),
    plot     = p_partial_cor,
    width    = 7,
    height   = 5
  )
  
  cat("Saved partial correlation plot | score_method:", score_method, "\n")
  
}

#### Plot Figure

## Set Parameters
score_method <- "ssgsea"    # zscore or ssgsea
weight       <- "function"  # function (standard, WLS via weights argument, abs weights) or prior (direct multiply wx transformation, signed weights. not the standard)
comp <- "comp2"

# Lineage filtering — uses exact OncotreeLineage values from Model.csv (one or both must = NULL)
keep <- NULL # only retain these lineages  e.g. c("Skin") or NULL, or epithelial c("Biliary Tract", "Bladder/Urinary Tract", "Bowel", "Breast", "Cervix", "Esophagus/Stomach", "Head and Neck", "Kidney", "Liver", "Lung", "Ovary/Fallopian Tube", "Pancreas", "Prostate", "Skin", "Thyroid", "Uterus", "Vulva/Vagina", "Ampulla of Vater", "Pleura")
remove <- NULL # drop these lineages  e.g. c("Lymphoid", "Myeloid") or NULL

# Adrenal Gland, Ampulla of Vater, Biliary Tract, Bladder/Urinary Tract, Bone, Bowel, Breast, Cervix, CNS/Brain, Embryonal, Esophagus/Stomach, Eye, Fibroblast, Hair, Head and Neck, Kidney, Liver, Lung, Lymphoid, Muscle, Myeloid , Normal, Other, Ovary/Fallopian Tube, Pancreas Peripheral Nervous System, Pleura, Prostate, Skin, Soft Tissue, Testis, Thyroid, Uterus, Vulva/Vagina

if (1) {
  
  ## Read in loadings to see where MED12 is pulled out the most
  PLS_loadings <- read.delim(
    file = paste0(path.pls, "PLS_Mode.canonical_X.CRISPR_Y.CTRP_X.loadings.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F
  ) %>%
    dplyr::filter(Loading == "GPX4") %>%
    tibble::column_to_rownames(var = "Loading")
  
  # MED12 strongest on comp 3 (negative side)
  # GPX4 strongest on comp 3 (psotive side)
  
  ## Read in PLS-C variates
  PLS_variates <- read.delim(
    file = paste0(path.pls, "PLS_Mode.canonical_X.CRISPR_Y.CTRP_X.variates.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F
  )
  
  weights_comp <- PLS_variates[[comp]]
  names(weights_comp) <- PLS_variates$Score
  
  # Read in signature score
  EMT_scores <- read.table(
    file = paste0(path.mel, "CCLE_Exhaustion_HALLMARK_EMT_", score_method, ".txt"),
    sep = "\t"
  )
  
  ## Read in CRISPR Data
  CRISPR <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  ## Read in Model file
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  ## Find common cell lines (all, before filtering)
  common_ids <- Reduce(intersect, list(
    rownames(EMT_scores),
    rownames(CRISPR),
    names(weights_comp)
  ))
  cat("Number of common cell lines (pre-filter):", length(common_ids))
  
  ## Build merged data frame
  df <- data.frame(
    ModelID  = common_ids,
    EMT      = EMT_scores[common_ids, "EMT"],
    GPX4     = CRISPR[common_ids, "GPX4"],
    MED12    = CRISPR[common_ids, "MED12"],
    W        = weights_comp[common_ids]
  )
  df <- merge(df, model[, c("ModelID", "OncotreeLineage", "OncotreePrimaryDisease")], by = "ModelID")
  
  ## Apply lineage filter
  if (!is.null(keep) && length(keep) > 0) {
    
    keep_pattern <- paste(keep, collapse = "|")
    df <- df %>%
      dplyr::filter(grepl(keep_pattern, OncotreeLineage, ignore.case = TRUE))
    cat("Keeping lineages matching:", paste(keep, collapse = ", "))
    
  } else if (!is.null(remove) && length(remove) > 0) {
    
    df <- df %>%
      dplyr::filter(!OncotreeLineage %in% remove)
    cat("Removing lineages:", paste(remove, collapse = ", "))
    
  } else {
    cat("No lineage filter applied — using all cell lines")
  }
  cat("Number of cell lines after filtering:", nrow(df))
  
  ## Absolute weights (required for WLS; also used for point sizing)
  df$W_abs        <- abs(df$W)
  df$variate_sign <- ifelse(df$W >= 0, "Positive Variate", "Negative Variate")
  
  ## Cancer type grouping
  df$cancer_group <- dplyr::case_when(
    # grepl("Brain",           df$OncotreeLineage, ignore.case = TRUE) ~ "Brain Cancer",
    # grepl("Stomach", df$OncotreeLineage, ignore.case = TRUE) ~ "Gastric Cancer",
    # grepl("Head|Neck",       df$OncotreeLineage, ignore.case = TRUE) ~ "Head and Neck Cancer",
    # grepl("Kidney|Renal",    df$OncotreeLineage, ignore.case = TRUE) ~ "Kidney Cancer",
    # grepl("Skin|Melanoma",   df$OncotreeLineage, ignore.case = TRUE) ~ "Skin Cancer",
    # TRUE ~ "Other"
    
    grepl("CNS/Brain",           df$OncotreeLineage, ignore.case = TRUE) ~ "CNS/Brain",
    grepl("Head|Neck",       df$OncotreeLineage, ignore.case = TRUE) ~ "Head and Neck Cancer",
    
    grepl("Myeloid|Lymphoid", df$OncotreeLineage, ignore.case = TRUE) ~ "Hematological",
    grepl("Bowel",    df$OncotreeLineage, ignore.case = TRUE) ~ "Bowel",
    grepl("Skin|Melanoma",   df$OncotreeLineage, ignore.case = TRUE) ~ "Skin Cancer",
    TRUE ~ "Other"
    
  )
  
  ## Reorder so "Other" plots first (behind)
  df <- df[order(df$cancer_group == "Other", decreasing = TRUE), ]
  
  ## Update common_ids to match filtered df
  common_ids <- df$ModelID
  
  ## Build save name from parameters
  weight_tag <- ifelse(weight == "function", "WLS", "PRIOR")
  
  filter_tag <- dplyr::case_when(
    !is.null(keep)   & length(keep)   > 0 ~ paste0("KEEP.", paste(keep,   collapse = ".")),
    !is.null(remove) & length(remove) > 0 ~ paste0("REM.",  paste(remove, collapse = ".")),
    TRUE ~ "ALL"
  )
  
  ## Genes to be highlighted in plots
  genes_gpx4 <- c("GPX4", "SEPSECS", "PSTK", "EEFSEC", "SECISBP2", "SEPHS2")
  genes_med  <- c("MED24", "MED16", "MED25", "MED10", "MED19", "MED1", "MED23", "MED9", "MED15", "MED12")
  
  #### Genome-wide CRISPR correlation with EMT signature
  if (0) {
    
    ## Subset CRISPR to filtered cell lines
    CRISPR_filtered <- CRISPR[common_ids, ]
    EMT_filtered    <- EMT_scores[common_ids, "EMT"]
    
    ## Correlate every gene with EMT score
    genome_wide_cor <- sapply(colnames(CRISPR_filtered), function(g) {
      cor(CRISPR_filtered[, g], EMT_filtered, use = "pairwise.complete.obs")
    })
    
    ## Format results
    genome_wide_df <- data.frame(
      gene = names(genome_wide_cor),
      r    = genome_wide_cor,
      row.names = NULL,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::arrange(dplyr::desc(r))
    
    ## Add rank and highlight columns
    genome_wide_df$rank       <- seq_len(nrow(genome_wide_df))
    genome_wide_df$gene_group <- dplyr::case_when(
      genome_wide_df$gene %in% genes_gpx4 ~ "GPX4 module",
      genome_wide_df$gene %in% genes_med  ~ "MED12 module",
      TRUE                                ~ "Other"
    )
    
    cat("Top 10 positively correlated genes:\n")
    print(head(genome_wide_df, 10))
    cat("Top 10 negatively correlated genes:\n")
    print(tail(genome_wide_df, 10))
    
    ## Save full results
    write.table(
      x         = genome_wide_df,
      file      = paste0(path.mel, "GenomeWide_CRISPR_EMT_Correlation_", filter_tag, "_", score_method, ".txt"),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
    
    ## Plot: ranked correlation (waterfall), highlighting gene modules
    p_genome <- ggplot2::ggplot(genome_wide_df, ggplot2::aes(x = rank, y = r, color = gene_group)) +
      ggplot2::geom_point(
        data = genome_wide_df %>% dplyr::filter(gene_group == "Other"),
        size = 0.3, alpha = 0.4
      ) +
      ggplot2::geom_point(
        data = genome_wide_df %>% dplyr::filter(gene_group != "Other"),
        size = 2.5, alpha = 1
      ) +
      ggrepel::geom_text_repel(
        data        = genome_wide_df %>% dplyr::filter(gene_group != "Other"),
        ggplot2::aes(label = gene),
        size        = 3,
        max.overlaps = Inf,
        box.padding  = 0.3,
        point.padding = 0.2
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "GPX4 module"  = "#1a3f6f",
          "MED12 module" = "#4b0070",
          "Other"        = "#C0C0C0"
        ),
        name = ""
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      ggplot2::labs(
        x = "Gene Rank",
        y = "Pearson r (CRISPR Dependency vs EMT Score)",
        title = paste0("Genome-wide CRISPR ~ EMT Correlation | ", filter_tag, " | ", score_method)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
    
    ggplot2::ggsave(
      filename = paste0(path.plots, "GenomeWide_CRISPR_EMT_Correlation_", filter_tag, "_", score_method, ".pdf"),
      plot     = p_genome,
      width    = 10,
      height   = 6
    )
    
    cat("Saved genome-wide correlation plot and table.")
    
  }
  
  ## Define weighted_r based on weight switch
  if (weight == "function") {
    
    weighted_r <- function(x, y, w) {
      fit       <- lm(y ~ x, weights = abs(w)) # lm (response ~ predictor)
      r2        <- summary(fit)$r.squared
      sign_coef <- sign(coef(fit)["x"])
      r         <- sign_coef * sqrt(r2)
      return(list(r = r, r2 = r2, fit = fit, beta = coef(fit)))
    }
    
    smooth_weight <- df$W_abs
    bar_w         <- abs(weights_comp[common_ids])
    
  } else {
    
    weighted_r <- function(x, y, w) {
      wx        <- w * x
      fit       <- lm(y ~ wx) # lm (response ~ predictor)
      r2        <- summary(fit)$r.squared
      sign_coef <- sign(coef(fit)["wx"])
      r         <- sign_coef * sqrt(r2)
      return(list(r = r, r2 = r2, fit = fit, beta = coef(fit)))
    }
    
    smooth_weight <- df$W_abs
    bar_w         <- weights_comp[common_ids]
    
  }
  
  ## Compute weighted r for GPX4 and MED12
  gpx4_wr  <- weighted_r(df$GPX4,  df$EMT, df$W)
  med12_wr <- weighted_r(df$MED12, df$EMT, df$W)
  
  ## Unweighted correlations for comparison
  gpx4_unw  <- cor(df$GPX4,  df$EMT)
  med12_unw <- cor(df$MED12, df$EMT)
  
  cat("GPX4  unweighted r:", round(gpx4_unw,  4), " | weighted r:", round(gpx4_wr$r,  4))
  cat("MED12 unweighted r:", round(med12_unw, 4), " | weighted r:", round(med12_wr$r, 4))
  
  ## Color palette
  cancer_colors <- c(
    # "Brain Cancer"         = "#87CEEB",
    # "Gastric Cancer"       = "#00BFFF",
    # "Head and Neck Cancer" = "#FFA500",
    # "Kidney Cancer"        = "#90EE90",
    # "Skin Cancer"          = "#FFB6C1",
    # "Other"                = "#C0C0C0"
    
    "CNS/Brain"         = "#850014",
    "Head and Neck Cancer"       = "#00BFFF",
    "Hematological" = "#FFA500",
    "Bowel"        = "#90EE90",
    "Skin Cancer"          = "#FFB6C1",
    "Other"                = "#C0C0C0"
  )
  
  ## Panel A: Bar chart (side-by-side raw unweighted vs weighted)
  all_genes  <- c(genes_gpx4, genes_med)
  all_genes  <- all_genes[all_genes %in% colnames(CRISPR)]
  
  bar_df <- data.frame(
    gene = all_genes,
    unweighted = sapply(all_genes, function(g) {
      cor(CRISPR[common_ids, g], EMT_scores[common_ids, "EMT"])
    }),
    weighted = sapply(all_genes, function(g) {
      weighted_r(
        x = CRISPR[common_ids, g],
        y = EMT_scores[common_ids, "EMT"],
        w = bar_w
      )$r
    }),
    stringsAsFactors = FALSE
  )
  
  bar_df$gene_group <- ifelse(bar_df$gene %in% genes_gpx4, "GPX4 module", "MED12 module")
  bar_df$gene       <- factor(bar_df$gene, levels = all_genes)
  
  bar_long <- bar_df %>%
    tidyr::pivot_longer(cols = c(unweighted, weighted), names_to = "type", values_to = "value") %>%
    dplyr::mutate(
      fill_group = dplyr::case_when(
        type == "unweighted" & gene_group == "GPX4 module"  ~ "GPX4 unweighted",
        type == "weighted"   & gene_group == "GPX4 module"  ~ "GPX4 weighted",
        type == "unweighted" & gene_group == "MED12 module" ~ "Mediator unweighted",
        type == "weighted"   & gene_group == "MED12 module" ~ "Mediator weighted"
      ),
      fill_group = factor(fill_group, levels = c(
        "GPX4 unweighted", "GPX4 weighted",
        "Mediator unweighted", "Mediator weighted"
      ))
    )
  
  bar_colors <- c(
    "GPX4 unweighted"     = "#1a3f6f",  # dark blue
    "GPX4 weighted"       = "#6aaed6",  # light blue
    "Mediator unweighted" = "#4b0070",  # dark purple
    "Mediator weighted"   = "#b07cc6"   # light purple
  )
  
  p_bar <- ggplot2::ggplot(bar_long, ggplot2::aes(x = gene, y = value, fill = fill_group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = bar_colors, name = "") +
    ggplot2::labs(x = "", y = "Correlation Coefficients") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  ## Panel B: GPX4 scatter
  p_gpx4 <- ggplot2::ggplot(df, ggplot2::aes(x = EMT, y = GPX4)) +
    ggplot2::geom_point(
      ggplot2::aes(size = W_abs, shape = variate_sign, color = cancer_group),
      alpha = 0.6
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(weight = smooth_weight), method = "lm", se = FALSE,
      color = "black", linewidth = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0(
        "r (unweighted) = ", round(gpx4_unw, 3),
        "\nr (weighted) = ",  round(gpx4_wr$r, 3)
      ),
      size = 3, fontface = "italic"
    ) +
    ggplot2::scale_color_manual(values = cancer_colors, name = "Cancer Types") +
    ggplot2::scale_shape_manual(
      values = c("Positive Variate" = 16, "Negative Variate" = 17),
      name = "Variate Score"
    ) +
    ggplot2::scale_size_continuous(range = c(0.5, 6), name = "Weight") +
    ggplot2::labs(
      # x = "Differentiation Score\n
      # Differentiated \u2190  \u2192 Un-differentiated
      # \n Epithelial \u2190  \u2192 Mesenchymal",
      # y = "GPX4 CRISPR Dependency\nMore Dependent \u2190\u2192 Less Dependent"
      
      x = "Differentiation Score",
      y = "GPX4 CRISPR Dependency"
      
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  ## Panel C: MED12 scatter
  p_med12 <- ggplot2::ggplot(df, ggplot2::aes(x = EMT, y = MED12)) +
    ggplot2::geom_point(
      ggplot2::aes(size = W_abs, shape = variate_sign, color = cancer_group),
      alpha = 0.6
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(weight = smooth_weight), method = "lm", se = FALSE,
      color = "black", linewidth = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0(
        "r (unweighted) = ", round(med12_unw, 3),
        "\nr (weighted) = ",  round(med12_wr$r, 3)
      ),
      size = 3, fontface = "italic"
    ) +
    ggplot2::scale_color_manual(values = cancer_colors, name = "Cancer Types") +
    ggplot2::scale_shape_manual(
      values = c("Positive Variate" = 16, "Negative Variate" = 17),
      name = "Variate Score"
    ) +
    ggplot2::scale_size_continuous(range = c(0.5, 6), name = "Weight") +
    ggplot2::labs(
    # x = "Differentiation Score\n
    #   Differentiated \u2190  \u2192 Un-differentiated
    #   \n Epithelial \u2190  \u2192 Mesenchymal",
    # y = "MED12 CRISPR Dependency"
      
      x = "Differentiation Score",
      y = "MED12 CRISPR Dependency"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  ## Shared legend
  p_legend <- ggplot2::ggplot(df, ggplot2::aes(x = EMT, y = GPX4)) +
    ggplot2::geom_point(
      ggplot2::aes(size = W_abs, shape = variate_sign, color = cancer_group),
      alpha = 0.6
    ) +
    ggplot2::scale_color_manual(values = cancer_colors, name = "Cancer Types") +
    ggplot2::scale_shape_manual(
      values = c("Positive Variate" = 16, "Negative Variate" = 17),
      name = "Variate Score"
    ) +
    ggplot2::scale_size_continuous(range = c(0.5, 6), name = "Weight") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
  
  shared_legend <- cowplot::get_legend(p_legend)
  
  ## Combine panels
  bottom_row <- cowplot::plot_grid(
    p_gpx4, p_med12, shared_legend,
    ncol = 3, rel_widths = c(1, 1, 0.4),
    labels = c("B", "C", "")
  )
  
  full_fig <- cowplot::plot_grid(
    p_bar, bottom_row,
    nrow = 2, rel_heights = c(1, 1),
    labels = c("A", "")
  )
  
  save_name <- paste0(
    "EMT_Differentiation_GPX4_MED12_Figure3",
    "_", filter_tag,
    "_", weight_tag,
    "_", score_method,
    "_", comp,
    ".pdf"
  )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, save_name),
    plot = full_fig, width = 10, height = 8
  )
  
  cat("Saved:", save_name)
  
}

#### Export Individual Panels
if (1) {
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "Panel_A_Bar_", filter_tag, "_", weight_tag, "_", score_method, "_", comp, ".pdf"),
    plot     = p_bar,
    width    = 11.5,
    height   = 4
  )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "Panel_B_GPX4_", filter_tag, "_", weight_tag, "_", score_method, "_", comp, ".pdf"),
    plot     = p_gpx4,
    width    = 4.5,
    height   = 3.5
  )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "Panel_C_MED12_", filter_tag, "_", weight_tag, "_", score_method, "_", comp, ".pdf"),
    plot     = p_med12,
    width    = 4.5,
    height   = 3.5
  )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "Panel_C_MED12_", filter_tag, "_", weight_tag, "_", score_method, "_", comp, "_FORLEGEND.pdf"),
    plot     = p_med12 + ggplot2::theme(legend.position = "right"),
    width    = 6,
    height   = 5
  )
  
  cat("Saved individual panels: A (bar), B (GPX4), C (MED12)\n")
  
}


#### KS Test: Cancer Type Enrichment Along Differentiation Score (All Lineages)
if (1) {
  
  ## Rank all cell lines by EMT score (ascending = most epithelial first)
  df_ks <- df %>%
    dplyr::arrange(EMT) %>%
    dplyr::mutate(rank_EMT = seq_len(nrow(.)))
  
  ## Filter to lineages with enough cell lines to test
  min_n <- 5
  lineage_counts <- table(df_ks$OncotreeLineage)
  lineages_test  <- names(lineage_counts[lineage_counts >= min_n])
  cat("Testing", length(lineages_test), "lineages with n >=", min_n, "\n")
  
  ks_results_lineage <- purrr::map_dfr(lineages_test, function(lin) {
    
    in_group  <- df_ks$EMT[df_ks$OncotreeLineage == lin]
    out_group <- df_ks$EMT[df_ks$OncotreeLineage != lin]
    
    ks     <- ks.test(in_group, out_group)
    med_in  <- median(in_group,  na.rm = TRUE)
    med_out <- median(out_group, na.rm = TRUE)
    
    data.frame(
      lineage        = lin,
      n              = length(in_group),
      median_EMT     = round(med_in,  4),
      median_other   = round(med_out, 4),
      direction      = ifelse(med_in > med_out, "Mesenchymal-enriched", "Epithelial-enriched"),
      D_statistic    = round(ks$statistic, 4),
      p_value        = ks$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  ks_results_lineage$p_adj <- p.adjust(ks_results_lineage$p_value, method = "BH")
  ks_results_lineage$sig   <- dplyr::case_when(
    ks_results_lineage$p_adj < 0.001 ~ "***",
    ks_results_lineage$p_adj < 0.01  ~ "**",
    ks_results_lineage$p_adj < 0.05  ~ "*",
    TRUE                             ~ "ns"
  )
  ks_results_lineage <- ks_results_lineage %>% dplyr::arrange(p_adj)
  
  cat("KS test results (all lineages, ranked by p_adj):\n")
  print(ks_results_lineage)
  
  write.table(
    x         = ks_results_lineage,
    file      = paste0(path.mel, "KS_AllLineages_EMT_Enrichment_", filter_tag, "_", score_method, ".txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  ## Plot: lollipop of D statistic, signed by direction, colored by significance
  ks_results_lineage$D_signed <- ifelse(
    ks_results_lineage$direction == "Mesenchymal-enriched",
    ks_results_lineage$D_statistic,
    -ks_results_lineage$D_statistic
  )
  ks_results_lineage$lineage <- factor(
    ks_results_lineage$lineage,
    levels = ks_results_lineage$lineage[order(ks_results_lineage$D_signed)]
  )
  
  p_ks_lineage <- ggplot2::ggplot(
    ks_results_lineage,
    ggplot2::aes(x = D_signed, y = lineage, color = sig)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = D_signed, y = lineage, yend = lineage),
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::scale_color_manual(
      values = c("***" = "#b2182b", "**" = "#ef8a62", "*" = "#fdbf6f", "ns" = "grey70"),
      name   = "BH-adjusted p"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0("n=", n)),
      hjust  = -0.3,
      size   = 2.5,
      color  = "grey40"
    ) +
    ggplot2::labs(
      x     = expression("KS D statistic (signed: negative" %->% "Epithelial, positive" %->% "Mesenchymal)"),
      y     = "",
      title = paste0("KS Test: All Lineages Along Differentiation Axis | ", score_method)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.y     = ggplot2::element_text(size = 8)
    )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "KS_AllLineages_EMT_Enrichment_", filter_tag, "_", score_method, ".pdf"),
    plot     = p_ks_lineage,
    width    = 9,
    height   = max(4, length(lineages_test) * 0.28)
  )
  
  cat("Saved all-lineage KS enrichment plot and table\n")
  
}


##### Drug Potentiality Score: PLS-C (CRISPR×CTRP) + PCA #####

### Set OS ###
OS <- "Mac" # Linux or Mac

if (OS == "Mac") {
  path.OS <- "/Users/jack/Library/CloudStorage/Box-Box/"
} else {
  path.OS <- "/media/testuser/SSD_4/jfreeland/Freeland/Github/"
}

## Set paths
path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.pls     <- paste0(path.wd, "DataSets/PLS/")
path.rcca    <- paste0(path.wd, "DataSets/rCCA/")
path.pca     <- paste0(path.wd, "DataSets/PCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.dp      <- paste0(path.wd, "DataSets/DrugPotentiality/")

## Method switch: "PLS" or "rCCA"
method   <- "PLS"  # "PLS" or "rCCA"

## PLS parameters (used when method == "PLS")
pls_mode <- "canonical"
X_source <- "CRISPR"
Y_source <- "CTRP"

## rCCA parameters (used when method == "rCCA")
rcca_mode      <- "shrinkage"  # "shrinkage" or "ridge"
lambda1_manual <- 0.20
lambda2_manual <- 0.10

X_source <- "CRISPR"
Y_source <- "CTRP"

## Filtered for all three data sets shared lines? (mirrors PLS script)
FilteredAll3 <- TRUE # TRUE or FALSE

## Number of components to summarize over (paper uses top 10)
n_comp <- 5

## Build Filtered_Tag (appended to file names, mirrors PLS script)
Filtered_Tag <- if (FilteredAll3) "_Filtered3" else ""

## Build file_tag, path, and component column prefix based on method
## Build file_tag (input — matches how PLS/rCCA files were saved, no n_comp)
## Build output_tag (output — adds n_comp to distinguish results)
if (method == "PLS") {
  file_tag    <- paste0("PLS_Mode.", pls_mode, "_X.", X_source, "_Y.", Y_source)
  output_tag  <- paste0("PLS_Mode.", pls_mode, "_X.", X_source, "_Y.", Y_source, "_ncomp", n_comp)
  path_drug   <- path.pls
  gene_side   <- "X"
  comp_prefix <- "comp"
} else {
  if (rcca_mode == "shrinkage") {
    file_tag   <- paste0("RCCA_shrinkage_X.", X_source, "_Y.", Y_source)
    output_tag <- paste0("RCCA_shrinkage_X.", X_source, "_Y.", Y_source, "_ncomp", n_comp)
  } else {
    file_tag   <- paste0("RCCA_ridge_lambda1.", format(lambda1_manual, digits = 3),
                         "_lambda2.", format(lambda2_manual, digits = 3),
                         "_X.", X_source, "_Y.", Y_source)
    output_tag <- paste0("RCCA_ridge_lambda1.", format(lambda1_manual, digits = 3),
                         "_lambda2.", format(lambda2_manual, digits = 3),
                         "_X.", X_source, "_Y.", Y_source, "_ncomp", n_comp)
  }
  path_drug   <- path.rcca
  gene_side   <- "Y"
  comp_prefix <- "X"
}

### Step 1: Load drug-gene co-embedding loadings (PLS-C or rCCA)
if (1) {
  
  drug_loadings_raw <- read.delim(
    file = paste0(path_drug, file_tag, Filtered_Tag, "_", gene_side, ".loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading))
  
  ## Column names differ: PLS uses comp1/comp2, rCCA uses X1/X2
  comp_cols_drug <- paste0(comp_prefix, 1:n_comp)
  comp_cols_drug <- comp_cols_drug[comp_cols_drug %in% names(drug_loadings_raw)]
  
  drug_loadings_sub <- drug_loadings_raw %>%
    dplyr::select(Loading, all_of(comp_cols_drug))
  
  cat(method, gene_side, "-loadings (gene side):", nrow(drug_loadings_sub), "genes x", length(comp_cols_drug), "components\n")
}

### Step 2: Load PCA loadings (gene-only signal)
if (1) {
  
  ## These are the loadings from PCA_from_file() run on CRISPR_common
  ## File name follows the pattern: <input>_prcomp_loadings.txt
  PCA_loadings <- read.delim(
    file = paste0(path.pca, "CRISPR_common_PCA_prcomp_loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading))
  
  ## Keep only PC1:PC{n_comp}
  comp_cols_pca <- paste0("PC", 1:n_comp)
  comp_cols_pca <- comp_cols_pca[comp_cols_pca %in% names(PCA_loadings)]
  
  PCA_loadings_sub <- PCA_loadings %>%
    dplyr::select(Loading, all_of(comp_cols_pca))
  
  cat("PCA loadings: ", nrow(PCA_loadings_sub), "genes x", length(comp_cols_pca), "components\n")
}

### Step 3: Summarize — max absolute loading per gene across top 10 components
if (1) {
  
  summarize_max_abs <- function(df, gene_col, comp_cols, value_suffix) {
    df %>%
      tidyr::pivot_longer(
        cols      = all_of(comp_cols),
        names_to  = "component",
        values_to = "loading"
      ) %>%
      dplyr::mutate(abs_loading = abs(loading)) %>%
      dplyr::group_by(.data[[gene_col]]) %>%
      dplyr::slice_max(abs_loading, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::rename(
        !!paste0("component_", value_suffix) := component,
        !!paste0("loading_", value_suffix)    := loading,
        !!paste0("abs_loading_", value_suffix) := abs_loading
      )
  }
  
  drug_max <- summarize_max_abs(drug_loadings_sub, "Loading", comp_cols_drug, "drug")
  PCA_max  <- summarize_max_abs(PCA_loadings_sub,  "Loading", comp_cols_pca,  "PCA")
  
  cat("Genes summarized —", method, ":", nrow(drug_max), " PCA:", nrow(PCA_max), "\n")
}

### Step 4: Fit log-normal distribution and convert to quantiles
if (1) {
  
  ## The paper fits each vector (PLS summarized loadings, PCA summarized loadings)
  ## separately to a log-normal distribution and uses the fitted CDF to get quantiles.
  
  fit_lognormal_quantile <- function(x, label = "") {
    
    x_pos <- x[x > 0 & !is.na(x)]
    
    ## Fit log-normal (and competitors) to check AIC
    fit_ln  <- fitdistrplus::fitdist(x_pos, "lnorm")
    fit_wb  <- fitdistrplus::fitdist(x_pos, "weibull")
    fit_gm  <- fitdistrplus::fitdist(x_pos, "gamma")
    
    aic_table <- data.frame(
      distribution = c("lognormal", "weibull", "gamma"),
      AIC          = c(fit_ln$aic, fit_wb$aic, fit_gm$aic)
    )
    cat("\nAIC comparison for", label, "loadings:\n")
    print(aic_table)
    cat("Best fit:", aic_table$distribution[which.min(aic_table$AIC)], "\n")
    
    ## Use log-normal as per the paper regardless (paper states lognormal best)
    meanlog <- fit_ln$estimate["meanlog"]
    sdlog   <- fit_ln$estimate["sdlog"]
    
    ## Convert each gene's abs loading to a quantile [0,1] via fitted log-normal CDF
    ## Genes with abs_loading = 0 get quantile = 0
    q <- plnorm(x, meanlog = meanlog, sdlog = sdlog)
    q[is.na(x)] <- NA
    
    list(
      quantile  = q,
      fit       = fit_ln,
      aic_table = aic_table
    )
  }
  
  drug_fit <- fit_lognormal_quantile(drug_max$abs_loading_drug, label = method)
  PCA_fit  <- fit_lognormal_quantile(PCA_max$abs_loading_PCA,  label = "PCA")
  
  drug_max$quantile_drug <- drug_fit$quantile
  PCA_max$quantile_PCA   <- PCA_fit$quantile
}

### Step 5: Merge and compute Drug Potentiality Score
if (1) {
  
  DP <- merge(
    drug_max %>% dplyr::select(Loading, abs_loading_drug, quantile_drug, component_drug),
    PCA_max  %>% dplyr::select(Loading, abs_loading_PCA,  quantile_PCA,  component_PCA),
    by = "Loading"
  )
  
  DP <- DP %>%
    dplyr::mutate(
      ## Drug Potentiality Score: drug co-signal quantile minus PCA quantile
      drug_potentiality = quantile_drug - quantile_PCA,
      
      ## Magnitude: Euclidean distance to origin using summarized abs loadings
      magnitude = sqrt(abs_loading_drug^2 + abs_loading_PCA^2),
      
      ## Categorical label
      category = dplyr::case_when(
        drug_potentiality >  0.05  ~ "Drug Rainforest",
        drug_potentiality < -0.05  ~ "Drug Desert",
        TRUE                       ~ "Neutral"
      )
    ) %>%
    dplyr::arrange(desc(drug_potentiality))
  
  cat("\nDrug Potentiality Score summary:\n")
  print(summary(DP$drug_potentiality))
  cat("\nCategory counts:\n")
  print(table(DP$category))
}

### Step 6: Add gene annotations for labeling
if (1) {
  
  ## Auto: top 10 and bottom 10 by drug potentiality score
  top10    <- DP %>% dplyr::arrange(desc(drug_potentiality)) %>% dplyr::slice_head(n = 5) %>% dplyr::pull(Loading)
  bottom10 <- DP %>% dplyr::arrange(drug_potentiality)      %>% dplyr::slice_head(n = 5) %>% dplyr::pull(Loading)
  
  ## Manual curation list (from paper figure + your existing annotation groups)
  genes_manual <- c(
    ## Desert
    "KRAS", "PTEN", "MED12",
    ## Forest
    "VEGFA", "EPCAM", "ERG",
    
    # "NGF", "ITGA5", "ITGA6", "ITGAB4", "ERG", "SOX5", "ELF4", "KDR", "EFG", "TAP1", "ICAM2", "SLAMF1", "CDC42BPA", "NCAM1", "VEGFA", "GABARAP", "VRK2",
    ## Neutral
    "PIK3CA", "GPX4", "BRAF", "MED12"
  )
  
  all_labels <- unique(c(top10, bottom10, genes_manual))
  DP$label   <- ifelse(DP$Loading %in% all_labels, DP$Loading, NA_character_)
  
  ## Tag label source — drives separate text colors in the plot
  DP$label_group <- dplyr::case_when(
    !is.na(DP$label) & DP$Loading %in% genes_manual ~ "manual",  # manual takes priority
    DP$Loading %in% top10                            ~ "top10",
    DP$Loading %in% bottom10                         ~ "bottom10",
    TRUE                                             ~ NA_character_
  )
}

### Step 7: Save output table
if (1) {
  
  write.table(
    x         = DP,
    file      = paste0(path.dp, "DrugPotentiality_", output_tag, Filtered_Tag, ".txt"),
    sep       = "\t", quote = FALSE, row.names = FALSE
  )
  cat("Saved drug potentiality table.\n")
}

### Step 8: Plot — replicating Figure 4B from the paper
if (1) {
  
  ## Percentile cutoffs — 25th and 75th define neutral band boundaries
  ## 10th/90th, 5th/95th, 1st/99th define the forest/desert gradient steps
  cutoffs <- quantile(DP$drug_potentiality, probs = c(0.01, 0.05, 0.10, 0.25, 0.75, 0.90, 0.95, 0.99))
  cat("\nPercentile cutoffs:\n"); print(round(cutoffs, 4))
  
  ## 9 bands total: 4 forest + 1 neutral (middle 50%) + 4 desert
  ## Neutral = genes between 25th and 75th percentile of drug_potentiality
  ## 50th Forest = genes between 75th percentile and 90th percentile
  ## 50th Desert = genes between 10th percentile and 25th percentile
  DP <- DP %>%
    dplyr::mutate(
      color_band = dplyr::case_when(
        drug_potentiality >= cutoffs["99%"] ~ "01_forest_1st",   # top 1%       — darkest green
        drug_potentiality >= cutoffs["95%"] ~ "02_forest_5th",   # top 1–5%     — dark green
        drug_potentiality >= cutoffs["90%"] ~ "03_forest_10th",  # top 5–10%    — mid green
        drug_potentiality >= cutoffs["75%"] ~ "04_forest_50th",  # top 10–25%   — light green
        drug_potentiality >= cutoffs["25%"] ~ "05_neutral",      # middle 50%   — purple
        drug_potentiality >= cutoffs["10%"] ~ "06_desert_50th",  # bottom 10–25% — light yellow
        drug_potentiality >= cutoffs["5%"]  ~ "07_desert_10th",  # bottom 5–10% — orange
        drug_potentiality >= cutoffs["1%"]  ~ "08_desert_5th",   # bottom 1–5%  — dark orange
        TRUE                                ~ "09_desert_1st"    # bottom 1%    — darkest brown
      )
    )
  
  band_colors <- c(
    "01_forest_1st"  = "#004d00",
    "02_forest_5th"  = "#1a7a1a",
    "03_forest_10th" = "#4caf50",
    "04_forest_50th" = "#b7e4b7",
    "05_neutral"     = "#9b72b0",
    "06_desert_50th" = "#f5c542",
    "07_desert_10th" = "#e07b00",
    "08_desert_5th"  = "#b85000",
    "09_desert_1st"  = "#8B3a00"
  )
  
  band_order <- c(
    "01_forest_1st", "02_forest_5th", "03_forest_10th", "04_forest_50th",
    "05_neutral",
    "06_desert_50th", "07_desert_10th", "08_desert_5th", "09_desert_1st"
  )
  band_midpoints <- c(0.995, 0.97, 0.925, 0.50, 0.0, -0.50, -0.70, -0.92, -0.995)
  
  DP <- DP %>%
    dplyr::mutate(color_band_num = band_midpoints[match(color_band, band_order)])
  
  ## Order so labeled points plot on top
  DP_plot <- DP %>%
    dplyr::arrange(!is.na(label))
  
  p_dp <- ggplot(DP_plot, aes(x = magnitude, y = drug_potentiality)) +
    
    geom_point(aes(color = color_band_num), size = 1.2, alpha = 0.7) +
    
    scale_color_gradientn(
      colors = c("#8B3a00", "#b85000", "#e07b00", "#f5c542",
                 "#9b72b0",
                 "#b7e4b7", "#4caf50", "#1a7a1a", "#004d00"),
      values = scales::rescale(c(-1, -0.92, -0.70, -0.50, 0.0, 0.50, 0.925, 0.97, 1)),
      limits = c(-1, 1),
      # name   = paste0("Drug Potentiality\n(", method, " Quant. \u2212 PCA Quant.)"),
      name   = "Score",
      guide  = guide_colorbar(barheight = unit(6, "cm"), barwidth = unit(0.4, "cm"))
    ) +
    
    ## Reference lines — now anchored to the actual percentile boundaries
    geom_hline(yintercept = cutoffs["75%"], linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["90%"], linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["95%"], linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["99%"], linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["25%"], linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["10%"], linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["5%"],  linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = cutoffs["1%"],  linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = 0,              linetype = "solid",  color = "black",  linewidth = 0.5) +
    
    ## Gene labels
 geom_text_repel(
      data          = DP_plot %>% dplyr::filter(label_group == "top10"),
      aes(label     = label),
      size          = 2.8,
      color         = "#004d00",
      fontface      = "bold",
      box.padding   = 0.8, # was 0.5
      point.padding = 0.5, # 0.5
      max.overlaps  = Inf, # Inf
      segment.size  = 0.4,
      segment.color = "#004d00",
      nudge_y       = 0.02,
      nudge_x       = -0.01
    ) +

    geom_text_repel(
      data          = DP_plot %>% dplyr::filter(label_group == "bottom10"),
      aes(label     = label),
      size          = 2.8,
      color         = "#8B0000",
      fontface      = "bold",
      box.padding   = 0.8, # was 0.5
      point.padding = 0.5, # 0.5
      max.overlaps  = Inf, # Inf
      segment.size  = 0.4,
      segment.color = "#8B0000"
    ) +

    geom_label_repel(
      data          = DP_plot %>% dplyr::filter(label_group == "manual"),
      aes(label     = label),
      size          = 3.2,
      color         = "black",
      fontface      = "bold",
      fill          = "white",
      label.size    = 0.3,
      label.padding = unit(0.15, "lines"),
      label.r       = unit(0.05, "lines"),
      box.padding   = 0.4,
      point.padding = 0.2,
      max.overlaps  = Inf,
      segment.size  = 0.3,
      segment.color = "grey50",
      min.segment.length = 0,
      nudge_x       = 0.03,
      nudge_y       = -0.02,
      direction     = "y",
      hjust         = 0
    ) +
    
    # geom_label_repel(
    #   data          = DP_plot %>% dplyr::filter(label_group == "manual"),
    #   aes(label     = label),
    #   size          = 3.2,
    #   color         = "black",
    #   fontface      = "bold",
    #   fill          = "white",
    #   label.size    = 0.3,
    #   label.padding = unit(0.15, "lines"),
    #   label.r       = unit(0.05, "lines"),
    #   box.padding   = 0.4,
    #   point.padding = 0.2,
    #   max.overlaps  = 30,
    #   segment.size  = 0.3,
    #   segment.color = "grey50"
    # ) +
    
    ## Axis annotations — labels now reference the correct percentile cutoffs
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["99%"] + 0.03,
             label = "1st Forest",  size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["95%"] + 0.03,
             label = "5th Forest",  size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["90%"] + 0.03,
             label = "10th Forest", size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["75%"] + 0.03,
             label = "50th Forest", size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["25%"] - 0.06,
             label = "50th Desert", size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["10%"] - 0.06,
             label = "10th Desert", size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["5%"]  - 0.06,
             label = "5th Desert",  size = 2.5, color = "grey30", hjust = 1) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 0.85, y = cutoffs["1%"]  - 0.06,
             label = "1st Desert",  size = 2.5, color = "grey30", hjust = 1) +
    
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 1.0, y =  0.85,
             label = "Drug\nRainforest", size = 4, color = "#004d00",
             fontface = "italic", hjust = 0.5) +
    annotate("text", x = max(DP$magnitude, na.rm = TRUE) * 1.0, y = -0.85,
             label = "Drug\nDesert", size = 4, color = "#8B0000",
             fontface = "italic", hjust = 0.5) +
    
    labs(
      # x     = "Distance to origin\n(Euclidean magnitude of drug + PCA loadings)",
      x = "Distance to origin",
      y     = paste0("Drug Potentiality Score (", method, " Quantile \u2212 PCA Quantile)"),
      # y = "Drug Potentiality Score",
      # title = paste0("Drug Potentiality Score\n(", file_tag, Filtered_Tag, ")")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.12))) +
    theme_classic(base_size = 11) +
    theme(
      legend.position   = "right",
      legend.key.height = unit(1.5, "cm"),
      plot.title        = element_text(size = 10)
    )
  
  ggsave(
    filename = paste0(path.plots, "DrugPotentiality_", method, "_", output_tag, Filtered_Tag, ".pdf"),
    plot     = p_dp,
    width    = 6,
    height   = 7,
    units    = "in",
    device   = cairo_pdf
  )
  
  cat("Plot saved.\n")
}

### Step 9: Top/Bottom gene tables for inspection
if (1) {
  
  cat("\n--- Top 30 Drug RAINFOREST genes (highest potentiality) ---\n")
  print(DP %>% dplyr::arrange(desc(drug_potentiality)) %>%
          dplyr::select(Loading, abs_loading_drug, quantile_drug, abs_loading_PCA, quantile_PCA, drug_potentiality, magnitude) %>%
          head(500))
  
  cat("\n--- Top 30 Drug DESERT genes (lowest potentiality) ---\n")
  print(DP %>% dplyr::arrange(drug_potentiality) %>%
          dplyr::select(Loading, abs_loading_drug, quantile_drug, abs_loading_PCA, quantile_PCA, drug_potentiality, magnitude) %>%
          head(30))
  
  cat("\n--- Top 30 HIGH MAGNITUDE genes near zero potentiality (neutrally targeted) ---\n")
  print(DP %>%
          dplyr::filter(abs(drug_potentiality) < 0.05) %>%
          dplyr::arrange(desc(magnitude)) %>%
          dplyr::select(Loading, abs_loading_drug, quantile_drug, abs_loading_PCA, quantile_PCA, drug_potentiality, magnitude) %>%
          head(30))
}

### Step 10: Figure 4A-style — fitted distributions + comparison scatters (6-panel)
if (1) {
  
  ## Uses drug_max, PCA_max, DP, band_colors, method, path.plots already in memory
  
  fit_panels <- function(abs_vals, label, dist_color, cdf_color) {
    
    x_pos  <- abs_vals[abs_vals > 0 & !is.na(abs_vals)]
    fit    <- fitdistrplus::fitdist(x_pos, "lnorm")
    ml     <- fit$estimate["meanlog"]
    sl     <- fit$estimate["sdlog"]
    x_grid <- seq(0, max(x_pos) * 1.1, length.out = 500)
    
    df_dist <- data.frame(x = x_grid, y = dlnorm(x_grid, meanlog = ml, sdlog = sl))
    df_cdf  <- data.frame(x = x_grid, y = plnorm(x_grid, meanlog = ml, sdlog = sl))
    
    ex_low  <- quantile(x_pos, 0.01)
    ex_neut <- quantile(x_pos, 0.50)
    ex_high <- quantile(x_pos, 0.99)
    
    example_pts <- data.frame(
      x     = c(ex_low, ex_neut, ex_high),
      y_cdf = plnorm(c(ex_low, ex_neut, ex_high), meanlog = ml, sdlog = sl),
      lab   = c(paste0("Low ", label), paste0("Neutral ", label), paste0("High ", label)),
      col   = c("#004d00", "#9b72b0", "#004d00")
    )
    
    p_dist <- ggplot() +
      geom_histogram(
        data = data.frame(x = x_pos), aes(x = x, y = after_stat(density)),
        bins = 60, fill = scales::alpha(dist_color, 0.4),
        color = scales::alpha(dist_color, 0.6), linewidth = 0.2
      ) +
      geom_line(data = df_dist, aes(x = x, y = y), color = dist_color, linewidth = 0.8) +
      labs(x = "Maximal Absolute Loading", y = "Density",
           title = paste0("Distribution of\n", label, " Maximal Loading")) +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(face = "italic", size = 8, hjust = 0.5))
    
    p_cdf <- ggplot() +
      geom_line(data = df_cdf, aes(x = x, y = y), color = cdf_color, linewidth = 0.8) +
      geom_segment(data = example_pts,
                   aes(x = x, xend = x, y = 0, yend = y_cdf),
                   linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_segment(data = example_pts,
                   aes(x = 0, xend = x, y = y_cdf, yend = y_cdf),
                   linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_point(data = example_pts, aes(x = x, y = y_cdf), size = 3, color = example_pts$col) +
      geom_text(data = example_pts, aes(x = x, y = y_cdf, label = lab),
                hjust = -0.15, size = 2.8, color = example_pts$col, fontface = "italic") +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.25))) +
      labs(x = "Maximal Absolute Loading", y = "F(x)",
           title = "Fitted Log Normal\nCumulative Distribution") +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(face = "italic", size = 8, hjust = 0.5))
    
    list(p_dist = p_dist, p_cdf = p_cdf)
  }
  
  panels_drug <- fit_panels(drug_max$abs_loading_drug, label = method,  dist_color = "#2166ac", cdf_color = "#2166ac")
  panels_pca  <- fit_panels(PCA_max$abs_loading_PCA,   label = "PCA",   dist_color = "#d95f02", cdf_color = "#d95f02")
  
  ## Scatter plots
  p_raw <- ggplot(DP, aes(x = abs_loading_PCA, y = abs_loading_drug, color = color_band)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_manual(values = band_colors, guide = "none") +
    labs(x = "PCA Maximal Loading", y = paste0(method, " Maximal Loading")) +
    theme_classic(base_size = 9)
  
  p_quant <- ggplot(DP, aes(x = quantile_PCA, y = quantile_drug, color = color_band)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_manual(values = band_colors, guide = "none") +
    labs(x = "Transformed PCA Maximal Loading", y = paste0("Transformed ", method, " Maximal Loading")) +
    theme_classic(base_size = 9)
  
  library(patchwork)
  
  p_all <- (panels_drug$p_dist | panels_drug$p_cdf) /
    (panels_pca$p_dist  | panels_pca$p_cdf) /
    (p_raw               | p_quant) +
    patchwork::plot_annotation(
      # title = paste0("Drug Potentiality Score = (Drug + Gene Co-signal) - Gene Signal\n",
      #                method, " vs PCA"),
      title = "Signal Normalization",
      theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
    )
  
  print(p_all)
  
  ggsave(
    filename = paste0(path.plots, "FigureA_", method, "_vs_PCA_", output_tag, Filtered_Tag, ".pdf"),
    plot     = p_all,
    width    = 6.5, height = 9, units = "in", device = cairo_pdf
  )
  
  cat("6-panel figure saved")
}

## Step 11: Extra table

write.table(
  x = DP %>%
    dplyr::select(Loading, quantile_drug, quantile_PCA,
                  drug_potentiality, magnitude, color_band,
                  component_drug, component_PCA) %>%
    dplyr::arrange(desc(drug_potentiality)),
  file = paste0(path.dp, "DP_inspection_", output_tag, Filtered_Tag, ".txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

##### Rank-Rank and RRHO plots (NOT USEFUL SO FAR) #####
library(RRHO)
library(lattice)
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(stringr)

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
path.plots   <- paste0(path.wd, "Plots/")

## Set PLS parameters
X_source <- "CRISPR" # CRISPR or CTRP
Y_source <- "CTRP"   # CRISPR or CTRP

mode  <- "canonical" # default = regression, symmetric = canonical

## Derived label for plot titles
mode_label <- if (mode == "canonical") "PLS-C" else "PLS-R"

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- character(0)  # e.g. c("Myeloid", "Lymphoid") or character(0)

## Filtered for all three data sets shared lines?
FilteredAll3 <- TRUE # TRUE or FALSE

## Load in files
excl_tag <- if (length(exclude_lineages) > 0) {
  paste0("_excl.", paste(exclude_lineages, collapse = "."))
} else {
  ""
}
file_tag <- paste0("PLS_Mode.", mode, "_X.", X_source, "_Y.", Y_source, excl_tag)

if (FilteredAll3 == TRUE) {
  Filtered_Tag <- "_Filtered3"
} else {
  Filtered_Tag <- character(0)
}

path1 <- paste0(path.pls, file_tag, Filtered_Tag, "_X.variates.txt")
path2 <- paste0(path.pls, file_tag, Filtered_Tag, "_Y.variates.txt")

path1_context <- "CRISPR"
path2_context <- "Drug"

Rank_1 <- read.delim(path1, sep = "\t", header = T)

Rank_2 <- read.delim(path2, sep = "\t", header = T)

Rank_1_filt <- Rank_1 %>%
  dplyr::arrange(desc(comp4)) %>% #desc
  dplyr::mutate(CRISPR = c(1:nrow(.))) %>%
  dplyr::select(Score, CRISPR)

Rank_2_filt <- Rank_2 %>%
  dplyr::arrange(desc(comp4)) %>% #desc
  dplyr::mutate(Drug = c(1:nrow(.))) %>%
  dplyr::select(Score, Drug)

rank_merged <- merge(Rank_1_filt, Rank_2_filt, by = "Score")

# correlation and plot
correlation <- cor(rank_merged$CRISPR, rank_merged$Drug, method = "spearman")

# lm_model <- lm(a ~ b, data = rank_merged)
# slope <- coef(lm_model)[2] #very similiar as just doing the correlation

length <- dim(rank_merged)[1]

### Corank
pdf(file = paste0(path.plots, "Corank_Variates_PLSC_CRISPRDrug_Comp1.pdf"), width = 7.25, height = 6.5)
ggplot(rank_merged, aes(x = CRISPR, y = Drug)) +
  theme(axis.ticks = element_blank()) +
  stat_density_2d(
    aes(fill = ..density..),
    geom = "raster", 
    contour = FALSE
  ) +
  scale_fill_distiller(
    palette= "Spectral",
    name = "Density"
  ) +
  scale_x_continuous(
    breaks = c(length*.25, length*.75),
    labels = c("Positive Loading", "Negative Loading"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(length*.25, length*.75),
    labels = c("Positive Loading", "Negative Loading"),
    expand = c(0, 0)
  ) +
  theme(
    legend.position = "right", 
    panel.border = element_rect(colour = "black", fill=NA, size=.75), 
    text = element_text(size = 15),
    axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size = 15)
    # axis.title = element_text(size = 12)
  ) +
  xlab("CRISPR Sample Order") + ylab("Drug Sample Order") +
  geom_point(size = .25, alpha = 0.2) +
  labs(title = paste0("CRISPR x CTRP PLS-C : Component 1 : Cor = ", signif(correlation, 3))) 
# geom_smooth(method = 'lm', color = 'blue', se = FALSE, size = 1)
dev.off()

### RRHO Plot
Rank_1_sign <- Rank_1 %>%
  dplyr::filter(Rank_1$gene %in% Shared_Pathways) %>%
  dplyr::arrange(desc(sign_log_padj)) %>%
  dplyr::select(gene, sign_log_padj) %>%
  dplyr::rename(a = sign_log_padj)

Rank_2_sign <- Rank_2 %>%
  dplyr::filter(Rank_2$gene %in% Shared_Pathways) %>%
  dplyr::arrange(desc(sign_log_padj)) %>%
  dplyr::select(gene, sign_log_padj) %>%
  dplyr::rename(b = sign_log_padj)

RRHO <- RRHO(Rank_1_filt, Rank_2_filt, BY = TRUE, alternative = "enrichment")
max <- max(RRHO$hypermat)

RRHO_sign <- RRHO(Rank_1_sign, Rank_2_sign, BY = TRUE, alternative = "enrichment")
max_sign <- max(RRHO_sign$hypermat)

rainbow <- colorRampPalette(brewer.pal(11, "RdYlBu"))(10000)

title_2 <- paste0("Plots/RRHO_rank_DESeq_", gsub(" ", "", Label1), gsub(" ", "", Label2), "_", gsub(" ", "", path1_context), "_vs_", gsub(" ", "", Label3), gsub(" ", "", Label4),  "_", gsub(" ", "", path2_context), ".png")

# png(file = title_2, width = 600, height = 600)
png(file = title_2, width = 600, height = 600)
levelplot(RRHO$hypermat,
          xlab = "",
          ylab = "",
          col.regions = rev(rainbow),
          colorkey = list(labels = list(cex = 2)), # adjust font size on legend
          scales = list(tck = c(0, 0), 
                        x = list(draw = FALSE), 
                        y = list(draw = FALSE)),
          main = paste0("Max Signal: ", sprintf("%.2f", max)))
dev.off()

title_3 <- paste0("Plots/RRHO_signlogp_DESeq_", gsub(" ", "", Label1), gsub(" ", "", Label2), "_", gsub(" ", "", path1_context), "_vs_", gsub(" ", "", Label3), gsub(" ", "", Label4),  "_", gsub(" ", "", path2_context), ".png")

png(file = title_3, width = 600, height = 600)
levelplot(RRHO_sign$hypermat,
          xlab = "",
          ylab = "",
          col.regions = rev(rainbow),
          colorkey = list(labels = list(cex = 2)), # adjust font size on legend
          scales = list(tck = c(0, 0), 
                        x = list(draw = FALSE), 
                        y = list(draw = FALSE)),
          main = paste0("Max Signal: ", sprintf("%.2f", max_sign)))
dev.off()

### new RRHO 
library(RedRibbon)

Rank_sign_merge <- merge(Rank_1_sign, Rank_2_sign, by = "gene")

rr <- RedRibbon(Rank_sign_merge, enrichment_mode="hyper")
quad <- quadrants(rr, algorithm="ea", permutation=TRUE, whole=FALSE)
gg <- ggRedRibbon(rr, quadrants=quad,
                  repel.force = 250) +
  coord_fixed(ratio = 1, clip = "off")

title_4 <- paste0("Plots/RRHO_signlogp_RedRibbon_DESeq_", gsub(" ", "", Label1), gsub(" ", "", Label2), "_", gsub(" ", "", path1_context), "_vs_", gsub(" ", "", Label3), gsub(" ", "", Label4),  "_", gsub(" ", "", path2_context), ".pdf")

pdf(file = title_4, width = 7.25, height = 6.5)
gg
dev.off()

## Doesnt work
Rank_merge <- merge(Rank_1_filt, Rank_2_filt, by = "gene")

rr_1 <- RedRibbon(Rank_merge, enrichment_mode="hyper")
quad_1 <- quadrants(rr_1, algorithm="ea", permutation=TRUE, whole=FALSE)
gg_1 <- ggRedRibbon(rr_1, quadrants=quad_1,
                    repel.force = 250) +
  coord_fixed(ratio = 1, clip = "off")

title_5 <- paste0("RNA_DDIT3/plots/RRHO_rank_RedRibbon_DESeq_", Label1, Label2, "vs", Label3, Label4, ".pdf")
pdf(file = title_4, width = 7.25, height = 6.5)
gg_1
dev.off()



##### Export Tables #####
library(flextable)
library(officer)

path.wd      <- paste0(path.OS, "WD_FDB_Freeland/")
path.pls     <- paste0(path.wd, "DataSets/PLS/")
path.stat    <- paste0(path.wd, "DataSets/Stats/")


df <- read.table(
  file = paste0(path.pls, "PLS_Mode.canonical_X.CRISPR_Y.CTRP_Filtered3_canonical_correlations.txt"),
  sep = "\t",
  header = T
)

colnames(df)[2] <- "r\n\ncor(CRISPR x CTRP)"

ft <- flextable(df) |>
  set_header_labels(values = setNames(colnames(df), colnames(df))) |>
  compose(
    j = 2, part = "header",
    value = as_paragraph("r", as_sub("cor(CRISPR x CTRP)"))
  ) |>
  autofit() |>
  theme_booktabs() |>
  align(align = "center", part = "all") |>
  bold(part = "header")

save_as_image(ft, path = paste0(path.stat, "canonical_correlations_table.png"))
