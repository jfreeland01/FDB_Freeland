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
        labs(title = paste0(mode_label, " | ", source_label, " loadings: ", comp1, " vs ", comp2)) +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".loadings_", comp1, "vs", comp2, Special_string, ".pdf"
        ),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "Group", "Loading")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "Group", "Loading")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "Group", "Loading")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "Group", "Loading")
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
        scale_color_manual(values = lineage_colors, name = "Lineage") +
        labs(
          title = paste0(mode_label, " | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
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
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, paste0("X.", X_source))
  plot_scores_side(y.variates.plot, paste0("Y.", Y_source))
  
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
        title = paste0(mode_label, " | ", source_label, " scores by cancer lineage (comps 1–10)"),
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
          title = paste0(mode_label, " | ", source_label, " scores — ", comp, " by cancer lineage"),
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
  
  plot_scores_boxplot(x.variates.plot, paste0("X.", X_source))
  plot_scores_boxplot(y.variates.plot, paste0("Y.", Y_source))
  
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
    Filtered_Tag <- ""
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
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP")   Y_data <- CTRP
  
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
  
}

#### 2. Execute to run PLS and save output files (requires Step 1)
if (1) {
  
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
  
  ## Helper function for NA-safe pattern detection
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
          # stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "p53 Pathway\nModulators",
          TRUE ~ "Other"
        ),
        group.na     = dplyr::if_else(is.na(Group), 1L, 0L),
        label.not.na = dplyr::if_else(!is.na(Group), Loading, NA_character_),
        mix.flag     = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
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
  
  ## Plotting colors
  group_colors <- c(
    "BRAF Signalling"        = "#F8766D",
    "BRAF & MEK\nInhibitors" = "#F8766D",
    "EGFR Signalling"        = "#DE8C00",
    "EGFR & HER2\nInhibitors"= "#DE8C00",
    "Ferroptosis"            = "#B79F00",
    "Ferroptosis\nInducers"  = "#B79F00",
    "MED12"                  = "#00BA38",
    "Other"                  = "grey80"
    # "p53 Pathway\nModulators" = "#619CFF",
    # "DNA Damage\nResponse"    = "#619CFF"
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
        labs(title = paste0(mode_label, " | ", source_label, " loadings: ", comp1, " vs ", comp2)) +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, Filtered_Tag, "_", source_label, ".loadings_", comp1, "vs", comp2, Special_string, ".pdf"
        ),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_loadings_side(X_plot, paste0("X.", X_source), "Group", "Loading")
  plot_loadings_side(Y_plot, paste0("Y.", Y_source), "Group", "Loading")
  
}

#### 4. Execute to plot PLS scores colored by cancer type (requires Step 1 + saved variates)
if (1) {
  
  ## Load model metadata
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  ## Add path.stat if not already defined
  path.stat <- paste0(path.wd, "DataSets/Stats/")
  
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
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "comp1"
      comp2_col <- paste0("comp", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1_col,
          y     = comp2_col,
          color = "lineage_label"
        )
      ) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_colors, name = "Lineage") +
        labs(
          title = paste0(mode_label, " | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
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
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, paste0("X.", X_source))
  plot_scores_side(y.variates.plot, paste0("Y.", Y_source))
  
  ## Helper: boxplots of scores per lineage for comps 1-10
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^comp([1-9]|10)$", names(df), value = TRUE)
    
    ## Order lineages by comp1 median for multi-facet overview
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
        title = paste0(mode_label, " | ", source_label, " scores by cancer lineage (comps 1–10)"),
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
          title = paste0(mode_label, " | ", source_label, " scores — ", comp, " by cancer lineage"),
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
  
  plot_scores_boxplot(x.variates.plot, paste0("X.", X_source))
  plot_scores_boxplot(y.variates.plot, paste0("Y.", Y_source))
  
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
  # NEURO                  = "Neuro",
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
                               title = "CRISPR GO:BP Results (FDR < 0.05)") {
  
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

pdf(file = paste0(path.plots, "GO_AllModules_CRISPR_V1.pdf"), height = 6, width = 8)
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

pdf(file = paste0(path.plots, "GO_AllModules_RNAi_V1.pdf"), height = 6, width = 8)
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
    title = "CRISPR Common",
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
    title = "RNAi Common",
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
path.pls     <- paste0(path.wd, "DataSets/PLS/")

#### Get Hallmark EMT gene set from MSigDB
if (0) {

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
score_method <- "zscore" # zscore or ssgsea

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

#### Plot Figure

## Set Parameters
score_method <- "ssgsea"    # zscore or ssgsea
weight       <- "function"  # function (standard, WLS via weights argument, abs weights) or prior (direct multiply wx transformation, signed weights. not the standard)

# Lineage filtering — uses exact OncotreeLineage values from Model.csv (one or both must = NULL)
keep <- c("Biliary Tract", "Bladder", "Bowel", "Breast", "Cervix", "Esophagus", "Head and Neck", "Kidney", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Skin", "Thyroid", "Uterus", "Vulva", "Ampulla of Vater", "Pleura") # only retain these lineages  e.g. c("Skin") or NULL, or epithelial c("Biliary Tract", "Bladder/Urinary Tract", "Bowel", "Breast", "Cervix", "Esophagus/Stomach", "Head and Neck", "Kidney", "Liver", "Lung", "Ovary/Fallopian Tube", "Pancreas", "Prostate", "Skin", "Thyroid", "Uterus", "Vulva/Vagina", "Ampulla of Vater", "Pleura")
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
  
  weights_comp <- PLS_variates$comp3
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
    grepl("Brain",           df$OncotreeLineage, ignore.case = TRUE) ~ "Brain Cancer",
    grepl("Stomach", df$OncotreeLineage, ignore.case = TRUE) ~ "Gastric Cancer",
    grepl("Head|Neck",       df$OncotreeLineage, ignore.case = TRUE) ~ "Head and Neck Cancer",
    grepl("Kidney|Renal",    df$OncotreeLineage, ignore.case = TRUE) ~ "Kidney Cancer",
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
  if (1) {
    
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
    "Brain Cancer"         = "#87CEEB",
    "Gastric Cancer"       = "#00BFFF",
    "Head and Neck Cancer" = "#FFA500",
    "Kidney Cancer"        = "#90EE90",
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
      x = "Differentiation Score\n
      Differentiated \u2190  \u2192 Un-differentiated
      \n Epithelial \u2190  \u2192 Mesenchymal",
      y = "GPX4 CRISPR Dependency\nMore Dependent \u2190\u2192 Less Dependent"
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
    ggplot2::labs(x = "Differentiation Score\n
      Differentiated \u2190  \u2192 Un-differentiated
      \n Epithelial \u2190  \u2192 Mesenchymal",
                  y = "MED12 CRISPR Dependency") +
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
    ".pdf"
  )
  
  ggplot2::ggsave(
    filename = paste0(path.plots, save_name),
    plot = full_fig, width = 10, height = 8
  )
  
  cat("Saved:", save_name)
  
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
  top10    <- DP %>% dplyr::arrange(desc(drug_potentiality)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(Loading)
  bottom10 <- DP %>% dplyr::arrange(drug_potentiality)      %>% dplyr::slice_head(n = 10) %>% dplyr::pull(Loading)
  
  ## Manual curation list (from paper figure + your existing annotation groups)
  genes_manual <- c(
    ## Desert
    "KRAS", "PTEN", "SOD2", "MED12",
    ## Forest
    "VEGFA", "EPCAM", "ERG",
    
    # "NGF", "ITGA5", "ITGA6", "ITGAB4", "ERG", "SOX5", "ELF4", "KDR", "EFG", "TAP1", "ICAM2", "SLAMF1", "CDC42BPA", "NCAM1", "VEGFA", "GABARAP", "VRK2",
    ## Neutral
    "PIK3CA", "GPX4", "BRAF"
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
      name   = paste0("Drug Potentiality\n(", method, " Quant. \u2212 PCA Quant.)"),
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
 # geom_text_repel(
 #      data          = DP_plot %>% dplyr::filter(label_group == "top10"),
 #      aes(label     = label),
 #      size          = 2.8,
 #      color         = "#004d00",
 #      fontface      = "bold",
 #      box.padding   = 0.5,
 #      point.padding = 0.3,
 #      max.overlaps  = 30,
 #      segment.size  = 0.4,
 #      segment.color = "#004d00"
 #    ) +
    
    # geom_text_repel(
    #   data          = DP_plot %>% dplyr::filter(label_group == "bottom10"),
    #   aes(label     = label),
    #   size          = 2.8,
    #   color         = "#8B0000",
    #   fontface      = "bold",
    #   box.padding   = 0.5,
    #   point.padding = 0.3,
    #   max.overlaps  = 30,
    #   segment.size  = 0.4,
    #   segment.color = "#8B0000"
    # ) +
    # 
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
      nudge_x       = 0.02,
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
      x     = "Distance to origin\n(Euclidean magnitude of drug + PCA loadings)",
      y     = paste0(method, " Quantile \u2212 PCA Quantile\n(Drug Potentiality Score)"),
      title = paste0("Drug Potentiality Score\n(", file_tag, Filtered_Tag, ")")
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
    width    = 7,
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
      title = paste0("Drug Potentiality Score = (Drug + Gene Co-signal) - Gene Signal\n",
                     method, " vs PCA"),
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

##### SCRAP BELOW !!!!!!!!!!!!!!!!!! #####
##### (OLD) WGCNA: CRISPR (BACKUP BEFORE DELETED TOP 5 SECTIONS) #####

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
merge_CutHeight <- 0.20   # after modules are built, any two modules who correlate above this threshold get merged. Higher = more aggressive merging = fewer final modules. 0.25 = modules >75% similar get collapsed.

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
      "HEATMAP_WGCNA_CRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_AllClusters.png"
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

if (0) {
  
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
if (0) {
  
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
      path.plots,"HEATMAP_WGCNA_RNAi_OrderedByCRISPR_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_AllClusters.png"
    ),
    width  = 3000,
    height = 3000,
    res    = 200
  )
  ComplexHeatmap::draw(p_RNAi)
  grDevices::dev.off()
  
}

#### Look at top CRISPR modules now in RNAi
if (0) {
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
                        soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, ".rds"))
  
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
               file = paste0(path.wd, "DataSets/WGCNA/GO_Enrichment_CRISPR_AllModules_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, ".xlsx"),
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
          file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", ".rds"))
  
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
    file = paste0(path.wd, "DataSets/WGCNA/ModulePreservation_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", ".txt"),
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
  
  ggsave(paste0(path.plots, "ModulePreservation_Zsummary_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", ".png"), p_preservation, width = 8, height = 7)
  
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
  
  ggsave(paste0(path.plots, "ModulePreservation_MedianRank_CRISPR_in_RNAi_SoftPower_", soft_power, "_MinModuleSize_", min_module_sz, "_mergeCutHeight_", merge_CutHeight, "_deepSplit_", deep_Split, "_", ".png"), p_rank, width = 8, height = 7)
  
}

##### (OLD) rCCA: CRISPR & CTRP #####

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

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- c("Myeloid", "Lymphoid")  # e.g. c("Myeloid", "Lymphoid") or character(0)

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
  
  ## Choose lambdas: either tuned or manual
  if (mode_rcca == "ridge") {
    
    if (tune_lambda) {
      
      grid1      <- c(0.10, 0.20, 0.30)
      grid2      <- c(0.05, 0.10, 0.20)
      ncomp_tune <- min(5L, ncomp)
      
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
      
      lambda1 <- lambda1_manual
      lambda2 <- lambda2_manual
    }
    
    excl_tag <- if (length(exclude_lineages) > 0) {
      paste0("_excl.", paste(exclude_lineages, collapse = "."))
    } else {
      ""
    }
    file_tag <- paste0(
      "RCCA_ridge",
      "_lambda1.", format(lambda1, digits = 3),
      "_lambda2.", format(lambda2, digits = 3),
      "_X.", X_source, "_Y.", Y_source, excl_tag
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    excl_tag <- if (length(exclude_lineages) > 0) {
      paste0("_excl.", paste(exclude_lineages, collapse = "."))
    } else {
      ""
    }
    file_tag <- paste0(
      "RCCA_shrinkage",
      "_X.", X_source, "_Y.", Y_source, excl_tag
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
  
  cancor.df <- data.frame(
    comp                  = seq_along(rcca_fit$cor),
    canonical_correlation = rcca_fit$cor
  )
  
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

#### 3. Execute to plot rCCA loadings (requires Step 1)
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
  
  ## Helper function for NA-safe pattern detection
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for CTRP loadings file (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
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
        group.atp5    = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na      = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na  = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag      = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
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
    
    percent.nas <- as.data.frame(colMeans(is.na(CRISPR_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ## Annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_crispr(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_crispr(Y_loadings, "Y")
  
  ## Plotting colors
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot(
        df,
        aes_string(x = comp1, y = comp2, color = color_col)
      ) +
        geom_point(size = 2.5) +
        
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col),
          size = 2
        ) +
        
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        labs(title = paste0("rCCA | ", source_label, " loadings: ", comp1, " vs ", comp2)) +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
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

#### 4. Execute to plot rCCA scores colored by cancer type (requires Step 1 + saved variates)
if(1) {
  
  ## Load model metadata
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  ## Load saved variates files
  x.variates.plot <- read.delim(
    file = paste0(path.rcca, file_tag, "_X.variates.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  y.variates.plot <- read.delim(
    file = paste0(path.rcca, file_tag, "_Y.variates.txt"),
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
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "X1"
      comp2_col <- paste0("X", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1_col,
          y     = comp2_col,
          color = "lineage_label"
        )
      ) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_colors, name = "Lineage") +
        labs(
          title = paste0("rCCA | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
          x     = comp1_col,
          y     = comp2_col
        ) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".scores_",
          comp1_col, "vs", comp2_col, ".pdf"
        ),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, paste0("X.", X_source))
  plot_scores_side(y.variates.plot, paste0("Y.", Y_source))
  
  
  ## Helper: boxplots of scores per lineage for comps 1-10
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^X([1-9]|10)$", names(df), value = TRUE)
    
    df_long <- df %>%
      dplyr::select(Score, OncotreeLineage, dplyr::all_of(comp_cols)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(comp_cols),
        names_to  = "Component",
        values_to = "Score_value"
      ) %>%
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
    
    ## Multi-facet overview PDF
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
        title = paste0("rCCA | ", source_label, " scores by cancer lineage (comps 1–10)"),
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
        path.plots, "Plot_", file_tag, "_", source_label, ".scores_boxplot_comps1to10.pdf"
      ),
      plot = p, width = 12, height = 18, units = "in", device = cairo_pdf
    )
    
    ## Individual per-comp PDFs
    for (comp in comp_cols) {
      
      df_comp <- df_long %>% dplyr::filter(Component == comp)
      
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
          title = paste0("rCCA | ", source_label, " scores — ", comp, " by cancer lineage"),
          x     = NULL,
          y     = "Score"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".scores_boxplot_", comp, ".pdf"
        ),
        plot = p_ind, width = 10, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_boxplot(x.variates.plot, paste0("X.", X_source))
  plot_scores_boxplot(y.variates.plot, paste0("Y.", Y_source))
  
}

##### (OLD) rCCA: RNAi & CTRP #####

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
mode_rcca      <- "shrinkage" # ridge (default) requires parameters or tuning, shrinkage

tune_lambda    <- FALSE   # set TRUE to run automatic tuning
lambda1_manual <- 0.20    # penalty on X (RNAi side if X_source == "RNAi")
lambda2_manual <- 0.10    # penalty on Y

## Cell lines to exclude by OncotreeLineage (set to character(0) to skip filtering)
exclude_lineages <- character(0)  # e.g. c("Myeloid", "Lymphoid") or character(0)

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
    file = paste0(path.dm, "Model.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::select(ModelID, CCLEName, OncotreeLineage)
  
  RNAi_t <- RNAi %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName")
  
  RNAi_t_ModelID <- merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName, -OncotreeLineage) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## Filter for shared cell lines, make matrix (for mixOmics), ensure numeric
  if (X_source == "RNAi") X_data <- RNAi_t_ModelID
  if (X_source == "CTRP") X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP") Y_data <- CTRP
  
  ids <- intersect(rownames(X_data), rownames(Y_data))
  
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
  
  if (mode_rcca == "ridge") {
    
    if (tune_lambda) {
      
      grid1      <- c(0.10, 0.20, 0.30)
      grid2      <- c(0.05, 0.10, 0.20)
      ncomp_tune <- min(5L, ncomp)
      
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
      
      lambda1 <- lambda1_manual
      lambda2 <- lambda2_manual
    }
    
    excl_tag <- if (length(exclude_lineages) > 0) {
      paste0("_excl.", paste(exclude_lineages, collapse = "."))
    } else {
      ""
    }
    file_tag <- paste0(
      "RCCA_ridge",
      "_lambda1.", format(lambda1, digits = 3),
      "_lambda2.", format(lambda2, digits = 3),
      "_X.", X_source, "_Y.", Y_source, excl_tag
    )
    
  } else if (mode_rcca == "shrinkage") {
    
    excl_tag <- if (length(exclude_lineages) > 0) {
      paste0("_excl.", paste(exclude_lineages, collapse = "."))
    } else {
      ""
    }
    file_tag <- paste0(
      "RCCA_shrinkage",
      "_X.", X_source, "_Y.", Y_source, excl_tag
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
  
  cancor.df <- data.frame(
    comp                  = seq_along(rcca_fit$cor),
    canonical_correlation = rcca_fit$cor
  )
  
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

#### 3. Execute to plot rCCA loadings (requires Step 1)
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
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
  }
  
  ## Helper function for NA-safe pattern detection
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for CTRP loadings file (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    
    ctrp.inform <- read.delim(
      file = paste0(path.ctrp, "CTRPv2.0._INFORMER_SET.txt"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    lk <- match(df$Loading, ctrp.inform$cpd_name)
    df$drug.target <- ctrp.inform$target_or_activity_of_compound[lk]
    
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
        group.atp5    = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na      = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na  = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag      = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
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
    
    percent.nas <- as.data.frame(colMeans(is.na(RNAi_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ## Annotate X- and Y- loadings based on actual sources
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_rnai(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_rnai(Y_loadings, "Y")
  
  ## Plotting colors
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot(
        df,
        aes_string(x = comp1, y = comp2, color = color_col)
      ) +
        geom_point(size = 2.5) +
        
        geom_text_repel(
          data = df %>% dplyr::filter(!is.na(.data[[color_col]])),
          aes_string(label = label_col),
          size = 2
        ) +
        
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        scale_color_manual(values = my_colors, na.value = "grey80") +
        labs(title = paste0("rCCA | ", source_label, " loadings: ", comp1, " vs ", comp2)) +
        theme_bw(base_size = 10)
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
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

#### 4. Execute to plot rCCA scores colored by cancer type (requires Step 1 + saved variates)
if(1) {
  
  ## Load model metadata
  model <- read.csv(paste0(path.dm, "Model.csv"))
  
  ## Load saved variates files
  x.variates.plot <- read.delim(
    file = paste0(path.rcca, file_tag, "_X.variates.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  y.variates.plot <- read.delim(
    file = paste0(path.rcca, file_tag, "_Y.variates.txt"),
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
    
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1_col <- "X1"
      comp2_col <- paste0("X", i)
      
      p <- ggplot(
        df,
        aes_string(
          x     = comp1_col,
          y     = comp2_col,
          color = "lineage_label"
        )
      ) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
        scale_color_manual(values = lineage_colors, name = "Lineage") +
        labs(
          title = paste0("rCCA | ", source_label, " scores: ", comp1_col, " vs ", comp2_col),
          x     = comp1_col,
          y     = comp2_col
        ) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".scores_",
          comp1_col, "vs", comp2_col, ".pdf"
        ),
        plot = p, width = 7, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_side(x.variates.plot, paste0("X.", X_source))
  plot_scores_side(y.variates.plot, paste0("Y.", Y_source))
  
  
  ## Helper: boxplots of scores per lineage for comps 1-10
  plot_scores_boxplot <- function(df, source_label) {
    
    comp_cols <- grep("^X([1-9]|10)$", names(df), value = TRUE)
    
    df_long <- df %>%
      dplyr::select(Score, OncotreeLineage, dplyr::all_of(comp_cols)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(comp_cols),
        names_to  = "Component",
        values_to = "Score_value"
      ) %>%
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
    
    ## Multi-facet overview PDF
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
        title = paste0("rCCA | ", source_label, " scores by cancer lineage (comps 1–10)"),
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
        path.plots, "Plot_", file_tag, "_", source_label, ".scores_boxplot_comps1to10.pdf"
      ),
      plot = p, width = 12, height = 18, units = "in", device = cairo_pdf
    )
    
    ## Individual per-comp PDFs
    for (comp in comp_cols) {
      
      df_comp <- df_long %>% dplyr::filter(Component == comp)
      
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
          title = paste0("rCCA | ", source_label, " scores — ", comp, " by cancer lineage"),
          x     = NULL,
          y     = "Score"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      
      ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".scores_boxplot_", comp, ".pdf"
        ),
        plot = p_ind, width = 10, height = 5, units = "in", device = cairo_pdf
      )
    }
  }
  
  plot_scores_boxplot(x.variates.plot, paste0("X.", X_source))
  plot_scores_boxplot(y.variates.plot, paste0("Y.", Y_source))
  
}

##### (OLD) Max loading: Scatter & GSEA #####

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
exclude_lineages_1 <- character(0) # for file1 (X1/Y1) c("Myeloid", "Lymphoid")
exclude_lineages_2 <- character(0)  # for file2 (X2/Y2) c("Myeloid", "Lymphoid")

### Create scatter plot and generate table of distances and theta
if(1){
  
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
  
  if (DimRedTec == "PLS") {
    file1_tag <- paste0("PLS_Mode.", mode, "_X.", X1_source, "_Y.", Y1_source, excl_tag_1)
    file2_tag <- paste0("PLS_Mode.", mode, "_X.", X2_source, "_Y.", Y2_source, excl_tag_2)
  }
  
  if (DimRedTec == "rCCA") {
    file1_tag <- paste0("rCCA_X.", X1_source, "_Y.", Y1_source, excl_tag_1)
    file2_tag <- paste0("rCCA_X.", X2_source, "_Y.", Y2_source, excl_tag_2)
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
  
  ## Build output tag for file/plot names
  out_tag <- paste0(
    mode,
    "_X1_", X1_source, "_vs_", Y1_source, excl_tag_1,
    "_X2_", X2_source, "_vs_", Y2_source, excl_tag_2
  )
  
  write.table(
    x    = Max,
    file = paste0(path.max, "MaxLoadingsDF_", out_tag, ".txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )
  
  ## plot
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
      data        = top5,
      aes(label   = gene),
      size        = 2.5,
      fontface    = "italic",
      show.legend = FALSE,
      max.overlaps = 20,
      segment.size = 0.3,
      segment.color = "grey50"
    ) +
    scale_color_manual(
      # values = c("RNAi" = "#7a1515", "Neutral" = "#BDBDBD", "CRISPR" = "#1a3a6b"),
      values = c("RNAi" = "#5E2F80", "Neutral" = "#BDBDBD", "CRISPR" = "#D47D37"),
      name   = "Bias"
    ) +
    labs(
      x     = "Maximal PLS-C RNAi Loading", # names(Max)[7]
      y     = "Maximal PLS-C CRISPR Loading", # names(Max)[4]
      # title = "Max Absolute Loadings per Gene comp 1-10"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic(base_size = 10) +
    theme(
      # legend.position   = c(0.18, 0.82),
      legend.background = element_blank(),
      legend.key        = element_blank()
      # legend.title      = element_text(size = 8),
      # legend.text       = element_text(size = 7) 
    )
  
  ggsave(
    filename = paste0(path.plots, "MaxLoadingsDF_", out_tag, "_Scatter.pdf"),
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
  filename = paste0(path.plots, "MaxLoadingsDF_", out_tag, "_GSEA_V2.pdf"),
  plot = p, width = 13, height = 9, units = "in", device = cairo_pdf
)
##### (OLD) WGCNA: RNAi #####

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
if (0) {
  
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
