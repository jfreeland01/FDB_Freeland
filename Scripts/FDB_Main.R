##### Set up #####
if(1) {
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
}

#### Imputation: CRISPR (NA's in Data) ####

## Set paths
path.wd <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
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

#### Imputation: RNAi (NA's in Data) ####

## Set paths
path.wd <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
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
  file = gsub(".csv$","_MFImputed.txt",file.rnai), 
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

## Set paths
path.wd <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.dm   <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp <- paste0(path.wd, "DataSets/CTRPv2/")

## Load in cell line info from depamp
models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F)

## lLoad in CTRP Data
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
  file = paste0(path.ctrp,"ctrp.curves.txt"),
  quote = F, col.names=T, row.names = F, sep = "\t")

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
  direction = "wide")

names(ctrpv2.ave.wide) <- gsub("^avg\\.","",names(ctrpv2.ave.wide))

file.ctrpv2.wide = paste0(path.ctrp,"ctrpv2.wide.txt")

write.table(
  x = ctrpv2.ave.wide, 
  file = file.ctrpv2.wide,
  quote = F, col.names=T, row.names = F, sep = "\t")

ctrpv2 <- ctrpv2.ave.wide

## Create a culled version
ctrpv3 <-  ctrpv2

## Rename one entry with no DepMap_ID
ctrpv3$DepMap_ID <-  ifelse(is.na(ctrpv3$DepMap_ID),"no.depmap.match",ctrpv3$DepMap_ID)
row.names(ctrpv3) <- ctrpv3$DepMap_ID
ctrpv3 <-  ctrpv3[,-1]

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

## Set paths
path.wd   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.dm   <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp <- paste0(path.wd, "DataSets/CTRPv2/")
path.pls  <- paste0(path.wd, "DataSets/PLS/")
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
  CRISPR_mat <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub(" .*", "", .))
  
  CTRP_mat <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE #, row.names = 1
  )
  
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
    comp_cols <- grep("^X\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "X1"
      comp2 <- paste0("X", i)
      
      p <- ggplot2::ggplot(
        df, ggplot2::aes_string(x = comp1, y = comp2, color = color_col, label = label_col)
      ) +
        ggplot2::geom_point(size = 2.5) +
        ggrepel::geom_text_repel(size = 2) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::scale_color_manual(values = my_colors, na.value = "grey80") +
        ggplot2::theme_bw(base_size = 10)
      
      ggplot2::ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## If side is CTRP: color by target.category, label by target.category
  ## If side is CRISPR: color by group, label by group
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "target.category", "target.category")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "group")
  }
    
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "target.category", "target.category")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "group")
  }
}

##### PLS: RNAi & CTRP #####

## Set paths
path.wd   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.dm   <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp <- paste0(path.wd, "DataSets/CTRPv2/")
path.pls  <- paste0(path.wd, "DataSets/PLS/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Set PLS parameters
X_source <- "RNAi"
Y_source <- "CTRP"

ncomp <- 15
mode  <- "canonical" # default = regression, symmetric = canonical

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
    tibble::rownames_to_column(var = "CCLEName")
  
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
  RNAi_mat <- read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub(" .*", "", .))
  
  CTRP_mat <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE #, row.names = 1
  )
  
  ## Helper function for NA-safe pattern detection (useful when labeling for plotting)
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for a CTRP loadings frame (drug metadata buckets)
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
    
    ## %NA per compound using non-imputed CTRP matrix
    percent.nas <- as.data.frame(colMeans(is.na(RNAi_mat)) * 100)
    names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ## annotate X- and Y- loadings based on actual sources
  debugonce(annotate_rnai)
  debugonce(annotate_ctrp)
  
  X_plot <- if (X_source == "CTRP") annotate_ctrp(X_loadings, "X") else annotate_rnai(X_loadings, "X")
  Y_plot <- if (Y_source == "CTRP") annotate_ctrp(Y_loadings, "Y") else annotate_rnai(Y_loadings, "Y")
  
  ## Plotting colors (always plot both sides)
  my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D",
                 "#00BFC4","#00B4F0","#619CFF","hotpink","purple","cyan")
  
  plot_loadings_side <- function(df, source_label, color_col, label_col) {
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot2::ggplot(
        df, ggplot2::aes_string(x = comp1, y = comp2, color = color_col, label = label_col)
      ) +
        ggplot2::geom_point(size = 2.5) +
        ggrepel::geom_text_repel(size = 2) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::scale_color_manual(values = my_colors, na.value = "grey80") +
        ggplot2::theme_bw(base_size = 10)
      
      ggplot2::ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## - If side is CTRP: color by target.category, label by target.category
  ## - If side is CRISPR: color by group, label by group
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "target.category", "target.category")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "group")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "target.category", "target.category")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "group")
  }
}

##### rCCA: CRISPR & CTRP #####

## Set paths
path.wd      <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.dm      <- paste0(path.wd, "DataSets/DepMap_25Q3/")
path.ctrp    <- paste0(path.wd, "DataSets/CTRPv2/")
path.rcca     <- paste0(path.wd, "DataSets/rCCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.general <- paste0(path.wd, "DataSets/General/")

## Set RCCA parameters
X_source <- "CTRP"    # "CRISPR" or "CTRP"
Y_source <- "CRISPR"  # "CRISPR" or "CTRP"

ncomp <- 15

## Regularization controls
tune_lambda    <- FALSE   # set TRUE to run automatic tuning
lambda1_manual <- 0.20    # default penalty on X (CRISPR side if X_source == "CRISPR")
lambda2_manual <- 0.10    # default penalty on Y

#### 1. Execute to prep for RCCA
if (1) {
  
  ## Read in data and set row names
  CRISPR <- utils::read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub("\\.\\..*", "", .))
  
  CTRP <- utils::read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  ## Filter for shared cell lines, make matrix (for mixOmics), ensure numeric
  if (X_source == "CRISPR") X_data <- CRISPR
  if (X_source == "CTRP")   X_data <- CTRP
  
  if (Y_source == "CRISPR") Y_data <- CRISPR
  if (Y_source == "CTRP")   Y_data <- CTRP
  
  ids <- base::intersect(base::rownames(X_data), base::rownames(Y_data))
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- base::lapply(X, base::as.numeric)
  Y[] <- base::lapply(Y, base::as.numeric)
  
  X <- base::as.matrix(X)
  Y <- base::as.matrix(Y)
}

#### 2. Execute to run RCCA and save output files (requires Step 1)
if (1) {
  
  ## Choose lambdas: either tuned or manual
  if (tune_lambda) {
    
    grid1 <- c(0.10, 0.20, 0.30)  # candidate lambdas for X
    grid2 <- c(0.05, 0.10, 0.20)  # candidate lambdas for Y
    ncomp_tune <- base::min(5L, ncomp)  # tune only first few components
    
    base::set.seed(1L)
    tune_time <- base::system.time({
      tune.out <- mixOmics::tune.rcc(
        X          = X,
        Y          = Y,
        grid1      = grid1,
        grid2      = grid2,
        ncomp      = ncomp_tune,
        validation = "loo"
      )
    })
    
    base::print(tune_time)                # show how long tuning took
    base::print(tune.out$opt.lambda1)     # chosen lambda1
    base::print(tune.out$opt.lambda2)     # chosen lambda2
    
    lambda1 <- tune.out$opt.lambda1
    lambda2 <- tune.out$opt.lambda2
    
  } else {
    ## Use manual defaults
    lambda1 <- lambda1_manual
    lambda2 <- lambda2_manual
  }
  
  ## For saving files
  file_tag <- paste0(
    "RCCA_lambda1.", base::format(lambda1, digits = 3),
    "_lambda2.", base::format(lambda2, digits = 3),
    "_X.", X_source, "_Y.", Y_source
  )
  
  ## Run RCCA
  rcca_fit <- mixOmics::rcc(
    X       = X,
    Y       = Y,
    ncomp   = ncomp,
    lambda1 = lambda1,
    lambda2 = lambda2
  )
  
  ## Canonical correlations per component
  print(rcca_fit$cor)
  
  ## Extract from rcca_fit object
  x.variates <- data.frame(rcca_fit$variates$X) %>%
    tibble::rownames_to_column(var = "Score")
  y.variates <- data.frame(rcca_fit$variates$Y) %>%
    tibble::rownames_to_column(var = "Score")
  
  x.loadings <- base::data.frame(rcca_fit$loadings$X) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(X1)
  y.loadings <- base::data.frame(rcca_fit$loadings$Y) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(X1)
  
  variates.X.Y <- base::merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = c(paste0(".", X_source), paste0(".", Y_source))
  )
  
  ## Canonical correlations data.frame
  cancor.df <- base::data.frame(
    comp                   = base::seq_along(rcca_fit$cor),
    canonical_correlation  = rcca_fit$cor
  )
  
  ## Save files
  utils::write.table(
    x = x.variates,
    file = paste0(path.rcca, file_tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = y.variates,
    file = paste0(path.rcca, file_tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = variates.X.Y,
    file = paste0(path.rcca, file_tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = x.loadings,
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = y.loadings,
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = cancor.df,
    file = paste0(path.rcca, file_tag, "_canonical_correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

#### 3. Execute to plot RCCA (requires Step 2)
if(1) {
  
  ## Load saved loading files
  X_loadings <- utils::read.delim(
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- utils::read.delim(
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Bring in raw matrices to compute %NA later
  CRISPR_mat <- read.delim(
    file = paste0(path.dm, "CRISPRGeneEffect.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub(" .*", "", .))
  
  CTRP_mat <- read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE #, row.names = 1
  )
  
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
    comp_cols <- grep("^comp\\d+$", names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot2::ggplot(
        df, ggplot2::aes_string(x = comp1, y = comp2, color = color_col, label = label_col)
      ) +
        ggplot2::geom_point(size = 2.5) +
        ggrepel::geom_text_repel(size = 2) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::scale_color_manual(values = my_colors, na.value = "grey80") +
        ggplot2::theme_bw(base_size = 10)
      
      ggplot2::ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## If side is CTRP: color by target.category, label by target.category
  ## If side is CRISPR: color by group, label by group
  
  debugonce(plot_loadings_side)
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "target.category", "target.category")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "group")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "target.category", "target.category")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "group")
  }
}

##### rCCA: RNAi & CTRP #####

## Set paths
path.wd      <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
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
tune_lambda    <- FALSE   # set TRUE to run automatic tuning
lambda1_manual <- 0.20    # penalty on X (RNAi side if X_source == "RNAi")
lambda2_manual <- 0.10    # penalty on Y

#### 1. Execute to prep for RCCA
if (1) {
  
  ## Read in data
  RNAi <- utils::read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  CTRP <- utils::read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  )
  
  ## Convert RNAi sample nomenclature (CCLEName -> ModelID via Model.csv)
  models <- utils::read.delim(
    file = paste0(path.dm,"Model.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE
  ) %>%
    dplyr::select(ModelID, CCLEName)
  
  RNAi_t <- RNAi %>%
    t() %>%
    base::data.frame() %>%
    tibble::rownames_to_column(var = "CCLEName")
  
  RNAi_t_ModelID <- base::merge(models, RNAi_t, by = "CCLEName") %>%
    dplyr::select(-CCLEName) %>%
    tibble::column_to_rownames(var = "ModelID")
  
  ## Filter for shared cell lines, make matrix (for mixOmics), ensure numeric
  if (X_source == "RNAi") X_data <- RNAi_t_ModelID
  if (X_source == "CTRP") X_data <- CTRP
  
  if (Y_source == "RNAi") Y_data <- RNAi_t_ModelID
  if (Y_source == "CTRP") Y_data <- CTRP
  
  ids <- base::intersect(base::rownames(X_data), base::rownames(Y_data))
  
  X <- X_data[ids, , drop = FALSE]
  Y <- Y_data[ids, , drop = FALSE]
  
  X[] <- base::lapply(X, base::as.numeric)
  Y[] <- base::lapply(Y, base::as.numeric)
  
  X <- base::as.matrix(X)
  Y <- base::as.matrix(Y)
  
  ## Enforce reasonable number of canonical components
  max_possible <- base::min(
    base::nrow(X) - 1L,
    base::ncol(X),
    base::ncol(Y)
  )
  ncomp <- base::min(ncomp, max_possible)
  base::cat("Using", ncomp, "canonical components (max possible =", max_possible, ")\n")
}

#### 2. Execute to run RCCA and save output files (requires Step 1)
if (1) {
  
  ## Choose lambdas: either tuned or manual
  if (tune_lambda) {
    
    grid1 <- c(0.10, 0.20, 0.30)  # candidate lambdas for X (RNAi)
    grid2 <- c(0.05, 0.10, 0.20)  # candidate lambdas for Y (CTRP)
    ncomp_tune <- base::min(5L, ncomp)  # tune only first few components
    
    base::set.seed(1L)
    tune_time <- base::system.time({
      tune.out <- mixOmics::tune.rcc(
        X          = X,
        Y          = Y,
        grid1      = grid1,
        grid2      = grid2,
        ncomp      = ncomp_tune,
        validation = "loo"
      )
    })
    
    base::print(tune_time)                # show how long tuning took
    base::print(tune.out$opt.lambda1)     # chosen lambda1
    base::print(tune.out$opt.lambda2)     # chosen lambda2
    
    lambda1 <- tune.out$opt.lambda1
    lambda2 <- tune.out$opt.lambda2
    
  } else {
    ## Use manual defaults
    lambda1 <- lambda1_manual
    lambda2 <- lambda2_manual
  }
  
  ## For saving files
  file_tag <- paste0(
    "RCCA_lambda1.", base::format(lambda1, digits = 3),
    "_lambda2.", base::format(lambda2, digits = 3),
    "_X.", X_source, "_Y.", Y_source
  )
  
  ## Run RCCA
  rcca_fit <- mixOmics::rcc(
    X       = X,
    Y       = Y,
    ncomp   = ncomp,
    lambda1 = lambda1,
    lambda2 = lambda2
  )
  
  ## Canonical correlations per component (full spectrum; first ncomp used)
  base::print(rcca_fit$cor[1:ncomp])
  
  ## Extract from rcca_fit object
  x.variates <- base::data.frame(rcca_fit$variates$X) %>%
    tibble::rownames_to_column(var = "Score")
  y.variates <- base::data.frame(rcca_fit$variates$Y) %>%
    tibble::rownames_to_column(var = "Score")
  
  x.loadings <- base::data.frame(rcca_fit$loadings$X) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(comp1)
  y.loadings <- base::data.frame(rcca_fit$loadings$Y) %>%
    tibble::rownames_to_column(var = "Loading") %>%
    dplyr::arrange(comp1)
  
  variates.X.Y <- base::merge(
    x = x.variates, y = y.variates, by = "Score",
    suffixes = c(paste0(".", X_source), paste0(".", Y_source))
  )
  
  ## Canonical correlations data.frame
  cancor.df <- base::data.frame(
    comp                  = base::seq_along(rcca_fit$cor),
    canonical_correlation = rcca_fit$cor
  )
  
  ## Save files
  if (!dir.exists(path.rcca)) dir.create(path.rcca, recursive = TRUE)
  
  utils::write.table(
    x = x.variates,
    file = paste0(path.rcca, file_tag, "_X.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = y.variates,
    file = paste0(path.rcca, file_tag, "_Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = variates.X.Y,
    file = paste0(path.rcca, file_tag, "_X.Y.variates.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = x.loadings,
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = y.loadings,
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  utils::write.table(
    x = cancor.df,
    file = paste0(path.rcca, file_tag, "_canonical_correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

#### 3. Execute to plot RCCA (requires Step 2)
if (1) {
  
  ## Load saved loading files
  X_loadings <- utils::read.delim(
    file = paste0(path.rcca, file_tag, "_X.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  Y_loadings <- utils::read.delim(
    file = paste0(path.rcca, file_tag, "_Y.loadings.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Bring in raw matrices to compute %NA later
  RNAi_mat <- utils::read.delim(
    file = paste0(path.dm, "D2_combined_gene_dep_scores.csv"),
    sep = ",", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1
  ) %>%
    dplyr::rename_with(~ sub(" .*", "", .))
  
  CTRP_mat <- utils::read.delim(
    file = paste0(path.ctrp, "ctrpv2.wide.txt"),
    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  ## Helper function for NA-safe pattern detection (useful when labeling for plotting)
  detect <- function(x, pattern) {
    stringr::str_detect(ifelse(is.na(x), "", x), stringr::regex(pattern, ignore_case = TRUE))
  }
  
  ### Annotation for a CTRP loadings frame (drug metadata buckets)
  annotate_ctrp <- function(df, side_label) {
    ctrp.inform <- utils::read.delim(
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
        group.atp5       = dplyr::if_else(stringr::str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
        group.na         = dplyr::if_else(is.na(group), 1L, 0L),
        group.atp5.na    = dplyr::if_else(is.na(group.atp5), 1L, 0L),
        label.not.na     = dplyr::if_else(!is.na(group), Loading, NA_character_),
        label.not.na.atp5 = dplyr::if_else(!is.na(group.atp5), Loading, NA_character_),
        mix.flag         = dplyr::if_else(stringr::str_detect(Loading, ":"), "dual drug", "single drug")
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
    percent.nas <- base::as.data.frame(base::colMeans(base::is.na(CTRP_mat)) * 100)
    base::names(percent.nas) <- "percent.nas"
    percent.nas <- tibble::rownames_to_column(percent.nas, var = "Loading")
    df <- dplyr::left_join(df, percent.nas, by = "Loading")
    
    df
  }
  
  ### Annotation for a RNAi loadings frame (gene metadata buckets)
  annotate_rnai <- function(df, side_label) {
    
    gene.info.all <- utils::read.delim(
      file = paste0(path.general, "Homo_sapiens.gene_info.20251028"),
      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
    )
    
    gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority != "-", ]
    gene.info.abr <- dplyr::select(gene.info, Symbol, description)
    
    df$Loading <- base::sub("\\.\\..*$", "", df$Loading)
    
    df <- base::merge(df, gene.info.abr, by.x = "Loading", by.y = "Symbol", all.x = TRUE)
    
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
    percent.nas <- base::as.data.frame(base::colMeans(base::is.na(RNAi_mat)) * 100)
    base::names(percent.nas) <- "percent.nas"
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
    comp_cols <- base::grep("^comp\\d+$", base::names(df), value = TRUE)
    if (length(comp_cols) < 2) return(invisible(NULL))
    
    for (i in 2:length(comp_cols)) {
      
      comp1 <- "comp1"
      comp2 <- paste0("comp", i)
      
      p <- ggplot2::ggplot(
        df, ggplot2::aes_string(x = comp1, y = comp2, color = color_col, label = label_col)
      ) +
        ggplot2::geom_point(size = 2.5) +
        ggrepel::geom_text_repel(size = 2) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
        ggplot2::scale_color_manual(values = my_colors, na.value = "grey80") +
        ggplot2::theme_bw(base_size = 10)
      
      ggplot2::ggsave(
        filename = paste0(
          path.plots, "Plot_", file_tag, "_", source_label, ".loadings_", comp1, "vs", comp2, ".pdf"
        ),
        plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
      )
    }
  }
  
  ## If side is CTRP: color by target.category, label by target.category
  ## If side is RNAi: color by group, label by group
  
  if (X_source == "CTRP") {
    plot_loadings_side(X_plot, paste0("X.", X_source), "target.category", "target.category")
  } else {
    plot_loadings_side(X_plot, paste0("X.", X_source), "group", "group")
  }
  
  if (Y_source == "CTRP") {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "target.category", "target.category")
  } else {
    plot_loadings_side(Y_plot, paste0("Y.", Y_source), "group", "group")
  }
}


##### Max loading scatter plot #####

## Set Paths (same as above section)
path.wd   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.pls  <- paste0(path.wd, "DataSets/PLS/")
path.rcca <- paste0(path.wd, "DataSets/rCCA/")
path.plots   <- paste0(path.wd, "Plots/")
path.max <- paste0(path.wd, "DataSets/MaxLoading/")

## set dim red technique. PLS, rCCA
DimRedTec <- "PLS"

## set parameters. CRISPR, RNAi, CTRP
X1_source <- "CRISPR"
Y1_source <- "CTRP"

X2_source <- "RNAi"
Y2_source <- "CTRP"

mode  <- "canonical" # default = regression, symmetric = canonical

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
  
  xcol <- names(Max)[4]   # 4th column name
  ycol <- names(Max)[7]   # 7th column name
  
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
    x      = Max[[4]],   # 4th column = abs_loading_CRISPR
    y      = Max[[7]],   # 7th column = abs_loading_RNAi
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
    theme_bw(base_size = 14) +
    labs(
      x = names(Max)[4],
      y = names(Max)[7],
      title = "Max Absolute Loadings per Gene comp 1-10"
    ) +
    scale_x_continuous(expand = expansion(mult = 0), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, NA))
  
  ggplot2::ggsave(
    filename = paste0(path.plots, "MaxLoadingsDF_", mode, "_X1_", X1_source, "_vs_", Y1_source, "_X2_", X2_source, "_vs_", Y2_source, ".pdf"),
    plot = p, width = 6, height = 4, units = "in", device = cairo_pdf
  )
  
}

##### GSEA #####

### Set paths and read functions
path.wd   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.max <- paste0(path.wd, "DataSets/MaxLoading/")
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

## Set Paths
path.wd    <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
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
plt <- ggplot2::ggplot(data_df, ggplot2::aes(x = Value, y = Category)) +
  ggplot2::geom_jitter(
    height = 0.2,
    width  = 0,
    ggplot2::aes(color = Category),
    size   = 1,
    shape  = 16
  ) +
  ggplot2::scale_y_discrete(
    labels   = custom_labels,
    sec.axis = ggplot2::dup_axis(
      labels = right_labels[levels(data_df$Category)],
      name   = ""
    )
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(max_rank * 0.15, max_rank * 0.85),
    labels = c("Enriched in\nRNAi", "Enriched in\nCRISPR")
  ) +
  ggplot2::labs(x = "Rank", y = "") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y.left  = ggplot2::element_text(size = 7),
    axis.text.y.right = ggplot2::element_text(size = 7, hjust = 0),
    axis.text.x       = ggplot2::element_text(size = 7),
    legend.position   = "none",
    panel.border      = ggplot2::element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor  = ggplot2::element_blank()
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

## Set paths
path.wd   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/"
path.ctrp <- paste0(path.wd, "DataSets/CTRPv2/")
path.general <- paste0(path.wd, "DataSets/General/")

ctrp.inform <- read.delim(paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F, check.names = F)

View(table(ctrp.inform$target_or_activity_of_compound))




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

