##### Set up #####
setwd("~/Documents/GitHub/FDB_Freeland/")

library(doParallel)
library(doRNG)
library(missForest)
library(dplyr)
library(mixOmics)
library(stringr)
library(ggplot2)
library(ggrepel)

#### Imputation: CRISPR (NA's in Data) ####

## pull in data
file.crispr <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/CRISPRGeneEffect.csv"

CRISPR <- read.delim(
  file = file.crispr,
  row.names = 1,
  stringsAsFactors = F,
  sep = ",",
  check.names = F
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
models <- read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F, check.names = F)

## load in CTRP Data
ctrp.expt   <- read.delim(paste0(path.ctrp,"v20.meta.per_experiment.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
ctrp.cell   <- read.delim(paste0(path.ctrp,"v20.meta.per_cell_line.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
ctrp.inform <- read.delim(paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
ctrp.curves <- read.delim(paste0(path.ctrp,"v20.data.curves_post_qc.txt"), sep = "\t", stringsAsFactors = F, check.names = F)

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

## create a culled version - rename one entry with no DepMap_ID - keep only drugs with data for > 20% of the cell lines
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

##### PLSR Run: CRISPR & CTRP #####

## set paths
path.dm   <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/"
path.ctrp <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/CTRPv2/"

## read in data and set row names
CRISPR <- read.delim(file = paste0(path.dm, "CRISPRGeneEffect_MFImputed.txt"), sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1)
CTRP <- read.delim(file = paste0(path.ctrp, "ctrpv2.wide_culled80_MFImputed.txt"), sep = "\t", stringsAsFactors = F, check.names = F, row.names = 1)

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

print(pls_fit$prop_expl_var$X)
print(pls_fit$prop_expl_var$Y)

## extract from pls_fit object
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

# rownames(x.exp_variance) = paste0("comp.",seq(1,nrow(x.exp_variance)))
# rownames(y.exp_variance) = paste0("comp.",seq(1,nrow(y.exp_variance)))

variates.X.Y = merge(x.variates,y.variates,by="Score",suffixes = c(".geneexp",".crispr"))

## save files
path.pls <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/PLS/"

write.table(
  x = x.variates,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_X.variates.txt"),
  sep = "\t", quote = F, row.names = F)

write.table(
  x = y.variates,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_Y.variates.txt"),
  sep = "\t", quote = F, row.names = F)

write.table(
  x = variates.X.Y,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_X.Y.variates.txt"),
  sep = "\t", quote = F, row.names = F)

write.table(
  x = x.loadings,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_X.loadings.txt"),
  sep = "\t", quote = F, row.names = F)

write.table(
  x = y.loadings,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_Y.loadings.txt"),
  sep = "\t", quote = F, row.names = F)

write.table(
  x = x.exp_variance,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_X.expvar.txt"),
  sep = "\t", quote = F, row.names = F)

write.table(
  x = y.exp_variance,
  file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_Y.expvar.txt"),
  sep = "\t", quote = F, row.names = F)

##### PLSR Plot: CRISPR & CTRP #####

#### load in drug data
path.ctrp <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/CTRPv2/"
path.pls <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/PLS/"

ctrp.inform <- read.delim(file = paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F, check.names = F)
Y_loadings <- read.delim(file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_Y.loadings.txt"), sep = "\t", stringsAsFactors = F, check.names = F)

## Prep for plotting
lookup_indices <- match(Y_loadings$Loading, ctrp.inform$cpd_name)
Y_loadings$drug.target <- ctrp.inform$target_or_activity_of_compound[lookup_indices]

## adding extra collumns
Y_loadings <- Y_loadings %>%
  dplyr::mutate(
    group = NA,
    group = case_when(
      stringr::str_detect(Loading, "^(selumetinib|PD318088|RAF265|regorafenib|PLX\\-4720|dabrafenib|GDC\\-0879)$") ~ "01 BRAFi.MEKi",
      str_detect(Loading, "^(erlotinib|afatinib|lapatinib|neratinib|canertinib|vandetanib|gefitinib)$") ~ "02 EGFRi.HER2i",
      stringr::str_detect(Loading, "^(1S\\,3R\\-RSL\\-3|ML210|erastin|ML162)") ~ "03 ferropt",
      stringr::str_detect(Loading, "^(nutlin\\-3|HBX\\-41108|KU\\-60019)$") ~ "04 MDM2i",
      stringr::str_detect(Loading, "^oligomycin[\\ .]?A$") ~ "05 oligomycinA",
      stringr::str_detect(Loading, "^dasatinib") ~ "06 SRC",

      stringr::str_detect(drug.target %||% "", "BCL2") & !str_detect(Loading, ":") ~ "07 BCL2+i",
      
      TRUE ~ NA_character_
    ),
    
    # group.atp5 only for oligomycin A
    group.atp5 = NA,
    group.atp5 = if_else(str_detect(Loading, "^oligomycin[\\ .]?A$"), "05 oligomycinA", NA_character_),
    
    # NA flags & labels
    group.na = if_else(is.na(group), 1L, 0L),
    group.atp5.na = if_else(is.na(group.atp5), 1L, 0L),
    label.not.na = if_else(!is.na(group), Loading, NA_character_),
    label.not.na.atp5 = if_else(!is.na(group.atp5), Loading, NA_character_),
    
    # mix flag (dual vs single drug) based on presence of ':'
    mix.flag = if_else(str_detect(Loading, ":"), "dual drug", "single drug")
  ) %>%
  arrange(desc(group.na))

## Target category bucketing
detect <- function(x, pattern) str_detect(x %||% "", regex(pattern, ignore_case = TRUE))

Y_loadings <- Y_loadings %>%
  mutate(target.category = NA_character_) %>%
  mutate(target.category = if_else(detect(drug.target, "DNA damage"), "DNA.damage", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "(micro|mi)rotubule"), "microtubule", target.category)) %>%  # handles typos
  mutate(target.category = if_else(detect(drug.target, "polo\\-like kinase 1|\\bPLK1\\b"), "PLK1", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "polo\\-like kinase 2|\\bPLK2\\b"), "PLK2", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "aurora kinase"), "aurora", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "DNA methyltransferase"), "DNA meth", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "DNA replication"), "DNA rep", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "nicotinamide phosphoribosyltransferase|\\bNAMPT\\b"), "NAMPT", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "dihydrofolate reductase|\\bDHFR\\b"), "DHFR", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "BCL2"), "BCL2.", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "XXXX"), "XXXX", target.category)) %>%
  mutate(target.category = if_else(detect(drug.target, "XXXX"), "XXXX", target.category))


## add in % NA's column ( Y from above section)
percent.nas <- as.data.frame(colMeans(is.na(Y)) * 100)
names(percent.nas) <- "percent.nas"
percent.nas <- percent.nas %>%
  tibble::rownames_to_column(var = "Loading")

Y_loadings <- Y_loadings %>%
  left_join(percent.nas, by = "Loading")

Y_loadings <- Y_loadings %>%
  dplyr::select(Loading, group, percent.nas, contains("target"), contains("flag"), everything())


#### load in CRISPR data
gene.info.all <- read.delim(
  file = "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/General/NCBI Gene Homo_sapiens.gene_info 20220719.txt", sep = "\t", stringsAsFactors = F, check.names = F)
gene.info <- gene.info.all[gene.info.all$Symbol_from_nomenclature_authority !="-",]
gene.info.abr <- gene.info %>% dplyr::select(Symbol, description)

X_loadings <- read.delim(file = paste0(path.pls, "PLS_Mode.Regression_X.CRISPR_Y.CTRP_X.loadings.txt"), sep = "\t", stringsAsFactors = F, check.names = F) %>%
  dplyr::mutate(Loading = sub("\\.\\..*$", "", Loading))

X_loadings <- merge(X_loadings, gene.info.abr, by.x="Loading", by.y="Symbol", all.x=T)

X_loadings <- X_loadings %>%
  dplyr::mutate(
    group = case_when(
      stringr::str_detect(Loading, "^(BRAF|MITF|MAPK1|SOX9|SOX10|PEA15|DUSP4)")     ~ "01 BRAF sig",
      stringr::str_detect(Loading, "^(EGFR|KLF5|STX4|GRHL2|PIK3CA|ERBB2)$")         ~ "02 EGFR sig",
      stringr::str_detect(Loading, "^(GPX4|SEPSECS|PSTK|EEFSEO|SEPHS2|SECISBP2)$")  ~ "03 ferropt",
      stringr::str_detect(Loading, "^MDM[24]$")                                     ~ "04 MDM2.MDM4",
      stringr::str_detect(Loading, "^ATP5")                                         ~ "05 ATP5",
      stringr::str_detect(Loading, "^(ABL|SRC|LCK|LYN)")                            ~ "06 dasa targets",
      stringr::str_detect(Loading, "^(BCL2|BCL2L1|BCL2L2|MCL1)$")                   ~ "07 BCL2+",
      stringr::str_detect(Loading, "^MYC(|N|L)")                                    ~ "08 MYC.",
      stringr::str_detect(Loading, "^(GRB2|CRKL)")                                  ~ "09 SRC-related",
      stringr::str_detect(Loading, "^TP53$")                                        ~ "10 TP53",
      stringr::str_detect(Loading, "^MED12$")                                       ~ "11 MED12",
      TRUE ~ NA_character_
    ),
    group.atp5 = if_else(str_detect(Loading, "^ATP5"), "05 ATP5", NA_character_)
  )

X_loadings <- X_loadings %>%
  dplyr::select(Loading, description, group, everything()) %>%
  dplyr::mutate(
    group.na         = if_else(is.na(group), 1L, 0L),
    group.atp5.na    = if_else(is.na(group.atp5), 1L, 0L),
    label.not.na     = if_else(!is.na(group), Loading, NA_character_),
    label.not.na.atp5= if_else(!is.na(group.atp5), Loading, NA_character_)
  ) %>%
  dplyr::arrange(desc(group.na))

# variates.X.Y.c.c.plotting = merge(samples %>% dplyr::select(stripped_cell_line_name,lineage_subtype,OncotreeLineage),variates.X.Y.c.c,by.x="stripped_cell_line_name",by.y="Score")

## plotting
my_colors <- c("#F8766D","#DE8C00","#B79F00","#00BA38","#00BF7D","#00BFC4","#00B4F0","#619CFF"
               ,"hotpink","purple","cyan")
my_colors_main <- my_colors

ggplot2::ggplot(Y_loadings, aes_string(x = "comp1", y = "comp2", color = "target.category", label = "target.category"))  + 
  geom_point() + 
  geom_text_repel() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) + 
  scale_color_manual(values = my_colors, na.value = "grey80") +
  theme_minimal()










