##### Set up #####
setwd("~/Documents/GitHub/FDB_Freeland/")

library(doParallel)
library(doRNG)
library(missForest)

#### Imputation (NA's in Data) ####

### CRISPR

# pull in data
file.crispr <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/CRISPRGeneEffect.csv"

CRISPR <- read.delim(
  file = file.crispr,
  row.names = 1,
  stringsAsFactors = F,
  sep = ","
)

table(colSums(is.na(CRISPR)))

# miss forest  
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

# save 
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

### CTRP

# pull in data

file.ctrp <- "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/CTRP drug screen/ctrp_gene_matrix.csv" ##??????
ctrp <- read.csv(file.ctrp)

path.dm <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/DepMap_25Q3/"

path.ctrp <- "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/CTRPv2/"

ctrp.cells <- read.csv("/Users/tgraeber/Dropbox/glab/data/functional databases FDB/CTRP drug screen/ctrp_gene_matrix_cellineinfo.csv")

ctrp.c.s <- as.data.frame(scale(ctrp, center = TRUE, scale = TRUE))  #mean centered columns with unit variance

ctrp$DepMap_ID <- ctrp.cells$DepMap_ID
ctrp.c.s$DepMap_ID <- ctrp$DepMap_ID
ctrp <- ctrp %>% dplyr::select(DepMap_ID, everything()) #move DepMap_ID to first column
ctrp.c.s <- ctrp.c.s %>% dplyr::select(DepMap_ID, everything()) #move DepMap_ID to first column
ctrp[1:7,1:5]
ctrp.c.s[1:7,1:5]

file_ctrp_targets = "/Users/jack/Library/CloudStorage/Box-Box/WD_FDB_Freeland/DataSets/CTRPv2/CTRPv2.0._INFORMER_SET.txt"
ctrp.targets = read.delim(file_ctrp_targets, sep = "\t", stringsAsFactors = F)

models = read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F)

ctrp.expt = read.delim(paste0(path.ctrp,"v20.meta.per_experiment.txt"), sep = "\t", stringsAsFactors = F)

temp = as.data.frame(table(ctrp.expt$master_ccl_id))
#table(temp$Freq)





ctrp.cell = read.delim(paste0(path.ctrp,"v20.meta.per_cell_line.txt"), sep = "\t", stringsAsFactors = F)
#master_ccl_id	ccl_name

ctrp.inform = read.delim(paste0(path.ctrp,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F)
#master_cpd_id	cpd_name	broad_cpd_id

ctrp.curves = read.delim(paste0(path.ctrp,"v20.data.curves_post_qc.txt"), sep = "\t", stringsAsFactors = F)
# #cancerous.only
# ctrp.curves = ctrp.curves[ctrp.curves$DepMap_ID %in% samples$DepMap_ID,]

#experiment_id	...	area_under_curve
ctrp.curves$master_ccl_id = ctrp.expt$master_ccl_id[match(ctrp.curves$experiment_id, ctrp.expt$experiment_id)]
ctrp.curves$ccl_name = ctrp.cell$ccl_name[match(ctrp.curves$master_ccl_id, ctrp.cell$master_ccl_id)]
ctrp.curves$DepMap_ID = models$ModelID[match(ctrp.curves$ccl_name, models$StrippedCellLineName)]
ctrp.curves$cpd_name = ctrp.inform$cpd_name[match(ctrp.curves$master_cpd_id, ctrp.inform$master_cpd_id)]
ctrp.curves <- ctrp.curves %>% dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, area_under_curve, experiment_id, everything())

not.mapped.celllines = ctrp.curves[is.na(ctrp.curves$DepMap_ID),]
not.mapped.celllines = not.mapped.celllines[,1:3] %>% dplyr::distinct()

#write.table(ctrp.curves, paste0(path.ctrp,"ctrp.curves.txt"), quote = F, col.names=T, row.names = F, sep = "\t")

ctrp.rsl3 = ctrp.curves[ctrp.curves$cpd_name == "1S,3R-RSL-3",]
table(table(ctrp.rsl3$DepMap_ID))

#write.table(ctrp.rsl3, paste0(path.ctrp,"ctrp.rsl3.txt"), quote = F, col.names=T, row.names = F, sep = "\t")

ctrp.curves.abr = ctrp.curves %>% dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, master_cpd_id, area_under_curve)

ctrpv2.ave = ctrp.curves.abr %>% group_by(DepMap_ID,ccl_name,master_ccl_id,cpd_name,master_cpd_id) %>% summarise(avg = mean(area_under_curve)) %>% ungroup()
ctrp.ave.rsl3 = ctrpv2.ave[ctrpv2.ave$cpd_name == "1S,3R-RSL-3",]

names(ctrpv2.ave)

ctrpv2.ave.wide = ctrpv2.ave %>% dplyr::select(DepMap_ID, ccl_name, cpd_name, avg)
ctrpv2.ave.wide.pre = as.data.frame(ctrpv2.ave %>% dplyr::select(DepMap_ID, cpd_name, avg))

ctrpv2.ave.wide = reshape(ctrpv2.ave.wide.pre, idvar = "DepMap_ID", v.names= c("avg"), timevar = "cpd_name", direction = "wide")
names(ctrpv2.ave.wide) = gsub("^avg\\.","",names(ctrpv2.ave.wide))

ctrpv2 = ctrpv2.ave.wide
table(colSums(is.na(ctrpv2)))
#fairly random summary distribution of NAs

#remove drugs with > 80% NAs----

ctrpv3 = ctrpv2
ctrpv3$DepMap_ID = ifelse(is.na(ctrpv3$DepMap_ID),"no.depmap.match",ctrpv3$DepMap_ID)
row.names(ctrpv3) = ctrpv3$DepMap_ID
ctrpv3 = ctrpv3[,-1]

percent.nas <- as.data.frame(colMeans(is.na(ctrpv3)) * 100)
names(percent.nas) = "percent.nas"

#keep only drugs with data for > 20% of the lines
percent.nas$eighty.percent.keep.flag = ifelse(percent.nas$percent.nas>80,0,1)
print("cut drugs")
percent.nas[percent.nas$eighty.percent.keep.flag==0,]

ctrpv3.culled = ctrpv3[,names(ctrpv3) %in% row.names(percent.nas[percent.nas$eighty.percent.keep.flag==1,])]

#
#missForest----
doParallel::registerDoParallel(cores = detectCores() - 2)
doRNG::registerDoRNG(seed = 999)
set.seed(999)

if (!file.exists(gsub(".csv$","-as.matrix_culled80_MFImputed.txt",file.ctrp))) {
  if (1) {
    print(Sys.time())
    ctrpv3.culled_mf = missForest(
      xmis = ctrpv3.culled,
      parallelize = "variables",
      verbose = T
    )
    print(Sys.time())
  }
  
  #ctrpv3.culled_mf$OOBerror
  ctrpv3.culled_mf.imp = ctrpv3.culled_mf$ximp
  Sys.time()
  ctrpv3.culled_mf.imp_t = as.data.frame(
    t(ctrpv3.culled_mf.imp),
    stringsAsFactors = F
  )
  Sys.time()
  
  write.table(
    x = ctrpv3, 
    file = gsub(".csv$","-as.matrix.txt",file.ctrp), 
    quote = F, 
    sep = "\t",
    col.names = NA
  )
  write.table(
    x = ctrpv3.culled, 
    file = gsub(".csv$","-as.matrix_culled80.txt",file.ctrp), 
    quote = F, 
    sep = "\t",
    col.names = NA
  )
  
  write.table(
    x = ctrpv3.culled_mf.imp, 
    file = gsub(".csv$","-as.matrix_culled80_MFImputed.txt",file.ctrp), 
    quote = F, 
    sep = "\t",
    col.names = NA
  )
  
  write.table(
    x = ctrpv3.culled_mf.imp_t, 
    file = gsub(".csv$","-as.matrix_culled80_MFImputed_sg.txt",file.ctrp), 
    quote = F, 
    sep = "\t",
    col.names = NA
  )
  
}

}







