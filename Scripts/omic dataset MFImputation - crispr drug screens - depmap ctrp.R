library(doParallel)
library(doRNG)
library(missForest)
library(dplyr)

# via Favour

#notes on time usage below

# Two sections
# 1-crispr 
# 2-ctrp drug scren further below


read_crispr_flag = 0
if (read_crispr_flag) {

#file.crispr = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 23Q4/CRISPRGeneEffect.csv"
  #file.crispr = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q2/data/CRISPRGeneEffect.csv"
  file.crispr = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q3/data/CRISPRGeneEffect.csv"
  # file.crispr2 = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q2/data/CRISPRGeneDependency.csv"
# file.crispr3 = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q2/data/ScreenGeneEffect.csv"

#data input----
ceres.q4 = read.delim(
  file = file.crispr,
  #file = "Achilles_gene_effect_20Q4.txt",
  row.names = 1,
  stringsAsFactors = F,
  sep = ","
)

# ceres.q4.dependency = read.delim(file = file.crispr2, row.names = 1, stringsAsFactors = F, sep = ",")
# ceres.q4.screen = read.delim(file = file.crispr3, row.names = 1, stringsAsFactors = F, sep = ",")

table(colSums(is.na(ceres.q4)))
# 0     1     2     3   123   125   776   785  1085  1086 
# 17107     2    95     1   426   156   129    15   511     1 

# 0     1     2     3    18    19    23   119   120   121   137   138   139   141   146   166   832   833   834 
# 17079    20     1     1   102     1     1   421     3     1   146     7     1     1     1     1   124     3     2 

# table(colSums(is.na(ceres.q4.dependency)))
# table(colSums(is.na(ceres.q4.screen)))

#
#missForest----
doParallel::registerDoParallel(cores = detectCores() - 2)
doRNG::registerDoRNG(seed = 999)
set.seed(999)

Sys.time()
ceres.q4_mf = missForest(
  xmis = ceres.q4,
  parallelize = "variables",
  verbose = T
)
Sys.time()

#ceres.q4_mf$OOBerror
ceres.q4_mf.imp = ceres.q4_mf$ximp
Sys.time()
ceres.q4_mf.imp_t = as.data.frame(
  t(ceres.q4_mf.imp),
  stringsAsFactors = F
)
Sys.time()


write.table(
  x = ceres.q4_mf.imp, 
  file = gsub(".(csv|txt)$","_MFImputed.txt",file.crispr), 
  quote = F, 
  sep = "\t",
  col.names = NA
)
#file = "Achilles_gene_effect_20Q4_MFImputed.txt", 

write.table(
  x = ceres.q4_mf.imp_t, 
  file = gsub(".(csv|txt)$","_MFImputed_sg.txt",file.crispr), 
  quote = F, 
  sep = "\t",
  col.names = NA
)
#file = "Achilles_gene_effect_20Q4_MFImputed_sg.txt", 

#


#example run on iMac Pro
# 2.3 GHz 18-Core Intel Xeon W
# 128 GB 2666 MHz DDR4

# 1 hr 20 min
# file.crispr = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q2/data/CRISPRGeneEffect.csv"
# > Sys.time()
# [1] "2025-08-30 09:46:01 PDT"
# > ceres.q4_mf = missForest(
#   +   xmis = ceres.q4,
#   +   parallelize = "variables",
#   +   verbose = T
#   + )
# parallelizing over the variables of the input data matrix 'xmis'
# missForest iteration 1 in progress...
# randomForest 4.7-1.2
# Type rfNews() to see new features/changes/bug fixes.
# 
# Attaching package: ‘randomForest’
# 
# The following object is masked from ‘package:ggplot2’:
#   
#   margin
# done!
#   estimated error(s): 0.08083031 
# difference(s): 6.638301e-05 
# time: 1251.926 seconds
# 
# missForest iteration 2 in progress...done!
#   estimated error(s): 0.08096675 
# difference(s): 4.503133e-05 
# time: 1224.723 seconds
# 
# missForest iteration 3 in progress...done!
#   estimated error(s): 0.08086907 
# difference(s): 4.464658e-05 
# time: 1183.191 seconds
# 
# missForest iteration 4 in progress...done!
#   estimated error(s): 0.08090863 
# difference(s): 4.481841e-05 
# time: 1179.247 seconds
# 
# > Sys.time()
# [1] "2025-08-30 11:06:54 PDT"

}


read_ctrp_flag = 1
version = "CTRPv1"
#version = "CTRPv2"

#v2 on mac book air
#"2025-10-18 17:35:47 PDT"
#"2025-10-18 17:51:09 PDT"
# on blue imac


if (read_ctrp_flag) {
  print(paste("read ctrp",format(Sys.time(), "%X")))
  #time.stamps$start.read.ctrp = Sys.time()  
  #colnames(prism)
  
  # NOTES:
  # ctrp dataframe is CTRP v1:              227 cell lines, 481 drugs
  # ctrpv2 dataframe is CTRPv2 (version 2): 840 cell lines, 545 drugs
  
  if (version == "CTRPv1") {
    
    # prep steps
    if (1) { 
      
      #data input----
      file.ctrp = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/CTRP drug screen/ctrp_gene_matrix.csv"
      ctrp <- read.csv(file.ctrp)
      
      # #cancerous.only
      # ctrp = ctrp[ctrp$DepMap_ID %in% samples$DepMap_ID,]
      
      ctrp.cells <- read.csv("/Users/tgraeber/Dropbox/glab/data/functional databases FDB/CTRP drug screen/ctrp_gene_matrix_cellineinfo.csv")
     
      # #cancerous.only
      # ctrp.cells = ctrp.cells[ctrp$DepMap_ID %in% samples$DepMap_ID,]
    
      ctrp.for.mfi = ctrp
      rownames(ctrp.for.mfi) <- ctrp.cells$DepMap_ID
    
      ctrp$DepMap_ID <- ctrp.cells$DepMap_ID
      ctrp <- ctrp %>% dplyr::select(DepMap_ID, everything()) #move DepMap_ID to first column

      
      #culled version - keep only drugs with data for > 20% of the cell lines
      
      # no entries without a DepMap ID
      table(is.na(ctrp$DepMap_ID))

      #remove drugs with > 80% NAs----
      percent.nas.v1 <- as.data.frame(colMeans(is.na(ctrp.for.mfi)) * 100)
      names(percent.nas.v1) = "percent.nas"
      
      #keep only drugs with data for > 20% of the cell lines
      percent.nas.v1$eighty.percent.keep.flag = ifelse(percent.nas.v1$percent.nas>80,0,1)
      print("cut drugs")
      print(percent.nas.v1[percent.nas.v1$eighty.percent.keep.flag==0,])
      
      ctrp.for.mfi.culled = ctrp.for.mfi[,names(ctrp.for.mfi) %in% row.names(percent.nas.v1[percent.nas.v1$eighty.percent.keep.flag==1,])]
      
    }
    
    #run imputation (missForrest) and save new files
    if (1) { 
      #
      #missForest----
      doParallel::registerDoParallel(cores = detectCores() - 2)
      doRNG::registerDoRNG(seed = 999)
      set.seed(999)
      
      if (!file.exists(gsub(".(csv|txt)$","_culled80_MFImputed.txt",file.ctrp))) {
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
          x = ctrpv3.culled_mf.imp, 
          file = gsub(".(csv|txt)$","_culled80_MFImputed.txt",file.ctrp), 
          quote = F, 
          sep = "\t",
          col.names = NA
        )
        
        write.table(
          x = ctrpv3.culled_mf.imp_t, 
          file = gsub(".(csv|txt)$","_culled80_MFImputed_sg.txt",file.ctrp), 
          quote = F, 
          sep = "\t",
          col.names = NA
        )
        
      }
    }
    
    
    
  }
  
  if (version == "CTRPv2") {
    
    # prep steps
    if (1) { 
      
      print(Sys.time())
      
      #data input----
      #file to be created:
      
      path.ctrpv2 = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/CTRPv2/"
      file.ctrpv2.wide = paste0(path.ctrpv2,"ctrpv2.wide.txt")
      
      #data compilation----
      
      #cell line info comes from DepMap
      version = "25Q3"
      if (version == "25Q3") {
        working.dir = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q3/"
        path.dm = "/Users/tgraeber/Dropbox/glab/data/functional databases FDB/DepMap 25Q3/data/"
      }
      

      # file_ctrp_targets = "/Users/tgraeber/Dropbox/glab/functional databases FDB/CTRPv2/CTRPv2.0._INFORMER_SET.txt"
      # ctrp.targets = read.delim(file_ctrp_targets, sep = "\t", stringsAsFactors = F)
      
      #drug.info <- read.csv("/Users/tgraeber/Dropbox/glab/data/functional databases FDB/primary-screen-replicate-treatment-info.csv")
      
      models = read.delim(paste0(path.dm,"Model.csv"), sep = ",", stringsAsFactors = F)
      #still need master_ccl_id to DepMap_ID
      # #cancerous.only
      # models = models[models$DepMap_ID %in% samples$DepMap_ID,]
      
      ctrp.expt = read.delim(paste0(path.ctrpv2,"v20.meta.per_experiment.txt"), sep = "\t", stringsAsFactors = F)
      #experiment_id master_ccl_id
      temp = as.data.frame(table(ctrp.expt$master_ccl_id))
      #table(temp$Freq)
      
      ctrp.cell = read.delim(paste0(path.ctrpv2,"v20.meta.per_cell_line.txt"), sep = "\t", stringsAsFactors = F)
      #master_ccl_id	ccl_name
      
      ctrp.inform = read.delim(paste0(path.ctrpv2,"CTRPv2.0._INFORMER_SET.txt"), sep = "\t", stringsAsFactors = F)
      #master_cpd_id	cpd_name	broad_cpd_id
      
      ctrp.curves = read.delim(paste0(path.ctrpv2,"v20.data.curves_post_qc.txt"), sep = "\t", stringsAsFactors = F)
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
      
      #write.table(ctrp.curves, paste0(path.ctrpv2,"ctrp.curves.txt"), quote = F, col.names=T, row.names = F, sep = "\t")
      
      ctrp.rsl3 = ctrp.curves[ctrp.curves$cpd_name == "1S,3R-RSL-3",]
      table(table(ctrp.rsl3$DepMap_ID))
      
      #write.table(ctrp.rsl3, paste0(path.ctrpv2,"ctrp.rsl3.txt"), quote = F, col.names=T, row.names = F, sep = "\t")
      
      ctrp.curves.abr = ctrp.curves %>% dplyr::select(DepMap_ID, ccl_name, master_ccl_id, cpd_name, master_cpd_id, area_under_curve)
      
      ctrpv2.ave = ctrp.curves.abr %>% group_by(DepMap_ID,ccl_name,master_ccl_id,cpd_name,master_cpd_id) %>% summarise(avg = mean(area_under_curve)) %>% ungroup()
      ctrp.ave.rsl3 = ctrpv2.ave[ctrpv2.ave$cpd_name == "1S,3R-RSL-3",]
      
      names(ctrpv2.ave)
      
      ctrpv2.ave.wide = ctrpv2.ave %>% dplyr::select(DepMap_ID, ccl_name, cpd_name, avg)
      ctrpv2.ave.wide.pre = as.data.frame(ctrpv2.ave %>% dplyr::select(DepMap_ID, cpd_name, avg))
      
      ctrpv2.ave.wide = reshape(ctrpv2.ave.wide.pre, idvar = "DepMap_ID", v.names= c("avg"), timevar = "cpd_name", direction = "wide")
      names(ctrpv2.ave.wide) = gsub("^avg\\.","",names(ctrpv2.ave.wide))
      
      ctrpv2 = ctrpv2.ave.wide
      
      if (0) {
        table(colSums(is.na(ctrpv2)))
        #fairly random summary distribution of NAs
      }
      
      
      #culled version - rename one entry with no DepMap_ID - keep only drugs with data for > 20% of the cell lines
      ctrpv3 = ctrpv2
      
      #rename one entry with no DepMap_ID
      ctrpv3$DepMap_ID = ifelse(is.na(ctrpv3$DepMap_ID),"no.depmap.match",ctrpv3$DepMap_ID)
      row.names(ctrpv3) = ctrpv3$DepMap_ID
      ctrpv3 = ctrpv3[,-1]
      
      #remove drugs with > 80% NAs----
      percent.nas <- as.data.frame(colMeans(is.na(ctrpv3)) * 100)
      names(percent.nas) = "percent.nas"
      
      #keep only drugs with data for > 20% of the cell lines
      percent.nas$eighty.percent.keep.flag = ifelse(percent.nas$percent.nas>80,0,1)
      print("cut drugs")
      print(percent.nas[percent.nas$eighty.percent.keep.flag==0,])
      
      ctrpv3.culled = ctrpv3[,names(ctrpv3) %in% row.names(percent.nas[percent.nas$eighty.percent.keep.flag==1,])]
      
      
      if (!file.exists(file.ctrpv2.wide)) {
        write.table(ctrpv2, file.ctrpv2.wide, quote = F, col.names=T, row.names = F, sep = "\t")
        write.table(ctrpv3.culled, gsub(".(csv|txt)$","_culled80.\\1",file.ctrpv2.wide), quote = F, col.names=T, row.names = F, sep = "\t")
      } else {
        print(paste0("WARNING: do you want to overwrite ",file.ctrpv2.wide," and the culled version?  If yes, delete these files and re-run"))
        if (0) { 
          file.remove(file.ctrpv2.wide) 
          file.remove(gsub(".(csv|txt)$","_culled80.\\1",file.ctrpv2.wide)) 
        }
        #read.delim(paste0("ERROR: do you want to overwrite file ",file.ctrpv2.wide,"?  If yes, delete this file and re-run"))
      }
    
      print(Sys.time())
      
    }
    
    #run imputation (missForrest) and save new files
    if (1) { 
      #
      #missForest----
      doParallel::registerDoParallel(cores = detectCores() - 2)
      doRNG::registerDoRNG(seed = 999)
      set.seed(999)
      
      if (!file.exists(gsub(".(csv|txt)$","_culled80_MFImputed.txt",file.ctrpv2.wide))) {
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
          x = ctrpv3.culled_mf.imp, 
          file = gsub(".(csv|txt)$","_culled80_MFImputed.txt",file.ctrpv2.wide), 
          quote = F, 
          sep = "\t",
          col.names = NA
        )
        
        write.table(
          x = ctrpv3.culled_mf.imp_t, 
          file = gsub(".(csv|txt)$","_culled80_MFImputed_sg.txt",file.ctrpv2.wide), 
          quote = F, 
          sep = "\t",
          col.names = NA
        )
        
      }
    }
    
  }
  
}
