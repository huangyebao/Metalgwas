#' @title 使用METAL软件对GWAS进行meta分析
#'
#' @description windows系统下使用METAL软件对全基因组关联研究(GWAS)数据进行META分析，
#' 注意：输入的数据格式必须为完整的TwoSampleMR格式文件，电脑内存小容易运行崩溃，
#' 且必须包含SNP、samplesize、effect_allele、other_allele、eaf、beta、se、pval列，
#' 且数据不能有缺失，若需要MAF过滤，需提前将eaf转换，并包含MAF列。
#'
#' @param GWASfile 字符向量，包含需要meta分析的所有TwoSampleMR格式的完整GWAS文件路径，注意必须是相同type类型。
#' @param Metal_path METAL软件的可执行文件路径，例如"./metal.exe"。
#' @param GWAS_name 输出文件的文件名。
#' @param save_path 文件保存路径。
#' @param type TwoSampleMR格式数据类型，outcome或exposure。
#' @param Analysis_Scheme 指定分析加权方法，samplesize或者SE。samplesize表示基于权重（通常为样本量）和效应方向，跨研究合并p值。
#' SE表示使用相应标准误差的倒数对效应大小估计值进行加权。
#' @param heterogeneity TRUE或FALSE，是否进行异质性分析。
#' @param Average_FREQ TRUE或FALSE，是否计算平均频率。
#' @param MinMax_FREQ TRUE或FALSE，是否计算最小和最大频率。
#' @param Filter_N 数值或FALSE，如果为数值，则添加一个基于样本数量的过滤。如果为FALSE，则不添加此过滤。
#' @param Filter_MAF 数值或FALSE，如果为数值，则添加一个基于次要等位基因频率(MAF)的过滤值，一般设置为0.01。
#' 如果为FALSE，则不添加此过滤。
#' @param Genomic_Control_Correction TRUE或FALSE，是否进行基因组控制校正。
#' @param Sample_Overlap_Correction TRUE或FALSE，是否进行样本重叠校正。
#' @import dplyr
#' @export

METAL_gwas <- function(GWASfile = c("./dat1.txt","./dat1.txt"),
                       Metal_path = "./metal.exe",
                       GWAS_name ="METAL_gwas",
                       save_path ="./",
                       type = "outcome",
                       Analysis_Scheme = "samplesize",
                       heterogeneity = TRUE,
                       Average_FREQ = TRUE,
                       MinMax_FREQ = FALSE,
                       Filter_N  = FALSE,
                       Filter_MAF = FALSE,
                       Genomic_Control_Correction = FALSE,
                       Sample_Overlap_Correction = FALSE){

  message("详见请参照官网：https://genome.sph.umich.edu/wiki/METAL_Documentation")
  Sys.sleep(2)
  cat("\t")
  stopifnot(type %in% c("exposure","outcome"))
  stopifnot(Analysis_Scheme %in% c("samplesize","SE"))
  stopifnot(heterogeneity %in% c(TRUE,FALSE))
  stopifnot(Average_FREQ %in% c(TRUE,FALSE))
  stopifnot(MinMax_FREQ %in% c(TRUE,FALSE))
  stopifnot(Genomic_Control_Correction %in% c(TRUE,FALSE))
  stopifnot(Sample_Overlap_Correction %in% c(TRUE,FALSE))
  outcome_metal<- outcome_metal
  metal <- c()
  if (Analysis_Scheme == "SE") {
    metal_se <- outcome_metal[c(1:2),]
    metal_se <- data.frame(a=metal_se)
    metal <- rbind(metal,metal_se)
  }

  if (Genomic_Control_Correction == TRUE) {
    metal_Genomic_Control_Correction <- outcome_metal[3,]
    metal_Genomic_Control_Correction  <- data.frame(a=metal_Genomic_Control_Correction)
    metal <- rbind(metal,metal_Genomic_Control_Correction)
  }
  if (Sample_Overlap_Correction == TRUE) {
    metal_Sample_Overlap_Correction <- outcome_metal[4,]
    metal_Sample_Overlap_Correction  <- data.frame(a=metal_Sample_Overlap_Correction)
    metal <- rbind(metal,metal_Sample_Overlap_Correction)
  }

  if (Average_FREQ == TRUE) {
    metal_Average_FREQ <- outcome_metal[5,]
    metal_Average_FREQ  <- data.frame(a=metal_Average_FREQ )
    metal <- rbind(metal,metal_Average_FREQ)
  }

  if (MinMax_FREQ == TRUE) {
    metal_MinMax_FREQ <- outcome_metal[6,]
    metal_MinMax_FREQ  <- data.frame(a=metal_MinMax_FREQ)
    metal <- rbind(metal,metal_MinMax_FREQ)
  }
  gwas_N <- "samplesize.outcome"
  if (type == "exposure") {
    gwas_N = gsub(pattern = "outcome",replacement = "exposure",x = gwas_N)
  }
  if (is.numeric(Filter_N)) {
    metal_Filter_N <- paste("ADDFILTER ", gwas_N,">",Filter_N)
    metal_Filter_N  <- data.frame(a=metal_Filter_N)
    metal <- rbind(metal,metal_Filter_N)
  }

  if (is.numeric(Filter_MAF)) {
    metal_Filter_MAF <- paste("ADDFILTER MAF >", Filter_N)
    metal_Filter_MAF  <- data.frame(a=metal_Filter_MAF)
    metal <- rbind(metal,metal_Filter_MAF)
  }
  gwas_col <- outcome_metal[c(9:14),]
  if (type == "exposure") {
    gwas_col = gsub(pattern = "outcome",replacement = "exposure",x = gwas_col)
  }
  gwas_col <- data.frame(a=gwas_col)
  metal <- rbind(metal,gwas_col)

  metal_filenames <- paste0("PROCESS  ", GWASfile)
  metal_filenames <- data.frame(a=metal_filenames)
  metal <- rbind(metal,metal_filenames)

  metal_OUTFILE  <- data.frame(a=paste0("OUTFILE  ",save_path,GWAS_name," ",".txt"))
  metal <- rbind(metal,metal_OUTFILE)


  if (heterogeneity == TRUE) {
    metal_het  <- data.frame(a="ANALYZE HETEROGENEITY")
    metal <- rbind(metal,metal_het)
  } else {
    metal_het  <- data.frame(a="ANALYZE")
    metal <- rbind(metal,metal_het)
  }
  write.table(metal,file=paste0(save_path,"metal.txt"),
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = F)
  system2(command = Metal_path,
          args = paste0(save_path,"metal.txt"))

}
