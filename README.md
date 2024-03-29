# Metalgwas
使用METAL软件对GWAS进行meta分析

软件和示例数据下载:https://csg.sph.umich.edu/abecasis/Metal/download/

METAL使用说明文档:https://genome.sph.umich.edu/wiki/METAL_Documentation


示例代码:

install.packages("remotes")

remotes::install_github("huangyebao/Metalgwas")

file1 <- system.file("data", "data1.txt", package = "Metalgwas")

file2 <- system.file("data", "data2.txt", package = "Metalgwas")

data1 <- data.table::fread(file1)

head(data1)

data2 <- data.table::fread(file2)

head(data2)

Metalgwas::METAL_gwas(GWASfile =c(file1,file2), Metal_path = "./metal/metal.exe",save_path = "./", GWAS_name = "METAL_gwas", type = "outcome",Analysis_Scheme = "SE", heterogeneity = TRUE, Average_FREQ = TRUE, MinMax_FREQ = TRUE, Filter_N = FALSE, Filter_MAF = FALSE, Genomic_Control_Correction = FALSE, Sample_Overlap_Correction = FALSE)
