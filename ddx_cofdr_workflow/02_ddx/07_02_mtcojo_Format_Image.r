# mtcojo organize image from FUMA
library(openxlsx)
source("02_ddx/06_02_DDx_Format_Image.r")

# directly use DDx_Format_Image()
# covidhgi = c("A2", "B2", "C2")
# t2phenos = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# for(t1 in covidhgi)
# {
	# DDx_Format_Image(t1, t2 = t2phenos, save_path = paste0("02_ddx/13_mtCOJO_Excel_Summary/Figures/", t1, "_respiratory.xlsx"), ddx_name = paste0(t1, "_respiratory"))
# }

