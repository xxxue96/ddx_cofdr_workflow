# Cofdr organize image from FUMA
library(openxlsx)

Cofdr_Format_Image <- function(t1, t2, save_path)
{
	phenos = paste0(t1, "_", t2)	
	gene2func = paste0("03_cofdr/05_02_Cofdr_FUMA_Result/GENE2FUNC/",phenos, ".png")
	
	# remove unexisted figure paths
	rem = which(unlist(lapply(gene2func, file.exists)))
	phenos = phenos[rem]
	gene2func = gene2func[rem]
	
	# @grids: rows and cols to be formatted
	# @heights: row heights
	format_header <- function(wb, sheet, grids, header, header_style, heights=13)
	{
		mergeCells(wb, sheet, cols = 1:grids[2], rows = grids[1])
		setRowHeights(wb, sheet, rows = grids[1], heights)
		writeData(wb, sheet, header, startCol = 1, startRow = grids[1], colNames = FALSE, rowNames = FALSE, sep = "\t")
		addStyle(wb, sheet, header_style, rows = grids[1], cols = 1, gridExpand = TRUE)
	}
	
	format_title <- function(wb, sheet, grids, header, title_style)
	{
		writeData(wb, sheet, header, startCol = grids[2], startRow = grids[1], colNames = FALSE, rowNames = FALSE, sep = "\t")
		addStyle(wb, sheet, title_style, rows = grids[1], cols = grids[2], gridExpand = TRUE)
	}
	
	format_img <- function(wb, sheet, figures, grids_header, grids_img, img_width, img_height, startRow_step, grids_evenCol, grids_oddCol)
	{
		for(i in 1:length(figures)){	
			figure = figures[i]
			pheno = phenos[i]
			
			format_title(wb, sheet, grids_header, pheno, title_style)
			insertImage(wb, sheet, figure, startRow = grids_img[1], startCol = grids_img[2], width = img_width, height = img_height)
			
			if(i%%2 != 0){
				grids_header[2] = grids_evenCol[1]
				grids_img[2] = grids_evenCol[2] # only change column number
			} else{
				grids_header[1] = grids_header[1] + startRow_step
				grids_header[2] = grids_oddCol[1] # change both column number and row number
				
				grids_img[1] = grids_img[1] + startRow_step
				grids_img[2] = grids_oddCol[2]
			}
		}
	}
	
	if(identical(gene2func, character(0))){
		return()	# no figures
	} else{
		# create workbook
		wb <- createWorkbook()
		addWorksheet(wb, "F1.sig_gene_DEG")

		# header style
		header_style <- createStyle(fontName = "Arial", fontSize = 11, halign = "left", textDecoration = "bold")
		subheader_style <- createStyle(fontName = "Arial", fontSize = 10, halign = "left", wrapText = TRUE)
		title_style <- createStyle(fontName = "Arial", fontSize = 9, halign = "center", textDecoration = "bold")

		# insert image
		h1 = "F1. Enrichment of prioritized genes in DEG Sets for all respiratory diseases"
		subh1 = "For colocalized genes detected by cofdr, gwas-pw and hyprcoloc, significantly enriched DEG sets (Pbon < 0.05) are highlighted in red."
		format_header(wb, "F1.sig_gene_DEG", c(1, 13), h1, header_style)
		format_header(wb, "F1.sig_gene_DEG", c(2, 13), subh1, subheader_style)
		format_img(wb, "F1.sig_gene_DEG", gene2func, grids_header=c(4, 3), grids_img=c(5, 1), img_width=4.63, img_height=2, startRow_step=13, grids_evenCol=c(10, 8), grids_oddCol=c(3, 1))	
		
		saveWorkbook(wb, file = save_path, overwrite = TRUE)
		
		return(save_path)
	}
}

# covidhgi = c("A2", "B1", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# for(t1 in covidhgi)
# {
	# Cofdr_Format_Image(t1, t2 = respiratory, save_path = paste0("03_cofdr/10_Excel_Summary/Figures/", t1, "_respiratory.xlsx"))
# }

