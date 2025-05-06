# DDx organize image from FUMA
library(openxlsx)

# @ddx_name: set a pseudo-name for the DDx-derived trait sumstats 
DDx_Format_Image <- function(t1, t2, save_path, ddx_name=NULL)
{
	if(is.null(ddx_name)){ phenos = paste0(c(t1,t2), collapse="_") } else{ phenos = ddx_name }
	SNP_Manhattan = paste0("02_ddx/04_FUMA_Results/Image/01_SNP_Manhattan/", phenos, ".png")
	Gene_Manhattan = paste0("02_ddx/04_FUMA_Results/Image/02_Gene_Manhattan/", phenos, ".png")
	MAGMA_Tissue = paste0("02_ddx/04_FUMA_Results/Image/03_MAGMA_Tissue/", phenos, ".png")
	SNP_Annot = paste0("02_ddx/04_FUMA_Results/Image/04_SNP_Annot/", phenos, ".png")
		
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
		
	# create workbook
	wb <- createWorkbook()
	addWorksheet(wb, "F1.SNP_Manhattan")
	addWorksheet(wb, "F2.Gene_Manhattan")
	addWorksheet(wb, "F3.MAGMA_Tissue_Enrichment")
	addWorksheet(wb, "F4.SNP_Annot")

	# header style
	header_style <- createStyle(fontName = "Arial", fontSize = 11, halign = "left", textDecoration = "bold")
	subheader_style <- createStyle(fontName = "Arial", fontSize = 10, halign = "left", wrapText = TRUE)
	title_style <- createStyle(fontName = "Arial", fontSize = 9, halign = "center", textDecoration = "bold")

	# insert image
	h1 = "F1. Manhattan plot of the input GWAS summary statistics"
	subh1 = "Genome wide significance (red dashed line in the plot) was defined at P = 5e-8"
	format_header(wb, "F1.SNP_Manhattan", c(1, 13), h1, header_style)
	format_header(wb, "F1.SNP_Manhattan", c(2, 13), subh1, subheader_style)
	format_img(wb, "F1.SNP_Manhattan", SNP_Manhattan, grids_header=c(4, 3), grids_img=c(5, 1), img_width=4.63, img_height=2, startRow_step=13, grids_evenCol=c(10, 8), grids_oddCol=c(3, 1))	
	
	h2 = "F2. Manhattan plot of the gene-based test as computed by MAGMA based on input GWAS summary statistics"
	subh2 = "Genome wide significance (red dashed line in the plot) was defined at P = 0.05/number of tested genes"
	format_header(wb, "F2.Gene_Manhattan", c(1, 13), h2, header_style)
	format_header(wb, "F2.Gene_Manhattan", c(2, 13), subh2, subheader_style)	
	format_img(wb, "F2.Gene_Manhattan", Gene_Manhattan, grids_header=c(4, 3), grids_img=c(5, 1), img_width=4.63, img_height=2.5, startRow_step=15, grids_evenCol=c(10, 8), grids_oddCol=c(3, 1))
	
	h3 = "F3. MAGMA tissue enrichment analysis across 53 tissue types"
	subh3 = "Significantly enriched tissues (highlighted in red) were defined at P = 0.05/53 = 9.43e-4"
	format_header(wb, "F3.MAGMA_Tissue_Enrichment", c(1, 16), h3, header_style)
	format_header(wb, "F3.MAGMA_Tissue_Enrichment", c(2, 16), subh3, subheader_style)	
	format_img(wb, "F3.MAGMA_Tissue_Enrichment", MAGMA_Tissue, grids_header=c(4, 4), grids_img=c(5, 1), img_width=5.71, img_height=3, startRow_step=18, grids_evenCol=c(13, 10), grids_oddCol=c(4, 1))
	
	h4 = "F4. Functional consequences of SNPs on genes"
	format_header(wb, "F4.SNP_Annot", c(1, 14), h4, header_style)
	format_img(wb, "F4.SNP_Annot", SNP_Annot, grids_header=c(3, 3), grids_img=c(4, 1), img_width=5, img_height=2.5, startRow_step=15, grids_evenCol=c(10, 8), grids_oddCol=c(3, 1))

	saveWorkbook(wb, file = save_path, overwrite = TRUE)
	
	return(save_path)
}

# covidhgi = c("A2", "B1", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")


# for(t1 in covidhgi)
# {
	# DDx_Format_Image(t1, t2 = respiratory, save_path = paste0("02_ddx/11_Excel_Summary/Figures/", t1, "_respiratory.xlsx"))
# }

