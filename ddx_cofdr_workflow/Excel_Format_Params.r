read_all_sheets = function(xlsxFile, ...) {
	sheet_names = openxlsx::getSheetNames(xlsxFile)
	sheet_list = as.list(rep(NA, length(sheet_names)))
	names(sheet_list) = sheet_names
	for (sn in sheet_names) {
		sheet_list[[sn]] = openxlsx::read.xlsx(xlsxFile, sheet=sn, ...)
	}
	return(sheet_list)
}

# @grids: rows and cols to be formatted
# @heights: row heights
format_header <- function(wb, sheet, grids, header, header_style, heights=13)
{
	mergeCells(wb, sheet, cols = 1:grids[2], rows = grids[1])
	setRowHeights(wb, sheet, rows = grids[1], heights)
	writeData(wb, sheet, header, startCol = 1, startRow = grids[1], colNames = FALSE, rowNames = FALSE, sep = "\t")
	addStyle(wb, sheet, header_style, rows = grids[1], cols = 1, gridExpand = TRUE)
}
	
format_main <- function(wb, sheet, dat, title_style, main_style, startRow = 3)
{
	writeData(wb, sheet, dat, startCol = 1, startRow, colNames = TRUE, rowNames = FALSE, sep = "\t", headerStyle = title_style, keepNA = TRUE, na.string = "NA")
	addStyle(wb, sheet, main_style, rows = (startRow+1):(nrow(dat)+4), cols = 1:ncol(dat), gridExpand = TRUE)
}
	
# header style
header_style <- createStyle(fontName = "Arial", fontSize = 11, halign = "left", textDecoration = "bold")
subheader_style <- createStyle(fontName = "Arial", fontSize = 10, halign = "left", wrapText = TRUE)
title_style <- createStyle(fontName = "Arial", fontSize = 10, halign = "left", fgFill = "#8DB4E2")
main_style <- createStyle(fontName = "Arial", fontSize = 10, halign = "left")
other_style <- createStyle(fontName = "Arial", fontSize = 10, halign = "left", textDecoration = "bold")
