load("cytof_test.aug13.RData")
library(Giotto)
dimGenePlot(cytof_test, expression_values="norm", genes=c("Cell_152Sm_CD45", "Cell_173Yb_CD45RO", "Cell_166Er_CD45RA", "Cell_151Eu_CD31", "Cell_159Tb_CD68", "Cell_170Er_CD3", "Cell_143Nd_Vimentin", "Cell_148Nd_PanCK", "Cell_146Nd_CD16", "Cell_175Lu_KIFC1", "Cell_158Gd_ECadherin", "Cell_156Gd_CD4", "Cell_163Dy_VEGF", "Cell_164Dy_HIF1a", "Cell_147Sm_CD163"), cow_n_col=3, point_size=1, dim_reduction_to_use="tsne", dim_reduction_name="tsne", point_border_stroke=0)
