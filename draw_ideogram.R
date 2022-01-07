library(RIdeogram)

setwd('/home/juan/jj/dataAnalysis/andy/andys_integration')

# Get human karyotype from package RIdeogram
data(human_karyotype, package ="RIdeogram")

# Import integration data
gene_density <- read.table("integrations_invitro_1M.tsv", sep = "\t", header = T, stringsAsFactors = F)

# Basic usage
ideogram(karyotype = human_karyotype, overlaid = gene_density)
#, colorset1 = c( "navy", "white","violet"), colorset2 = c("#b35806", "#f7f7f7", "#542788"))
convertSVG("chromosome.svg", device = "png")
#convertSVG("all_HBRV_int_inVitro.svg", device = "png")
