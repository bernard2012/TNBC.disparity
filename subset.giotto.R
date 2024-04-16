#sample<-commandArgs(trailingOnly=T)[1]
#roi<-as.integer(commandArgs(trailingOnly=T)[2])

do_one_subset <- function(samplex){
subset1=subsetGiotto(cytof_test, cell_ids=cell_metadata[sample==samplex][["cell_ID"]])
subset1=createSpatialNetwork(subset1, name="delauney_network")
#kmeans_annot
cell_proximities=cellProximityEnrichment(subset1, cluster_column="k20.10000_annot", spatial_network_name="delauney_network", number_of_simulations=1000, set_seed=FALSE)
write.table(file=paste0("mar23/prox.", samplex, ".txt"), cell_proximities$enrichm_res, quot=F, sep="\t")
}
