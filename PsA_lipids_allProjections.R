# Set the working directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd( "/home/joern/.Datenplatte/Joerns Dateien/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/" )

# Function to load the main functions file
loadMainFunctions <- function( ) {
  mainFunctionsFile <- "ProjectionsBiomed_MainFunctions_6.R"

  if ( file.exists( mainFunctionsFile ) ) {
    message( "Loading main functions from ", mainFunctionsFile )
    source( mainFunctionsFile )
  } else {
    stop( "The file ", mainFunctionsFile, " does not exist in the current directory." )
  }
}

# Call the function to load the main functions
loadMainFunctions( )

######################### Load data ##############################################

PsA_lipids <- read.csv("PsA_controls_lipid_markers_reduced.csv", row.names = 1)
rownames(PsA_lipids)


# Create a data frame

PsA_lipids <- prepare_dataset( PsA_lipids )
names(PsA_lipids)
table(PsA_lipids$Target)
PsA_lipids$Target
PsA_lipids$Label


DatasetNames <- c( "PsA_lipids" )
MethodsList <- c( "PCA", "ICA", "IPCA", "LDA", "PLSDA", "PLSLDA", "MDS", "isomap", "tSNE", "Umap", "LLE", "NMF" )
ClusterAlgsList <- c( "none" ) #, "kmeans", "kmedoids", "ward.D2", "single", "average", "complete", "median", "centroid")
ClusterNumberMethodsList <- "orig"

# Ensuring PsA_lipids is available in the global environment
assign( "PsA_lipids", PsA_lipids, envir = .GlobalEnv )

# Create projections and plots
set.seed( 42 )
projectionsAndPlots <- perform_analysis( datasets = DatasetNames, projection_methods = MethodsList,
                                         clustering_methods = ClusterAlgsList, cluster_number_methods = ClusterNumberMethodsList,
                                         label_points = FALSE, highlight_misclassified = FALSE,
                                         #selected_cluster_metrics =  "cluster_accuracy",
                                         highlight_best_clustering = T,
                                         seed = 42, nProc = parallel::detectcores()-1)

# # Extract all plots into a single list
cPlotsPsA_lipids <- combine_all_plots( datasets = DatasetNames, projection_plots = projectionsAndPlots$projections_plots,
                                     projection_methods = MethodsList, clustering_methods = ClusterAlgsList,
                                     cluster_number_methods = ClusterNumberMethodsList )
print( cPlotsPsA_lipids )

cPlotsPsA_lipids_orig <- cowplot::plot_grid(
  s3d_grob,
  projectionsAndPlots$projections_plots$PsA_lipids$PCA.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$ICA.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$IPCA.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$LDA.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$PLSDA.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$MDS.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$isomap.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$tSNE.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$Umap.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$LLE.none[[1]],
  projectionsAndPlots$projections_plots$PsA_lipids$NMF.none[[1]],
  labels = "AUTO",
  ncol = 3
)

print( cPlotsPsA_lipids_orig )

ggsave(
  filename = paste0( "Combined_plot_orig_", DatasetNames[i], ".svg" ),
  plot = cPlotsPsA_lipids_orig, width = 12, height = 16
)


