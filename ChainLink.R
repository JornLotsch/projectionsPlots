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

######################### Chainlink equation ##############################################

set.seed(412)
# Parameters for the rings
n_points <- 200  # Number of points per ring
radius <- 5      # Radius of the rings
noise_level <- 0.1 # Add small random noise to points

# Shift parameters (half the radius)
x_shift <- 0 #radius / 2
y_shift <- 0 #radius / 2
z_shift <- radius

# Generate the horizontal ring (shifted in the X-Z plane)
theta1 <- seq( 0, 2 * pi, length.out = n_points )
x1 <- radius * cos( theta1 ) + x_shift
y1 <- rep( 0, n_points ) + y_shift
z1 <- radius * sin( theta1 ) + z_shift

# Add some noise
x1 <- x1 + rnorm( n_points, sd = noise_level )
y1 <- y1 + rnorm( n_points, sd = noise_level )
z1 <- z1 + rnorm( n_points, sd = noise_level )

# Generate the vertical ring (in the Y-Z plane)
theta2 <- seq( 0, 2 * pi, length.out = n_points )
x2 <- rep( 0, n_points )
y2 <- radius * cos( theta2 )
z2 <- radius * sin( theta2 )

# Add some noise
x2 <- x2 + rnorm( n_points, sd = noise_level )
y2 <- y2 + rnorm( n_points, sd = noise_level )
z2 <- z2 + rnorm( n_points, sd = noise_level )

# Combine both rings into a single dataset
Chainlink2 <- data.frame(
  x = c( x1, x2 ),
  y = c( y1, y2 ),
  z = c( z1, z2 ),
  Target = c( rep( 1, n_points ), rep( 2, n_points ) )
)


s3d_plot_chain <- function( ) {
  scatterplot3d::scatterplot3d( Chainlink2[, 1:3], color = ggthemes::colorblind_pal( )( 8 )[Chainlink2$Target],
                                pch = c( 15, 19 )[Chainlink2$Target], xlab = "x", ylab = "y", zlab = "z", main = "Original" )
}

s3d_grob_chain <- ggplotify::as.grob( s3d_plot_chain )

# Create a data frame

ChainLink2 <- prepare_dataset( input_X = Chainlink2[, 1:3], Target = Chainlink2$Target )
dim( ChainLink2 )

DatasetNames <- c( "ChainLink2" )
MethodsList <- c( "PCA", "ICA", "IPCA", "LDA", "PLSDA", "PLSLDA", "MDS", "isomap", "tSNE", "Umap", "LLE", "NMF" )
ClusterAlgsList <- c( "none" ) #, "kmeans", "kmedoids", "ward.D2", "single", "average", "complete", "median", "centroid")
ClusterNumberMethodsList <- "orig"

# Ensuring ChainLink is available in the global environment
assign( "ChainLink2", ChainLink2, envir = .GlobalEnv )

# Create projections and plots
set.seed( 42 )
projectionsAndPlots <- perform_analysis( datasets = DatasetNames, projection_methods = MethodsList,
                                         clustering_methods = ClusterAlgsList, cluster_number_methods = ClusterNumberMethodsList,
                                         label_points = FALSE, highlight_misclassified = FALSE,
                                         #selected_cluster_metrics =  "cluster_accuracy",
                                         highlight_best_clustering = T,
                                         seed = 42 )

# # Extract all plots into a single list
cPlotsChainLink <- combine_all_plots( datasets = DatasetNames, projection_plots = projectionsAndPlots$projections_plots,
                                      projection_methods = MethodsList, clustering_methods = ClusterAlgsList,
                                      cluster_number_methods = ClusterNumberMethodsList )
print( cPlotsChainLink )

cPlotsChainLink_orig<-cowplot::plot_grid(
  s3d_grob_chain,
  projectionsAndPlots$projections_plots$ChainLink2$PCA.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$ICA.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$IPCA.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$LDA.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$PLSDA.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$MDS.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$isomap.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$tSNE.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$Umap.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$LLE.none[[1]],
  projectionsAndPlots$projections_plots$ChainLink2$NMF.none[[1]],
  labels="AUTO",
  ncol=3
)

print( cPlotsChainLink_orig )

ggsave(
  filename = paste0( "Combined_plot_orig_", DatasetNames[i], ".svg" ),
  plot = cPlotsChainLink_orig, width = 12, height = 16
)


