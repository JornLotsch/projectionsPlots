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

######################### Golfball equation ##############################################
sphere <- function( n, r = 1.0, surface_only = FALSE, center = cbind( 0.0, 0.0, 0.0 ) ) {
  phi <- runif( n, 0.0, 2.0 * pi )
  cos_theta <- runif( n, -1.0, 1.0 )
  sin_theta <- sqrt( ( 1.0 - cos_theta ) * ( 1.0 + cos_theta ) )
  radius <- r
  if ( surface_only == FALSE ) {
    radius <- r * runif( n, 0.0, 1.0 )^( 1.0 / 3.0 )
  }

  x <- radius * sin_theta * cos( phi )
  y <- radius * sin_theta * sin( phi )
  z <- radius * cos_theta

  # if radius is fixed, we could check it
  # rr = sqrt(x^2+y^2+z^2)
  # print(rr)

  ma <- cbind( x + center[1], y + center[2], z + center[3] )
  colnames( ma ) <- paste0( "X", 1:3 )
  return( ma )
}


#################################### Produce original data  ########################################################################

nPoints = 30
nPointsOcta = 30 #nPoints
nSpheres = 2
equalDensity = T
Cornerspheres = T


set.seed( 42 )

if ( equalDensity ) {
  if ( !Cornerspheres ) {
    assign( paste0( "nGolfball", nSpheres ), do.call( rbind, lapply( 1:nSpheres,
                                                                     function( k ) rsphere( n = nPoints * pi * k^2,
                                                                                            r = k,
                                                                                            surface_only = T ) ) ) )
    ClassesSpheres <- rep( 1:nSpheres, unlist( lapply( 1:nSpheres,
                                                       function( k ) nPoints * pi * k^2 ) ) )
  } else {

    Xcenters <- cbind(
      X1 = c( 0, 0, 2, -2, 2, -2, 2, -2, 2, -2 ),
      X2 = c( 0, 0, 2, 2, -2, -2, 2, 2, -2, -2 ),
      X3 = c( 0, 0, -2, -2, -2, -2, 2, 2, 2, 2 )
    )
    set.seed( 42 )
    assign( paste0( "nGolfball", nSpheres ), do.call( rbind, lapply( 1:10,
                                                                     function( k ) rsphere( n = ifelse( k < 3, nPoints * pi * k^2, nPointsOcta ),
                                                                                            r = ifelse( k < 3, k, 0.2 ),
                                                                                            surface_only = ifelse(k < 3, TRUE, TRUE),
                                                                                            center = Xcenters[k,] ) ) ) )
    ClassesSpheres <- rep( 1:10, unlist( lapply( 1:10,
                                                 function( k ) ifelse( k < 3, nPoints * pi * k^2, nPointsOcta ) ) ) )

  }
} else {
  assign( paste0( "nGolfball", nSpheres ), do.call( rbind, lapply( 1:nSpheres,
                                                                   function( k ) rsphere( n = nPoints,
                                                                                          r = k,
                                                                                          surface_only = T ) ) ) )
  ClassesSpheres <- rep( 1 + nSpheres, each = nPoints )
}


scatterplot3d::scatterplot3d( get( paste0( "nGolfball", nSpheres ) ), color = cb_palette[ClassesSpheres],
                              pch = c( 15:26 )[ClassesSpheres], xlab = "x", ylab = "y", zlab = "z", main = "Original" )

s3d_plot <- function( ) {
  scatterplot3d::scatterplot3d( get( paste0( "nGolfball", nSpheres ) ), color = cb_palette[ClassesSpheres],
                                pch = c( 15:26 )[ClassesSpheres], xlab = "x", ylab = "y", zlab = "z", main = "Original" )
}

s3d_grob <- ggplotify::as.grob( s3d_plot )

# Create a data frame

Golfball2 <- prepare_dataset( input_X = get( paste0( "nGolfball", nSpheres ) ), Target = ClassesSpheres )
#Golfball2 <- prepare_dataset(input_X = FCPS::Chainlink$Data, Target = FCPS::Chainlink$Cls)
dim( Golfball2 )

DatasetNames <- c( "Golfball2" )
MethodsList <- c( "PCA", "ICA", "IPCA", "LDA", "PLSDA", "PLSLDA", "MDS", "isomap", "tSNE", "Umap", "LLE", "NMF" )
ClusterAlgsList <- c( "none" ) #, "kmeans", "kmedoids", "ward.D2", "single", "average", "complete", "median", "centroid")
ClusterNumberMethodsList <- "orig"

# Ensuring Golfball2 is available in the global environment
assign( "Golfball2", Golfball2, envir = .GlobalEnv )

# Create projections and plots
set.seed( 42 )
projectionsAndPlots <- perform_analysis( datasets = DatasetNames, projection_methods = MethodsList,
                                         clustering_methods = ClusterAlgsList, cluster_number_methods = ClusterNumberMethodsList,
                                         label_points = FALSE, highlight_misclassified = FALSE,
                                         #selected_cluster_metrics =  "cluster_accuracy",
                                         highlight_best_clustering = F,
                                         seed = 42, nProc = 4 )

# # Extract all plots into a single list
cPlotsGolfball <- combine_all_plots( datasets = DatasetNames, projection_plots = projectionsAndPlots$projections_plots,
                                     projection_methods = MethodsList, clustering_methods = ClusterAlgsList,
                                     cluster_number_methods = ClusterNumberMethodsList )
print( cPlotsGolfball )

cPlotsGolfball_orig <- cowplot::plot_grid(
  s3d_grob,
  projectionsAndPlots$projections_plots$Golfball2$PCA.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$ICA.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$IPCA.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$LDA.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$PLSDA.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$MDS.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$isomap.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$tSNE.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$Umap.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$LLE.none[[1]],
  projectionsAndPlots$projections_plots$Golfball2$NMF.none[[1]],
  labels = "AUTO",
  ncol = 3
)

print( cPlotsGolfball_orig )

ggsave(
  filename = paste0( "Combined_plot_orig_Octa_non_dense_", DatasetNames[i], ".svg" ),
  plot = cPlotsGolfball_orig, width = 12, height = 16
)


