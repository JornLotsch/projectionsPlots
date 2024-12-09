# dfData = cbind.data.frame(Target = as.integer(as.factor(iris$Species)), iris[, 1:4])
# 
# dfData_projected <- performProjection(X = within(dfData, rm(Target)), method = "PCA", Target = dfData$Target)
# 
# plotVoronoiTargetProjection(X = dfData_projected$Projected, targets = dfData_projected$UniqueData$Target, LabelPoints = T)

# Set the working directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("/home/joern/.Datenplatte/Joerns Dateien/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/")

# Function to load the main functions file
loadMainFunctions <- function() {
  mainFunctionsFile <- "ProjectionsBiomed_MainFunctions_5.R"

  if (file.exists(mainFunctionsFile)) {
    message("Loading main functions from ", mainFunctionsFile)
    source(mainFunctionsFile)
  } else {
    stop("The file ", mainFunctionsFile, " does not exist in the current directory.")
  }
}

# Call the function to load the main functions
loadMainFunctions()

# Example call
# Create a data frame

dfData <- prepare_dataset(input_X = iris[, 1:4], Target = iris$Species)
# dfData <- prepare_data(X = within(mtcars, rm(gear)), Y = mtcars$gear)
# dfData <- WineData
#dfData <- cbind.data.frame(Target = as.integer(as.factor(iris$Species)), iris[, 1:4])
DatasetNames <- c("dfData")
MethodsList <- c("PCA", "ICA", "IPCA", "PLSDA", "MDS", "isomap", "tSNE", "Umap", "LLE", "NMF")
ClusterAlgsList <- c("none")#, "kmeans", "kmedoids", "ward.D2", "single", "average", "complete", "median", "centroid")
ClusterNumberMethodsList <- "orig"

# Ensuring dfData is available in the global environment
assign("dfData", dfData, envir = .GlobalEnv)

# Create projections and plots
set.seed(42)
projectionsAndPlots <- perform_analysis(datasets = DatasetNames, projection_methods = MethodsList, 
                                        clustering_methods =  ClusterAlgsList, cluster_number_methods = ClusterNumberMethodsList,
                                        label_points = FALSE, highlight_misclassified = FALSE,
                                        #selected_cluster_metrics =  "cluster_accuracy",
                                        highlight_best_clustering = T, 
                                        seed = 42)

# # Extract all plots into a single list
# allPlots <- unlist(lapply(projectionsAndPlots, function(projectionList) {
#   if (is.null(projectionList)) return(NULL)
#   unlist(projectionList, recursive = FALSE, use.names = FALSE)
# }), recursive = FALSE, use.names = FALSE)
# 
# # Combine all plots into a single grid
# combinedPlot <- cowplot::plot_grid(plotlist = allPlots, ncol = length(ClusterAlgsList) * length(ClusterNumberMethodsList), nrow = length(MethodsList))
# 
# # Print the combined plot
# print(combinedPlot)


# Example usage:
combine_all_plots(datasets = DatasetNames, projection_plots = projectionsAndPlots$projections_plots, 
                  projection_methods = MethodsList, clustering_methods = ClusterAlgsList, 
                  cluster_number_methods = ClusterNumberMethodsList)
projectionsAndPlots$cluster_quality_results$best_combination


# plot_default_ggplot_colors <- function(n = 6) {
# # Generate default ggplot2 colors
# default_colors <- scales::hue_pal()(n)
# # Create a data frame for plotting
# color_df <- data.frame(
# x = 1:n,
# y = rep(1, n),
# fill = factor(1:n)
# )
# # Create the plot
# ggplot(color_df, aes(x = x, y = y, fill = fill)) +
# geom_tile(color = "black", size = 1) +
# scale_fill_manual(values = default_colors) +
# geom_text(aes(label = default_colors), vjust = 3) +
# theme_minimal() +
# theme(axis.text = element_blank(),
# axis.title = element_blank(),
# panel.grid = element_blank(),
# legend.position = "none") +
# labs(title = "Default ggplot2 colors")
# }
# 
# plot_default_ggplot_colors(11)
# plot_default_ggplot_colors(111)
