# IrisData = cbind.data.frame(Target = as.integer(as.factor(iris$Species)), iris[, 1:4])
# 
# IrisData_projected <- performProjection(X = within(IrisData, rm(Target)), method = "PCA", Target = IrisData$Target)
# 
# plotVoronoiTargetProjection(X = IrisData_projected$Projected, targets = IrisData_projected$UniqueData$Target, LabelPoints = T)

# Set the working directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("/home/joern/.Datenplatte/Joerns Dateien/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/")

# Function to load the main functions file
loadMainFunctions <- function() {
  mainFunctionsFile <- "ProjectionsBiomed_MainFunctions_6.R"

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

IrisData <- prepare_dataset(input_X = iris[, 1:4], Target = iris$Species)
# IrisData <- prepare_data(X = within(mtcars, rm(gear)), Y = mtcars$gear)
# IrisData <- WineData
#IrisData <- cbind.data.frame(Target = as.integer(as.factor(iris$Species)), iris[, 1:4])
DatasetNames <- c("IrisData")
MethodsList <- c("PCA", "ICA", "IPCA", "LDA", "PLSDA", "MDS", "isomap", "tSNE", "Umap", "LLE", "NMF", "autoencoder")
ClusterAlgsList <- c("none")#, "kmeans", "kmedoids", "ward.D2", "single", "average", "complete", "median", "centroid")
ClusterNumberMethodsList <- "orig"

# Ensuring IrisData is available in the global environment
assign("IrisData", IrisData, envir = .GlobalEnv)

# Create projections and plots
set.seed(42)
projectionsAndPlots <- perform_analysis(datasets = DatasetNames, projection_methods = MethodsList, 
                                        clustering_methods =  ClusterAlgsList, cluster_number_methods = ClusterNumberMethodsList,
                                        label_points = FALSE, highlight_misclassified = FALSE,
                                        #selected_cluster_metrics =  "cluster_accuracy",
                                        highlight_best_clustering = T, 
                                        seed = 42, 
                                        nProc = parallel::detectcores()-1)

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
cPlots <- combine_all_plots(datasets = DatasetNames, projection_plots = projectionsAndPlots$projections_plots, 
                  projection_methods = MethodsList, clustering_methods = ClusterAlgsList, 
                  cluster_number_methods = ClusterNumberMethodsList)
print(cPlots)
projectionsAndPlots$cluster_quality_results$best_combination


for (i in seq_along(DatasetNames)) {
ggsave(
  filename = paste0("Combined_plot_", DatasetNames[i], ".svg"),
  plot = cPlots[[i]], width = 12, height = 16
)
}




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
