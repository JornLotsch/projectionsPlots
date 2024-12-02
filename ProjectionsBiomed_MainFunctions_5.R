############### Functions ProjectionsBiomed new #########################################################
# Package loading function, avoiding installation during the function execution
# Constants for package names
REQUIRED_PACKAGES <- c(
  "pbmcapply", "deldir", "ggplot2", "ggrepel", "ggplotify",
  "grid", "cowplot", "cluster", "FactoMineR",
  "fastICA", "MASS", "Rtsne", "RDRToolbox", "mixOmics",
  "umap", "caret", "pracma", "combinat", "NbClust", "parallel"
)

# Load required R packages
load_required_packages <- function(packages) {
  missing_pkgs <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
  invisible(lapply(packages, library, character.only = TRUE))
}

# Load required packages at the start
load_required_packages(REQUIRED_PACKAGES)

# Main analysis function
perform_analysis <- function(dataset_names, method_list, clustering_algorithms = "none",
                             cluster_number_methods = "orig", label_points = TRUE, color_misclassified = FALSE,
                             dots_color_is = "Target", cells_color_is = "Cluster") {
  projections <- list()

  for (dataset_name in dataset_names) {
    message("Processing dataset: ", dataset_name)

    dataset <- retrieve_dataset(dataset_name)
    if (is.null(dataset)) next

    projection_results <- apply_projection_methods(dataset, method_list)
    plots_per_dataset <- process_and_create_plots(projections = projection_results, methods_list = method_list, cluster_algs = clustering_algorithms,
                                                  cluster_number_methods =  cluster_number_methods, dataset_name = dataset_name,
                                                  label_points = label_points, color_misclassified = color_misclassified, 
                                                  dots_color_is = dots_color_is, cells_color_is = cells_color_is)
    projections[[dataset_name]] <- unlist(plots_per_dataset, recursive = FALSE)
  }

  return(projections)
}

prepare_data <- function(X, Y = NULL) {
  # Convert input to a data frame if it is a matrix
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }

  # Validate input
  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("Input must be a data frame or matrix.")
  }

  # If Y is not provided, check for "Target" in input
  if (is.null(Y)) {
    if (!"Target" %in% colnames(X)) {
      stop("Column 'Target' is missing from the input.")
    }
    if (ncol(X) < 3) {
      stop("Input must have at least three columns when 'Target' is included.")
    }
    data <- X

  } else {
    # If Y is provided separately
    if (ncol(X) < 2) {
      stop("Input matrix must have at least two columns.")
    }

    if (length(Y) != nrow(X)) {
      stop("The length of Y must be equal to the number of rows in X.")
    }

    data <- as.data.frame(X)
    data$Target <- Y
  }

  # Convert Target to numeric and provide a mapping message
  if (!is.numeric(data$Target)) {
    class_mapping <- as.integer(as.factor(data$Target))
    levels <- levels(as.factor(data$Target))
    data$Target <- class_mapping
    class_map <- setNames(seq_along(levels), levels)

    mapping_message <- paste(sprintf("'%s' is now class %d", names(class_map), class_map), collapse = ", ")

    message("Class mapping: ", mapping_message)
  }

  return(data)
}
# Retrieve datasets
retrieve_dataset <- function(dataset_name) {
  if (!exists(dataset_name, envir = .GlobalEnv)) {
    warning("Dataset ", dataset_name, " not found in global environment. Skipping.")
    return(NULL)
  }

  ds <- get(dataset_name, envir = .GlobalEnv)
  if (!"Target" %in% names(ds)) {
    warning("Error: 'Target' column not found in dataset ", dataset_name, ". Skipping.")
    return(NULL)
  }

  return(ds)
}

# Apply projection methods
apply_projection_methods <- function(dataset, methods_list) {
  data_clean <- dataset[, !(names(dataset) %in% "Target")]
  target <- dataset$Target

  projections <- pbmcapply::pbmclapply(methods_list, function(method) {
    message("Applying projection method: ", method)
    result <- tryCatch({
      proj_result <- performProjection(X = data_clean, method = method, Target = target)
      if (!is.list(proj_result) || !"Projected" %in% names(proj_result))
        stop("Error: Invalid projection result structure.")
      return(proj_result)
    }, error = function(err) {
      message("Error in performProjection with method ", method, ": ", err$message)
      return(NULL)
    })
    return(result)
  }, mc.cores = min(length(methods_list), parallel::detectCores() - 1))

  names(projections) <- methods_list
  return(projections)
}

# Process projections and create plots
process_and_create_plots <- function(projections, methods_list, cluster_algs, cluster_number_methods,
                                     dataset_name, label_points, color_misclassified, dots_color_is, cells_color_is) {
  lapply(methods_list, function(method) {
    if (is.null(projections[[method]])) {
      warning("Projection for method ", method, " returned NULL. Skipping.")
      return(NULL)
    }

    proj_data <- as.data.frame(projections[[method]]$Projected)
    if (!is.data.frame(proj_data)) {
      warning("Error: Result for method ", method, " is not a data frame. Skipping.")
      return(NULL)
    }

    clustering_plots <- apply_clustering(projection_result = projections[[method]], projected_data = proj_data, cluster_algs = cluster_algs, 
                                         cluster_number_methods = cluster_number_methods,dataset_name =  dataset_name, 
                                         method_name =method, label_points = label_points, color_misclassified = color_misclassified, 
                                         dots_color_is = dots_color_is, cells_color_is = cells_color_is)
    names(clustering_plots) <- cluster_algs
    return(clustering_plots)
  })
}

# Apply clustering algorithms and plot results
apply_clustering <- function(projection_result, projected_data, cluster_algs, cluster_number_methods,
                             dataset_name, method_name, label_points, color_misclassified, 
                             dots_color_is, cells_color_is) {
  lapply(cluster_algs, function(cluster_alg) {
    tryCatch({
      if (is.null(projection_result$UniqueData$Target)) {
        stop("Required data is missing for clustering.")
      }
      clusters_results <- lapply(cluster_number_methods, function(cluster_number_method) {
        clusters <- performClustering(X = projected_data, Target = projection_result$UniqueData$Target,
                                      method = cluster_alg, ClusterNumberMethod = cluster_number_method)
        clusters <- renameClusters(trueCls = projection_result$UniqueData$Target, clusters)
        combined_result <- create_combined_result(projected_data, clusters, projection_result)
        plot <- create_projection_plot(combined_result = combined_result, dataset_name = dataset_name, method_name = method_name, 
                                         cluster_alg = cluster_alg, cluster_number_method = cluster_number_method,
                                       label_points = label_points, color_misclassified = color_misclassified, 
                                       dots_color_is = dots_color_is, cells_color_is = cells_color_is)
        return(plot)
      })

      names(clusters_results) <- cluster_number_methods
      return(clusters_results)
    }, error = function(err) {
      message("Error in clustering with algorithm ", cluster_alg, ": ", err$message)
      return(NULL)
    })
  })
}

# Combine results for plotting
create_combined_result <- function(projected_data, clusters, projection_result) {
  combined_result <- cbind.data.frame(
    projected_data,
    Cluster = clusters,
    Target = projection_result$UniqueData$Target,
    Label = projection_result$UniqueData$Label
  )
  combined_result$Misclassified <- ifelse(combined_result$Cluster == combined_result$Target, 0, 1)
  return(combined_result)
}

# Create projection plot
create_projection_plot <- function(combined_result, dataset_name, method_name, cluster_alg,
                                   cluster_number_method, label_points, color_misclassified, dots_color_is, cells_color_is) {
  plot <- plotVoronoiTargetProjection(
    X = combined_result[, 1:2],
    targets = combined_result$Target,
    clusters = combined_result$Cluster,
    labels = combined_result$Label,
    misclassified = combined_result$Misclassified,
    LabelPoints = label_points,
    ColorMisClassified = color_misclassified,
    dots_color_is = dots_color_is, 
    cells_color_is = cells_color_is
  )
  plot_title <- paste(dataset_name, ": ", method_name, "- ", cluster_alg)
  if (cluster_alg != "none") {
    plot_title <- paste(plot_title, "- ", cluster_number_method)
  }
  plot <- plot + labs(title = plot_title)

  return(plot)
}

# Function to perform projections
performProjection <- function(X, method = "PCA", scaleX = TRUE, Target = NULL, switchDims = FALSE, seed = 42, labels = NULL) {
  
  # Check and remove duplicates 
  if (!is.null(Target) & method != "PCA") {
    dupl <- which(duplicated(X))
    Target <- Target[-dupl]
  }
  
  if (!is.null(labels) & method != "PCA") {
    dupl <- which(duplicated(X))
    labels <- labels[-dupl]
  }
  
  if (method != "PCA") X <- X[!duplicated(X),]
  
  # Scale data
  if (scaleX & method != "IPCA") {
    X <- scale(X)
  }
  
  # Check if row names are provided or use labels argument
  if (is.null(labels)) {
    if (!is.null(rownames(X))) {
      labels <- rownames(X)
    } else {
      labels <- seq_len(nrow(X))
    }
  }
  
  # Combine X and Target into a single data frame
  if (!is.null(Target)) {
    df_unique <- data.frame(X, Target = Target, Label = labels)
  } else {
    df_unique <- data.frame(X, Label = labels)
  }
  
  # Extract X and Target separately after removing duplicates
  X_unique <- df_unique[, setdiff(names(df_unique), c("Target", "Label")), drop = FALSE]
  if (!is.null(Target)) {
    Target_unique <- df_unique$Target
  }
  
  # If PLSDA is selected, Target must be provided
  if (method == "PLSDA" && is.null(Target)) {
    stop("PLSDA selected but no Target provided. Projection cannot be performed.")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Process data based on the selected projection method
  switch(method,
         "none" = {
           # No projection, use the data as is
           Proj <- X_unique
         },
         "PCA" = {
           # Perform Principal Component Analysis
           pca1 <- FactoMineR::PCA(X_unique, graph = FALSE)
           ncomp <- max(2, sum(pca1$eig[, 1] > 1))
           
           Proj.pca <- FactoMineR::PCA(X_unique, graph = FALSE, ncp = ncomp)
           Proj <- Proj.pca$ind$coord
         },
         "ICA" = {
           # Independent Component Analysis
           Proj.ica <- fastICA::fastICA(X_unique,
                                        n.comp = 2,
                                        alg.typ = "parallel", fun = "logcosh", alpha = 1,
                                        method = "C", row.norm = FALSE, maxit = 200,
                                        tol = 0.0001, verbose = FALSE
           )
           Proj <- Proj.ica$S
         },
         "MDS" = {
           # Multidimensional Scaling
           d <- dist(X_unique)
           fit <- MASS::isoMDS(d, k = 2)
           Proj <- fit$points
         },
         "tSNE" = {
           # t-Distributed Stochastic Neighbor Embedding
           perpl <- ifelse(nrow(X_unique) / 3 - 1 / 3 < 30, nrow(X_unique) / 3 - 0.34, 30)
           Proj <- Rtsne::Rtsne(X_unique, perplexity = perpl)[[2]]
         },
         "isomap" = {
           # Perform Isomap projection
           Proj <- RDRToolbox::Isomap(as.matrix(X_unique), dims = 2, k = 5)$dim2
         },
         "PLSDA" = {
           # Partial Least Squares Discriminant Analysis
           res.plsda <- mixOmics::plsda(X_unique, factor(Target_unique))
           res.plsda.plot <- mixOmics::plotIndiv(res.plsda)
           Proj <- res.plsda.plot$df[c("x", "y")]
         },
         "IPCA" = {
           # Independent Principal Component Analysis with potential dimension switch
           res.ipca <- mixOmics::ipca(as.matrix(X_unique),
                                      ncomp = 2,
                                      mode = "deflation",
                                      fun = "logcosh",
                                      scale = F,
                                      w.init = NULL,
                                      max.iter = 500,
                                      tol = 1e-04
           )
           res.ipca.plot <- mixOmics::plotIndiv(res.ipca)
           Proj <- res.ipca.plot$df[c("x", "y")]
           # Conditionally switch dimensions if needed
           if (switchDims) {
             if (res.ipca$prop_expl_var$X[1] < res.ipca$prop_expl_var$X[2]) {
               Proj <- res.ipca.plot$df[c("y", "x")]
             }
           }
         },
         "Umap" = {
           # UMAP projection
           res.umap <- umap::umap(X_unique)
           Proj <- res.umap$layout
         }, {
           # No projection, use the data as is
           Proj <- X_unique
         }
  )
  
  # Return the projection and unique data frame
  return(list(Projected = Proj, UniqueData = df_unique))
}

# Function to determine the number of clusters
findOptimalClusters <- function(X, method, ClusterNumberMethod, seed, nProc = 1, Target = NULL) {
  # Load necessary libraries
  if (!requireNamespace("NbClust", quietly = TRUE)) {
    stop("Package 'NbClust' is required but not installed.")
  }
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    stop("Package 'FactoMineR' is required but not installed.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required but not installed.")
  }
  
  set.seed(seed)
  
  # Function to calculate mode
  calculateMode <- function(v) {
    if (length(v) == 0) return(NA)
    freqTable <- table(v)
    modes <- as.numeric(names(freqTable)[freqTable == max(freqTable)])
    if (length(modes) > 1) modes else modes[1]
  }
  
  nClusters <- switch(ClusterNumberMethod,
                      "NbClust" = {
                        # Some workarounds for methods lacking in NbClust
                        method <- switch(method,
                                         "kmedoids" = "kmeans",
                                         "diana" = "ward.D2",
                                         "HCPC" = "ward.D2",
                                         "hcpc" = "ward.D2",
                                         method)
                        
                        # Try if NbClust runs through
                        # Parallel computing as a workaround to silence NbClust outputs
                        # Only the first core's results will be read
                        NBres <- parallel::mclapply(c("all", "kl"), function(i) {
                          try(suppressMessages(NbClust::NbClust(data = X, diss = NULL, distance = "euclidean",
                                                                min.nc = 2, max.nc = 9, method = method, index = i)), silent = FALSE)
                        }, mc.cores = nProc)
                        
                        # Check if the first trail was successful, else run NbClust method by method                      
                        if (!inherits(NBres[[1]], "try-error")) {
                          length(unique(NBres[[1]]$Best.partition))
                        } else {
                          NbClustIndices <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot",
                                              "trcovw", "tracew", "friedman", "rubin", "cindex", "db",
                                              "silhouette", "duda", "pseudot2", "beale", "ratkowsky",
                                              "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma",
                                              "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")
                          nClustersTest <- parallel::mclapply(NbClustIndices, function(i) {
                            set.seed(seed)
                            nc <- try(suppressWarnings(NbClust::NbClust(data = X, diss = NULL, distance = "euclidean",
                                                                        min.nc = 2, max.nc = 9, method = method, index = i)), silent = FALSE)
                            if (!inherits(nc, "try-error")) {
                              length(unique(nc$Best.partition))
                            } else {
                              NA
                            }
                          }, mc.cores = nProc)
                          
                          nC <- na.omit(unlist(nClustersTest))
                          calculateMode(nC)
                        }
                      },
                      "HCPC" = {
                        res.hcpc <- FactoMineR::HCPC(X, consol = TRUE, method = "ward", nb.clust = -1, iter.max = 100, graph = FALSE, nstart = 10)
                        length(unique(res.hcpc$data.clust$clust))
                      },
                      if (!is.null(Target)) {
                        nClusters <- length(unique(Target))
                      } else {
                        nClusters <- 1
                      }
  )
 
   

if (is.na(nClusters)) {
  warning("Unable to determine the number of clusters.\nSetting it to the orginal classes (when given) or 1")
  if (!is.null(Target)) {
    nClusters <- length(unique(Target))
  } else {
    nClusters <- 1
  }
}

return(nClusters)
}

# Function to perform clustering
performClustering <- function(X, Target = NULL, method = "kmeans", seed = 42, ClusterNumberMethod = "orig") {
  
  # Check if input X is a non-empty data frame
  if (!is.data.frame(X) || nrow(X) == 0) stop("Input data must be a non-empty data frame.")
  
  # Define valid clustering methods
  valid_methods <- c("none", "kmeans", "kmedoids", "diana", "hcpc", "ward.D2", "single", "average", "median", "complete", "centroid")
  
  # Validate the selected method
  if (!(method %in% valid_methods)) stop(paste("Invalid method. Choose from:", paste(valid_methods, collapse = ", ")))
  
  # If Target is NULL, initialize it with default values
  if (is.null(Target)) Target <- rep(1, nrow(X))
  
  # Determine the number of clusters using the custom function
  
  if (method == "none") {
    nClusters <- length(unique(Target))
  } else {
    nClusters <- findOptimalClusters(X, method = method, ClusterNumberMethod = ClusterNumberMethod, seed = seed, nProc = nProc, Target = Target)
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Perform clustering based on the specified method
  cluster_result <- switch(method,
                           "none" = Target, # Return Target directly if method is "none"
                           "kmeans" = {
                             # Perform K-means clustering
                             kmeans(X, centers = nClusters, nstart = 100)$cluster
                           },
                           "kmedoids" = {
                             # Perform K-medoids clustering using PAM
                             cluster::pam(X, k = nClusters)$clustering
                           },
                           "diana" = {
                             # Perform divisive clustering
                             cluster::diana(X)$clustering
                           },
                           "hcpc" = {
                             # Perform hierarchical clustering and principal component analysis
                             res.pca <- FactoMineR::PCA(X, graph = FALSE)
                             FactoMineR::HCPC(res.pca, nb.clust = nClusters, graph = FALSE, nstart = 100, method = linkage)$data.clust$clust
                           }, {
                             # Perform other hierarchical clustering methods
                             clusterObject <- stats::hclust(dist(X), method = method)
                             as.numeric(cutree(clusterObject, k = nClusters))
                           })
  
  # Return the clustering result
  return(cluster_result)
}


# Function to align cluster names with liekely class labels of the prior classification
renameClusters <- function(trueCls, currentCls, K = 9) {
  
  # Helper function to reduce clusters
  ReduceClsToK <- function(Cls, K = 9) {
    uniqueCls <- unique(Cls)
    while (length(uniqueCls) > K) {
      counts <- table(Cls)
      to_merge <- names(sort(counts)[1:2])
      Cls[Cls %in% to_merge] <- to_merge[1]
      uniqueCls <- unique(Cls)
    }
    return(Cls)
  }
  
  # Preprocess input
  trueCls[!is.finite(trueCls)] <- 9999
  currentCls[!is.finite(currentCls)] <- 9999
  
  # Warnings
  if (length(unique(trueCls)) > 9) {
    warning("Too many clusters in PriorCls. Consider using cloud computing.")
  }
  if (length(unique(currentCls)) > K) {
    warning("Too many clusters in CurrentCls. Combining smallest clusters.")
    currentCls <- ReduceClsToK(currentCls, K)
  }
  
  # Normalize clusters
  trueCls <- as.numeric(factor(trueCls))
  currentCls <- as.numeric(factor(currentCls))
  
  # Get unique labels
  uniqueLabels <- sort(unique(c(trueCls, currentCls)))
  nLabels <- length(uniqueLabels)
  
  # Generate permutations
  permutations <- permn(nLabels)
  
  # Find best permutation
  bestAccuracy <- 0
  bestPermutation <- NULL
  
  for (perm in permutations) {
    newLabels <- perm[match(currentCls, seq_along(perm))]
    accuracy <- sum(trueCls == newLabels) / length(trueCls)
    if (accuracy > bestAccuracy) {
      bestAccuracy <- accuracy
      bestPermutation <- perm
    }
  }
  
  # Rename clusters
  renamedCls <- bestPermutation[match(currentCls, seq_along(bestPermutation))]
  
  return(renamedCls)
}

# Function to plot various versions of Voronoi cell presentation of projection and clustering results
plotVoronoiTargetProjection <- function(X, targets, clusters = NULL,
                                        labels = NULL, LabelPoints = FALSE, 
                                        misclassified = NULL, ColorMisClassified = FALSE, 
                                        dots_color_is = "Target", cells_color_is = "Cluster", 
                                        palette_target = NULL, palette_cluster = NULL,
                                        create_projection_plot) {
  
  # colorblind palette extended by random further colors
  cb_palette <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  
  shape_values <- c(17, 16, 15, 18, 2, 1, 0, 5, 7, 8, 9, 19, 11, 12)
  
  # Inital checks of arguments
  if (ncol(X) < 2) {
    stop("The input data `X` must have at least two columns.")
  }
  
  if (length(clusters) != nrow(X)) {
    stop("The length of `clusters` must match the number of rows in `X`.")
  }
  
  if (length(misclassified) != length(labels)) {
    stop("The length of `misclassified` must match the number of labels.")
  }
  
  if (!is.null(palette_target) && length(palette_target) < length(unique(targets))) {
    stop("The `palette_target` must provide enough color values for the number of target groups.")
  }
  
  if (!is.null(palette_cluster) && length(palette_cluster) < length(unique(clusters))) {
    stop("The `palette_cluster` must provide enough color values for the number of clusters.")
  }
  
  # Check if row names are provided or use labels argument
  if (is.null(labels)) {
    if (!is.null(rownames(X))) {
      labels <- rownames(X)
    } else {
      labels <- seq_len(nrow(X))
    }
  }
  
  # Check if the length of targets matches the number of rows in X
  if (nrow(X) != length(targets)) {
    stop("Die Länge von `targets` muss mit der Anzahl der Zeilen in `X` übereinstimmen.")
  }
  
  # Check if information on clusters is provided
  if (is.null(clusters)) clusters <- targets
  
  # Check if information on misclassified cases is provided
  if (is.null(misclassified)) misclassified <- rep(0, length(labels))
  
  # Data preparation
  plotData <- data.frame(Proj1 = X[, 1], Proj2 = X[, 2], Target = targets,
                         Clusters = clusters, Label = labels, Misclassified = misclassified)
  
  ifelse(dots_color_is == "Target", plotData$DotsInformation <- plotData$Target, plotData$DotsInformation <- plotData$Clusters)
  ifelse(cells_color_is == "Cluster", plotData$CellsInformation <- plotData$Clusters, plotData$CellsInformation <- plotData$Target)
  
  # Voronoi diagram computation
  voronoiData <- deldir::deldir(plotData$Proj1, plotData$Proj2)
  
  # Convert Voronoi tessellation to a data frame for plotting
  vor_polys <- deldir::tile.list(voronoiData)
  voronoi_df <- do.call(rbind, lapply(1:length(vor_polys), function(i) {
    data.frame(vor_polys[[i]]$x, vor_polys[[i]]$y, id = i)
  }))
  colnames(voronoi_df) <- c("x", "y", "id")
  
  # Define the default color palettes if none is provided
  if (is.null(palette_target)) {
    palette_target <- cb_palette
  }
  
  if (is.null(palette_cluster)) {
    palette_cluster <- cb_palette
  }
  
  # Create the plot with ggplot2
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = voronoi_df, ggplot2::aes(x = x, y = y, group = id, fill = as.factor(plotData$CellsInformation[id])),
                          alpha = 0.3, color = NA) +
    ggplot2::geom_point(data = plotData, ggplot2::aes(x = Proj1, y = Proj2,
                                                      color = as.factor(DotsInformation), shape = as.factor(DotsInformation))) +
    ggplot2::theme_light() +
    theme(legend.position = c(0.5, 0.08),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.background = element_rect(color = "transparent", fill = ggplot2::alpha("white", 0.2))) +
    ggplot2::labs(x = "Dim 1", y = "Dim 2", color = "Target", fill = "Cluster", shape = "Target") +
    scale_shape_manual(values = shape_values) +
    ggplot2::scale_fill_manual(values = palette_cluster) +
    ggplot2::scale_color_manual(values = palette_target) +
    guides(shape = guide_legend(override.aes = list(label = "")))
  
  
  # Conditional labeling of points
  if (LabelPoints) {
    plot <- plot +
      ggrepel::geom_text_repel(data = plotData,
                               ggplot2::aes(x = Proj1, y = Proj2, color = as.factor(DotsInformation),
                                            label = Label, fontface = 2), vjust = -1, size = 3, max.overlaps = Inf)
  }
  
  # Color misclassified cases in red
  if (sum(plotData$Target[plotData$Misclassified == 1]) > 0 & ColorMisClassified) {
    suppressMessages({
      plot <-
        plot +
        geom_text(
          aes(x = Inf, y = Inf,
              label = paste0(round(100 * sum(plotData$Misclassified) / length(plotData$Misclassified), 1), "% misclassified")),
          hjust = 1,
          vjust = 1
        )
    })
    
    q <- ggplot_build(plot)
    q$data[[2]]$colour[plotData$Misclassified == 1] <- "red"
    gtable <- ggplot2::ggplot_gtable(q)
    library(ggplotify)
    plot <- ggplotify::as.ggplot(function() grid::grid.draw(gtable))
    
  }
  
  return(plot)
}
