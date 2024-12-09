de ############### Functions ProjectionsBiomed new #########################################################

############### Libraries ###############

# Package loading function, avoiding installation during the function execution
# Constants for required package names
REQUIRED_PACKAGES <- c(
  "pbmcapply", "deldir", "ggplot2", "ggrepel", "ggplotify",
  "grid", "cowplot", "cluster", "FactoMineR", "ABCanalysis",
  "fastICA", "MASS", "Rtsne", "RDRToolbox", "mixOmics",
  "umap", "caret", "pracma", "combinat", "NbClust", "parallel",
  "dplyr", "ape", "combinat", "NMF", "RANN"
)

############### Main function ###############

# Function to load required R packages
load_required_packages <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  invisible(lapply(packages, library, character.only = TRUE))
}

# Load required packages at the start
load_required_packages(REQUIRED_PACKAGES)

# Function to perform analysis
perform_analysis <- function(datasets, projection_methods = "PCA",
                             clustering_methods = "none",
                             cluster_number_methods = "orig",
                             selected_cluster_metrics = c("cluster_accuracy", "Silhouette_index", "Dunn_index",
                                                 "Rand_index", "DaviesBouldin_index", "dbcv_index"),
                             distance_metric = "euclidean",
                             highlight_best_clustering = FALSE,
                             method_for_hcpc = "ward",
                             label_points = TRUE,
                             highlight_misclassified = FALSE,
                             points_colored_for = "Target",
                             cells_colored_for = "Cluster",
                             palette_target = NULL,
                             palette_cluster = NULL,
                             seed = 42) {

  projections_plots <- list() # Initialize the projections_plots variable
  all_projection_results <- list() # To store all projection results
  all_clustering_results <- list() # To store all clustering results
  cluster_quality_results <- list() # To store the cluster scores

  for (dataset_name in datasets) {
    message("Processing dataset: ", dataset_name)
    dataset <- retrieve_dataset(dataset_name)
    if (is.null(dataset)) next

    local_projection_methods <- projection_methods # Local copy of projection methods
    if (length(unique(dataset$Target)) == 1 && "PLSDA" %in% local_projection_methods) {
      message("'Target' has only one level. Removing PLSDA from projection methods list.")
      local_projection_methods <- setdiff(local_projection_methods, "PLSDA")
      if (length(local_projection_methods) < 1) {
        stop("No available projection methods for dataset: ", dataset_name)
      }
    }

    message("Projection data")
    projection_results <- generate_projections(data = dataset, projection_methods = local_projection_methods, seed = seed)
    all_projection_results[[dataset_name]] <- projection_results

    message("Clustering data")
    clustering_results <- apply_clustering(projection_results = projection_results,
                                           clustering_methods = clustering_methods,
                                           cluster_number_methods = cluster_number_methods,
                                           method_for_hcpc = method_for_hcpc,
                                           distance_metric = distance_metric,
                                           seed = seed)
    all_clustering_results[[dataset_name]] <- clustering_results

    message("Calculating cluster quality or stability indices")
    cluster_indices <- calculate_cluster_scores_df(clustering_results = clustering_results,
                                                   projection_methods = local_projection_methods,
                                                   clustering_methods = clustering_methods,
                                                   cluster_number_methods = cluster_number_methods,
                                                   distance_metric = distance_metric)

    cluster_quality_results <- process_cluster_scores(cluster_scores_df = cluster_indices, selected_cluster_metrics = selected_cluster_metrics)

    message("Creating plots")
    plots <- create_dataset_plots(clustering_results = clustering_results,
                                  projection_methods = local_projection_methods,
                                  clustering_methods = clustering_methods,
                                  cluster_number_methods = cluster_number_methods,
                                  dataset_name = dataset_name,
                                  label_points = label_points,
                                  highlight_misclassified = highlight_misclassified,
                                  points_colored_for = points_colored_for,
                                  cells_colored_for = cells_colored_for,
                                  palette_target = palette_target,
                                  palette_cluster = palette_cluster)

    if (highlight_best_clustering) {

      if (is.data.frame(cluster_quality_results$best_combination) &&
        nrow(cluster_quality_results$best_combination) > 0) {

        for (i in 1:nrow(cluster_quality_results$best_combination)) {
          # Library and color palette code is commented out
          # Modify the plot with suppressed messages and warnings
          plots[[cluster_quality_results$best_combination[i, "projection_method"]]][[cluster_quality_results$best_combination[i, "clustering_method"]]][[cluster_quality_results$best_combination[i, "cluster_number_method"]]] <-
            suppressMessages(
              suppressWarnings(
                plots[[cluster_quality_results$best_combination[i, "projection_method"]]][[cluster_quality_results$best_combination[i, "clustering_method"]]][[cluster_quality_results$best_combination[i, "cluster_number_method"]]] +
                  scale_color_hue() +
                  scale_fill_hue()
              )
            )
        }
      } else {
        message("The `best_combination` is not a valid data frame or it has no rows.")
      }
    }

    projections_plots[[dataset_name]] <- unlist(plots, recursive = FALSE)
  }
  return(list(projections_plots = projections_plots,
              projection_results = all_projection_results,
              clustering_results = all_clustering_results,
              cluster_quality_results = cluster_quality_results)
  )
}


############### Functions for data preparations ###############

# Function to prepare datasets
prepare_dataset <- function(input_X, Target = NULL, Label = NULL) {
  input_X <- if (is.matrix(input_X)) as.data.frame(input_X) else input_X
  if (!is.data.frame(input_X) && !is.matrix(input_X)) stop("Input must be a data frame or matrix.")

  # Use existing "Target" column if Target is NULL and "Target" exists in input_X
  output_Y <- if (is.null(Target)) {
    if (!"Target" %in% colnames(input_X)) {
      message("'Target' is missing, creating 'Target' = 1.")
      rep(1, nrow(input_X))
    } else {
      input_X$Target
    }
  } else {
    Target
  }

  # Use existing "Label" column if Labels is NULL and "Target" exists in input_X
  output_L <- if (is.null(Label)) {
    if (!"Label" %in% colnames(input_X)) {
      message("Taking row names as case labels.")
      rownames(input_X)
    } else {
      message("Taking 'Label' column as case labels.")
      input_X$Label
    }
  } else {
    if (length(Label) != nrow(input_X)) {
      message("Length of 'Label' does not match numbr of rows in input matrix\nTaking row names as case labels.")
      rownames(input_X)
    } else {
      Label
    }
  }


  data_frame <- as.data.frame(input_X)
  if (ncol(data_frame) < 3) stop("Matrix needs at least three columns including 'Target'.")

  # Ensure "Target" is numeric after checking availability
  data_frame$Target <- output_Y
  if (!is.numeric(data_frame$Target)) {
    data_frame$Target <- as.numeric(factor(data_frame$Target))
  }

  # Add "Label" is numeric after checking availability
  data_frame$Label <- output_L

  return(data_frame)
}


# Function to retrieve datasets
retrieve_dataset <- function(dataset_name) {
  if (!exists(dataset_name, envir = .GlobalEnv)) {
    warning("Dataset ", dataset_name, " not found. Skipping.")
    return(NULL)
  }

  data_set <- get(dataset_name, envir = .GlobalEnv)
  if (!is.data.frame(data_set)) {
    warning("Dataset ", dataset_name, " is not a data frame. Skipping.")
    return(NULL)
  }

  if (!"Target" %in% names(data_set)) {
    warning("No 'Target' column found in dataset ", dataset_name, ". Skipping.")
    return(NULL)
  }

  return(data_set)
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


############### Functions for data projection ###############

# Function to generate projections
generate_projections <- function(data, projection_methods, seed = 42) {
  data_clean <- data[, !(names(data) %in% c("Target", "Label"))]
  target <- data$Target
  labels <- data$Label

  projections <- pbmcapply::pbmclapply(projection_methods, function(projecion_method) {
    message("Applying projection: ", projecion_method)
    tryCatch({
      projection <- performProjection(X = data_clean, projecion_method = projecion_method, Target = target, labels = labels, seed = seed)
      if (!is.list(projection) || !"Projected" %in% names(projection))
        stop("Invalid projection result.")
      return(projection)
    }, error = function(err) {
      message("Error with method ", projecion_method, ": ", err$message)
      return(NULL)
    })
  }, mc.cores = min(length(projection_methods), parallel::detectCores() - 1))

  names(projections) <- projection_methods
  return(projections)
}


# Function to perform projections
performProjection <-
  function(X, projecion_method = "PCA", scaleX = TRUE, Target = NULL, switchDims = FALSE, labels = NULL, seed = 42) {

    # Define valid  methods
    valid_projection_methods <- c("none", "PCA", "ICA", "MDS", "tSNE", "isomap", "PLSDA", "IPCA", "Umap", "NMF", "LLE")

    # Validate the selected method
    if (!(projecion_method %in% valid_projection_methods))
      stop(paste("Invalid projection method. Choose from:", paste(valid_projection_methods, collapse = ", ")))

    # Helper function for capturing integer(0)
    is.integer0 <- function(x) {
      is.integer(x) && length(x) == 0L
    }

    # Check and remove duplicates 
    dupl <- which(duplicated(X))

    if (!is.null(Target) && projecion_method != "PCA" && !is.integer0(dupl)) {
      Target <- Target[-dupl]
    }

    if (!is.null(labels) && projecion_method != "PCA" && !is.integer0(dupl)) {
      labels <- labels[-dupl]
    }

    if (projecion_method != "PCA") X <- X[!duplicated(X),]

    # Scale data
    if (scaleX & !projecion_method %in% c("IPCA", "NMF")) {
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
    if (projecion_method == "PLSDA" && is.null(Target)) {
      stop("PLSDA selected but no Target provided. Projection cannot be performed.")
    }

    # Process data based on the selected projection method
    set.seed(seed)
    switch(projecion_method,
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
      # Helper fucntion to sort ICA componets for importance
      sort_ica_components <- function(S, A) {
        # S: Source matrix (ICA components)
        # A: Mixing matrix
        # Calculate the variance explained by each componentQ
        var_explained <- apply(S, 2, var)

        # Sort the components based on variance explained
        order_components <- order(var_explained, decreasing = TRUE)

        # Reorder the columns of S and A
        S_sorted <- S[, order_components]
        A_sorted <- A[, order_components]

        return(list(S = S_sorted, A = A_sorted))
      }

      # Independent Component Analysis
      # Determine a suitable n.com value from n PCs in a PCA
      pca1 <- FactoMineR::PCA(X_unique, graph = FALSE)
      ncomp <- max(2, sum(pca1$eig[, 1] > 1))

      # Run ICA
      Proj.ica <- fastICA::fastICA(X_unique,
                                          n.comp = ncomp,
                                          alg.typ = "parallel", fun = "logcosh", alpha = 1,
                                          method = "C", row.norm = FALSE, maxit = 1000,
                                          tol = 1e-04, verbose = FALSE
             )

      # Sort ICA components so that the first two are most relevant
      Proj.ica_sorted <- sort_ica_components(S = Proj.ica$S, A = Proj.ica$A)
      Proj <- Proj.ica_sorted$S
    },
           "MDS" = {

      k_estimate <- floor(nrow(X_unique) / 2)
      # Perform MDS with maximum dimensions
      # Multidimensional Scaling
      d <- dist(X_unique)
      mds1 <- MASS::isoMDS(d, k = k_estimate)

      # Calculate variance explained by each dimension
      var_explained <- apply(mds1$points, 2, var)

      # Extract the most profitable components
      k_important <- max(2, length(ABCanalysis::ABCanalysis(var_explained)$Aind))

      # Repeat MDS with a better estimate for k
      mds.fit <- MASS::isoMDS(d, k = k_important)
      Proj <- mds.fit$points
    },
           "tSNE" = {
      # t-Distributed Stochastic Neighbor Embedding
      perpl <- ifelse(nrow(X_unique) / 3 - 1 / 3 < 30, nrow(X_unique) / 3 - 0.34, 30)
      Proj <- Rtsne::Rtsne(X_unique, perplexity = perpl)[[2]]
    },
           "isomap" = {
      # Perform Isomap projection with 10 dimensions
      isomapfit1 <- RDRToolbox::Isomap(as.matrix(X_unique), dims = 1:10, k = 5)

      # Determine how many dimensions are needed for explaining 80% of the variance
      variance_explained <- apply(isomapfit1[[length(isomapfit1)]], 2, var)
      cumulative_variance <- cumsum(variance_explained) / sum(variance_explained)

      # Determine the dimensionality of the final Isomap
      good_dims <- max(2, which(cumulative_variance >= 0.8)[1])

      # Repeat Isomap with a better estimate for dims
      Proj <- RDRToolbox::Isomap(as.matrix(X_unique), dims = good_dims, k = 5)$dim2
    },
           "PLSDA" = {
      # Partial Least Squares Discriminant Analysis
      res.plsda <- mixOmics::plsda(X_unique, factor(Target_unique))
      res.plsda.plot <- mixOmics::plotIndiv(res.plsda)
      Proj <- res.plsda.plot$df[c("x", "y")]
    },
           "IPCA" = {
      # Independent Principal Component Analysis with potential dimension switch
      # Determine a suitable n.com value from n PCs in a PCA
      pca1 <- FactoMineR::PCA(X_unique, graph = FALSE)
      ncomp <- max(2, sum(pca1$eig[, 1] > 1))

      # Run IPCA with that ncomp value
      res.ipca <- mixOmics::ipca(as.matrix(X_unique),
                                        ncomp = ncomp,
                                        mode = "deflation",
                                        fun = "logcosh",
                                        scale = F,
                                        w.init = NULL,
                                        max.iter = 1000,
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
    },
    "NMF" = {
      # Function to determine optimal number of components
      determine_k <- function(data, max_k = 10, threshold = 0.9) {
        best_k <- 1
        best_explained_variance <- 0

        for (k in 2:max_k) {
          nmf_result <- NMF::nmf(data, k,) # Verbose output for progress
          # Calculate explained variance
          basis_matrix <- basis(nmf_result)
          variance_explained <- sum(basis_matrix ^ 2) / sum(data ^ 2)

          if (variance_explained > threshold) {
            best_k <- k
            break
          }

          if (variance_explained > best_explained_variance) {
            best_k <- k
            best_explained_variance <- variance_explained
          }
        }

        return(best_k)
      }

      # Function for Sequential NMF with non-negativity check
      sequential_nmf <- function(data, max_k = 10, threshold = 0.9) {
        # Check for non-positive values
        if (any(data <= 0)) {
          message("Input matrix contains non-positive values. Shifting all values to be positive.")
          min_value <- min(data)
          # Shift the data to make all values positive
          data <- data - min_value + 1e-10 # Small constant to ensure all values are positive
        }

        k <- determine_k(data, max_k, threshold)

        nmf_result <- NMF::nmf(data, k)

        # Use coefficient matrix directly; no need to transpose for this case.
        coef_matrix <- basis(nmf_result)

        # Each row corresponds to a sample and columns are the components.
        projected_data <- coef_matrix

        return(list(projected_data = projected_data, basis_matrix = basis(nmf_result)))
      }

      nmf.res <- suppressMessages(sequential_nmf(data = as.matrix(X_unique)))

      Proj <- nmf.res$projected_data

    },
    "LLE" = {
      # Local linear embedding
      # lle_optimized <- function(data, max_dim = 10, max_k = 30) {
      #   # Ensure data is a numeric matrix
      #   X <- as.matrix(data)
      #   X <- apply(X, 2, as.numeric)
      #   
      #   lle <- function(X, k, d) {
      #     n <- nrow(X)
      #     p <- ncol(X)
      #     
      #     # Find k nearest neighbors
      #     nn <- RANN::nn2(X, k = k + 1)
      #     neighbors <- nn$nn.idx[, -1]
      #     
      #     # Compute weights
      #     W <- matrix(0, n, n)
      #     for (i in 1:n) {
      #       z <- X[neighbors[i,], ] - matrix(X[i,], k, p, byrow = TRUE)
      #       C <- tcrossprod(z)
      #       C <- C + diag(k) * .Machine$double.eps * sum(abs(C))
      #       w <- solve(C, rep(1, k))
      #       w <- w / sum(w)
      #       W[i, neighbors[i,]] <- w
      #     }
      #     
      #     # Compute embedding
      #     M <- diag(n) - W - t(W) + t(W) %*% W
      #     eig <- eigen(M)
      #     Y <- eig$vectors[, (n-d):(n-1)]
      #     
      #     return(Y)
      #   }
      #   
      #   optimize_params <- function(X, max_k, max_d) {
      #     n <- nrow(X)
      #     errors <- matrix(Inf, max_k, max_d)
      #     
      #     for (k in 2:max_k) {
      #       for (d in 1:min(max_d, k-1)) {
      #         tryCatch({
      #           Y <- lle(X, k, d)
      #           errors[k,d] <- mean((X - Y %*% t(Y) %*% X)^2)
      #         }, error = function(e) {
      #           message("Error for k=", k, ", d=", d, ": ", e$message)
      #         })
      #       }
      #     }
      #     
      #     best <- which(errors == min(errors, na.rm = TRUE), arr.ind = TRUE)[1,]
      #     return(list(k = best[1], d = best[2]))
      #   }
      #   
      #   params <- optimize_params(X, max_k, max_dim)
      #   projection <- lle(X, params$k, params$d)
      #   
      #   var_order <- order(apply(projection, 2, var), decreasing = TRUE)
      #   projection <- projection[, var_order]
      #   
      #   return(list(projection = projection, k = params$k, d = params$d))
      # }
      # 
      library(RANN)
      library(Rtsne)
      
      lle_optimized <- function(data, max_dim = 10, max_k = 30, perplexity = 30) {
        # Ensure data is a numeric matrix
        X <- as.matrix(data)
        X <- scale(X)  # Standardize the data
        
        lle <- function(X, k, d, reg = 1e-3) {
          n <- nrow(X)
          p <- ncol(X)
          
          # Find k nearest neighbors
          nn <- RANN::nn2(X, k = k + 1)
          neighbors <- nn$nn.idx[, -1]
          
          # Compute weights
          W <- matrix(0, n, n)
          for (i in 1:n) {
            z <- X[neighbors[i,], , drop = FALSE] - matrix(X[i,], k, p, byrow = TRUE)
            C <- tcrossprod(z)
            C <- C + diag(k) * reg * sum(diag(C))  # Regularization
            w <- solve(C, rep(1, k))
            w <- w / sum(w)
            W[i, neighbors[i,]] <- w
          }
          
          # Compute embedding
          M <- diag(n) - W - t(W) + t(W) %*% W
          eig <- eigen(M, symmetric = TRUE)
          Y <- eig$vectors[, (n-d):(n-1), drop = FALSE]
          
          return(Y)
        }
        
        optimize_params <- function(X, max_k, max_d, perplexity) {
          n <- nrow(X)
          scores <- matrix(Inf, max_k, max_d)
          
          # Compute t-SNE embedding for reference
          tsne_emb <- Rtsne(X, dims = 2, perplexity = perplexity, verbose = FALSE)$Y
          
          for (k in seq(5, max_k, by = 5)) {  # Step by 5 for efficiency
            for (d in 1:min(max_d, k-1)) {
              tryCatch({
                Y <- lle(X, k, d)
                
                # Ensure Y has at least 2 columns for KL divergence calculation
                if (ncol(Y) < 2) {
                  Y <- cbind(Y, rep(0, nrow(Y)))
                }
                
                # Compute KL divergence between t-SNE and LLE embeddings
                kl_div <- KL_divergence(tsne_emb, Y[,1:2, drop = FALSE])
                
                # Compute trustworthiness
                trust <- trustworthiness(X, Y, k = min(10, k-1))
                
                # Combine scores (lower is better)
                scores[k,d] <- kl_div - trust
                
              }, error = function(e) {
                message("Error for k=", k, ", d=", d, ": ", e$message)
              })
            }
          }
          
          best <- which(scores == min(scores, na.rm = TRUE), arr.ind = TRUE)[1,]
          return(list(k = best[1], d = best[2]))
        }
        
        # Helper functions
        KL_divergence <- function(P, Q) {
          P <- as.matrix(dist(P))
          Q <- as.matrix(dist(Q))
          P <- P / sum(P)
          Q <- Q / sum(Q)
          return(sum(P * log((P + 1e-10) / (Q + 1e-10))))
        }
        
        trustworthiness <- function(X, Y, k) {
          n <- nrow(X)
          r_x <- apply(as.matrix(dist(X)), 2, rank)
          r_y <- apply(as.matrix(dist(Y)), 2, rank)
          sum_penalty <- sum(sapply(1:n, function(i) {
            sum(pmax(r_x[i, r_y[i,] <= k] - k, 0))
          }))
          return(1 - 2 / (n * k * (2*n - 3*k - 1)) * sum_penalty)
        }
        
        params <- optimize_params(X, max_k, max_dim, perplexity)
        params$d <- max(2,params$d)
        projection <- lle(X, params$k, params$d)
        
        return(list(projection = projection, k = params$k, d = params$d))
      }
      
      
      result.LLE <- lle_optimized(X_unique)
      Proj <- result.LLE$projection
    }, {
      # No projection, use the data as is
      Proj <- X_unique
    }
    )

    # Return the projection and unique data frame
    return(list(Projected = Proj, UniqueData = df_unique))
  }


############### Functions for clustering ###############

# Apply clustering algorithms
apply_clustering <- function(projection_results, clustering_methods, cluster_number_methods, method_for_hcpc, distance_metric = "euclidean", seed = 42) {

  cluster_list <- pbmcapply::pbmclapply(names(projection_results), function(projection) {
    projection_result <- projection_results[[projection]]
    tryCatch({
      if (is.null(projection_result$UniqueData$Target)) {
        stop("Required data is missing for clustering.")
      }
      projected_data <- as.data.frame(projection_result$Projected)
      clusters_per_projection <- lapply(clustering_methods, function(cluster_alg) {
        clusters_results <- lapply(cluster_number_methods, function(cluster_number_method) {
          set.seed(seed)
          clusters <- performClustering(X = projected_data, Target = projection_result$UniqueData$Target,
                                        clustering_method = cluster_alg, cluster_number_method = cluster_number_method,
                                        method_for_hcpc = method_for_hcpc,
                                        distance_metric = distance_metric)
          clusters <- renameClusters(trueCls = projection_result$UniqueData$Target, clusters)
          combined_result <- create_combined_result(projected_data, clusters, projection_result)
          return(combined_result)
        })
        names(clusters_results) <- cluster_number_methods
        return(clusters_results)
      })
      names(clusters_per_projection) <- clustering_methods
      return(clusters_per_projection)

    }, error = function(err) {
      message("Error in clustering with projection ", projection, ": ", err$message)
      return(NULL)
    })
  }, mc.cores = min(length(projection_methods), parallel::detectCores() - 1))

  names(cluster_list) <- names(projection_results)
  return(cluster_list)
}


# Function to determine the number of clusters
findOptimalClusters <- function(X, clustering_method, cluster_number_method, nProc = 1, Target = NULL,
                                method_for_hcpc = "ward", distance_metric = "euclidean") {

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

  # Function to calculate mode
  calculateMode <- function(v) {
    if (length(v) == 0) return(NA)
    freqTable <- table(v)
    modes <- as.numeric(names(freqTable)[freqTable == max(freqTable)])
    if (length(modes) > 1) modes else modes[1]
  }

  nClusters <- switch(cluster_number_method,
                      "NbClust" = {
    # Adjust clustering method for compatibility with NbClust
    clustering_method <- switch(clustering_method,
                                                    "kmedoids" = "kmeans",
                                                    "diana" = "ward.D2",
                                                    "hcpc" = "ward.D2",
                                                    clustering_method
                        )

    # Try running NbClust using parallel processing
    nbclustIndices <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew",
                                            "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
                                            "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
                                            "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")

    # Determine number of clusters by running NbClust method by method
    nClustersTest <- parallel::mclapply(nbclustIndices, function(i) {
      res <- try(suppressWarnings(NbClust::NbClust(data = X, diss = NULL, distance = distance_metric,
                                                                       min.nc = 2, max.nc = 9, method = clustering_method, index = i)), silent = TRUE)
      if (!inherits(res, "try-error")) {
        length(unique(res$Best.partition))
      } else {
        warning(paste("Error with index", i))
        NA
      }
    }, mc.cores = nProc)

    # Compute mode of results to determine optimal number of clusters
    nC <- na.omit(unlist(nClustersTest))
    result <- if (length(nC) > 0) calculateMode(nC) else NA
    result
  },
                      "hcpc" = {
    res.hcpc <- FactoMineR::HCPC(X, consol = TRUE, method = method_for_hcpc, metric = distance_metric, nb.clust = -1,
                                                     iter.max = 100, graph = FALSE, nstart = 100)
    length(unique(res.hcpc$data.clust$clust))
  },
  # Default case if Target is provided or none matches
                      if (!is.null(Target)) {
                        length(unique(Target))
                      } else {
    1
  }
  )

  if (is.na(nClusters)) {
    warning("Unable to determine the number of clusters. Setting to 1 or the original classes if provided.")
    if (!is.null(Target)) {
      nClusters <- length(unique(Target))
    } else {
      nClusters <- 1
    }
  }

  return(nClusters)
}


# Function to perform clustering
performClustering <- function(X, Target = NULL, clustering_method = "kmeans", cluster_number_method = "orig", method_for_hcpc = "ward",
                              distance_metric = "euclidean") {

  # Check if input X is a non-empty data frame
  if (!is.data.frame(X) || nrow(X) == 0) stop("Input data must be a non-empty data frame.")

  # Define valid clustering methods
  valid_clustering_methods <- c("none", "kmeans", "kmedoids",
                                "diana", "hcpc", "ward.D2", "single", "average", "median", "complete", "centroid")
  valid_hcpc_methods <- c("average", "single", "complete", "ward", "weighted", "flexible", "gaverage")
  valid_distances <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  valid_cluster_number_methods <- c("orig", "NbClust", "hcpc")

  # Validate the selected method
  if (!(clustering_method %in% valid_clustering_methods))
    stop(paste("Invalid clustering method. Choose from:", paste(valid_clustering_methods, collapse = ", ")))
  if (!(method_for_hcpc %in% valid_hcpc_methods))
    stop(paste("Invalid method for HCPC clustering. Choose from:", paste(valid_hcpc_methods, collapse = ", ")))
  if (!(distance_metric %in% valid_distances))
    stop(paste("Invalid distance metric. Choose from:", paste(valid_distances, collapse = ", ")))
  if (!(cluster_number_method %in% valid_cluster_number_methods))
    stop(paste("Invalid method for cluster number determination. Choose from:", paste(valid_cluster_number_methods, collapse = ", ")))

  # If Target is NULL, initialize it with default values
  if (is.null(Target)) Target <- rep(1, nrow(X))

  # Determine the number of clusters using the custom function

  if (clustering_method == "none") {
    nClusters <- length(unique(Target))
  } else {
    nClusters <- findOptimalClusters(X, clustering_method = clustering_method,
                                     cluster_number_method = cluster_number_method, nProc = nProc, Target = Target)
  }

  # Perform clustering based on the specified method
  cluster_result <- switch(clustering_method,
                           "none" = Target, # Return Target directly if method is "none"
                           "kmeans" = {
    # Perform K-means clustering
    kmeans(X, centers = nClusters, nstart = 100,)$cluster
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
    FactoMineR::HCPC(res.pca, nb.clust = nClusters, metric = distance_metric,
                                              graph = FALSE, nstart = 100, method = method_for_hcpc)$data.clust$clust
  }, {
    # Perform other hierarchical clustering methods
    clusterObject <- stats::hclust(dist(X, method = distance_metric), method = clustering_method)
    as.numeric(cutree(clusterObject, k = nClusters))
  })

  # Return the clustering result
  return(cluster_result)
}


# Function to align cluster names with likely class labels of the prior classification
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
  permutations <- combinat::permn(nLabels)

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


# Function to calculate Adjusted Rand Index (ARI) without using mclust
calculate_adjusted_rand_index <- function(Target, Clusters) {
  # Helper function to compute a contingency table
  contingency_table <- function(x, y) {
    table(factor(x, levels = unique(c(x, y))), factor(y, levels = unique(c(x, y))))
  }

  # Compute the contingency table
  cont_table <- contingency_table(Target, Clusters)

  # Sum over rows and columns
  sum_comb_n <- function(n) sum(choose(n, 2))

  # Total number of pairs
  N <- sum(cont_table)

  # Sum of combinations for each cell
  sum_comb_c_ij <- sum(mapply(choose, as.numeric(cont_table), 2))

  # Sum of combinations for each row
  sum_comb_a_i <- sum_comb_n(rowSums(cont_table))

  # Sum of combinations for each column
  sum_comb_b_j <- sum_comb_n(colSums(cont_table))

  # Calculate ARI
  expected_index <- sum_comb_a_i * sum_comb_b_j / choose(N, 2)
  max_index <- (sum_comb_a_i + sum_comb_b_j) / 2
  ari <- (sum_comb_c_ij - expected_index) / (max_index - expected_index)

  return(ari)
}


# Function to calculate Dunn's Index using a distance matrix
calculate_dunn_index <- function(distance_matrix, clusters) {
  # Convert the distance vector to a complete distance matrix
  distance_matrix <- as.matrix(distance_matrix)

  # Number of clusters
  unique_clusters <- unique(clusters)
  k <- length(unique_clusters)

  # Initialize variables for inter-cluster and intra-cluster distances
  min_inter_cluster_dist <- Inf
  max_intra_cluster_dist <- 0

  # Calculate inter-cluster distances
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      cluster_i_indices <- which(clusters == unique_clusters[i])
      cluster_j_indices <- which(clusters == unique_clusters[j])

      for (p in cluster_i_indices) {
        for (q in cluster_j_indices) {
          dist <- distance_matrix[p, q]
          if (!is.na(dist) && dist < min_inter_cluster_dist) {
            min_inter_cluster_dist <- dist
          }
        }
      }
    }
  }

  # Calculate intra-cluster distances
  for (i in 1:k) {
    cluster_indices <- which(clusters == unique_clusters[i])
    if (length(cluster_indices) > 1) {
      # Ensure there's more than one point in the cluster
      for (p in 1:(length(cluster_indices) - 1)) {
        for (q in (p + 1):length(cluster_indices)) {
          dist <- distance_matrix[cluster_indices[p], cluster_indices[q]]
          if (!is.na(dist) && dist > max_intra_cluster_dist) {
            max_intra_cluster_dist <- dist
          }
        }
      }
    }
  }

  # Ensure the max_intra_cluster_dist is valid for division
  if (max_intra_cluster_dist == 0) {
    stop("Intra-cluster distance is zero, which might indicate that each cluster has a single point.")
  }

  # Calculate Dunn's Index
  dunn_index <- min_inter_cluster_dist / max_intra_cluster_dist

  return(dunn_index)
}


# Calculate the Davies-Bold cluster index
calculate_davies_bouldin_index <- function(clusters, data) {

  # Number of clusters
  k <- length(unique(clusters))

  # Ensure data is numeric
  numeric_data <- data[, sapply(data, is.numeric)]

  # Calculate centroids for each cluster
  centroids <- aggregate(numeric_data, by = list(cluster = clusters), FUN = mean)
  rownames(centroids) <- centroids$cluster
  centroids <- centroids[, -1] # Remove the cluster column

  # Calculate average distance within clusters
  avg_within_dist <- tapply(1:nrow(data), clusters, function(idx) {
    cluster_data <- data[idx,]
    centroid <- centroids[as.character(clusters[idx[1]]),]
    mean(apply(cluster_data, 1, function(x) sqrt(sum((x - centroid) ^ 2))))
  })

  emicheps <- .Machine$double.eps ^ 0.5

  # Calculate Davies-Bouldin index
  db_index <- mean(sapply(1:k, function(i) {
    max_ratio <- 0
    for (j in 1:k) {
      if (i != j) {
        dist_between <- sqrt(sum((centroids[i,] - centroids[j,]) ^ 2))
        ratio <- (avg_within_dist[i] + avg_within_dist[j]) / (dist_between + emicheps)
        max_ratio <- max(max_ratio, ratio)
      }
    }
    max_ratio
  }))

  return(db_index)
}


# Calculate the Density-Based Clustering Validation (DBCV) index
calculate_dbcv <- function(data, clusters, dist_method = "euclidean") {

  # Validate input
  if (!dist_method %in% c("euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski")) {
    stop("Unsupported distance method.")
  }

  # Calculate the DBCV index (using a simplified heuristic approach)
  intra_density <- numeric(max(clusters))
  inter_density <- numeric(max(clusters))

  for (i in unique(clusters)) {
    cluster_points <- data[clusters == i,]
    if (nrow(cluster_points) > 1) {
      intra_cluster_dist_matrix <- dist(cluster_points, method = dist_method)
      intra_density[i] <- mean(intra_cluster_dist_matrix)
    }
  }

  for (i in unique(clusters)) {
    for (j in unique(clusters)) {
      if (i < j) {
        combined_points <- rbind(data[clusters == i,], data[clusters == j,])
        inter_cluster_dist_matrix <- dist(combined_points, method = dist_method)
        inter_density[i] <- mean(inter_cluster_dist_matrix)
      }
    }
  }

  dbcv_index <- sum(intra_density) / (sum(intra_density) + sum(inter_density))

  return(dbcv_index)
}


# Calculate various cluster indices and return as a data frame
calculate_cluster_scores_df <- function(clustering_results, projection_methods, clustering_methods, cluster_number_methods,
                                        distance_metric = "euclidean") {

  cluster_scores <- pbmcapply::pbmclapply(projection_methods, function(projection_method) {
    if (is.null(clustering_results[[projection_method]])) {
      warning("Projection for method ", projection_method, " returned NULL. Skipping.")
      return(NULL)
    }

    cluster_scores_per_alg <- lapply(clustering_methods, function(cluster_alg) {
      cluster_scores_per_number_method <- lapply(cluster_number_methods, function(cluster_number_method) {

        combined_result <- clustering_results[[projection_method]][[cluster_alg]][[cluster_number_method]]

        projected_data <- combined_result[, 1:2]
        Target <- combined_result$Target
        Clusters <- combined_result$Cluster

        cluster_accuracy <- sum(Clusters == Target) / length(Clusters)

        distance_matrix <- stats::dist(projected_data, method = distance_metric)
        Silhouette_index <- mean(cluster::silhouette(Clusters, distance_matrix)[, "sil_width"])
        Dunn_index <- calculate_dunn_index(distance_matrix = distance_matrix, clusters = Clusters)
        Rand_index <- calculate_adjusted_rand_index(Clusters, Target)
        DaviesBouldin_index <- calculate_davies_bouldin_index(clusters = Clusters, data = projected_data)
        dbcv_index <- calculate_dbcv(data = projected_data, clusters = Clusters, dist_method = distance_metric)

        return(data.frame(
          projection_method = projection_method,
          clustering_method = cluster_alg,
          cluster_number_method = cluster_number_method,
          cluster_accuracy = cluster_accuracy,
          Silhouette_index = Silhouette_index,
          Dunn_index = Dunn_index,
          Rand_index = Rand_index,
          DaviesBouldin_index = DaviesBouldin_index,
          dbcv_index = dbcv_index,
          stringsAsFactors = FALSE
        ))
      })
      do.call(rbind, cluster_scores_per_number_method)
    })
    do.call(rbind, cluster_scores_per_alg)
  }, mc.cores = min(length(projection_methods), parallel::detectCores() - 1))

  cluster_scores_df <- do.call(rbind, cluster_scores)
  return(cluster_scores_df)
}


# Rank the metrics and find the best combination of methods
process_cluster_scores <-
  function(cluster_scores_df, selected_cluster_metrics = c("cluster_accuracy", "Silhouette_index", "Dunn_index",
                                                   "Rand_index", "DaviesBouldin_index", "dbcv_index")) {
    # Validate selected metrics
    valid_cluster_metrics <- c("cluster_accuracy", "Silhouette_index", "Dunn_index", "Rand_index", "DaviesBouldin_index", "dbcv_index")
    selected_cluster_metrics <- intersect(selected_cluster_metrics, valid_cluster_metrics)

    if (length(selected_cluster_metrics) == 0) {
      message("No valid cluster metrics selected for ranking.\n Reverting to all impleneted metrics: ")
      message(valid_cluster_metrics)
      selected_cluster_metrics <- valid_cluster_metrics
    }

    # Ranking scores: higher is better for most, lower is better for Davies-Bouldin
    rank_columns <- sapply(selected_cluster_metrics, function(metric) {
      if (metric == "DaviesBouldin_index") {
        -cluster_scores_df[[metric]] # Lower is better, so negate values to rank
      } else {
        cluster_scores_df[[metric]] # Higher is better
      }
    })

    # Helper function for row medians
    rowMedians <- function(x, na.rm = TRUE) {
      if (!is.matrix(x)) {
        stop("Input must be a matrix")
      }

      apply(x, 1, median, na.rm = na.rm)
    }

    # Calculate ranks for each metric
    rank_df <- apply(rank_columns, 2, rank, ties.method = "min")

    # Create a combined rank by summing across selected metrics
    combined_rank <- rowMedians(rank_df)
    cluster_scores_df$combined_rank <- combined_rank

    # Find the row with the highest combined rank (best combination)
    best_index <- which(combined_rank == max(combined_rank))

    # Extract the best combination information
    best_combination <- as.data.frame(cluster_scores_df[best_index, c("projection_method", "clustering_method", "cluster_number_method")])

    # Add columns for individual ranks to the data frame
    rank_names <- paste0(selected_cluster_metrics, "_rank")
    cluster_scores_df[rank_names] <- rank_df

    # Return results as a list
    return(list(
      ranked_data = cluster_scores_df,
      best_combination = best_combination
    ))
  }

############### Functions for plotting ###############

# Function to plot the projected data and prior classes or clusters
create_dataset_plots <- function(clustering_results, projection_methods, clustering_methods, cluster_number_methods,
                                 dataset_name, label_points, highlight_misclassified, points_colored_for, cells_colored_for,
                                 palette_target = NULL, palette_cluster = NULL) {

  projection_plots <- pbmcapply::pbmclapply(projection_methods, function(projection_method) {
    if (is.null(clustering_results[[projection_method]])) {
      warning("Projection for method ", projection_method, " returned NULL. Skipping.")
      return(ggplot2::ggplot() + ggplot2::theme_void())
      #      return(NULL)
    }

    cluster_algs_plots <- lapply(clustering_methods, function(cluster_alg) {
      cluster_number_method_plots <- lapply(cluster_number_methods, function(cluster_number_method) {

        combined_result <- clustering_results[[projection_method]][[cluster_alg]][[cluster_number_method]]

        plot <- plotVoronoiTargetProjection(
          X = combined_result[, 1:2],
          targets = combined_result$Target,
          clusters = combined_result$Cluster,
          labels = combined_result$Label,
          misclassified = combined_result$Misclassified,
          LabelPoints = label_points,
          ColorMisClassified = highlight_misclassified,
          points_colored_for = points_colored_for,
          cells_colored_for = cells_colored_for,
          palette_target = palette_target,
          palette_cluster = palette_cluster
        )
        plot_title <- paste(dataset_name, ": ", projection_method, "- ", cluster_alg)
        if (cluster_alg != "none") {
          plot_title <- paste(plot_title, "- ", cluster_number_method)
        } else {
          plot <- plot + labs(fill = "Target")
        }
        plot <- plot + labs(title = plot_title)
        return(plot) # Ensuring plot is returned
      })
      names(cluster_number_method_plots) <- cluster_number_methods
      return(cluster_number_method_plots)
    })
    names(cluster_algs_plots) <- clustering_methods
    return(cluster_algs_plots)
  }, mc.cores = min(length(projection_methods), parallel::detectCores() - 1))
  names(projection_plots) <- projection_methods
  return(projection_plots)
}


# Function to plot Voronoi cells with projection and clustering results
plotVoronoiTargetProjection <- function(X, targets, clusters = NULL,
                                        labels = NULL, LabelPoints = FALSE,
                                        misclassified = NULL, ColorMisClassified = FALSE,
                                        points_colored_for = "Target", cells_colored_for = "Cluster",
                                        palette_target = NULL, palette_cluster = NULL) {
  # Extended colorblind palette
  cb_palette <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  shape_values <- c(17, 16, 15, 18, 2, 1, 0, 5, 7, 8, 9, 19, 11, 12)

  # Initial checks of arguments
  if (ncol(X) < 2) stop("The input data `X` must have at least two columns.")
  if (nrow(X) != length(targets)) stop("The length of `targets` must match the number of rows in `X`.")

  if (!is.null(clusters) && length(clusters) != nrow(X))
    stop("The length of `clusters` must match the number of rows in `X`.")
  if (!is.null(misclassified) && length(misclassified) != length(labels))
    stop("The length of `misclassified` must match the number of labels.")

  if (!is.null(palette_target) && length(palette_target) < length(unique(targets))) {
    stop("The `palette_target` must provide enough color values for the number of target groups.")
  }
  if (!is.null(palette_cluster) && length(palette_cluster) < length(unique(clusters))) {
    stop("The `palette_cluster` must provide enough color values for the number of clusters.")
  }

  # Define default palettes if none are provided
  if (is.null(palette_target)) palette_target <- rep(cb_palette, length.out = length(unique(targets)))
  if (is.null(palette_cluster)) palette_cluster <- rep(cb_palette, length.out = length(unique(clusters)))

  # Assign labels if not present
  if (is.null(labels)) {
    if (!is.null(rownames(X))) {
      labels <- rownames(X)
    } else {
      labels <- seq_len(nrow(X))
    }
  }

  # Default values for clusters and misclassified if NULL
  if (is.null(clusters)) clusters <- targets
  if (is.null(misclassified)) misclassified <- rep(0, nrow(X))

  # Data preparation
  plotData <- data.frame(Proj1 = X[, 1], Proj2 = X[, 2], Target = targets,
                         Clusters = clusters, Label = labels, Misclassified = misclassified)

  plotData$DotsInformation <- if (points_colored_for == "Target") plotData$Target else plotData$Clusters
  plotData$CellsInformation <- if (cells_colored_for == "Cluster") plotData$Clusters else plotData$Target

  # Voronoi diagram computation
  voronoiData <- deldir::deldir(plotData$Proj1, plotData$Proj2)

  # Convert Voronoi tessellation to a data frame for plotting
  vor_polys <- deldir::tile.list(voronoiData)
  voronoi_df <- do.call(rbind, lapply(1:length(vor_polys), function(i) {
    data.frame(x = vor_polys[[i]]$x, y = vor_polys[[i]]$y, id = i)
  }))

  # Create plot with ggplot2
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = voronoi_df, ggplot2::aes(x = x, y = y, group = id, fill = as.factor(plotData$CellsInformation[id])),
                          alpha = 0.3, color = NA) +
    ggplot2::geom_point(data = plotData, ggplot2::aes(x = Proj1, y = Proj2,
                                                      color = as.factor(DotsInformation), shape = as.factor(DotsInformation))) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.5, 0.08), legend.direction = "horizontal",
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
  if (sum(plotData$Misclassified) > 0 & ColorMisClassified) {
    suppressMessages({
      plot <- plot +
        geom_text(
          aes(x = Inf, y = Inf,
              label = paste0(round(100 * sum(plotData$Misclassified) / nrow(X), 1), "% misclassified")),
          hjust = 1,
          vjust = 1
        )
    })
    q <- ggplot_build(plot)
    q$data[[2]]$colour[plotData$Misclassified == 1] <- "red"
    gtable <- ggplot2::ggplot_gtable(q)
    plot <- ggplotify::as.ggplot(function() grid::grid.draw(gtable))
  }

  return(plot)
}


# Function to extract all plots from the generated lists and combine them into a figure
combine_all_plots <- function(datasets, projection_plots, projection_methods, clustering_methods, cluster_number_methods) {
  if (any(sapply(list(datasets, projection_plots, projection_methods, clustering_methods, cluster_number_methods), is.null))) {
    stop("All input lists must be non-null and contain elements.")
  }

  if (any(sapply(list(datasets, projection_plots, projection_methods, clustering_methods, cluster_number_methods), length) == 0)) {
    stop("All input lists must contain non-empty elements.")
  }

  calculate_plots_per_dataset <- function()
    length(projection_methods) * length(clustering_methods) * length(cluster_number_methods)

  plots_per_dataset <- calculate_plots_per_dataset()

  if (plots_per_dataset == 0) {
    stop("No valid combinations found in input lists.")
  }

  all_plots <- unlist(lapply(projection_plots, function(projection) {
    if (!is.null(projection)) {
      return(unlist(projection, recursive = FALSE, use.names = FALSE))
    }
  }), recursive = FALSE, use.names = FALSE)

  if (length(all_plots) == 0) {
    stop("No plots available to combine.")
  }

  figure_count <- length(all_plots) %/% plots_per_dataset

  if (figure_count == 0) {
    stop("Insufficient plots to create figures based on input methods.")
  }

  if (length(clustering_methods) == 1 & clustering_methods == "none") {
    columns <- calculate_compact_matrix_dims(n_items = length(projection_methods))$ncol
    rows <- calculate_compact_matrix_dims(n_items = length(projection_methods))$nrow
  } else {
    columns <- length(clustering_methods) * length(cluster_number_methods)
    rows <- length(projection_methods)
  }

  create_combined_plot <- function(idx) {
    start_idx <- 1 + (idx - 1) * plots_per_dataset
    end_idx <- plots_per_dataset + (idx - 1) * plots_per_dataset

    if (end_idx > length(all_plots)) {
      stop(sprintf("Index out of range when creating plot: %s", idx))
    }

    return(cowplot::plot_grid(plotlist = all_plots[start_idx:end_idx], ncol = columns, nrow = rows))
  }

  if (figure_count > 1) {
    figures <- pbmcapply::pbmclapply(
      seq_along(datasets),
      create_combined_plot,
      mc.cores = min(figure_count, parallel::detectCores())
    )
  } else {
    figures <- lapply(seq_along(datasets), create_combined_plot)
  }

  names(figures) <- datasets
  return(figures)
}


# Function to claulate a compact matrix for a combined results figure 
calculate_compact_matrix_dims <- function(n_items) {
  if (n_items < 1) {
    stop("Number of items must be at least 1.")
  }

  # Calculate the number of rows (ceiling of square root)
  nrow <- ceiling(sqrt(n_items))

  # Calculate the number of columns (ceiling of items divided by rows)
  ncol <- ceiling(n_items / nrow)

  return(list(nrow = nrow, ncol = ncol))
}




