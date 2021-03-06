#' Function to plot a PCA of the data
#'
#' Function was originally part of the OlinkAnalyze package, description below:
#' Generates a PCA projection of all samples from NPX data along two principal components (default PC2 vs. PC1) including the explained variance and dots colored by QC_Warning using stats::prcomp and ggplot2::ggplot. 
#' The values are by default scaled and centered in the PCA and proteins with missing NPX values are by default removed from the corresponding assay. 
#' Unique sample names are required. 
#' Imputation by the median is done for assays with missingness <10\% for multi-plate projects and <5\% for single plate projects.
#' 
#' @param df data frame in long format with Sample Id, NPX and column of choice for colors
#' @param color_g Character value indicating which column to use for colors (default QC_Warning)
#' @param x_val Integer indicating which principal component to plot along the x-axis (default 1)
#' @param y_val Integer indicating which principal component to plot along the y-axis (default 2)
#' @param label_samples Logical. If TRUE, points are replaced with SampleID (default FALSE) 
#' @param drop_assays Logical. All assays with any missing values will be dropped. Takes precedence over sample drop.
#' @param drop_samples Logical. All samples with any missing values will be dropped.
#' @param n_loadings Integer. Will plot the top n_loadings based on size.
#' @param loadings_list Character vector indicating for which GeneID's to plot as loadings. It is possible to use n_loadings and loadings_list simultaneously.
#' @param verbose Logical. Whether warnings about the number of samples and/or assays dropped or imputed should be printed to the console.
#' 
#' @return An object of class "ggplot"
#' @keywords NPX, PCA
#' 
#' @export
#' @import dplyr stringr tidyr ggfortify ggrepel gridExtra

olink_pca_plot <- function(df, 
                           color_g = "QC_Warning", 
                           x_val = 1, 
                           y_val = 2, 
                           label_samples = FALSE, 
                           drop_assays = FALSE,
                           drop_samples = FALSE, 
                           n_loadings = 0, 
                           loadings_list = NULL,
                           verbose = TRUE,
                           return_prcomp=FALSE,
                           return_df_wide_matrix=FALSE,
                           predicted_prcomp=NULL,
                           plot_scores = TRUE){ 
  
  # assigns colours to points based on grouping
  if (color_g == "QC_Warning"){
    
    df_temp <- df %>% 
      group_by(SampleID) %>% 
      mutate(QC_Warning = if_else(any(QC_Warning == "Warning"), "Warning", "Pass")) %>% 
      ungroup()
    
    colors_for_pca <- df_temp %>%
      group_by(SampleID) %>% 
      summarise(pca_colors = unique(!!rlang::ensym(color_g))) %>%
      ungroup()
    
    
  } else {
    
    number_of_sample_w_more_than_one_color <- df %>% 
      group_by(SampleID) %>% 
      summarise(n_colors = n_distinct(!!rlang::ensym(color_g), na.rm = T)) %>%
      ungroup() %>%
      filter(n_colors > 1) %>%
      nrow(.)
    
    if(number_of_sample_w_more_than_one_color > 0) {
      
      stop(paste0("There are ", number_of_sample_w_more_than_one_color, " samples that do not have a unique color. Only one color per sample is allowed."))
      
    }else{
      
      df_temp <- df
      
      colors_for_pca <- df_temp %>%
        group_by(SampleID) %>% 
        summarise(pca_colors = unique(!!rlang::ensym(color_g))) %>%
        ungroup()
      
    }
    
  }
  
  # Checking if there are any proteins with 0 variance, they are filtered out
  df_temp <- df_temp %>% 
    group_by(GeneID) %>%
    mutate(assay_var = var(NPX, na.rm = T)) %>%
    ungroup() %>%
    filter(!(assay_var == 0 | is.na(assay_var))) %>%
    dplyr::select(-assay_var)

  # wide format
  df_wide <- df_temp %>% 
    dplyr::select(SampleID, GeneID, NPX) %>% 
    filter(!is.na(NPX)) %>% 
    pivot_wider(names_from=GeneID, values_from=NPX)
  
  percent_missingness <- colSums(is.na(df_wide[, -c(1:2)]))/nrow(df_wide)
  
  # assays with missingness > 10% are dropped from the PCA
  PERCENT_CUTOFF <- 0.25
  
  if(any(percent_missingness > PERCENT_CUTOFF)){
    
    removed_assays_index <- which(percent_missingness > PERCENT_CUTOFF)
    percent_missingness <- percent_missingness[-removed_assays_index]
    
    removed_assays_index <- removed_assays_index + 2
    removed_assays <- colnames(df_wide)[removed_assays_index]
    
    df_wide <- df_wide[, -removed_assays_index]
    
    if(verbose){
      warning(paste0("There are ",
                     paste0(length(removed_assays)), 
                     " assay(s) dropped due to high missingness (>",
                     round(PERCENT_CUTOFF*100),
                     "%)."))
    }
    
    if(!is.null(loadings_list)){
      
      dropped_loadings <- intersect(removed_assays, 
                                    loadings_list)
      
      
      if(length(dropped_loadings) > 0){
        
        if(verbose){
          warning(paste0("The loading(s) ",
                         paste0(dropped_loadings, collapse=", "),
                         " from the loadings_list are dropped due to high missingness. "))
        }
        
        loadings_list <- setdiff(loadings_list, dropped_loadings)
        
        if(length(loadings_list) == 0){
          
          loadings_list <- NULL
          
        }
      }
      
    }
    
  }
  
  # convert long format to matrix
  df_wide <- df_wide %>% 
    left_join(colors_for_pca,
              by = c('SampleID')) %>%
    dplyr::select(SampleID, pca_colors, everything()) 
  
  df_wide_matrix <- df_wide %>% 
    dplyr::select(-pca_colors) %>%
    column_to_rownames('SampleID')
  
  # scales and centres data (imputes if there are any missing values)
  df_wide_matrix <- preProcess(df_wide_matrix, method = c("knnImpute", "scale", "center")) %>%
    predict(df_wide_matrix) %>%
    as.matrix()
  
  if (return_df_wide_matrix) {
    return(df_wide_matrix)
  }
  
  if(!all(colSums(is.na(df_wide_matrix[, -c(1:2)])) == 0)){
    stop('Missingness imputation failed.')
  }
  
  # calculate pca via svd
  if (is.null(predicted_prcomp)) {  # applies prcomp to calculate PCA
    pca_fit <- prcomp(df_wide_matrix, scale. = FALSE, center = FALSE)  # already scaled and centered data
    
  } else {  # if the prcomp is already calculated on another dataset, apply it to the new dataset
    
    pca_fit <- predicted_prcomp  # loadings etc. are going to stay the same so set pca_fit to be the previous value of prcomp
    pca_fit$x <- predict(predicted_prcomp, newdata=df_wide_matrix)  # but the PCs will change so set this for the new dataset via predict
  }
  
  if (return_prcomp) {
    return(pca_fit)
  }
  
  #Standardizing and selecting components
  
  scaling_factor_lambda <- pca_fit$sdev*sqrt(nrow(df_wide_matrix))
  
  PCX <- pca_fit$x[,x_val]/scaling_factor_lambda[x_val]
  PCY <- pca_fit$x[,y_val]/scaling_factor_lambda[y_val]
  PoV <- pca_fit$sdev^2/sum(pca_fit$sdev^2)
  LX <- pca_fit$rotation[, x_val]
  LY <- pca_fit$rotation[, y_val]
  
  observation_names <- df_wide$SampleID
  observation_colors <- df_wide$pca_colors
  
  scores <- cbind(PCX, PCY)
  loadings <- data.frame(variables = rownames(pca_fit$rotation), LX, LY)
  
  range_PX <- c(-abs(min(PCX, na.rm = TRUE)), abs(max(PCX, na.rm = TRUE)))
  range_PY <- c(-abs(min(PCY, na.rm = TRUE)), abs(max(PCY, na.rm = TRUE)))
  range_LX <- c(-abs(min(LX, na.rm = TRUE)), abs(max(LX, na.rm = TRUE)))
  range_LY <- c(-abs(min(LY, na.rm = TRUE)), abs(max(LY, na.rm = TRUE)))
  
  loadings_scaling_factor <- 0.8/max(range_LX/range_PX, range_LY/range_PY)
  
  #Plotting
  
  pca_plot <- ggplot(scores, aes(x = PCX, y = PCY)) +
    xlab(paste0("PC", x_val,  " (", round(PoV[x_val]*100, digits = 2), "%)")) +
    ylab(paste0("PC", y_val, " (", round(PoV[y_val]*100, digits = 2), "%)")) 
  
  #Drawing scores
  
  if (plot_scores) {
    plot_alpha <- 0.85
  } else {
    plot_alpha <- 0
  }
  
  if(label_samples){
    
    pca_plot <- pca_plot +
      geom_text(aes(label = observation_names, color = observation_colors), size = 3, alpha=plot_alpha) +
      labs(color = color_g) +
      guides(size = FALSE)
    
  } else {
    
    pca_plot <- pca_plot +
      geom_point(aes(color = observation_colors), size=1.7, alpha=plot_alpha) +
      guides(alpha=FALSE) +
      labs(color = color_g) +
      guides(size = FALSE)
    
  }
  
  #Drawing loadings
  
  if(n_loadings > 0 | !is.null(loadings_list)) {
    
    N_loadings <- data.frame(matrix(vector(), 0, ncol(loadings)),
                             stringsAsFactors=F)
    colnames(N_loadings) <- colnames(loadings)
    
    L_loadings <- N_loadings
    
    if (length(n_loadings) > 1) {
      
      # n_loadings is a vector of protein names
      N_loadings <- loadings[loadings$variables %in% n_loadings,]
      
    } else if (n_loadings > 0) {
      
      #Largest loadings based on Pythagoras
      
      N_loadings <- loadings %>%
        mutate(abs_loading = sqrt(LX^2 + LY^2)) %>%
        arrange(desc(abs_loading)) %>%
        head(n_loadings)
    }
    
    if(!is.null(loadings_list)){
      
      #Selected loadings
      
      L_loadings <- loadings %>%
        filter(variables %in% loadings_list)
    }
    
    loadings <- rbind(N_loadings, 
                      L_loadings) %>%
      distinct()
    
    pca_plot <- pca_plot +
      geom_segment(data = loadings,
                   aes(x = 0, 
                       y = 0,
                       xend = LX*loadings_scaling_factor, 
                       yend = LY*loadings_scaling_factor), 
                   arrow = arrow(length = unit(1/2, "picas")),
                   color = "black") +
      geom_label_repel(data = loadings, 
                       aes(x = LX*loadings_scaling_factor,
                           y = LY*loadings_scaling_factor,
                           label = variables), 
                       box.padding = 1, 
                       show.legend = F,
                       segment.colour = 'gray',
                       size=2.25)
    
    if (!plot_scores) {
      pca_plot <- pca_plot +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        geom_path(data=circleFun(max(loadings$abs_loading) * loadings_scaling_factor), aes(x, y))
    }
  }
  
  return(pca_plot)
}

circleFun <- function(radius, center = c(0,0), npoints=1000) {
  
  tt <- seq(0, 2*pi, length.out = npoints)
  xx <- center[1] + radius * cos(tt)
  yy <- center[2] + radius * sin(tt)
  
  return(data.frame(x = xx, y = yy))
}
