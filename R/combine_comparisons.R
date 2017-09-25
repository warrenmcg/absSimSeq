combine_comparisons <- function(comparisons, which_comparison = "alr", 
                                group_name = "denom", group_labels = NULL,
                                split = NULL, file_prefix = NULL) {
  cbbPalette <- c("#000000", "#D55E00", "#56B4E9", "#E69F00", "#009E73", "#F0E442",
                  "#0072B2", "#CC79A7")
  cols2include <- c("target_id", "tpr", "fdr", "fpr")
  if (!is.null(split)) {
    pattern <- split[1]
    col <- as.integer(split[2])
    splits <- strsplit(group_labels, pattern, fixed = TRUE)
    s_vs_k <- as.factor(sapply(splits, function(x) x[col]))
    tool <- as.factor(sapply(splits, function(x) x[-col]))
    group_df <- data.frame(comp = which_comparison, group = group_labels,
                           linetype = as.integer(s_vs_k),
                           tool = as.integer(tool))
  } else {
    group_df <- data.frame(comp = which_comparison, group = group_labels)
  }
  
  result <- lapply(seq_along(comparisons), function(index){
    comparison <- comparisons[[index]]
    data_type <- as.character(group_df$comp[index])
    result <- comparison[[data_type]][, cols2include]
    result$Method <- group_df$group[index]
    result
  })
  
  result <- do.call(rbind, result)
  
  g <- ggplot2::ggplot(result[which(!is.na(result$tpr)), ],
                       ggplot2::aes(x = fpr, y = tpr, color = Method,
                                    group = Method, linetype = Method)) +
    ggplot2::geom_line() + ggplot2::xlim(c(0,0.25)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 3,
                         alpha = 0.3)
  if (split) {
    g <- g + ggplot2::scale_linetype_manual(values = group_df$linetype) +
      ggplot2::scale_color_manual(values = cbbPalette[group_df$tool])
  } else {
    g <- g + ggplot2::scale_color_manual(values = cbbPalette[as.integer(
      factor(group_df$group))])
  }
  g <- g +  ggplot2::theme_bw() +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.5,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 16,
                                                        face = "bold"),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 20),
                   axis.line.x = ggplot2::element_line(color="black",
                                                       size = 1),
                   axis.line.y = ggplot2::element_line(color="black",
                                                       size = 1)) +
    ggplot2::labs(group = eval(group_name))
  
  g2 <- ggplot2::ggplot(result[which(!is.na(result$tpr)), ],
                        ggplot2::aes(x = fdr, y = tpr, color = Method,
                                     linetype = Method)) +
    ggplot2::geom_path(size = 0.8, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = 3,
                        alpha = 0.3) + 
    ggplot2::geom_vline(xintercept = 0.1, linetype = 3,
                        alpha = 0.3)
  if (split) {
    g2 <- g2 + ggplot2::scale_linetype_manual(values = group_df$linetype) +
      ggplot2::scale_color_manual(values = cbbPalette[group_df$tool])
  }
  g2 <- g2 + ggplot2::theme_bw() +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.5,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 16,
                                                        face = "bold"),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 20),
                   axis.line.x = ggplot2::element_line(color="black",
                                                       size = 1),
                   axis.line.y = ggplot2::element_line(color="black",
                                                       size = 1)) +
    ggplot2::labs(group = eval(group_name))
  
  pdf(paste0(file_prefix, "_roc.pdf"))
  print(g)
  invisible(dev.off())
  
  pdf(paste0(file_prefix, "_sensVsFDR.pdf"))
  print(g2)
  invisible(dev.off())
  
  result
}
