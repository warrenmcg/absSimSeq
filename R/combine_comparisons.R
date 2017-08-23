combine_comparisons <- function(comparisons, which_comparison = "alr", 
                                group_name = "denom",
                                group_labels = NULL, file_prefix = NULL) {
  cols2include <- c("target_id", "tpr", "fdr", "fpr")
  group_df <- data.frame(comp = which_comparison, label = group_labels)
  
  result <- lapply(seq_along(comparisons), function(index){
    comparison <- comparisons[[index]]
    data_type <- as.character(group_df$comp[index])
    result <- comparison[[data_type]][, cols2include]
    result$group <- group_df$label[index]
    result
  })
  
  result <- do.call(rbind, result)
  
  g <- ggplot2::ggplot(result[which(!is.na(result$tpr)), ],
                       ggplot2::aes(x = fpr, y = tpr, color = group,
                                    group = group)) + 
    ggplot2::geom_line() + ggplot2::xlim(c(0,0.25)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha = 0.3) +
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
                   legend.title = ggplot2::element_text(size = 20)) +
    ggplot2::labs(color = eval(group_name))
  
  g2 <- ggplot2::ggplot(result[which(!is.na(result$tpr)), ],
                        ggplot2::aes(x = fdr, y = tpr, color = group)) + 
    ggplot2::geom_path(size = 0.8, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, linetype="dashed",
                        alpha = 0.3) + 
    ggplot2::geom_vline(xintercept = 0.1, linetype="dashed",
                        alpha = 0.3) +
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
                   legend.title = ggplot2::element_text(size = 20)) +
    ggplot2::labs(color = eval(group_name))
  
  pdf(paste0(file_prefix, "_roc.pdf"))
  print(g)
  invisible(dev.off())
  
  pdf(paste0(file_prefix, "_sensVsFDR.pdf"))
  print(g2)
  invisible(dev.off())
  
  result
}
