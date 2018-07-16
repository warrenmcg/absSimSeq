#' @importFrom dplyr select mutate
#' @export
extract_oracles <- function(final_results, min_tpm = 1) {
  reg_results <- lapply(final_results$results, function(run) {
    run <- as.data.frame(run[c(1:5)])
    res <- dplyr::select(run, .data$abs_fold_changes,
                         ctr_tpm = .data$copy_numbers_per_cell.ctr_copy_numbers)
    res <- dplyr::mutate(res, target_id = rownames(res),
                         log_fc = log(.data$abs_fold_changes),
                         is_de = (.data$abs_fold_changes != 1))
    res <- dplyr::select(res, target_id, ctr_tpm, log_fc, is_de)
  })
  alr_results <- lapply(final_results$alr_data, function(run) {
    run <- as.data.frame(run)
    res <- dplyr::select(run, .data$alr_diff,
                         ctr_ratio =  .data$alr_means.ctr_copy_numbers)
    res <- dplyr::mutate(res, target_id = rownames(res),
                         log_fc = ifelse(abs(.data$alr_diff) < 1e-5,
                                         0, .data$alr_diff),
                         is_de = abs(.data$alr_diff) >= 1e-5)
    res <- dplyr::select(res, .data$target_id, .data$ctr_ratio, .data$log_fc,
                         .data$is_de)
  })
  list(main_oracle = reg_results,
       alr_oracle = alr_results)
}

#' @importFrom methods is
#' @importFrom dplyr select
#' @export
rename_fc_list <- function(de_list) {
  method_names <- names(de_list)
  new_list <- lapply(method_names, function(method) {
    method_list <- de_list[[method]]
    if(is(method_list, "list")) {
      lapply(method_list, function(x) {
        if (grepl("wt$", method)) {
          dplyr::select(x, .data$target_id, .data$pval, .data$qval,
                        log_fc = .data$b)
        } else if (grepl("lrt$", method)) {
          dplyr::select(x, .data$target_id, .data$pval, .data$qval)
        } else if (grepl("ALDEx2", method)) {
          dplyr::select(x, .data$target_id, .data$pval, .data$qval,
                        log_fc = .data$effect)
        } else {
          dplyr::select(x, .data$target_id, .data$pval, .data$qval,
                        log_fc = .data$beta)
        }
      })
    } else {
      x <- method_list
      if (grepl("wt$", method)) {
        dplyr::select(x, .data$target_id, .data$pval, .data$qval,
                      log_fc = .data$b)
      } else if (grepl("lrt$", method)) {
        dplyr::select(x, .data$target_id, .data$pval, .data$qval)
      } else if (grepl("ALDEx2", method)) {
        dplyr::select(x, .data$target_id, .data$pval, .data$qval,
                      log_fc = .data$effect)
      } else {
        dplyr::select(x, .data$target_id, .data$pval, .data$qval,
                      log_fc = .data$beta)
      }
    }
  })
  names(new_list) <- method_names
  new_list
}
