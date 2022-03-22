#' Calculates Proactive Conservation Prioritization Index
#'
#' \code{pcpi} Calculates the Proactive Conservation Prioritization Index (ref), a new tool to prioritize species for conservation, which can incorporates information about future threats.
#'
#' @importFrom magrittr "%>%"
#'
#' @param sp character. Names of the taxa being evaluated.
#' @param var_out numeric. Threat variables. higher values must indicate increased threat.
#' @param var_in numeric. Interacting variables. Will modulate the effect of threat variables.
#' @param weight_out numeric. Weights for threat variables
#' @param weight_in numeric. Matrix of weights for the combination of interacting variables and threat variables.
#'
#' @return Data frame with PCPI and rank.
#'
#' @examples
#'
#' data(reptile_df)
#'
#' # No interactions or weighting
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = reptile_df[5:7]
#'   )
#'
#' # No interactions but weighting for variables
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = reptile_df[5:7],
#'     weight_out = weight_out
#'   )
#'
#' # With interactions
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = reptile_df[5:7],
#'     var_in = reptile_df[3:4]
#'   )
#'
#' # With interactions and weighting for both variables and interactions
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = reptile_df[5:7],
#'     var_in = reptile_df[3:4],
#'     weight_out = weight_out,
#'     weight_in = weight_in
#'   )
#'
#' @export

pcpi <-
  function(sp,
           var_out,
           var_in = NULL,
           weight_out = NULL,
           weight_in = NULL) {

    if (is.null(var_in)) {
      var_in <- matrix(rep(1, nrow(var_out)))
      colnames(var_in) <- c("var_in")
    }

    if (all(class(var_in) == "numeric")) {
      var_in <- matrix(var_in)
      colnames(var_in) <- c("var_in")
    }

    if (is.null(weight_out)) weight_out <- rep(1, ncol(var_out))
    if (is.null(weight_in)) weight_in <- matrix(1, ncol(var_out), ncol(var_in))
    if (all(class(weight_in) == "numeric")) weight_in <- matrix(weight_in)
    f <- 1 / 10^100

    sv <- caret::preProcess(var_out, method = c("range"))
    scaled_var_out <- predict(sv, var_out) + f

    if (max(var_in) != min(var_in)) {
      siv <- caret::preProcess(var_in, method = c("range"))
      scaled_var_in <- predict(siv, var_in) + f
    } else {
      scaled_var_in <- var_in
    }

    weight_in_list <-
      lapply(1:nrow(weight_in), function(i) {
        t(weight_in[i, ] * t(scaled_var_in))
      })

    weight_out_vec <- weight_out / sum(weight_out)

    weighted_vars <-
      lapply(1:ncol(scaled_var_out), function(i) {
        weight_mat <- weight_in_list[[i]]

        lapply(1:ncol(weight_mat), function(j) {
          scaled_var_out[, i] * weight_mat[, j]
        }) %>%
          do.call(cbind, .) %>%
          rowSums() * weight_out_vec[i]
      }) %>%
      do.call(cbind, .)

    scaled_weighted_vars <-
      apply(weighted_vars, 2, function(z) {
        ws <- caret::preProcess(as.data.frame(z), method = c("range"))
        predict(ws, as.data.frame(z)) + f
      }) %>%
      do.call(cbind, .)

    pcpi <-
      apply(scaled_weighted_vars, 1, prod, na.rm = TRUE)

    pcpi <- pcpi / max(pcpi)

    rank <- rep(NA, length(pcpi))
    rank[order(pcpi, decreasing = TRUE)] <- 1:length(pcpi)

    data.frame(sp, pcpi, rank)
  }
