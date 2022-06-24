#' Calculates Proactive Conservation Prioritization Index
#'
#' \code{pcpi} Calculates the Proactive Conservation Prioritization Index (ref), a new tool to prioritize species for conservation, which can incorporates information about future threats.
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
#' # Load data
#' data(reptile_df)
#'
#' # Calculate inverse range area
#' reptile_df$inv_range_area <- 1/(reptile_df$range_area)
#'
#' var_out <- reptile_df[,6:8]
#' var_in <- reptile_df[,3:4]
#' weight_out <- c(1, 1.5, 3)
#' weight_in <-
#'   matrix(1:6, ncol(var_out), ncol(var_in))
#'
#' # No interactions or weighting
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = var_out
#'   )
#'
#' # No interactions but weighting for variables
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = var_out,
#'     weight_out = weight_out
#'   )
#'
#' # With interactions
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = var_out,
#'     var_in = var_in
#'   )
#'
#' # With interactions and weighting for both variables and interactions
#'
#' reptile_pcpi <-
#'   pcpi(
#'     sp = reptile_df$sp,
#'     var_out = var_out,
#'     var_in = var_in,
#'     weight_out = weight_out,
#'     weight_in = weight_in
#'   )
#'
#' @importFrom stats predict
#'
#' @export

pcpi <-
  function(sp,
           var_out,
           var_in = NULL,
           weight_out = NULL,
           weight_in = NULL) {

    gm_mean <- function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

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

    sv <- caret::preProcess(log(var_out + 1), method = c("range"))
    scaled_var_out <- predict(sv, log(var_out + 1)) + f

    if (max(var_in) != min(var_in)) {
      siv <- caret::preProcess(log(var_in + 1), method = c("range"))
      scaled_var_in <- predict(siv, log(var_in + 1)) + f
    } else {
      scaled_var_in <- var_in
    }

    weight_in_list <-
      lapply(1:nrow(weight_in), function(i) {
        t(weight_in[i, ] * t(scaled_var_in))
      })

    weight_out_vec <- weight_out / sum(weight_out)

    wv_list <-
      lapply(1:ncol(scaled_var_out), function(i) {
        weight_mat <- weight_in_list[[i]]

        wm_list <-
          lapply(1:ncol(weight_mat), function(j) {
            scaled_var_out[, i] * weight_mat[, j]
          })

        wm_mat <-
          do.call(cbind, wm_list)

        rowSums(wm_mat)
      })

      weighted_vars <- do.call(cbind, wv_list)

    scaled_weighted_list <-
      sapply(1:ncol(weighted_vars), function(i) {
        z <- weighted_vars[,i]
        ws <- caret::preProcess(as.data.frame(z), method = c("range"))
        xx <- (predict(ws, as.data.frame(z)))^(1/weight_out_vec[i]) + f
      })

    scaled_weighted_vars <-
      do.call(cbind, scaled_weighted_list)



    pcpi <-
      apply(scaled_weighted_vars, 1, gm_mean, na.rm = TRUE)

    pcpi <- pcpi / max(pcpi)

    rank <- rep(NA, length(pcpi))
    rank[order(pcpi, decreasing = TRUE)] <- 1:length(pcpi)

    data.frame(sp, pcpi, rank)
  }
