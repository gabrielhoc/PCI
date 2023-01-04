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
#' @importFrom stats predict
#'
#' @export

pcpi <- function (sp, var_out, var_in = NULL, weight_out = NULL, weight_in = NULL) {

    # geometric mean function
    am_mean <- function(x, weights, na.rm = TRUE) {

      sum(x*weights)/sum(weights)

    }

    # create matrix with 1 if "var_in" is not supplied

    if (is.null(var_in)) {

      var_in <- matrix(rep(1, nrow(var_out)))

      colnames(var_in) <- c("var_in")

    }

    if (all(class(var_in) == "numeric")) {

      var_in <- matrix(var_in)

      colnames(var_in) <- c("var_in")

    }

    # if weights are not supplied, assign 1
    if (is.null(weight_out)){

      weight_out <- rep(1, ncol(var_out))
    }

    if (is.null(weight_in)){

      weight_in <- matrix(1, ncol(var_out), ncol(var_in))
    }

    if (all(class(weight_in) == "numeric")) {

      weight_in <- matrix(weight_in)
    }

    # small number to add to zeros
    f <- 0

    # sv_out<- caret::preProcess(weight_out, method = c("range"), rangeBounds = c(f, 1))
    # sv_in <- caret::preProcess(weight_in, method = c("range"), rangeBounds = c(f, 1))
    #
    # weight_out <- predict(sv_out, weight_out)
    # weight_in <- predict(sv_in, weight_in)

    # calculate in weights
    weight_in_list <-
      lapply(1:nrow(weight_in), function(i) {

        t(weight_in[i, ] * t(var_in))

      })

    # apply in weights
    wv_list <- lapply(1:ncol(var_out), function(i) {

      weight_mat <- weight_in_list[[i]]

      wm_list <- lapply(1:ncol(weight_mat), function(j) {

        var_out[, i] * weight_mat[, j]

      })

      wm_mat <- do.call(cbind, wm_list)

      lapply(1:nrow(wm_mat), function(k) {

        sum(wm_mat[k,])

      }) %>%
        do.call(rbind, .)

    })

    weighted_vars <- do.call(cbind, wv_list)

    # scale weighted variables
    # scaled_weighted_list <- sapply(1:ncol(weighted_vars), function(i) {
    #
    #   z <- weighted_vars[, i]
    #
    #   ws <- caret::preProcess(as.data.frame(z), method = c("range"), rangeBounds = c(f, 1))
    #
    #   xx <- predict(ws, as.data.frame(z))
    #
    # })
    #
    # scaled_weighted_vars <- do.call(cbind, scaled_weighted_list)

    # weight and average variables

    scaled_vars <-
      apply(weighted_vars, 2, function(x){

        ps <- caret::preProcess(data.frame(x), method = c("range"), rangeBounds = c(f, 1))

        predict(ps, data.frame(x))[,1]

    })

    pcpi <- apply(scaled_vars, 1, am_mean, weight_out, na.rm = TRUE)

    rank <- rank(-pcpi)

    out <- data.frame(sp, pcpi, rank)
  }
