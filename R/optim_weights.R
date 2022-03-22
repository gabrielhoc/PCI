#' Optimizes weighting for the calculation of Proactive Conservation Prioritization Index
#'
#' \code{optim_weights} Optimizes weights for calculating Proactive Conservation Prioritization Index (ref),
#'
#' @param sp character. Names of the taxa being evaluated.
#' @param var_out numeric. Threat variables. higher values must indicate increased threat.
#' @param var_in numeric. Interacting variables. Will modulate the effect of threat variables.
#' @param weight_out numeric. Weights for threat variables
#' @param weight_in numeric. Matrix of weights for the combination of interacting variables and threat variables.
#' @param compare numeric. Threat rank towards which weights will be optimized.
#' @param type character. Optimize weights for threat variables ("out"), for interacting variables ("in") or for both ("both").
#'
#' @details The Pearson correlation between the calculated rank and 'compare' is displayed as the weights are optimized.
#'
#' @return Vector ("out"), matrix ("in") or list ("both") with optimal weights.
#'
#' @examples
#'
#' # Load data
#' data(reptile_df)
#'
#' # Remove "DD" and "NE" species
#' reptile_df_iucn <-
#'   reptile_df[!reptile_df$iucn_cat %in% c("DD", "NE")]
#'
#' # Convert IUCN categories to numeric
#' reptile_df_iucn$iucn_cat <- as.numeric(reptile_df_iucn$iucn_cat)
#'
#' # Optimize weights
#' optim_weights <-
#'   optim_weights(sp = reptile_df_iucn$sp,
#'                 var_out = reptile_df_iucn[5:7],
#'                 var_in = reptile_df_iucn[3:4],
#'                 compare = reptile_df_iucn$iucn_cat)
#'
#' @export

optim_weights <-
  function(sp,
           var_out,
           var_in  = NULL,
           weight_out = NULL,
           weight_in = NULL,
           compare,
           type = "both",
           ...) {
    optim_args <- list(...)

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

    switch(type,
           `out` = {
             xx <-
               do.call(optim,
                       args = append(
                         optim_args,
                         list(
                           par = weight_out,
                           fn = function(x) {
                             ranks <-
                               pcpi(
                                 sp = sp,
                                 var_out = var_out,
                                 var_in = var_in,
                                 weight_out = x,
                                 weight_in = weight_in
                               )

                             out <-
                               1 - cor(
                                 compare,
                                 ranks$rank
                               )

                             cat(paste("correlation = ", 1 - out))
                             flush.console()

                             out
                           }
                         )
                       )
               )

             xx$par
           },
           `in` = {
             xx <-
               do.call(optim,
                       args = append(
                         optim_args,
                         list(
                           par = c(weight_in),
                           fn = function(x) {
                             weight_in_mat <-
                               matrix(x, ncol(var_out), ncol(var_in))

                             ranks <-
                               pcpi(
                                 sp = sp,
                                 var_out = var_out,
                                 var_in = var_in,
                                 weight_out = weight_out,
                                 weight_in = weight_in_mat
                               )

                             out <-
                               1 - cor(
                                 compare,
                                 ranks$rank
                               )

                             cat(paste("correlation = ", 1 - out))
                             flush.console()

                             out
                           }
                         )
                       )
               )

             weight_in_mat_xx <-
               matrix(xx$par, nrow(weight_in), ncol(weight_in))
           },
           `both` = {
             xx <-
               do.call(
                 what = optim,
                 args = append(
                   optim_args,
                   list(
                     par = c(weight_out, weight_in),
                     fn = function(x) {
                       part_out <- x[1:length(weight_out)]
                       part_in <- x[(1:length(weight_in)) + length(weight_out)]

                       weight_in_mat <-
                         matrix(part_in, nrow(weight_in), ncol(weight_in))

                       ranks <-
                         pcpi(
                           sp = sp,
                           var_out = var_out,
                           var_in = var_in,
                           weight_out = part_out,
                           weight_in = weight_in_mat
                         )

                       out <-
                         1 - cor(
                           compare,
                           ranks$rank
                         )

                       message(paste("correlation = ", 1 - out), "\r", appendLF=FALSE)
                       flush.console()

                       out
                     }
                   )
                 )
               )

             part_out_xx <- xx$par[1:length(weight_out)]
             part_in_xx <- xx$par[(1:length(weight_in)) + length(weight_out)]

             weight_in_mat_xx <-
               matrix(part_in_xx, nrow(weight_in), ncol(weight_in))

             list(
               weight_out = part_out_xx,
               weight_in = weight_in_mat_xx
             )
           }
    )
  }
