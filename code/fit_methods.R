
get_auc <- function(predictor, response) {
    pROC::roc(response = response, predictor = predictor)$auc
}

#' Fit all competitors and return betahat, lfdr, and pi0hat.
#'
#'
#' This function won't work yet if ofinterest is not equal to ncol(X).
#'
#' @param Y An n by p matrix of numerics. The normalized count data.
#' @param X An n by k matrix. The covariates.
#' @param num_sv A positive integer. The number of hidden confounders.
#' @param ofinterest A positive integer. The column of X whose
#'     betahats to return.
#'
#'
#'
fit_all_competitors <- function(Y, X, num_sv, ofinterest = ncol(X), control_genes = NULL) {

    pi0hat_vec <- c()


    trash <- tryCatch({
        ## OLS-----------------------------------------------------------------
        ols_fit <- get_ols(log_counts = t(Y), condition = X[, 2])
        ols_out <- fit_freq_methods(out_obj = ols_fit)

        betahat_df <- data.frame(ols = ols_fit$betahat)
        lfdr_df    <- data.frame(ols = ols_out$q_storey$lfdr)
        pi0hat_vec <- c(pi0hat_vec, ols_out$q_storey$pi0)
        TRUE
    },
    error = function(e){NULL})
    if (is.null(trash)) {
        betahat_df <- data.frame(ols = rep(NA, length = ncol(Y)))
        lfdr_df    <- data.frame(ols = rep(NA, length = ncol(Y)))
        pi0hat_vec <- c(pi0hat_vec, NA)
    }
    trash <- TRUE

    trash <- tryCatch({
        ## OLS + ASH ---------------------------------------------------------
        ash_ols <- ashr::ash(betahat = ols_fit$betahat, sebetahat = ols_fit$sebetahat)
        betahat_df$ols_ash <- ash_ols$PosteriorMean
        lfdr_df$ols_ash <- ash_ols$lfdr
        pi0hat_vec <- c(pi0hat_vec, ash_ols$fitted.g$pi[1])
        TRUE
    },
    error = function(e){NULL})
    if (is.null(trash)) {
        betahat_df$ols_ash <- rep(NA, length = ncol(Y))
        lfdr_df$ols_ash    <- rep(NA, length = ncol(Y))
        pi0hat_vec         <- c(pi0hat_vec, NA)
    }
    trash <- TRUE



    trash <- tryCatch({
        ## succotash---------------------------------------------------------
        succ_out <- succotashr::succotash(Y = Y, X = X, k = num_sv,
                                          fa_method = "pca", num_em_runs = 3,
                                          optmethod = "em")

        betahat_df$succotash <- succ_out$betahat
        lfdr_df$succotash    <- succ_out$lfdr
        pi0hat_vec           <- c(pi0hat_vec, succ_out$pi0)
        succ_lfsr            <- data.frame(normal = succ_out$lfsr)
        TRUE
    },
    error = function(e){NULL})
    if (is.null(trash)) {
        betahat_df$succotash <- rep(NA, length = ncol(Y))
        lfdr_df$succotash    <- rep(NA, length = ncol(Y))
        pi0hat_vec           <- c(pi0hat_vec, NA)
        succ_lfsr            <- data.frame(normal = rep(NA, length = ncol(Y)))
    }
    trash <- TRUE

    ## succotash-t---------------------------------------------------------
    ## succ_out_t <- succotashr::succotash(Y = Y, X = X, k = num_sv,
    ##                                   fa_method = "pca", num_em_runs = 3,
    ##                                   likelihood = "t", mix_type = "uniform")

    ## betahat_df$succotash_t <- succ_out_t$betahat
    ## lfdr_df$succotash_t    <- succ_out_t$lfdr
    ## pi0hat_vec             <- c(pi0hat_vec, succ_out_t$pi0)

    ## lfsr data frame for succotash -------------------------------------
    ## succ_lfsr <- data.frame(normal = succ_out$lfsr, t = succ_out_t$lfsr)


    trash <- tryCatch({
        ## RR CATE-----------------------------------------------------------
        ## error if only use one surrogate variable.
        num_sv_cate <- max(num_sv, 2)
        if (nrow(Y) - num_sv_cate - ncol(X) <= 0) {
            num_sv_cate <- 0
        }
        cate_rr <- cate::cate(~Treatment, Y = Y,
                              X.data = data.frame(Treatment = X[, 2]),
                              r = num_sv_cate, fa.method = "pc", adj.method = "rr",
                              calibrate = FALSE)
        cate_rr_out           <- list()
        cate_rr_out$betahat   <- cate_rr$beta
        cate_rr_out$sebetahat <- sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(nrow(X))
        cate_rr_out$pvalue    <- cate_rr$beta.p.value
        cate_qv               <- fit_freq_methods(out_obj = cate_rr_out)

        betahat_df$cate_rr <- cate_rr_out$betahat
        lfdr_df$cate_rr    <- cate_qv$q_storey$lfdr
        pi0hat_vec         <- c(pi0hat_vec, cate_qv$q_storey$pi0)
        TRUE
    },
    error = function(e){NULL})
    if (is.null(trash)) {
        betahat_df$cate_rr <- rep(NA, length = ncol(Y))
        lfdr_df$cate_rr    <- rep(NA, length = ncol(Y))
        pi0hat_vec         <- c(pi0hat_vec, NA)
    }
    trash <- TRUE
    ## RR CATE + ASH ----------------------------------------------------
    ## ash_rr_cate <- fit_ash(cate_rr_out)
    ## betahat_df$cate_rr_ash <- ash_rr_cate$PosteriorMean
    ## lfdr_df$cate_rr_ash <- ash_rr_cate$lfdr
    ## pi0hat_vec <- c(pi0hat_vec, ash_rr_cate$fitted.g$pi[1])

    ## LEAPP Sparse-------------------------------------------------------
    ## Sometimes stalls when sparse = TRUE
    ## trash2 <- tryCatch({
    ##     trash <- R.utils::evalWithTimeout({
    ##         leapp_sparse <- leapp::leapp(dat = t(Y), pred.prim = X[, ofinterest],
    ##                                      pred.covar = X[, -ofinterest],
    ##                                      num.fac = num_sv)
    ##         leapp_sparse_out         <- list()
    ##         leapp_sparse_out$betahat <- leapp_sparse$gamma
    ##         leapp_sparse_out$pvalue  <- leapp_sparse$p
    ##         leapp_sparse_qv          <- fit_freq_methods(out_obj = leapp_sparse_out)

    ##         betahat_df$leapp_sparse <- leapp_sparse_out$betahat
    ##         lfdr_df$leapp_sparse    <- leapp_sparse_qv$q_storey$lfdr
    ##         ##pi0hat_vec            <- c(pi0hat_vec, leapp_sparse_qv$q_storey$pi0)
    ##         pi0hat_vec              <- c(pi0hat_vec,
    ##                                      mean(abs(leapp_sparse_out$betahat) <= 10 ^ -10))
    ##         TRUE
    ##     },
    ##     timeout = 120,
    ##     onTimeout = "silent")
    ## },
    ## error = function(e){NULL})
    ## if(is.null(trash2)) {
    ##     betahat_df$leapp_sparse <- rep(NA, length = ncol(Y))
    ##     lfdr_df$leapp_sparse    <- rep(NA, length = ncol(Y))
    ##     pi0hat_vec              <- c(pi0hat_vec, NA)
    ## } else if (is.null(trash)) {
    ##     betahat_df$leapp_sparse <- rep(NA, length = ncol(Y))
    ##     lfdr_df$leapp_sparse    <- rep(NA, length = ncol(Y))
    ##     pi0hat_vec              <- c(pi0hat_vec, NA)
    ## }


    ## LEAPP Ridge--------------------------------------------------------
    ## Sometimes gives p-values that are all exactly zero and throws an error.
    trash <- tryCatch({
        leapp_ridge <- leapp::leapp(dat = t(Y), pred.prim = X[, ofinterest],
                                    pred.covar = X[, -ofinterest],
                                    num.fac = num_sv, method = "hard", sparse = FALSE)
        leapp_ridge_out         <- list()
        leapp_ridge_out$betahat <- leapp_ridge$gamma
        leapp_ridge_out$pvalue  <- leapp_ridge$p
        leapp_ridge_qv          <- fit_freq_methods(out_obj = leapp_ridge_out)

        betahat_df$leapp_ridge <- leapp_ridge_out$betahat
        lfdr_df$leapp_ridge    <- leapp_ridge_qv$q_storey$lfdr
        pi0hat_vec             <- c(pi0hat_vec, leapp_ridge_qv$q_storey$pi0)
        TRUE
    },
    error = function(e){NULL})
    if(is.null(trash)) {
        betahat_df$leapp_ridge <- rep(NA, length = ncol(Y))
        lfdr_df$leapp_ridge    <- rep(NA, length = ncol(Y))
        pi0hat_vec             <- c(pi0hat_vec, NA)
    }


    trash <- tryCatch({
        ## SVA----------------------------------------------------------------
        trash      <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
        X.sv       <- cbind(X, sva_out$sv)
        sva_ols    <- get_ols(log_counts = t(Y), condition = X.sv[, -1])
        sva_ols_qv <- fit_freq_methods(out_obj = sva_ols)

        betahat_df$sva_ols <- sva_ols$betahat
        lfdr_df$sva_ols    <- sva_ols_qv$q_storey$lfdr
        pi0hat_vec         <- c(pi0hat_vec, sva_ols_qv$q_storey$pi0)
        TRUE
    },
    error = function(e){NULL})
    if (is.null(trash)) {
        betahat_df$sva_ols <- rep(NA, length = ncol(Y))
        lfdr_df$sva_ols    <- rep(NA, length = ncol(Y))
        pi0hat_vec         <- c(pi0hat_vec, NA)
    }
    trash <- TRUE

    ## Negative control methods--------------------------------------------
    if (!is.null(control_genes)) {

        ruvash_out <- ashr::ash_ruv(Y = Y, X = X, k = num_sv, ctl = as.logical(control_genes))
        betahat_df$ruvash <- ruvash_out$PosteriorMean
        lfdr_df$ruvash <- ruvash_out$lfdr
        pi0hat_vec <- c(pi0hat_vec, ruvash_out$fitted.g$pi[1])


        trash <- tryCatch({
            ## Negative Control CATE --------------------------------------
            ## error if only use one surrogate variable.
            cate_nc <- cate::cate(~Treatment, Y = Y,
                                  X.data = data.frame(Treatment = X[, 2]),
                                  r = num_sv_cate, fa.method = "pc", adj.method = "nc",
                                  nc = as.logical(control_genes), calibrate = FALSE)
            cate_nc_out           <- list()
            cate_nc_out$betahat   <- cate_nc$beta
            cate_nc_out$sebetahat <- sqrt(cate_nc$beta.cov.row * cate_nc$beta.cov.col) /
                sqrt(nrow(X))
            cate_nc_out$pvalue    <- cate_nc$beta.p.value
            cate_nc_qv            <- fit_freq_methods(out_obj = cate_nc_out)

            betahat_df$cate_nc <- cate_nc_out$betahat
            lfdr_df$cate_nc    <- cate_nc_qv$q_storey$lfdr
            pi0hat_vec         <- c(pi0hat_vec, cate_nc_qv$q_storey$pi0)
            TRUE
        },
        error = function(e){NULL})
        if (is.null(trash)) {
            betahat_df$cate_nc <- rep(NA, length = ncol(Y))
            lfdr_df$cate_nc    <- rep(NA, length = ncol(Y))
            pi0hat_vec         <- c(pi0hat_vec, NA)
        }
        trash <- TRUE

        ## NC CATE + ASH
        ## ash_nc_cate <- fit_ash(cate_nc_out)
        ## betahat_df$cate_nc_ash <- ash_nc_cate$PosteriorMean
        ## lfdr_df$cate_nc_ash <- ash_nc_cate$lfdr
        ## pi0hat_vec <- c(pi0hat_vec, ash_nc_cate$fitted.g$pi[1])

        ## Exact same thing as RUV2
        ## ## Negative Control SVA ---------------------------------------------
        ## trash <- capture.output(sva_nc_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv,
        ##                                             controls = control_genes))
        ## X.sv_nc <- cbind(X, sva_nc_out$sv)
        ## sva_nc_ols <- get_ols(log_counts = t(Y), condition = X.sv_nc[, -1])
        ## sva_nc_ols_qv <- fit_freq_methods(out_obj = sva_nc_ols)

        ## betahat_df$sva_nc_ols <- sva_nc_ols$betahat
        ## lfdr_df$sva_nc_ols <- sva_nc_ols_qv$q_storey$lfdr
        ## pi0hat_vec <- c(pi0hat_vec, sva_nc_ols_qv$q_storey$pi0)


        trash <- tryCatch({
            ## RUV2-------------------------------------------------------
            ruv_ruv2 <- ruv::RUV2(Y = Y, X = as.matrix(X[, ofinterest]),
                                  ctl = as.logical(control_genes),
                                  k = num_sv, Z = as.matrix(X[, -ofinterest]))
            ruv_ruv2_list         <- list()
            ruv_ruv2_list$betahat <- c(ruv_ruv2$betahat)
            ruv_ruv2_list$pvalue  <- c(ruv_ruv2$p)
            ruv_ruv2_list$df      <- rep(ruv_ruv2$df, length = ncol(Y))
            ruv2_qv               <- fit_freq_methods(out_obj = ruv_ruv2_list)

            betahat_df$ruv2 <- ruv_ruv2_list$betahat
            lfdr_df$ruv2    <- ruv2_qv$q_storey$lfdr
            pi0hat_vec      <- c(pi0hat_vec, ruv2_qv$q_storey$pi0)
            TRUE
        },
        error = function(e){NULL})
        if (is.null(trash)) {
            betahat_df$ruv2 <- rep(NA, length = ncol(Y))
            lfdr_df$ruv2    <- rep(NA, length = ncol(Y))
            pi0hat_vec      <- c(pi0hat_vec, NA)
        }
        trash <- TRUE

        trash <- tryCatch({
            ## RUV4 -------------------------------------------------------
            ruv_ruv4 <- ruv::RUV4(Y = Y, X = as.matrix(X[, ofinterest]),
                                  ctl = as.logical(control_genes),
                                  k = num_sv, Z = as.matrix(X[, -ofinterest]))
            ruv_ruv4_list         <- list()
            ruv_ruv4_list$betahat <- c(ruv_ruv4$betahat)
            ruv_ruv4_list$pvalue  <- c(ruv_ruv4$p)
            ruv_ruv4_list$df      <- rep(ruv_ruv4$df, length = ncol(Y))
            ruv4_qv               <- fit_freq_methods(out_obj = ruv_ruv4_list)

            betahat_df$ruv4 <- ruv_ruv4_list$betahat
            lfdr_df$ruv4    <- ruv4_qv$q_storey$lfdr
            pi0hat_vec      <- c(pi0hat_vec, ruv4_qv$q_storey$pi0)
            TRUE
        },
        error = function(e){NULL})
        if (is.null(trash)) {
            betahat_df$ruv4 <- rep(NA, length = ncol(Y))
            lfdr_df$ruv4    <- rep(NA, length = ncol(Y))
            pi0hat_vec      <- c(pi0hat_vec, NA)
        }
        trash <- TRUE
    }

    return(list(betahat_df = betahat_df, lfdr_df = lfdr_df, pi0hat_vec = pi0hat_vec, succ_lfsr = succ_lfsr))
}
