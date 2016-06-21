succ_then_ash <- function (Y, X, ctl = NULL, k = NULL,
                           cov_of_interest = ncol(X),
                           likelihood = c("normal", "t"),
                           ash_args = list(), limmashrink = FALSE,
                           include_intercept = TRUE, gls = TRUE,
                           fa_func = ashr:::pca_naive, fa_args = list()) {
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <=
        ncol(X))
    assertthat::assert_that(is.logical(gls))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.list(ash_args))
    assertthat::assert_that(is.null(ash_args$betahat))
    assertthat::assert_that(is.null(ash_args$sebetahat))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y))
    assertthat::assert_that(is.null(fa_args$r))
    assertthat::assert_that(is.function(fa_func))
    likelihood <- match.arg(likelihood)
    if (include_intercept) {
        X_scaled <- apply(X, 2, function(x) {
            x/sqrt(sum(x^2))
        })
        int_term <- rep(1, length = nrow(X))/sqrt(nrow(X))
        any_int <- any(colSums((int_term - X_scaled)^2) < 10^-14)
        if (!any_int) {
            X <- cbind(X, rep(1, length = nrow(X)))
        }
    }
    if (is.null(k)) {
        if (requireNamespace("sva", quietly = TRUE)) {
            message("Number of confounders not provided so being estimated with package sva.")
            k <- sva::num.sv(dat = t(Y), mod = X)
        } else {
            stop("If sva is not installed, then k needs to be provided. To install sva, run in R\n   source(\"https://bioconductor.org/biocLite.R\")\n   biocLite(\"sva\")")
        }
    }
    assertthat::assert_that(k + ncol(X) < nrow(X))
    if (k >= sum(ctl)) {
        stop("k is larger than the number of control genes so model not identified.\nReduce k or increase the number of control genes.\nYou can also try out succotashr. To install succotashr, run in R:\n    install.packages(\"devtools\")\n    devtools::install_github(\"dcgerard/succotashr\")")
    }
    X <- X[, c((1:ncol(X))[-cov_of_interest], cov_of_interest),
        drop = FALSE]
    cov_of_interest <- ncol(X)
    qr_x <- qr(X)
    Q <- qr.Q(qr_x, complete = TRUE) * sign(qr.R(qr_x)[cov_of_interest,
        cov_of_interest])
    Y_tilde <- crossprod(Q, Y)[cov_of_interest:nrow(Y), , drop = FALSE]
    fa_args$Y <- Y_tilde[2:nrow(Y_tilde), , drop = FALSE]
    fa_args$r <- k
    fa_out <- do.call(what = fa_func, args = fa_args)
    alpha <- fa_out$alpha
    sig_diag <- fa_out$sig_diag
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::are_equal(length(sig_diag), ncol(Y))
    assertthat::assert_that(all(sig_diag > 0))
    if (k != 0) {
        assertthat::assert_that(is.matrix(alpha))
        assertthat::are_equal(nrow(alpha), k)
        assertthat::are_equal(ncol(alpha), ncol(Y))
    } else {
        assertthat::assert_that(is.null(alpha))
    }
    if (requireNamespace("limma", quietly = TRUE) & limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag, df = nrow(X) -
            ncol(X) - k)
        sig_diag <- limma_out$var.post
    } else if (!requireNamespace("limma", quietly = TRUE) & limmashrink) {
        stop("limmashrink = TRUE but limma not installed. To install limma, run in R:\n    source(\"https://bioconductor.org/biocLite.R\")\n    biocLite(\"limma\")")
    }
    fnorm_x <- abs(qr.R(qr_x)[cov_of_interest, cov_of_interest])
    betahat_ols <- t(Y_tilde[1, , drop = FALSE])/fnorm_x
    alpha_scaled <- alpha/fnorm_x
    sig_diag_scaled <- sig_diag/(fnorm_x^2)
    Yc <- betahat_ols[ctl, , drop = FALSE]
    alphac <- alpha_scaled[ctl, , drop = FALSE]
    sig_diagc <- sig_diag_scaled[ctl]

    succout <- succotashr::succotash_given_alpha(Y = Yc, alpha = alphac, sig_diag = sig_diagc)
    Z1 <- succout$Z
    errordist <- list()
    for (index in 1:length(betahat_ols)) {
        errordist[[index]] <- ashr::normalmix(pi = succout$pi_vals,
                                              mean = rep(0, length(succout$pi_vals)),
                                         sd = sqrt(succout$tau_seq ^ 2 + succout$scale_val *
                                                                     sig_diag_scaled[index]))
    }

    betahat <- betahat_ols - alpha_scaled %*% Z1

    ruvashZ <- crossprod(solve(crossprod(alphac, 1 / sig_diagc * alphac)),
                         crossprod(alphac, 1 / sig_diagc * Yc))

    ash_args$betahat <- c(betahat)
    ash_args$errordist <- errordist
    ash_out <- do.call(what = stramash::stramash.workhorse, args = ash_args)
    ash_out$ruv <- list()
    ash_out$ruv$betahat_ols <- betahat_ols
    ash_out$ruv$sebetahat_ols <- sqrt(sig_diag_scaled)
    ash_out$ruv$betahat <- betahat
    ash_out$ruv$alphahat <- alpha
    ash_out$ruv$input <- ash_args
    ash_out$ruv$sigma2 <- sig_diag
    ash_out$ruv$fnorm_x <- fnorm_x
    ash_out$ruv$Z1 <- Z1
    ash_out$ruv$ruvashZ <- ruvashZ
    ash_out$succ <- succout
    return(ash_out)
}
