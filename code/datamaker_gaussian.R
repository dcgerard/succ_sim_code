## generate gaussian data
datamaker_gaussian <- function(args) {

    assertthat::assert_that(is.list(args))

    dfargs <- default_gaussian_args(args)

    n <- dfargs$Nsamp
    p <- dfargs$Ngene
    k <- dfargs$Nconfounders

    X <- cbind(rep(1, length = n), c(rep(1, floor(n / 2)), rep(0, n - floor(n / 2))))

    beta <- matrix(NA, nrow = ncol(X), ncol = p)

    beta[1, ] <- rnorm(p)

    if (is.null(dfargs$nullpi)) {
        beta[2, ] <- rmixnorm(pi_vals = dfargs$pi_vals,
                              mu_seq = dfargs$mu_seq,
                              tau_seq = dfargs$tau_seq,
                              p = p)
    } else {
        which_null <- sort(sample(1:p, round(p * dfargs$nullpi)))
        isnullgene <- rep(FALSE, length = p)
        isnullgene[which_null] <- TRUE
        beta[2, isnullgene] <- 0
        beta[2, !isnullgene] <- rmixnorm(pi_vals = dfargs$pi_vals,
                                         mu_seq = dfargs$mu_seq,
                                         tau_seq = dfargs$tau_seq,
                                         p = sum(!isnullgene))
    }

    alpha <- matrix(rnorm(k * p, sd = sqrt(dfargs$sig2a)), nrow = k)
    Z <- matrix(rnorm(n * k, sd = sqrt(dfargs$sig2z)), nrow = n)
    E <- matrix(rnorm(n * p, sd = sqrt(dfargs$sig2e)), nrow = n)
    Y <- X %*% beta + Z %*% alpha + E

    out <- list()
    out$dfargs <- dfargs
    out$input <- list()
    out$input$Y <- Y
    out$input$X <- X
    out$input$beta <- beta
    out$input$Z <- Z
    out$input$alpha <- alpha
    out$input$E <- E
    out$input$which_null <- which_null
    out$input$isnullgene <- isnullgene

    return(out)
}


## give default arguments to datamaker_gaussian
default_gaussian_args <- function(args) {

    if (is.null(args$pi_vals)) {
        pi_vals <- 1
    }

    if (is.null(args$tau_seq)) {
        tau_seq <- 1
    }

    if (is.null(args$sig2e)) { ## error variance
        args$sig2e <- 1
    }

    if (is.null(args$sig2z)) {
        args$sig2z <- 1
    }

    if (is.null(args$sig2a)) {
        args$sig2a <- 1
    }

    return(args)
}
