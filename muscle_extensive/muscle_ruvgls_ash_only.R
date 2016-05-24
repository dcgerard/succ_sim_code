one_rep <- function(new_params, current_params) {
    start.time <- proc.time()
    source("../code/datamaker_only_counts.R")
    source("../code/fit_methods.R")
    args_val <- append(current_params, new_params)
    set.seed(new_params$current_seed)

    d_out <- datamaker_counts_only(args_val)

    which_null <- d_out$meta$null

    half_null <- which_null
    half_null[half_null == 1][sample(1:sum(which_null), size = sum(which_null) / 2)] <- 0

    beta_true <- rep(0, length = args_val$Ngene)
    beta_true[!which_null] <- d_out$meta$true_log2foldchange

    X <- as.matrix(model.matrix(~d_out$input$condition))
    colnames(X) <- c("Intercept", "Treatment")
    Y <- t(log2(as.matrix(d_out$input$counts + 1)))

    ## OLS + ASH ---------------------------------------------------------
    ash_ruv_out <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null))
    betahat_df <- data.frame(ash = ash_ruv_out$PosteriorMean)
    lfdr_df <- data.frame(ash = ash_ruv_out$lfdr)
    pi0hat_vec <- ash_ruv_out$fitted.g$pi[1]

    fit_out <- list()
    fit_out$betahat_df <- betahat_df
    fit_out$lfdr_df <- lfdr_df
    fit_out$pi0hat_vec <- pi0hat_vec

    fit_out$betahat_df$beta_true <- beta_true
    fit_out$lfdr_df$which_null   <- which_null

    tot.time <- proc.time() - start.time

    return(fit_out)
}

itermax <- 200

## these change
Nsamp_seq    <- c(5, 10, 20)
nullpi_seq   <- c(0.5, 0.9)
alt_type_seq <- c("spiky", "near_normal", "flattop", "skew", "big_normal", "bimodal")

par_vals <- expand.grid(list(1:itermax, Nsamp_seq, nullpi_seq, alt_type_seq))
colnames(par_vals) <- c("current_seed", "Nsamp", "nullpi", "alt_type")

null_df <- expand.grid(list(1:itermax, Nsamp_seq))
colnames(null_df) <- c("current_seed", "Nsamp")
null_df$nullpi <- 1
null_df$alt_type <- "all_null"

par_vals <- rbind(par_vals, null_df)

par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
    par_list[[list_index]] <- list()
    inner_adjust <- 0
    for (inner_list_index in 1:ncol(par_vals)) {
        if (colnames(par_vals)[inner_list_index] == "alt_type") {
            par_list[[list_index]]$alt_type = "mixnorm"
            current_alt_type <- par_vals[list_index, inner_list_index]
            if (current_alt_type == "spiky") {
                par_list[[list_index]]$pi_vals <- c(0.4, 0.2, 0.2, 0.2)
                par_list[[list_index]]$tau_seq <- c(0.25, 0.5, 1, 2)
                par_list[[list_index]]$mu_seq  <- c(0, 0, 0, 0)
            } else if (current_alt_type == "near_normal") {
                par_list[[list_index]]$pi_vals <- c(2/3, 1/3)
                par_list[[list_index]]$tau_seq <- c(1, 2)
                par_list[[list_index]]$mu_seq  <- c(0, 0)
            } else if (current_alt_type == "flattop") {
                par_list[[list_index]]$pi_vals <- rep(1/7, 7)
                par_list[[list_index]]$tau_seq <- rep(0.5, 7)
                par_list[[list_index]]$mu_seq  <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
            } else if (current_alt_type == "skew") {
                par_list[[list_index]]$pi_vals <- c(1/4, 1/4, 1/3, 1/6)
                par_list[[list_index]]$tau_seq <- c(2, 1.5, 1, 1)
                par_list[[list_index]]$mu_seq  <- c(-2, -1, 0, 1)
            } else if (current_alt_type == "big_normal") {
                par_list[[list_index]]$pi_vals <- c(1)
                par_list[[list_index]]$tau_seq <- c(4)
                par_list[[list_index]]$mu_seq  <- c(0)
            } else if (current_alt_type == "bimodal") {
                par_list[[list_index]]$pi_vals <- c(0.5, 0.5)
                par_list[[list_index]]$tau_seq <- c(1, 1)
                par_list[[list_index]]$mu_seq  <- c(-2, 2)
            }
            inner_adjust <- inner_adjust + 3
        } else {
            par_list[[list_index]][[inner_list_index + inner_adjust]] <-
                par_vals[list_index, inner_list_index]
            names(par_list[[list_index]])[inner_list_index + inner_adjust] <-
                colnames(par_vals)[inner_list_index]
        }
    }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 1
args_val$tissue       <- "muscle"
args_val$path         <- "../../../data/gtex_tissue_gene_reads/"
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores()-1)
sout <- parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val)
stopCluster(cl)


ruvash_betahat   <- sout[1, ]
ruvash_lfdr      <- sout[2, ]
ruvash_pi0hat    <- sout[3, ]

save(ruvash_betahat, file = "betahat_muscle_ruvash.Rd")
save(ruvash_lfdr, file = "lfdr_muscle_ruvash.Rd")
save(ruvash_pi0hat, file = "pi0hat_muscle_ruvash.Rd")
