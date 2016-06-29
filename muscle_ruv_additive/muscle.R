library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    source("../code/datamaker_only_counts.R")
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

    num_sv <- sva::num.sv(t(Y), mod = X, method = "be")

    start.time <- proc.time()
    ash_args <- list()
    ruv_adinf_NTF <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   posthoc_inflate = TRUE,
                                   limmashrink = FALSE,
                                   additive_inflate = TRUE)
    ruv_adinf_NTT <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   posthoc_inflate = TRUE,
                                   limmashrink = TRUE,
                                   additive_inflate = TRUE)
    ruv_adinf_NFF <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   posthoc_inflate = FALSE,
                                   limmashrink = FALSE,
                                   additive_inflate = TRUE)
    ruv_adinf_NFT <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   posthoc_inflate = FALSE,
                                   limmashrink = TRUE,
                                   additive_inflate = TRUE)
    tot.time <- proc.time() - start.time

    postmean_df <- data.frame(ash_TF = ruv_adinf_NTF$PosteriorMean,
                              ash_TT = ruv_adinf_NTT$PosteriorMean,
                              ash_FF = ruv_adinf_NFF$PosteriorMean,
                              ash_FT = ruv_adinf_NFT$PosteriorMean)
    pi0vec      <- c(ruv_adinf_NTF$fitted.g$pi[1],
                     ruv_adinf_NTT$fitted.g$pi[1],
                     ruv_adinf_NFF$fitted.g$pi[1],
                     ruv_adinf_NFT$fitted.g$pi[1])
    lfdr_df     <- data.frame(ash_TF = ruv_adinf_NTF$lfdr,
                              ash_TT = ruv_adinf_NTT$lfdr,
                              ash_FF = ruv_adinf_NFF$lfdr,
                              ash_FT = ruv_adinf_NFT$lfdr)

    nmeth <- length(pi0vec)

    mse <- colMeans((postmean_df - beta_true) ^ 2)
    auc <- rep(NA, nmeth)
    if(args_val$nullpi < 1) {
        for (index in 1:nmeth) {
            auc[index] <- pROC::roc(predictor = lfdr_df[, index], response = which_null)$auc
        }
    }

    mult_scale_vec <- c(ruv_adinf_NTF$ruv$multiplier,
                        ruv_adinf_NTT$ruv$multiplier,
                        ruv_adinf_NFF$ruv$multiplier,
                        ruv_adinf_NFT$ruv$multiplier)
    add_scale_vec  <- c(ruv_adinf_NTF$ruv$additiver,
                        ruv_adinf_NTT$ruv$additiver,
                        ruv_adinf_NFF$ruv$additiver,
                        ruv_adinf_NFT$ruv$additiver)

    return_vec <- c(pi0vec, mse, auc, mult_scale_vec, add_scale_vec)

    return(return_vec)
}

itermax <- 200

## these change
Nsamp_seq  <- c(5, 10, 20)
nullpi_seq <- c(0.5, 0.9, 1)

par_vals <- expand.grid(list(1:itermax, Nsamp_seq, nullpi_seq))
colnames(par_vals) <- c("current_seed", "Nsamp", "nullpi")
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
    par_list[[list_index]] <- list()
    for (inner_list_index in 1:ncol(par_vals)) {
        par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
        names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
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
sout <- t(parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)

## ## on RCC, use this
## np <- mpi.universe.size() - 1
## cluster <- makeMPIcluster(np)
## sout <- t(snow::parSapply(cl = cluster, X = par_list, FUN = one_rep, current_params = args_val))
## stopCluster(cluster)
## mpi.exit()


save(sout, file = "sout_ash_adinf.Rd")

colnames(sout) <- rep(c("ruvash_adinf_TF", "ruvash_adinf_TT",
                        "ruvash_adinf_FF", "ruvash_adinf_FT"), times = 5)

pi0_mat  <- cbind(par_vals, sout[, 1:4])
mse_mat  <- cbind(par_vals, sout[, 5:8])
auc_mat  <- cbind(par_vals, sout[, 9:12])
mult_mat <- cbind(par_vals, sout[, 13:16])
add_mat  <- cbind(par_vals, sout[, 17:20])

write.csv(pi0_mat, file = "pi0_ruvash_adinf.csv", row.names = FALSE)
write.csv(mse_mat, file = "mse_ruvash_adinf.csv", row.names = FALSE)
write.csv(auc_mat, file = "auc_ruvash_adinf.csv", row.names = FALSE)
write.csv(mult_mat, file = "mult_ruvash_adinf.csv", row.names = FALSE)
write.csv(add_mat, file = "add_ruvash_adinf.csv", row.names = FALSE)
