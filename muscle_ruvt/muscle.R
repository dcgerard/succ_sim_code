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
    ruv_out <- ruv::RUV4(Y = Y, X = X[, "Treatment", drop = FALSE],
                         ctl = as.logical(half_null), k = num_sv)
    betahat   <- c(ruv_out$betahat)
    sebetahat <- sqrt(ruv_out$sigma2 * ruv_out$multiplier)
    df <- nrow(X) - ncol(X) - num_sv
    ash_out <- ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, df = df)
    tot.time <- proc.time() - start.time



    mse <- mean((ash_out$PosteriorMean - beta_true) ^ 2)
    pi0 <- ash_out$fitted.g$pi[1]
    auc <- NA
    if(args_val$nullpi < 1) {
        auc <- pROC::roc(predictor = ash_out$lfdr, response = which_null)$auc
    }

    return_vec <- c(pi0, mse, auc)

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


write.csv(sout, file = "ruv_asht_mat.csv", row.names = FALSE)


colnames(sout) <- c("pi0hat", "mse", "auc")
ruv_asht_df <- as.data.frame(cbind(sout, par_vals))
save(ruv_asht_df, file = "ruv_asht_df.Rd")
