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

    ruvash_out <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null), k = num_sv)

    tstat_ruvash <- ruvash_out$ruv$betahat / ruvash_out$ruv$sebetahat

    sdmult <- sqrt(nrow(X) / (nrow(X) - ncol(X) - num_sv))
    tstat_ruvash_premult <- tstat_ruvash * sdmult

    ruv4out <- ruv::RUV4(Y = Y, X = X[, 2, drop = FALSE], ctl = as.logical(half_null),
                         k = num_sv)
    ruv4out_sebetahat <- sqrt(ruv4out$sigma2 * ruv4out$multiplier)

    ruv_se <- cbind(ruv4out_sebetahat, ruvash_out$ruv$sebetahat, ruvash_out$ruv$sebetahat_ols)


    return_vec <- c(mean(tstat_ruvash), sd(tstat_ruvash),
                    mean(tstat_ruvash_premult), sd(tstat_ruvash_premult),
                    mean(ruv4out$t), sd(ruv4out$t),
                    ruvash_out$ruv$multiplier, sdmult)

    return(return_vec)
}

itermax <- 200

## these change
Nsamp_seq  <- c(5, 10, 20)

par_vals <- expand.grid(list(1:itermax, Nsamp_seq))
colnames(par_vals) <- c("current_seed", "Nsamp")

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
args_val$nullpi       <- 1
args_val$poisthin     <- FALSE

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


write.csv(sout, file = "null_mat.csv", row.names = FALSE)
write.csv(sout, "/home/david/Dropbox/null_mat.csv", row.names = FALSE)
