library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
        start.time <- proc.time()
        num_methods <- 10 ## change this if number of methods in fit_method.R changes.

        trash2 <- R.utils::evalWithTimeout({
            trash <- tryCatch({
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

                num_sv <- sva::num.sv(t(Y), mod = X, method = "be")

                fit_out <- fit_all_competitors(Y = Y, X = X, num_sv = num_sv,
                                               control_genes = half_null)


                pi0hat <- fit_out$pi0hat

                tot.time <- proc.time() - start.time

                TRUE
            },
            error = function(e) { NULL })
            TRUE
        },
        timeout = 2700,
        onTimeout = "silent")

        if (is.null(trash)) {
            pi0hat <- rep(NA, num_methods)
        }

        return(pi0hat)
}


itermax <- 30

## these change
Nsamp_seq    <- c(2, 10, 50)
alt_type_seq <- c("spiky", "near_normal", "flattop", "skew", "big_normal", "bimodal")

par_vals <- expand.grid(list(1:itermax, Nsamp_seq, alt_type_seq))
par_vals$nullpi <- seq(0.01, 0.99, length = itermax)
colnames(par_vals) <- c("current_seed", "Nsamp", "alt_type", "nullpi")

par_vals$poisthin <- TRUE

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
                par_list[[list_index]]$tau_seq <- c(0.25, 0.5, 1, 2) /
                    sqrt(2 * par_vals$Nsamp[list_index] - 2)
                par_list[[list_index]]$mu_seq  <- c(0, 0, 0, 0)
            } else if (current_alt_type == "near_normal") {
                par_list[[list_index]]$pi_vals <- c(2/3, 1/3)
                par_list[[list_index]]$tau_seq <- c(1, 2) /
                    sqrt(2 * par_vals$Nsamp[list_index] - 2)
                par_list[[list_index]]$mu_seq  <- c(0, 0)
            } else if (current_alt_type == "flattop") {
                par_list[[list_index]]$pi_vals <- rep(1/7, 7)
                par_list[[list_index]]$tau_seq <- rep(0.5, 7) /
                    sqrt(2 * par_vals$Nsamp[list_index] - 2)
                par_list[[list_index]]$mu_seq  <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
            } else if (current_alt_type == "skew") {
                par_list[[list_index]]$pi_vals <- c(1/4, 1/4, 1/3, 1/6)
                par_list[[list_index]]$tau_seq <- c(2, 1.5, 1, 1) /
                    sqrt(2 * par_vals$Nsamp[list_index] - 2)
                par_list[[list_index]]$mu_seq  <- c(-2, -1, 0, 1)
            } else if (current_alt_type == "big_normal") {
                par_list[[list_index]]$pi_vals <- c(1)
                par_list[[list_index]]$tau_seq <- c(4) /
                    sqrt(2 * par_vals$Nsamp[list_index] - 2)
                par_list[[list_index]]$mu_seq  <- c(0)
            } else if (current_alt_type == "bimodal") {
                par_list[[list_index]]$pi_vals <- c(0.5, 0.5)
                par_list[[list_index]]$tau_seq <- c(1, 1) /
                    sqrt(2 * par_vals$Nsamp[list_index] - 2)
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
args_val$Ngene        <- 10000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

##one_out <- one_rep(new_params = par_list[[10]], current_params = args_val)

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores()-1)
sout <- parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val)
stopCluster(cl)

## ## on RCC, use this
## np <- mpi.universe.size() - 1
## cluster <- makeMPIcluster(np)
## sout <- snow::parSapply(cl = cluster, X = par_list, FUN = one_rep, current_params = args_val)
## stopCluster(cluster)
## mpi.exit()


save(sout, file = "piohat_out.Rd")

pi0hat_df <- as.data.frame(t(sout))

colnames(pi0hat_df) <- c("ols", "ash", "succotash", "cate_rr", "leapp_sparse", "leapp_ridge",
                         "sva", "cate_nc", "ruv2", "ruv4")
pi0hat_df <- cbind(pi0hat_df, par_vals)

write.csv(pi0hat_df, file = "pi0hat_df.csv", row.names = FALSE)

write.csv(par_vals, file = "sim_settings.csv", row.names = FALSE)
