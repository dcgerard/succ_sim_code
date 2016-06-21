library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    start.time <- proc.time()
    source("../code/datamaker_only_counts.R")
    source("../code/datamaker_gaussian.R")
    source("../code/fit_methods.R")
    args_val <- append(current_params, new_params)
    set.seed(new_params$current_seed)

    d_out <- datamaker_gaussian(args_val)

    which_null <- d_out$input$isnullgene

    half_null <- which_null
    half_null[half_null][sample(1:sum(which_null), size = sum(which_null) / 2)] <- FALSE

    beta_true <- d_out$input$beta[2, ]

    X <- d_out$input$X
    colnames(X) <- c("Intercept", "Treatment")
    Y <- d_out$input$Y

    num_sv <- args_val$Nconfounders

    fit_out <- fit_all_competitors(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = half_null)


    fit_out$betahat_df$beta_true <- beta_true
    fit_out$lfdr_df$which_null   <- which_null
    fit_out$succ_lfsr$which_null <- which_null

    tot.time <- proc.time() - start.time

    return(fit_out)
}


itermax <- 200

## these change
Nsamp_seq    <- c(10, 20, 40)
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
args_val$Nconfounders <- 5


## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores() - 2)
sout <- parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val)
stopCluster(cl)

## ## on RCC, use this
## np <- mpi.universe.size() - 1
## cluster <- makeMPIcluster(np)
## sout <- snow::parSapply(cl = cluster, X = par_list, FUN = one_rep, current_params = args_val)
## stopCluster(cluster)
## mpi.exit()


betahat   <- sout[1, ]
lfdr      <- sout[2, ]
pi0hat    <- sout[3, ]
succ_lfsr <- sout[4, ]

save(betahat, file = "betahat_muscle.Rd")
save(lfdr, file = "lfdr_muscle.Rd")
save(pi0hat, file = "pi0hat_muscle.Rd")
save(succ_lfsr, file = "succ_lfsr_muscle.Rd")

write.csv(par_vals, file = "sim_settings.csv", row.names = FALSE)


method_vec <- colnames(betahat[[1]])
method_vec <- method_vec[-length(method_vec)]
save(method_vec, file = "method_vec.Rd")
