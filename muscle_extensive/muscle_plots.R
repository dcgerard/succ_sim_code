library(ggplot2)
library(reshape2)
library(dplyr)
library(pROC)



## pi0hat plots ------------------------------------------------------------------------
load("pi0hat_muscle.Rd")
load("method_vec.Rd")
par_vals <- read.csv("sim_settings.csv")

nullpi_seq   <- unique(par_vals$nullpi)
nsamp_seq    <- unique(par_vals$Nsamp)
alt_type_seq <- unique(par_vals$alt_type)


pdf(file = "./fig/pi0hat_different_alts.pdf")
for (alt_type_index in 1:(length(alt_type_seq) - 1)) {
    which_alt <- par_vals$alt_type == alt_type_seq[alt_type_index]

    pi0mat <- as.data.frame(t(sapply(pi0hat[which_alt], FUN = c)))
    names(pi0mat) <- method_vec
    pi0mat$nullpi <- par_vals$nullpi[which_alt]
    pi0mat$Nsamp  <- par_vals$Nsamp[which_alt]

    ## pi0_plot
    dummy_dat <- expand.grid(nullpi_seq[-length(nullpi_seq)], nsamp_seq)
    colnames(dummy_dat) <- c("nullpi", "Nsamp")

    name_vec <- colnames(pi0mat)
    colnames(pi0mat) <- gsub("pi0hat_", "", x = name_vec)
    long_dat <- melt(pi0mat, id.vars = c("nullpi", "Nsamp"),
                     measure.vars = colnames(pi0mat)[1:(ncol(pi0mat) - 2)])
    p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
    p <- p + facet_grid(nullpi~Nsamp) + ylab(expression(hat(pi)[0]))
    p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
    p <- p + geom_hline(data = dummy_dat, aes(yintercept = nullpi), lty = 2, color = "red", lwd = 1)
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
    p <- p + ggtitle(paste("Estimates of pi0 When Using Muscle Tissue, Alternative =",
                           alt_type_seq[alt_type_index]))
    print(p)
}

which_alt <- par_vals$alt_type == "all_null"
pi0mat <- as.data.frame(t(sapply(pi0hat[which_alt], FUN = c)))
names(pi0mat) <- method_vec
pi0mat$nullpi <- par_vals$nullpi[which_alt]
pi0mat$Nsamp  <- par_vals$Nsamp[which_alt]
dummy_dat <- data.frame(Nsamp = nsamp_seq, nullpi = rep(1, length = length(nsamp_seq)))

name_vec <- colnames(pi0mat)
colnames(pi0mat) <- gsub("pi0hat_", "", x = name_vec)
long_dat <- melt(pi0mat, id.vars = c("Nsamp"),
                 measure.vars = colnames(pi0mat)[1:(ncol(pi0mat) - 2)])
p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
p <- p + facet_grid(.~Nsamp) + ylab(expression(hat(pi)[0]))
p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
p <- p + geom_hline(data = dummy_dat, aes(yintercept = nullpi), lty = 2, color = "red", lwd = 1)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
p <- p + ggtitle("Estimates of pi0 When Using Muscle Tissue and All Null")
print(p)
dev.off()



## lfsr plots ------------------------------------------------------------------------
rm(list = ls())
load("betahat_muscle.Rd")
get_beta_sign <- function(mat) {
    return(sign(mat$beta_true))
}
sign_beta <- sapply(betahat, FUN = get_beta_sign)

get_betahat_sign <- function(mat) {
    return(sign(mat$succotash))
}
sign_betahat <- sapply(betahat, FUN = get_betahat_sign)
rm("betahat")

load("succ_lfsr_muscle.Rd")
get_lfsr <- function(mat) {
    return(mat$normal)
}
lfsr <- sapply(succ_lfsr, FUN = get_lfsr)

par_vals <- read.csv("sim_settings.csv")


lfsr_curve <- function(lfsr_vec, betahat_vec, betatrue_vec, props = seq(0.01, 1, by  = 0.01)) {

    ## ordering_lfsr <- order(lfsr_vec)
    ## lfsr_vec     <- lfsr_vec[ordering_lfsr]
    ## betahat_vec  <- betahat_vec[ordering_lfsr]
    ## betatrue_vec <- betatrue_vec[ordering_lfsr]
    fsp <- rep(NA, length = length(props))
    for (pindex in 1:length(props)) {
        current_p <- props[pindex]
        which_small <- lfsr_vec <= current_p
        if (sum(which_small) == 0) {
            fsp[pindex] <- 0
        } else {
            betahat_small <- betahat_vec[which_small]
            betatrue_small <- betatrue_vec[which_small]
            fsp[pindex] <- 1 - mean(betahat_small == betatrue_small)
        }
    }
    return(fsp)
}

lfsr_curve_concat <- function(concat_vec, props = seq(0.01, 1, by  = 0.01)) {
    num_entries <- length(concat_vec) / 3
    lfsr_vec     <- concat_vec[1:num_entries]
    betahat_vec  <- concat_vec[(num_entries + 1):(2 * num_entries)]
    betatrue_vec <- concat_vec[(2 * num_entries + 1):(3 * num_entries)]
    fsp <- lfsr_curve(lfsr_vec = lfsr_vec, betahat_vec = betahat_vec,
                      betatrue_vec = betatrue_vec, props = props)
    return(fsp)
}


overall_mat <- rbind(lfsr, sign_betahat, sign_beta)

fsp <- as.data.frame(t(apply(overall_mat, 2, FUN = lfsr_curve_concat)))
props = seq(0.01, 1, by  = 0.01)
names(fsp)   <- props
fsp$nullpi   <- par_vals$nullpi
fsp$Nsamp    <- par_vals$Nsamp
fsp$alt_type <- par_vals$alt_type
alt_type_seq <- unique(par_vals$alt_type)

fsp_dat <- melt(fsp, id.vars = c("nullpi", "Nsamp", "alt_type"))
fsp_dat$variable <- as.numeric(levels(fsp_dat$variable))[fsp_dat$variable]

names(fsp_dat) <- c("nullpi", "Nsamp", "alt_type", "lfsr_max", "fsp")

summary_lfsr <- aggregate(fsp ~ nullpi + Nsamp + alt_type + lfsr_max, data = fsp_dat,
                          FUN = quantile, probs = c(0.025, 0.5, 0.975))

summ <- cbind(summary_lfsr[, 1:4], summary_lfsr[,5])
colnames(summ) <- c("nullpi", "Nsamp", "alt_type", "lfsr_max", "Lower", "Median", "Upper")

pdf(file = "./fig/lfsr_succ.pdf")
for (alt_type_index in 1:length(alt_type_seq)) {
    current_alt_type <- alt_type_seq[alt_type_index]
    summ_sub <- summ[summ$alt_type == current_alt_type, ]
    p <- ggplot(data = summ_sub, aes(x = lfsr_max, y = Median), color = "blue") +
        facet_grid(nullpi~Nsamp) +
        geom_path() + geom_abline(slope = 1, intercept = 0, color = I("black"), lty = 2) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1, fill = I("green")) +
        ylim(c(0, 1)) + ylab("False Sign Proportion") + xlab("Max lfsr") +
        ggtitle(paste("FSP Below Max LFSR when Alt Type =", current_alt_type))
    print(p)
}
dev.off()


## Same thing but with s-values -----------------------------------------------
get_svalue <- function(lfsr_vec) {
    order_lfsr <- order(lfsr_vec)
    lfsr_vec <- lfsr_vec[order_lfsr]
    svalue <- cumsum(lfsr_vec) / 1:length(lfsr_vec)
    return(svalue[order(order_lfsr)])
}

svalue <- apply(lfsr, 2, get_svalue)

overall_mat <- rbind(svalue, sign_betahat, sign_beta)

fsp <- as.data.frame(t(apply(overall_mat, 2, FUN = lfsr_curve_concat)))
props = seq(0.01, 1, by  = 0.01)
names(fsp)   <- props
fsp$nullpi   <- par_vals$nullpi
fsp$Nsamp    <- par_vals$Nsamp
fsp$alt_type <- par_vals$alt_type
alt_type_seq <- unique(par_vals$alt_type)

fsp_dat <- melt(fsp, id.vars = c("nullpi", "Nsamp", "alt_type"))
fsp_dat$variable <- as.numeric(levels(fsp_dat$variable))[fsp_dat$variable]

names(fsp_dat) <- c("nullpi", "Nsamp", "alt_type", "svalue_max", "fsp")

summary_svalue <- aggregate(fsp ~ nullpi + Nsamp + alt_type + svalue_max, data = fsp_dat,
                          FUN = quantile, probs = c(0.025, 0.5, 0.975))

summ <- cbind(summary_svalue[, 1:4], summary_svalue[,5])
colnames(summ) <- c("nullpi", "Nsamp", "alt_type", "svalue_max", "Lower", "Median", "Upper")

pdf(file = "./fig/svalue_succ.pdf")
for (alt_type_index in 1:length(alt_type_seq)) {
    current_alt_type <- alt_type_seq[alt_type_index]
    summ_sub <- summ[summ$alt_type == current_alt_type, ]
    p <- ggplot(data = summ_sub, aes(x = svalue_max, y = Median), color = "blue") +
        facet_grid(nullpi~Nsamp) +
        geom_path() + geom_abline(slope = 1, intercept = 0, color = I("black"), lty = 2) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1, fill = I("green")) +
        ylim(c(0, 1)) + ylab("False Sign Proportion") + xlab("Max svalue") +
        ggtitle(paste("FSP Below Max svalue when Alt Type =", current_alt_type))
    print(p)
}
dev.off()


## AUC Boxplots ---------------------------------------------------------------------
rm(list = ls())
load("lfdr_muscle.Rd")
par_vals <- read.csv("sim_settings.csv")

get_auc <- function(lfdr_mat) {
    null_col <- which(colnames(lfdr_mat) == "which_null")
    auc_vec <- rep(NA, length = ncol(lfdr_mat) - 1)
    for(index in (1:length(auc_vec))[-null_col]) {
        if (all(is.na(lfdr_mat[, index]))) {
            auc_vec[index] <- NA
        } else {
            auc_vec[index] <- roc(response = lfdr_mat$which_null,
                                  predictor = c(lfdr_mat[, index]))$auc
            }
    }
    return(auc_vec)
}

## which_allnull <- par_vals$alt_type == "all_null"
## auc_mat <- as.data.frame(t(sapply(lfdr[!which_allnull], FUN = get_auc)))
## names(auc_mat) <- colnames(lfdr[[1]])[1:(ncol(lfdr[[1]]) - 1)]
## auc_mat$nullpi   <- par_vals$nullpi[!which_allnull]
## auc_mat$Nsamp    <- par_vals$Nsamp[!which_allnull]
## auc_mat$alt_type <- par_vals$alt_type[!which_allnull]
## write.csv(auc_mat, file = "auc_mat.csv", row.names = FALSE)

auc_mat <- read.csv("auc_mat.csv")
alt_type_seq <- unique(auc_mat$alt_type)


pdf(file = "./fig/auc_different_alts.pdf")
for (alt_type_index in 1:length(alt_type_seq)) {
    current_alt_type <- alt_type_seq[alt_type_index]

    auc_mat_sub <- auc_mat[auc_mat$alt_type == current_alt_type,]
    auc_mat_sub <- dplyr::select(auc_mat_sub, -alt_type)


    ## create dummy dat to plot max median auc
    nullpi_seq <- unique(auc_mat_sub$nullpi)
    nsamp_seq <- unique(auc_mat_sub$Nsamp)
    dummy_dat <- expand.grid(nullpi_seq, nsamp_seq)
    colnames(dummy_dat) <- c("nullpi", "Nsamp")
    med_mat <- matrix(NA, nrow = length(nullpi_seq) * length(nsamp_seq), ncol = ncol(auc_mat_sub) - 3)
    for (index in 1:(ncol(auc_mat_sub) - 3)) {
        form1 <- as.formula(paste(colnames(auc_mat_sub)[index], "~ nullpi + Nsamp"))
        out1 <- aggregate(form1, FUN = median, na.rm = TRUE,
                          data = auc_mat_sub)
        med_mat[, index] <- out1[, 3]
    }
    dummy_dat2 <- cbind(expand.grid(nullpi_seq[nullpi_seq != 1], nsamp_seq), apply(med_mat, 1, max))
    colnames(dummy_dat2) <- c("nullpi", "Nsamp", "max_med")
    ## end dummy dat


    long_dat <- melt(auc_mat_sub, id.vars = c("nullpi", "Nsamp"))
    p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
    p <- p + facet_grid(nullpi~Nsamp) + ylab(expression(hat(pi)[0]))
    p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
    p <- p + geom_hline(data = dummy_dat2, aes(yintercept = max_med), lty = 2,
                        color = "red", lwd = 0.6)
    p <- p + ggtitle(paste("AUC When Using Muscle Tissue, Alternative =",
                           alt_type_seq[alt_type_index]))
    print(p)
}
dev.off()


## Mean ROC Curves ------------------------------------------------------------------------
rm(list = ls())
load("lfdr_muscle.Rd")
par_vals <- read.csv("sim_settings.csv")
nullpi_seq <- unique(par_vals$nullpi)
Nsamp_seq <- unique(par_vals$Nsamp)
alt_type_seq <- unique(par_vals$alt_type)

get_fd <- function(predictor_vec, response_vec) {
    pord <- order(predictor_vec)
    oresp <- response_vec[pord]
    fd <- cumsum(oresp)
    return(fd)
}

get_roc <- function(lfdr_mat) {
    roc_out <- apply(lfdr_mat[-10], 2, get_fd, response_vec = lfdr_mat$which_null)

    fpr <- roc_out / sum(lfdr_mat$which_null)
    tpr <- (1:nrow(roc_out) - roc_out) / sum(1 - lfdr_mat$which_null)

    ##matplot(fpr, tpr, type = "l", col = 1:9, lty = 1)
    ##legend("bottomright", colnames(lfdr_mat), col = 1:9, lty = 1)
    return(list(tpr = tpr, fpr = fpr))
}

roc_out <- lapply(lfdr, get_roc)


ft_prsum <- function(ftpr_list) {
    tpr_mean <- matrix(0, nrow = nrow(ftpr_list[[1]][[1]]), ncol = ncol(ftpr_list[[1]][[1]]))
    fpr_mean <- matrix(0, nrow = nrow(ftpr_list[[1]][[2]]), ncol = ncol(ftpr_list[[1]][[2]]))

    for (index in 1:length(ftpr_list)) {
        tpr_mean <- tpr_mean + ftpr_list[[index]][[1]]
        fpr_mean <- fpr_mean + ftpr_list[[index]][[2]]
    }
    tpr_mean <- tpr_mean / length(ftpr_list)
    fpr_mean <- fpr_mean / length(ftpr_list)
    return(list(tpr = tpr_mean, fpr = fpr_mean))
}


ave_ft_list <- list()
list_index <- 1
for (nullpi_index in 1:(length(nullpi_seq) - 1)) {
    nullpi_current <- nullpi_seq[nullpi_index]
    for(Nsamp_index in 1:length(Nsamp_seq)) {
        Nsamp_current <- Nsamp_seq[Nsamp_index]
        for(alt_type_index in 1:(length(alt_type_seq) - 1)) {
            alt_type_current <- alt_type_seq[alt_type_index]
            current_list <- par_vals$nullpi == nullpi_current & par_vals$Nsamp == Nsamp_current &
                par_vals$alt_type == alt_type_current
            ave_ft_list[[list_index]] <- ft_prsum(roc_out[current_list])
            list_index <- list_index + 1
        }
    }
}



summary_par_vals <- expand.grid(alt_type_seq[1:(length(alt_type_seq) - 1)], Nsamp_seq, nullpi_seq[1:(length(nullpi_seq) - 1)])
names(summary_par_vals) <- c("alt_type", "Nsamp", "nullpi")


list_index <- 1
fpr <- ave_ft_list[[list_index]][[2]]
tpr <- ave_ft_list[[list_index]][[1]]

temp_mat          <- cbind(melt(fpr)[, 2:3], melt(tpr)$value)
temp_mat$nullpi   <- summary_par_vals$nullpi[list_index]
temp_mat$Nsamp    <- summary_par_vals$Nsamp[list_index]
temp_mat$alt_type <- summary_par_vals$alt_type[list_index]
colnames(temp_mat) <- c("Method", "fpr", "tpr", "nullpi", "Nsamp", "alt_type")

roc_mat <- temp_mat
for (list_index in 2:nrow(summary_par_vals)) {
    fpr <- ave_ft_list[[list_index]][[2]]
    tpr <- ave_ft_list[[list_index]][[1]]

    temp_mat          <- cbind(melt(fpr)[, 2:3], melt(tpr)$value)
    temp_mat$nullpi   <- summary_par_vals$nullpi[list_index]
    temp_mat$Nsamp    <- summary_par_vals$Nsamp[list_index]
    temp_mat$alt_type <- summary_par_vals$alt_type[list_index]
    colnames(temp_mat) <- c("Method", "fpr", "tpr", "nullpi", "Nsamp", "alt_type")

    roc_mat <- rbind(roc_mat, temp_mat)
}


pdf(file = "./fig/ave_roc.pdf", height = 6, width = 9)
for (current_alt_type in alt_type_seq[1:(length(alt_type_seq) - 1)]) {
submat <- dplyr::filter(roc_mat, alt_type == current_alt_type)
p <- ggplot(data = submat, mapping = aes(x = fpr, y = tpr, group = Method, color = Method)) +
    facet_grid(nullpi~Nsamp) +
    geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
    geom_line(alpha = I(0.7), size = I(0.3)) + xlab("FPR") + ylab("TPR") +
    ggtitle(paste("Pointwise Mean ROC when Alt Type =", current_alt_type)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
print(p)
}
dev.off()
