library(ggplot2)
library(reshape2)
library(tidyr)

pi0_mat <- read.csv("../muscle/pi0_mat.csv")
mse_mat <- read.csv("../muscle/mse_mat.csv")
auc_mat <- read.csv("../muscle/auc_mat.csv")

pi0_ruvash_adinf <- read.csv("pi0_ruvash_adinf.csv")
mse_ruvash_adinf <- read.csv("mse_ruvash_adinf.csv")
auc_ruvash_adinf <- read.csv("auc_ruvash_adinf.csv")

mult_mat <- read.csv("mult_ruvash_adinf.csv")
add_mat <- read.csv("add_ruvash_adinf.csv")

pi0_mat <- cbind(pi0_mat, pi0_ruvash_adinf[, -(1:4)])
mse_mat <- cbind(mse_mat, mse_ruvash_adinf[, -(1:4)])
auc_mat <- cbind(auc_mat, auc_ruvash_adinf[, -(1:4)])
auc_mat <- auc_mat[auc_mat$nullpi != 1, ]

pdf(file = "muscle_tissue_adinf.pdf")
## pi0_plot
nullpi_seq <- unique(pi0_mat$nullpi)
nsamp_seq <- unique(pi0_mat$Nsamp)
dummy_dat <- expand.grid(nullpi_seq, nsamp_seq)
colnames(dummy_dat) <- c("nullpi", "Nsamp")

name_vec <- colnames(pi0_mat)
colnames(pi0_mat) <- gsub("pi0hat_", "", x = name_vec)
long_dat <- melt(pi0_mat, id.vars = c("nullpi", "Nsamp"),
                 measure.vars = colnames(pi0_mat)[5:ncol(pi0_mat)])
p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
p <- p + facet_grid(nullpi~Nsamp) + ylab(expression(hat(pi)[0]))
p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
p <- p + geom_hline(data = dummy_dat, aes(yintercept = nullpi), lty = 2, color = "red", lwd = 1)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
p <- p + ggtitle("Estimates of pi0 When Using Muscle Tissue")
print(p)


## auc_plot
nullpi_seq <- unique(pi0_mat$nullpi)
nsamp_seq <- unique(pi0_mat$Nsamp)
dummy_dat <- expand.grid(nullpi_seq, nsamp_seq)
colnames(dummy_dat) <- c("nullpi", "Nsamp")
med_mat <- matrix(NA, nrow = (length(nullpi_seq) - 1) * length(nsamp_seq), ncol = ncol(auc_mat) - 4)
for (index in 5:ncol(auc_mat)) {
    form1 <- as.formula(paste(colnames(auc_mat)[index], "~ nullpi + Nsamp"))
    out1 <- aggregate(form1, FUN = median, na.rm = TRUE,
                      data = auc_mat)
    med_mat[, index - 4] <- out1[, 3]
}
dummy_dat2 <- cbind(expand.grid(nullpi_seq[nullpi_seq != 1], nsamp_seq), apply(med_mat, 1, max))
colnames(dummy_dat2) <- c("nullpi", "Nsamp", "max_med")

name_vec <- colnames(auc_mat)
colnames(auc_mat) <- gsub("auc_", "", x = name_vec)
long_dat <- melt(auc_mat, id.vars = c("nullpi", "Nsamp"),
                 measure.vars = colnames(auc_mat)[5:ncol(auc_mat)])
p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
p <- p + facet_grid(nullpi~Nsamp) + ylab("AUC")
p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
p <- p + geom_hline(data = dummy_dat2, aes(yintercept = max_med), lty = 2, color = "red", lwd = 1)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
p <- p + ggtitle("AUC When Using Muscle Tissue")
print(p)

## mse_plot
name_vec <- colnames(mse_mat)
colnames(mse_mat) <- gsub("mse_", "", x = name_vec)
long_dat <- melt(mse_mat, id.vars = c("nullpi", "Nsamp"),
                 measure.vars = colnames(mse_mat)[5:ncol(mse_mat)])
p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
p <- p + facet_grid(nullpi~Nsamp) + ylab("MSE")
p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
p <- p + ylim(0, max(mse_mat[, -c(1:4, 9)], na.rm = TRUE))
p <- p + ggtitle("MSE When Using Muscle Tissue")
print(p)


long_dat <- melt(add_mat, id.vars = c("Nsamp", "nullpi"),
                 measure.vars = colnames(add_mat)[5:ncol(add_mat)])
p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
p <- p + facet_grid(nullpi~Nsamp) + ylab("Additive Inflation")
p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
p <- p + ggtitle("Additive Inflation When Using Muscle Tissue")
print(p)

long_dat <- melt(mult_mat, id.vars = c("Nsamp", "nullpi"),
                 measure.vars = colnames(mult_mat)[5:ncol(mult_mat)])
p <- ggplot(data = long_dat, mapping = aes(x = variable, y = value, fill = variable))
p <- p + facet_grid(nullpi~Nsamp) + ylab("Mult Inflation")
p <- p + geom_boxplot(size = 0.4, outlier.size = 0.4)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
p <- p + ggtitle("Mult Inflation When Using Muscle Tissue")
print(p)

dev.off()
