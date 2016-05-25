library(ggplot2)
library(reshape2)

pi0hat <- read.csv("pi0hat_df.csv")
pi0hat_ruvash <- read.csv("pi0hat_df_ruvash.csv")
pi0hat$ruvash <- pi0hat_ruvash$ruvash
pi0hat <- dplyr::select(pi0hat, -c(leapp_sparse, ols, ash, leapp_ridge))


longdat <- melt(pi0hat, measure.vars = c(1:6, 12))

pdf(file = "pi0hat_muscle.pdf", height = 20, width = 10)
ggplot(data = longdat, mapping = aes(x = nullpi, y = value, group = variable,
                                     color = variable)) +
    facet_grid(alt_type~Nsamp) +
    geom_point(pch = I(1)) +
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(hat(pi)[0])) +
    xlab(expression(pi[0]))
dev.off()


ggplot(data = pi0hat, mapping = aes(x = nullpi, y = ruvash, color = I("blue"))) +
    facet_grid(alt_type~Nsamp) +
    geom_point(pch = I(1)) +
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(hat(pi)[0])) +
    xlab(expression(pi[0])) +
    ggtitle("RUVASH pi0hat.")
