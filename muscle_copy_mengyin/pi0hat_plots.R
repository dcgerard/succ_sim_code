library(ggplot2)
library(reshape2)

pi0hat <- read.csv("pi0hat_df.csv")

longdat <- melt(pi0hat, measure.vars = 1:10)

pdf(file = "pi0hat_muscle.pdf", height = 20, width = 10)
ggplot(data = longdat, mapping = aes(x = nullpi, y = value, group = variable,
                                     color = variable)) +
    facet_grid(alt_type~Nsamp) +
    geom_point(pch = I(1)) +
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(hat(pi)[0])) +
    xlab(expression(pi[0]))
dev.off()


head(pi0hat)
