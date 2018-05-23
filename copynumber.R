library(jags)
purity <- 0.6
ploidy <- 2.18
fuzz = 0.5  # apply Gaussian noise ~ N(0, fuzz) to ideal logr

ideal_logr <- function(purity, cn, ploidy) {
    log2( (2 * (1-purity) + purity * cn) / (2 * (1 - purity) + purity * ploidy))
}

d <- data.frame(sizes = c(1000, 500, 1500, 1000, 3000, 1000, 1000, 2000, 300, 6000, 1000),
                cns   = c(   2,   0,    1,    4,    2,    5,    1,    0,   7,    2,    3))
logr <- unlist(apply(d, 1, function(row) rnorm(row[["sizes"]], ideal_logr(purity, row[["cns"]], ploidy), fuzz)))

plot(logr, pch = 20, cex = .75, col = "firebrick", ylim = c(-4.5, 2))
segment_borders <- c(0, cumsum(d$sizes))
abline(v = segment_borders, col = "grey90")
centres <- segment_borders[1:(length(segment_borders) - 1)] + diff(segment_borders) / 2
text(d$cns, x = centres, y = -3, col = "black")

jmod <- jags.model("~/code/stan/copynumber/copynumber.bug",
                   data = list(logr = logr,
                               nsegments = nrow(d),
                               starts = c(1, 1+cumsum(d$sizes)[1:(length(d$sizes-2))]),
                               ends = cumsum(d$sizes),
                               psi_a = 4,  # psi_a, psi_b: parameters
                               psi_b = 2,  # of gamma prior on psi (ploidy)
                               rho_a = 1,  # rho_a, rho_b: parameters
                               rho_b = 1,  # of beta prior on rho (purity)
                               sd = 0.1),  # standard deviation of normal distribution
                   n.chains = 1,
                   n.adapt = 100)

# burn-in
update(jmod, 1000)

jsamples <- jags.samples(jmod,
                         c("CNplus1", "rho", "psi"),
                         1000, thin = 1)

jsamples

text(as.character(round(summary(jsamples$CNplus1, mean)$stat - 1, 2)), x = centres, y = -4, col = "orange")
