library(rjags)
purity <- 0.6
ploidy <- 2.18
fuzz = 0.5  # apply Gaussian noise ~ N(0, fuzz) to ideal logr

ideal_logr <- function(purity, cn, ploidy) {
    log2( (2 * (1-purity) + purity * cn) / (2 * (1 - purity) + purity * ploidy))
}

ideal_vaf <- function(purity, cn_a, cn_b) {
    ((1 - purity + purity * cn_b) /
     (2 - 2 * purity + purity * (cn_a + cn_b)))
}

make_symmetric <- function(v) {
    s <- sample(1:length(v), length(v) / 2, replace = FALSE)
    v[s] <- 1 - v[s]
    v
}

mirror <- function(v) {
    0.5 - abs(0.5 - v)
}

d <- data.frame(sizes = c(1000, 500, 1500, 1000, 3000, 1000, 1000, 2000, 300, 6000, 1000),
                cnt   = c(   2,   0,    1,    4,    2,    5,    1,    0,   7,    2,    3),
                cna   = c(   1,   0,    1,    3,    2,    3,    0,    0,   7,    1,    1),
                cnb   = c(   1,   0,    0,    1,    0,    2,    1,    0,   0,    1,    2))
logr <- unlist(apply(d, 1, function(row) rnorm(row[["sizes"]], ideal_logr(purity, row[["cnt"]], ploidy), fuzz)))

vaf <- unlist(apply(d, 1, function(row) {
                        ivaf <- ideal_vaf(purity, row[["cna"]], row[["cnb"]])
                        rbeta(row[["sizes"]], ivaf * 25, (1 - ivaf) * 25)
                    }))

par(mfrow = c(2, 1))
plot(logr, pch = 20, cex = .75, col = "firebrick", ylim = c(-4.5, 2))
plot(make_symmetric(vaf), pch = 20, cex = .5, col = "firebrick", ylim = c(0, 0))
plot((vaf), pch = 20, cex = .5, col = "firebrick", ylim = c(0, 0))
segment_borders <- c(0, cumsum(d$sizes))
abline(v = segment_borders, col = "grey90")
centres <- segment_borders[1:(length(segment_borders) - 1)] + diff(segment_borders) / 2
text(d$cna, x = centres, y = 0.0, col = "black")
text(d$cnb, x = centres, y = 0.5, col = "black")
par(mfrow = c(1, 1))

jmod <- jags.model("~/code/jags/copynumber/copynumber.bugs",
                   data = list(logr = logr,
                               vaf = mirror(vaf),
                               nsegments = nrow(d),
                               starts = c(1, 1+cumsum(d$sizes)[1:(length(d$sizes-2))]),
                               ends = cumsum(d$sizes),
                               psi_a = 40,  # psi_a, psi_b: parameters
                               psi_b = 20,  # of gamma prior on psi (ploidy)
                               rho_a = 1,  # rho_a, rho_b: parameters
                               rho_b = 1), # of beta prior on rho (purity)
                   n.chains = 1,
                   n.adapt = 100)

# burn-in
update(jmod, 100)

jsamples <- jags.samples(jmod,
                         c("cn_a", "cn_b", "rho", "psi", "normal_sd", "beta_conc"),
                         100, thin = 1)

jsamples

text(as.character(round(summary(jsamples$cn_a, mean)$stat - 1, 2)), x = centres, y = -0.5, col = "orange")
text(as.character(round(summary(jsamples$cn_b, mean)$stat - 1, 2)), x = centres, y = -1.0, col = "orange")
