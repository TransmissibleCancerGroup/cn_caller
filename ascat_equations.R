na <- function(rho, ploidy, logr, vaf) {
    (rho - 1 + 2^logr * (1 - vaf) * (2 * (1 - rho) + rho * ploidy)) / rho
}

nb <- function(rho, ploidy, logr, vaf) {
    (rho - 1 + 2^logr * vaf * (2 * (1 - rho) + rho * ploidy)) / rho
}

