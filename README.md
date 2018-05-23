#Copy number estimation from logR using ASCAT formula
=====

Requires logR values and previously defined segments (doesn't do segmentation!)


##Setting up the model

You need to provide:

- a vector of logR values
- number of segments
- indices of start of each segment
- indices of end of each segment
- prior parameters



Build the model and initialise it with data:    
    
    jmod <- jags.model("copynumber.bug",
        data = list(
            logr = logr_vector,
            nsegments = number_of_segments,
            starts = left_hand_position_of_breakpoints,
            ends = right_hand_position_of_breakpoints,
            psi_a = 4,    # psi_a, psi_b: parameters
            psi_b = 2,    # of gamma prior on psi (ploidy)
            rho_a = 1,    # rho_a, rho_b: parameters
            rho_b = 1,    # of beta prior on rho (purity)
            sd = 0.1),    # standard dev of normal distribution
        n.chains = 1,   # number of MCMC chains to run
        n.adapt = 100)  # number of cycles to spend tuning the model before running
        

## Running the model

    # burn-in 1000 steps
    update(jmod, 1000)
    
    # sample 1000 iterations of the MCMC (n.iter / thin = 1000)
    jsamples <- jags.samples(
        jmod,
        c("CNplus1", "rho", "psi"),
        n.iter = 1000,
        thin = 1)


## Getting the results

    CN.est <- summary(jsamples$CNplus1, mean)$stat - 1
    rho.est <- summary(jsamples$rho, mean)$stat
    psi.est <- summary(jsamples$psi, mean)$stat
