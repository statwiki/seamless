
## =============================================================================
## Wang et al. (2023) - Rank-Based Approach for Inferential Seamless Phase 2/3
## R Simulation Code
## =============================================================================

## Install required packages if not already installed
## install.packages(c("MASS", "Mediana", "mvtnorm"))


## install.packages('mvtnorm', repos='http://cran.us.r-project.org', lib = "C:/apps/R/local-lib", dependencies = TRUE)

library(tzdb, lib.loc = "C:/apps/R/local-lib")
library(withr, lib.loc = "C:/apps/R/local-lib")
library(utf8, lib.loc = "C:/apps/R/local-lib")
library(readr, lib.loc = "C:/apps/R/local-lib")



library(dplyr, lib.loc = "C:/apps/R/local-lib")
library(tidyr, lib.loc = "C:/apps/R/local-lib")
library(haven, lib.loc = "C:/apps/R/local-lib")
library(survival, lib.loc = "C:/apps/R/local-lib")
library(rmarkdown, lib.loc = "C:/apps/R/local-lib")


library(MASS, lib.loc = "C:/apps/R/local-lib")  ## for mvrnorm()
library(Mediana, lib.loc = "C:/apps/R/local-lib") ## for AdjustPvalues() with DunnettAdj
library(mvtnorm, lib.loc = "C:/apps/R/local-lib") ## for pmvnorm() if needed
# =============================================================================
# SECTION 1: HELPER FUNCTIONS
# =============================================================================

#' Build the 2m x 2m Correlation Matrix (Table 1 in paper)
#' For m doses vs 1 control, balanced design (rho_e = rho_b = 0.5)
#'
#' @param m    Number of treatment doses
#' @param rho  Correlation between efficacy and biomarker endpoints (test statistics) 
#' @return A 2m x 2m correlation matrix
#' Order: E1, E2, ..., Em, B1, B2, ..., Bm

build_sigma <- function(m, rho) {
  total_dim <- 2 * m
  Sigma <- matrix(0, nrow = total_dim, ncol = total_dim)
  
  ## Correlation between efficacy test statistics (rho_e = 0.5 for balanced)
  rho_e <- 0.5
  ## Correlation between biomarker test statistics (rho_b = 0.5 for balanced)
  rho_b <- 0.5
  ## Cross-correlation between E_j and B_j (same dose): rho
  rho_eb1 <- rho
  ## Cross-correlation between E_j and B_p (different doses): 0.5 * rho
  rho_eb2 <- 0.5 * rho
  
  for (i in 1:total_dim) {
    for (j in 1:total_dim) {
      if (i == j) {
        Sigma[i, j] <- 1  # diagonal
      } else if (i <= m && j <= m) {
        # Efficacy vs Efficacy (different doses)
        Sigma[i, j] <- rho_e
      } else if (i > m && j > m) {
        # Biomarker vs Biomarker (different doses)
        Sigma[i, j] <- rho_b
      } else {
        # Cross: Efficacy vs Biomarker
        ei <- ifelse(i <= m, i, i - m)  # dose index for row
        bj <- ifelse(j <= m, j, j - m)  # dose index for col
        if (ei == bj) {
          Sigma[i, j] <- rho_eb1  # same dose
        } else {
          Sigma[i, j] <- rho_eb2  # different doses
        }
      }
    }
  }
  
  return(Sigma)
}


#' Rank-Based Dunnett Adjusted P-value
#' Uses Mediana::AdjustPvalues() with DunnettAdj procedure
#'
#' @param p01  The unadjusted p-value for the selected dose
#' @param r    Rank of the selected dose's biomarker response (1 = lowest, m = highest)
#' @return Adjusted p-value

rank_based_dunnett_pval <- function(p01, r) {
  # Construct pval vector: place p01 at position 1, fill remaining r-1 with 1
  # This reduces the multiplicity adjustment to dimension r
  pval <- c(p01, rep(1, r - 1))
  
  adj <- AdjustPvalues(
    pval,
    proc       = "DunnettAdj",
    par        = parameters(n = Inf)  # Normal approximation (df = Inf)
  )
  
  return(adj[1])  # Return adjusted p-value for the selected dose
}


#' Sidak Adjusted P-value
#' p_sidak = 1 - (1 - p)^m
#'
#' @param p01  Unadjusted p-value
#' @param m    Total number of doses
#' @return Sidak adjusted p-value

sidak_pval <- function(p01, m) {
  return(1 - (1 - p01)^m)
}


#' Traditional Dunnett Adjusted P-value
#' Equivalent to rank-based Dunnett with r = m
#'
#' @param p01  Unadjusted p-value
#' @param m    Total number of doses
#' @return Dunnett adjusted p-value

traditional_dunnett_pval <- function(p01, m) {
  rank_based_dunnett_pval(p01, m)
}


#' Weighted Inverse Normal Combination Test (Lehmacher & Wassmer, 1999)
#' Combines p1 (Phase 2) and p2 (Phase 3) into a single global p-value
#'
#' @param p1  Phase 2 adjusted p-value
#' @param p2  Phase 3 p-value
#' @param w1  Pre-specified weight for Phase 2 (default: n1/(n1+n2))
#' @return Combined global p-value

## function to return the final p-value
combine_pvalues <- function(p1, p2, w1) {
  w2 <- sqrt(1 - w1^2)
  z_combined <- w1 * qnorm(1 - p1) + w2 * qnorm(1 - p2)  ## NOTE: w2 = sqrt(1-w1)
  ## Per formula (2): w1*Phi^-1(1-p1) + sqrt(1-w1)*Phi^-1(1-p2)
  ## Re-implement exactly per paper:
  combined <- 1 - pnorm(sqrt(w1) * qnorm(1 - p1) + sqrt(1 - w1) * qnorm(1 - p2))
  return(combined)
}

# =============================================================================
# SECTION 2: SIMULATION PARAMETERS
# =============================================================================

set.seed(100)

m  <- 3       ## Number of treatment doses
B  <- 10000   ## Number of simulation runs (paper uses 50,000)
n1 <- 50      ## Phase 2 sample size per arm
n2 <- 300     ## Phase 3 sample size per arm (selected dose + control)
alpha <- 0.025  ## One-sided significance level

## Pre-specified weight for Phase 2 in the combination test
## w1 = n1 / (n1 + n2) = 50 / (50 + 300) = 1/7
w1 <- n1 / (n1 + n2)
cat(sprintf("Combination weight w1 = %.4f\n", w1))

## Correlation values to evaluate (rho from 0 to 1, step 0.1)
## correlation bewteen biomarker and efficacy endpoint 
rho_seq <- seq(0, 1, by = 0.1)

# =============================================================================
# SECTION 3: TYPE I ERROR SIMULATION (Table 2 in paper)
# =============================================================================

cat("\n=== TABLE 2: FWER Simulation Results ===\n")
cat(sprintf("%-6s | %6s %6s %6s | %6s %6s %6s | %6s %6s %6s\n",
            "rho",
            "Sid_r1", "Sid_r2", "Sid_r3",
            "Dun_r1", "Dun_r2", "Dun_r3",
            "RBD_r1", "RBD_r2", "RBD_r3"))
cat(strrep("-", 75), "\n")

ptm <- Sys.time()

fwer_results <- data.frame()

for (rho in rho_seq) {
  
  ## Build the 2m x 2m correlation matrix (null hypothesis: mu = 0)
  sig <- build_sigma(m, rho)
  
  ## Simulate Phase 2 test statistics (efficacy + biomarker for m doses)
  ## Columns 1:m = efficacy stats (E1, E2, E3)
  ## Columns (m+1):(2m) = biomarker stats (B1, B2, B3)
  z1 <- mvrnorm(n = B, mu = rep(0, 2 * m), Sigma = sig)
  ## Assuming test statistics follow normal distribution 
  E_phase2 <- z1[, 1:m]          ## Efficacy test statistics in Phase 2
  B_phase2 <- z1[, (m+1):(2*m)]  ## Biomarker test statistics in Phase 2
  
  # Simulate Phase 3 efficacy test statistics (independent, standard normal under H0)
  z2 <- rnorm(n = B)
  
  ## Compute p-values for Phase 3
  p2 <- 1 - pnorm(z2)
  
  ## Initialize rejection counters
  reject_sidak  <- rep(0, m)
  reject_dunnett <- rep(0, m)
  reject_rbd    <- rep(0, m)
    
    for (r in 1:m) {
    
    ## For each simulation run, assume the selected dose has biomarker rank r
    ## i.e., the dose whose biomarker statistic is the r-th order statistic
    
    ## Get the efficacy statistic of the dose with biomarker rank r
    ## Sort biomarker stats and find corresponding efficacy stats
    biomarker_ranks <- apply(B_phase2, 1, rank)  # rank within each row
    
    ## For rank r: find the column index whose biomarker rank == r
    selected_col <- apply(biomarker_ranks, 2, function(x) which(x == r))
    ## Note: biomarker_ranks is m x B (transposed from apply), need to handle:
    biomarker_ranks_T <- t(apply(B_phase2, 1, rank))  # B x m matrix
    
    ## Efficacy statistic for the dose with biomarker rank r (across B simulations)
    e_jr <- sapply(1:B, function(i) {
      col_idx <- which(biomarker_ranks_T[i, ] == r)
      E_phase2[i, col_idx[1]]
    })
    
    ## Unadjusted p-value for Phase 2 (one-sided)
    p1_raw <- 1 - pnorm(e_jr)
    
    ## ---- Sidak Adjustment ----
    p1_sidak <- sapply(p1_raw, function(p) sidak_pval(p, m))
    p_combined_sidak <- mapply(function(p1, p2_i) combine_pvalues(p1, p2_i, w1),
                               p1_sidak, p2)
    reject_sidak[r] <- mean(p_combined_sidak < alpha)
    
    ## ---- Traditional Dunnett Adjustment ----
    p1_dunnett <- sapply(p1_raw, function(p) traditional_dunnett_pval(p, m))
    p_combined_dunnett <- mapply(function(p1, p2_i) combine_pvalues(p1, p2_i, w1),
                                 p1_dunnett, p2)
    reject_dunnett[r] <- mean(p_combined_dunnett < alpha)
    
    ## ---- Rank-Based Dunnett Adjustment ----
    p1_rbd <- sapply(p1_raw, function(p) rank_based_dunnett_pval(p, r))
    p_combined_rbd <- mapply(function(p1, p2_i) combine_pvalues(p1, p2_i, w1),
                             p1_rbd, p2)
    reject_rbd[r] <- mean(p_combined_rbd < alpha)
  }
  
  cat(sprintf("%.1f   | %6.4f %6.4f %6.4f | %6.4f %6.4f %6.4f | %6.4f %6.4f %6.4f\n",
              rho,
              reject_sidak[1], reject_sidak[2], reject_sidak[3],
              reject_dunnett[1], reject_dunnett[2], reject_dunnett[3],
              reject_rbd[1], reject_rbd[2], reject_rbd[3]))
  
  fwer_results <- rbind(fwer_results, data.frame(
    rho = rho,
    Sidak_r1 = reject_sidak[1], Sidak_r2 = reject_sidak[2], Sidak_r3 = reject_sidak[3],
    Dunnett_r1 = reject_dunnett[1], Dunnett_r2 = reject_dunnett[2], Dunnett_r3 = reject_dunnett[3],
    RBD_r1 = reject_rbd[1], RBD_r2 = reject_rbd[2], RBD_r3 = reject_rbd[3]
  ))
}

Sys.time() - ptm

## =============================================================================
## SECTION 4: EMPIRICAL CRITICAL VALUES (Table 3 in paper)
## =============================================================================

ptm <- Sys.time()
cat("\n=== TABLE 3: Empirical Critical Values (one-sided alpha = 0.025) ===\n")
cat(sprintf("%-6s | %8s %8s %8s | %8s | %8s %8s %8s\n",
            "rho",
            "EJr_r1", "EJr_r2", "EJr_r3",
            "Sidak",
            "RBD_r1", "RBD_r2", "RBD_r3"))
cat(strrep("-", 80), "\n")

B_ecdf <- 10000  # Paper uses 50,000 for empirical cutoffs

## Sidak cutoffs (fixed, independent of rho)
sidak_cutoff_1 <- qnorm(1 - alpha / m)                 # m=1 equivalent
sidak_cutoff_2 <- qnorm(1 - (1 - (1-alpha)^(1/2)))    # m=2 equivalent -- simplified
sidak_cutoff_3 <- qnorm(1 - (1 - (1-alpha)^(1/3)))    # m=3

## Dunnett cutoffs using Mediana (r=1,2,3)
dunnett_cutoffs <- sapply(1:m, function(r) {
  ## Find z such that P(T(r) > z) = alpha under Dunnett
  ## Equivalent: find quantile of Dunnett adjusted distribution
  ## Use uniroot approach
  uniroot(function(z) {
    p_raw <- 1 - pnorm(z)
    p_adj <- rank_based_dunnett_pval(p_raw, r)
    p_adj - alpha
  }, interval = c(1, 4))$root
})

for (rho in rho_seq) {
  sig <- build_sigma(m, rho)
  z1 <- mvrnorm(n = B_ecdf, mu = rep(0, 2 * m), Sigma = sig)
  
  E_phase2 <- z1[, 1:m]
  B_phase2 <- z1[, (m+1):(2*m)]
  biomarker_ranks_T <- t(apply(B_phase2, 1, rank))
  
  ejr_cutoffs <- numeric(m)
  
  for (r in 1:m) {
    # Efficacy statistic for dose with biomarker rank r
    e_jr <- sapply(1:B_ecdf, function(i) {
      col_idx <- which(biomarker_ranks_T[i, ] == r)
      E_phase2[i, col_idx[1]]
    })
    # Empirical 97.5th percentile = critical value for one-sided alpha=0.025
    ejr_cutoffs[r] <- quantile(e_jr, 1 - alpha)
  }
  
  cat(sprintf("%.1f   | %8.4f %8.4f %8.4f | %8.4f | %8.4f %8.4f %8.4f\n",
              rho,
              ejr_cutoffs[1], ejr_cutoffs[2], ejr_cutoffs[3],
              ifelse(rho %in% c(0, 0.2, 1.0),
                     c(sidak_cutoff_1, sidak_cutoff_2, sidak_cutoff_3)[3], NA),
              dunnett_cutoffs[1], dunnett_cutoffs[2], dunnett_cutoffs[3]))
}

Sys.time() - ptm 

## =============================================================================
## SECTION 5: POWER SIMULATION (Table 4 in paper)
## =============================================================================

cat("\n=== TABLE 4: Power Simulation Results ===\n")

## Power scenario: binary biomarker, PFS-based efficacy
## Response rates: RRc=0.5, RRt1=0.6, RRt2=0.65, RRt3=0.7
## HR=0.7 under alternative; 331 events needed for 90% power
## Log-rank test stat ~ N(log(HR) * sqrt(E/4), 1)  :: normal approximation

HR   <- 0.7
E    <- 331   ## expected events for 90% power
RRc  <- 0.5
RRt  <- c(0.6, 0.65, 0.7)  # biomarker response rates for doses 1, 2, 3

B

ptm <- Sys.time()
## Non-centrality parameter for log-rank test (Phase 2 + Phase 3 combined)
## Phase 2: 50 subjects/arm => expected events depends on event rate (assume ~60%)
event_rate <- 0.60
e_phase2   <- n1 * event_rate  # events in Phase 2 per arm
e_phase3   <- n2 * event_rate  # events in Phase 3 per arm

## Non-centrality for Phase 2 log-rank: ncp = -log(HR) * sqrt(e/4)
ncp_phase2 <- -log(HR) * sqrt(e_phase2 / 4)  # positive = alternative

## Non-centrality for Phase 3 log-rank
ncp_phase3 <- -log(HR) * sqrt(e_phase3 / 4)

cat(sprintf("Phase 2 log-rank NCP: %.4f\n", ncp_phase2))
cat(sprintf("Phase 3 log-rank NCP: %.4f\n\n", ncp_phase3))

cat(sprintf("%-6s | %6s %6s %6s | %6s %6s %6s | %6s %6s %6s\n",
            "rho",
            "Sid_r1", "Sid_r2", "Sid_r3",
            "Dun_r1", "Dun_r2", "Dun_r3",
            "RBD_r1", "RBD_r2", "RBD_r3"))
cat(strrep("-", 75), "\n")

power_results <- data.frame()

for (rho in rho_seq) {
  
  sig <- build_sigma(m, rho)
  
  ## Simulate under alternative: mean vector is non-centrality parameters
  ## Efficacy means (log-rank NCP for each dose): same for all doses under alt
  ## Biomarker means: derived from normal approx of binary test statistic
  ## For binary endpoint: E[B_j] = (RRtj - RRc) / sqrt(RRc*(1-RRc)/n + RRtj*(1-RRtj)/n)
  
  biomarker_ncp <- sapply(RRt, function(RRtj) {
    se <- sqrt(RRc * (1 - RRc) / n1 + RRtj * (1 - RRtj) / n1)
    (RRtj - RRc) / se
  })
  
  # Mean vector: [E1_ncp, E2_ncp, E3_ncp, B1_ncp, B2_ncp, B3_ncp]
  mu_alt <- c(rep(ncp_phase2, m), biomarker_ncp)
  
  # Simulate test statistics under alternative
  z1 <- mvrnorm(n = B, mu = mu_alt, Sigma = sig)
  
  E_phase2 <- z1[, 1:m]
  B_phase2 <- z1[, (m+1):(2*m)]
  biomarker_ranks_T <- t(apply(B_phase2, 1, rank))
  
  # Phase 3: simulate efficacy under alternative
  z2 <- rnorm(n = B, mean = ncp_phase3, sd = 1)
  p2 <- 1 - pnorm(z2)
  
  reject_sidak   <- rep(0, m)
  reject_dunnett <- rep(0, m)
  reject_rbd     <- rep(0, m)
  
  for (r in 1:m) {
    e_jr <- sapply(1:B, function(i) {
      col_idx <- which(biomarker_ranks_T[i, ] == r)
      E_phase2[i, col_idx[1]]
    })
    
    p1_raw <- 1 - pnorm(e_jr)
    
    # Sidak
    p1_sidak <- sapply(p1_raw, function(p) sidak_pval(p, m))
    p_comb_sidak <- mapply(combine_pvalues, p1_sidak, p2, MoreArgs = list(w1 = w1))
    reject_sidak[r] <- mean(p_comb_sidak < alpha)
    
    # Traditional Dunnett
    p1_dunnett <- sapply(p1_raw, function(p) traditional_dunnett_pval(p, m))
    p_comb_dunnett <- mapply(combine_pvalues, p1_dunnett, p2, MoreArgs = list(w1 = w1))
    reject_dunnett[r] <- mean(p_comb_dunnett < alpha)
    
    # Rank-Based Dunnett
    p1_rbd <- sapply(p1_raw, function(p) rank_based_dunnett_pval(p, r))
    p_comb_rbd <- mapply(combine_pvalues, p1_rbd, p2, MoreArgs = list(w1 = w1))
    reject_rbd[r] <- mean(p_comb_rbd < alpha)
  }
  
  cat(sprintf("%.1f   | %6.4f %6.4f %6.4f | %6.4f %6.4f %6.4f | %6.4f %6.4f %6.4f\n",
              rho,
              reject_sidak[1], reject_sidak[2], reject_sidak[3],
              reject_dunnett[1], reject_dunnett[2], reject_dunnett[3],
              reject_rbd[1], reject_rbd[2], reject_rbd[3]))
  
  power_results <- rbind(power_results, data.frame(
    rho = rho,
    Sidak_r1 = reject_sidak[1], Sidak_r2 = reject_sidak[2], Sidak_r3 = reject_sidak[3],
    Dunnett_r1 = reject_dunnett[1], Dunnett_r2 = reject_dunnett[2], Dunnett_r3 = reject_dunnett[3],
    RBD_r1 = reject_rbd[1], RBD_r2 = reject_rbd[2], RBD_r3 = reject_rbd[3]
  ))
}

cat("\nSimulation complete.\n")

Sys.time() - ptm
