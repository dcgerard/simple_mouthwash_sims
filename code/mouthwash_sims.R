# devtools::install_github("dcgerard/seqgendiff")
library(seqgendiff)
# devtools::install_github("dcgerard/vicar")
library(vicar)
library(sva)
library(pROC)
library(limma)

one_rep <- function(new_params, current_params) {
  args_val <- append(current_params, new_params)
  set.seed(new_params$current_seed)

  ## Choose all of the genes because already got top expressed
  stopifnot(args_val$Ngene == ncol(args_val$mat))
  d_out <- seqgendiff::poisthin(mat = args_val$mat,
                                nsamp = args_val$Nsamp,
                                ngene = args_val$Ngene,
                                skip_gene = args_val$skip_gene,
                                signal_params = list(mean = 0, sd = args_val$log2foldsd),
                                gvec = rep(TRUE, length(args_val$Ngene)),
                                gselect = "custom",
                                prop_null = args_val$nullpi,
                                alpha = 0)

  which_null <- abs(d_out$beta) < 10 ^ -6
  nnull         <- sum(which_null)
  control_genes <- which_null
  control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE

  beta_true <- d_out$beta

  X <- d_out$X
  colnames(X) <- c("Intercept", "Treatment")
  Y <- log2(d_out$Y + 1)

  num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

  ## ASH on VOOM-LIMMA-EBAYES estimates
  voom_out   <- limma::voom(counts = t(Y), design = X)
  limma_out  <- limma::lmFit(object = voom_out)
  ebayes_out <- limma::ebayes(fit = limma_out)
  betahat    <- limma_out$coefficients[, 2]
  sebetahat  <- sqrt(ebayes_out$s2.post) * limma_out$stdev.unscaled[, 2]
  df         <- ebayes_out$df.total[1]
  pvalues    <- ebayes_out$p.value[, 2]
  ashout     <- ashr::ash(betahat = betahat, sebetahat = sebetahat, df = df)

  ## Mouthwash
  mout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, cov_of_interest = 2)

  ## Return summaries ------------------------------------------------
  return_vec <- rep(NA, 6)
  names(return_vec) <- c("mouth_pi0", "ash_pi0",
                         "mouth_auc", "ash_auc",
                         "mouth_mse", "ash_mse")

  ## pi0hat ----------------------------------------------------------
  return_vec[1] <- mout$pi0
  return_vec[2] <- ashr::get_pi0(ashout)

  ## auc ------------------------------------------------------------
  return_vec[3] <- pROC::roc(predictor = mout$result$lfdr, response = which_null)$auc
  return_vec[4] <- pROC::roc(predictor = ashr::get_lfdr(ashout), response = which_null)$auc

  ## mse ------------------------------------------------------------
  return_vec[5] <- mean((mout$result$PosteriorMean - beta_true) ^ 2)
  return_vec[6] <- mean((ashr::get_pm(ashout) - beta_true) ^ 2)

  return(return_vec)
}

itermax <- 100
seed_start <- 1

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(6, 10, 20, 40)
ncontrol_seq <- c(100)

par_vals <- expand.grid(list((1 + seed_start):(itermax + seed_start),
                             nullpi_seq, Nsamp_seq, ncontrol_seq))
colnames(par_vals) <- c("current_seed", "nullpi", "Nsamp", "ncontrols")
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
  par_list[[list_index]] <- list()
  for (inner_list_index in 1:ncol(par_vals)) {
    par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
    names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
  }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 0.8
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## Create muscle_mat with most expressed genes
mat <- t(as.matrix(read.csv("./data/muscle.csv",
                            header = TRUE)[, -c(1,2)]))
args_val$mat <- mat[, order(apply(mat, 2, median), decreasing = TRUE)[1:args_val$Ngene]]
rm(mat)

# Number of threads to use for multithreaded computing. This must be
# specified in the command-line shell; e.g., to use 8 threads, run
# command
#
#  R CMD BATCH '--args nc=8' mouthwash_sims.R
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}

## run once -------------------------------------------------------
one_rep(new_params = par_list[[1]], current_params = args_val)

# run in parallel -------------------------------------------------
library(parallel)
cl   <- makeCluster(nc)
cat("Running multithreaded computations with",nc,"threads.\n")
sout <- t(parallel::parSapply(cl = cl, X = par_list, FUN = one_rep,
                              current_params = args_val))
stopCluster(cl)

saveRDS(cbind(par_vals, sout), "./output/sims_out.Rds")
