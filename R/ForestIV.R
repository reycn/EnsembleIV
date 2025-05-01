#' ForestIV Main Function
#'
#' This function implements the main ForestIV approach.
#'
#' @param data_test Testing dataframe for random forest, must have a column named "actual" that contains the ground truth, and all trees' predictions.
#' @param data_unlabel Unlabel dataframe for random forest, must have all trees' predictions.
#' @param control A character vector of control variable names. Pass an empty vector if there are no control variables
#' @param method "Lasso" for ForestIV method and "IIV" for EnsembleIV method.
#' @param iterative Whether to perform iterative IV selection or not, default to TRUE. Only relevant when method = "Lasso"
#' @param ntree Number of trees in the random forest.
#' @param model_unbias Unbiased estimation.
#' @param family Model specification, same as in the family parameter in glm.
#' @param diagnostic Whether to output diagnostic correlations for instrument validity and strength, default to TRUE.
#' @param select_method method of IV selection. One of "optimal" (LASSO based), "top3", and "PCA".
#' @return ForestIV estimation results
#' @export


ForestIV = function(data_test, data_unlabel, control, method, iterative = TRUE, ntree, model_unbias, family, diagnostic, select_method="PCA") {
  run_per_tree <- function(data_test, data_unlabel, control, method, iterative = TRUE, ntree, model_unbias, family, diagnostic, select_method="PCA", i) {
    # use i-th tree as the endogenous covariate
    regressor = paste0("X", i)
    if (method == "Lasso") {
      output = LassoSelect(data_test, data_unlabel, iterative, ntree, regressor)
      IVs = output$IVs
      data_unlabel_new = data_unlabel
    }
    if (method == "IIV") {
      output = IIVSelect(data_test, data_unlabel, ntree, regressor, select_method)
      IVs = output$IVs
      data_unlabel_new = output$data_unlabel_new
    }

    # Estimate 2SLS
    # TODO: family parameter unnecessarily messy for now. Will clean up after deciding what to do with GLMs
    if (length(IVs) > 0) {
      if (family$family == "gaussian" & family$link == "identity") {
        f = make_formula(regressor, control, IVs, "all")
        model_IV = AER::ivreg(f, data = data_unlabel_new)
        beta_IV = stats::coef(model_IV)
        vcov_IV = stats::vcov(model_IV)
        se_IV = sqrt(diag(vcov_IV))
        # the convergence code is not useful for linear case - 2SLS always converges
        convergence = 0
      } else {
        link = switch(family$family,
          "binomial" = "logit",
          "poisson" = "logadd",
        )
        # implementation via OneSampleMR, only supports logit and poisson
        f = make_formula(regressor, control, IVs, "all")
        model_IV = OneSampleMR::tsri(f, data = data_unlabel_new, link = link)
        # tsri produces a lot of estimates, what we care about starts with (Intercept) and excludes "resres"
        allcoefs = stats::coef(model_IV$fit)
        start = which(names(allcoefs) == "(Intercept)")
        varlist = c(start, start + 1, (start + 3):length(allcoefs))
        beta_IV = stats::coef(model_IV$fit)[varlist]
        vcov_IV = stats::vcov(model_IV$fit)[varlist, varlist]
        se_IV = sqrt(diag(vcov_IV))
        # note that the TSRI estimation may not converge - we want to keep only the convergent ones when reporting results
        convergence = model_IV$fit$algoInfo$convergence

        # implementation via ivtools, support all glm families
        # f1 = make_formula(regressor, control, IVs, "XZ")
        # fitX.LZ = stats::glm(f1, data = data_unlabel_new)
        # f2 = make_formula(regressor, control, IVs, "YX")
        # fitY.LX = stats::glm(f2, data = data_unlabel_new, family = family)
        # model_IV = ivtools::ivglm(estmethod = "ts", fitX.LZ = fitX.LZ, fitY.LX = fitY.LX, data = data_unlabel_new, ctrl = TRUE)
        # # for control function estimation, a variable named "R" is attached to the end and need to be removed
        # nvar = length(model_IV$est)
        # beta_IV = model_IV$est[-nvar]
        # vcov_IV = model_IV$vcov[-nvar,-nvar]
        # se_IV = sqrt(diag(vcov_IV))
      }
      H_stats = hotelling(beta_IV, vcov_IV, model_unbias)
      corrs = output$correlations
      
      # print(i)
      return(c(beta_IV, se_IV, H_stats, convergence, corrs))
    }
  }
  result = list()
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' needed for parallel execution. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' needed for parallel execution. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' needed for parallel execution. Please install it.", call. = FALSE)
  }
  numCores <- parallel::detectCores() - 1
  doParallel::registerDoParallel(numCores) # use multicore, set to the number of our cores
  result <- foreach::foreach(i = 1:ntree, .packages = c("AER", "OneSampleMR", "stats")) %dopar% {
    run_per_tree(data_test, data_unlabel, control, method, iterative, ntree, model_unbias, family, diagnostic, select_method, i)
  }
  print(paste0("Parallel executing... Number of cores: ", length(numCores)))
  result = do.call(rbind.data.frame, result)
  print(paste0("Parallel execution complete. Number of results: ", length(result)))

  # Handle cases where no results were generated
  if (is.null(result) || nrow(result) == 0) {
    warning("No valid results obtained from any tree. Returning NULL.")
    # Return NULL or an appropriately structured empty data frame
    return(NULL) 
  }
  
  # Calculate the number of coefficients (k) from the result dimensions
  # Expecting 2k (betas+SEs) + 1 (Hotelling) + 1 (Convergence) + 4 (corrs) = 2k + 6 columns
  num_extra_cols = 6 
  if (ncol(result) < num_extra_cols || (ncol(result) - num_extra_cols) %% 2 != 0) {
      stop("Unexpected number of columns in the results matrix. Check run_per_tree output.")
  }
  num_coeffs = (ncol(result) - num_extra_cols) / 2

  if (diagnostic) {
    # Keep all 2k + 6 columns
    colnames(result) = c(paste0("beta_", c(1:num_coeffs)), paste0("se_", c(1:num_coeffs)), "Hotelling", "Convergence", "pp_abs_before", "pe_abs_before", "pp_abs_after", "pe_abs_after")
    return(result)
  }
  else {
    # Keep only the first 2k + 2 columns (betas, SEs, Hotelling, Convergence)
    num_cols_to_keep = 2 * num_coeffs + 2
    result = result[, 1:num_cols_to_keep, drop = FALSE] # Use drop=FALSE for robustness if num_coeffs=0 somehow
    colnames(result) = c(paste0("beta_", c(1:num_coeffs)), paste0("se_", c(1:num_coeffs)), "Hotelling", "Convergence")
    return(result)
  }
}


