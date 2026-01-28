#' ATE in Percentage Points Under Subgroup Heterogeneity (Post-Estimation)
#'
#' @description
#' \code{ate_pct} is a post-estimation function for computing average treatment effects (ATEs)
#' in percentage points when treatment effects are heterogeneous across observable subgroups.
#'
#' The function extracts subgroup-specific average log-point treatment effects \eqn{\tau_g}
#' (the coefficients on subgroup indicator variables) and their covariance matrix from a fitted
#' model, constructs subgroup weights \eqn{w_g}, and reports three ATE-in-percentage-points
#' functionals: \eqn{\bar{\tau}}, \eqn{\rho_a}, and \eqn{\rho_b}. Standard errors are computed
#' using delta-method formulas that account for estimation uncertainty in subgroup effects and,
#' when applicable, in subgroup weights.
#'
#' @details
#' Let \eqn{\tau = (\tau_1,\ldots,\tau_G)'} denote subgroup average log-point effects and
#' \eqn{w = (w_1,\ldots,w_G)'} subgroup weights with \eqn{\sum_g w_g = 1}. The reported targets are:
#' \describe{
#'   \item{\code{taubar}}{\eqn{\bar{\tau} = \sum_g w_g \tau_g}.}
#'   \item{\code{rho_a}}{\eqn{\rho_a = \exp(\bar{\tau}) - 1}.}
#'   \item{\code{rho_b}}{\eqn{\rho_b = \sum_g w_g \exp(\tau_g) - 1}.}
#' }
#'
#' \strong{Weights.}
#' If \code{groupsize} is supplied, weights are constructed as
#' \eqn{w_g = N_g / \sum_h N_h} where \eqn{N_g} are the supplied subgroup sizes.
#' Otherwise, \code{ate_pct} uses the model's estimation sample (via \code{model.frame(model)})
#' and counts treated observations in each subgroup indicator to form \eqn{N_g}.
#' When subgroup sizes are computed from the data, subgroup indicators must be 0/1 and
#' mutually exclusive in the estimation sample.
#'
#' \strong{Inference.}
#' Standard errors are computed by the delta method. Weights are treated as estimated unless
#' \code{truew = TRUE}, in which case the weight-variance component is set to zero
#' (i.e., \eqn{\Sigma_w = 0}).
#'
#' @param model A fitted model object with \code{coef()} and \code{vcov()} methods
#'   (e.g., from \code{lm()}). The estimation sample is taken to be \code{model.frame(model)}.
#' @param group_vars Character vector of subgroup indicator variable names. These must
#'   correspond to coefficient names in \code{coef(model)} and columns present in the model frame.
#' @param truew Logical. If \code{TRUE}, treats subgroup weights as fixed/known
#'   and sets the weight-variance component to zero.
#' @param groupsize Optional numeric vector of subgroup sizes \eqn{N_g}, one per element of
#'   \code{group_vars} and in the same order. All values must be positive.
#'   If supplied, weights are \eqn{w_g = N_g / \sum_h N_h} and \eqn{p_T = (\sum_g N_g)/N}.
#'   \strong{Important:} \code{groupsize} should correspond to the same estimation sample
#'   used to fit \code{model}.
#' @param VC_model Optional variance-covariance matrix for the fitted model coefficients.
#'   If \code{NULL}, uses \code{vcov(model)}. Must be square and have row/column names
#'   containing \code{group_vars}.
#'
#' @return
#' An object of class \code{"ate_pct"}: a list containing
#' \describe{
#'   \item{\code{coefficients}}{Named numeric vector \code{c(taubar, rho_a, rho_b)}.}
#'   \item{\code{vcov}}{3x3 variance-covariance matrix for \code{coefficients} (delta-method).}
#'   \item{\code{resvec}}{Concatenation \code{c(taubar, rho_a, rho_b, sd_taubar, sd_rho_a, sd_rho_b)}.}
#'   \item{\code{G}}{Number of subgroups.}
#'   \item{\code{w}}{Group weights \eqn{w_g} (sum to 1).}
#'   \item{\code{N}}{Number of observations in the estimation sample.}
#'   \item{\code{N_T}}{Total treated count used to construct weights (\eqn{N_T=\sum_g N_g}).}
#'   \item{\code{p_T}}{Target-sample share \eqn{p_T = N_T/N}.}
#'   \item{\code{tau}}{Subgroup log-point effects \eqn{\tau_g}.}
#'   \item{\code{Sigma_tau}}{Estimated covariance matrix of \eqn{\tau}.}
#'   \item{\code{Sigma_w}}{Estimated covariance matrix of \eqn{w} (zero if \code{truew=TRUE}).}
#'   \item{\code{delta}}{Stacked vector \eqn{(\tau', w')'}.}
#'   \item{\code{Sigma_delta}}{Block-diagonal covariance matrix of \eqn{(\tau,w)}.}
#' }
#'
#' @section Methods:
#' The returned object supports \code{print()} and \code{summary()} methods:
#' \itemize{
#'   \item \code{print.ate_pct}: prints estimates, standard errors, z-stats, and p-values.
#'   \item \code{summary.ate_pct}: adds normal-approximation confidence intervals.
#' }
#'
#' @seealso
#' \code{\link[stats]{lm}}, \code{\link[stats]{coef}}, \code{\link[stats]{vcov}}
#'
#' @references
#' Zeng, Y. \emph{Estimation and Inference on Average Treatment Effects in Percentage Points under Heterogeneity}.
#' Working paper.
#'
#' @author
#' Ying Zeng \email{zengying17@@gmail.com}. School of Economics and WISE, Xiamen University.
#' Code repository: \url{https://github.com/zengying17/ate_pct-r}.
#'
#' @examples
#' ## Suppose variables gr1, gr2, gr3 are indicators for three sub-treatment groups and we are interested in the ATE in percentage points of these three treatment groups. 
#' ## First run a semi-log regression with subgroup-specific treatment effects, and have the heteroscedasticity-robust standard errors (optional).
#' ## (y must be log-transformed)
#' # library("sandwich")
#' # load('ate_pct_example.Rda')
#' ## Run the regression and have the heteroscedasticity-robust standard errors.
#' # reg_res <- lm(lny~x+gr1+gr2+gr3,data=df)
#' # VChet <- vcovHC(reg_res, type = "HC1")
#' 
#' Example 1: ATE in percentage points for the three groups, weights calculated with sample size. 
#' # out1 <- ate_pct(reg_res, c("gr1","gr2","gr3"), VC_model=VChet)
#' # out1
#' # summary(out1)
#'
#' ## Example 2: weights are treated as fixed (no sampling variation), and uses the conventional VC matrix for log points estimator.
#' # out2 <- ate_pct(reg_res, c("gr1","gr2","gr3"), truew=TRUE)
#'
#' ## Example 3: ATE in percentage points for the first two groups, and assgin equal weights (so weight are fixed and known.)
#' # out3 <- ate_pct(reg_res, c("gr1","gr2"), truew=TRUE, groupsize=c(1,1))
#'
#' ## Example 4: Provide subgroup sizes externally. The command below is this is exactly the same as \code{ate_pct gr1 gr2 gr3} as 15, 24, 37 are the sample size for the three treatment groups.
#' # out4 <- ate_pct(reg_res, c("gr1","gr2","gr3"), groupsize=c(15,24,37), VC_model=VChet)
#'
#' @export
#' @keywords models
ate_pct <- function(model, group_vars, truew = FALSE, groupsize = NULL, VC_model = NULL) {

  # ---- Dependencies: block-diagonal helper (Matrix if available, else base-R fallback)
  bdiag <- if (requireNamespace("Matrix", quietly = TRUE)) {
    Matrix::bdiag
  } else {
    function(...) {
      mats <- list(...)
      if (length(mats) == 1L && is.list(mats[[1L]])) mats <- mats[[1L]]
      mats <- lapply(mats, as.matrix)
      nr <- vapply(mats, nrow, 0L)
      nc <- vapply(mats, ncol, 0L)
      out <- matrix(0, nrow = sum(nr), ncol = sum(nc))
      r0 <- 0L; c0 <- 0L
      for (k in seq_along(mats)) {
        r1 <- r0 + nr[k]; c1 <- c0 + nc[k]
        out[(r0 + 1L):r1, (c0 + 1L):c1] <- mats[[k]]
        r0 <- r1; c0 <- c1
      }
      out
    }
  }
  
  # ---- Basic checks
  if (missing(model) || is.null(model)) stop("'model' must be a fitted model object.")
  if (missing(group_vars) || length(group_vars) == 0) stop("'group_vars' must be a non-empty character vector.")
  group_vars <- as.character(group_vars)

  beta <- stats::coef(model)
  if (is.null(names(beta))) stop("Model coefficients are unnamed; cannot match 'group_vars'.")

  if (!all(group_vars %in% names(beta))) {
    missing_vars <- group_vars[!group_vars %in% names(beta)]
    stop("The following group variables are not in the model coefficients: ",
         paste(missing_vars, collapse = ", "))
  }
  
  vcov_matrix <- if (is.null(VC_model)) stats::vcov(model) else VC_model
  vcov_matrix <- as.matrix(vcov_matrix)

  if (nrow(vcov_matrix) != ncol(vcov_matrix)) stop("'VC_model' / vcov(model) must be square.")
  if (is.null(rownames(vcov_matrix)) || is.null(colnames(vcov_matrix))) {
    stop("'VC_model' / vcov(model) must have rownames and colnames.")
  }
    # Symmetrize to avoid pure-asymmetry numerical artifacts
  vcov_matrix <- (vcov_matrix + t(vcov_matrix)) / 2

  if (!all(group_vars %in% rownames(vcov_matrix)) || !all(group_vars %in% colnames(vcov_matrix))) {
    missing_r <- setdiff(group_vars, rownames(vcov_matrix))
    missing_c <- setdiff(group_vars, colnames(vcov_matrix))
    msg <- c()
    if (length(missing_r) > 0) msg <- c(msg, paste0("rownames missing: ", paste(missing_r, collapse = ", ")))
    if (length(missing_c) > 0) msg <- c(msg, paste0("colnames missing: ", paste(missing_c, collapse = ", ")))
    stop("The following group variables are missing in 'VC_model' / vcov(model): ", paste(msg, collapse = "; "))
  }
  
  mf <- stats::model.frame(model)   # exact estimation sample
  N  <- nrow(mf)
  if (N == 0) stop("Empty model frame after NA handling.")
  G <- length(group_vars)
  
  ## safely fetch a full data.frame to source extra columns (for control), if possible
   if (!all(group_vars %in% names(mf))) {
    missing_in_mf <- group_vars[!group_vars %in% names(mf)]
    stop("The following group variables are not in the model frame: ",
         paste(missing_in_mf, collapse = ", "),
         ". Ensure they are included as regressors.")
  }
group_mat <- mf[, group_vars, drop = FALSE]
  colnames(group_mat) <- group_vars
  tau <- as.numeric(beta[group_vars])
  names(tau) <- group_vars
  Sigma_tau <- as.matrix(vcov_matrix[group_vars, group_vars, drop = FALSE]) 
  
  ## ---- Counts and weights (within treated share)
  if (!is.null(groupsize)) {
    if (length(groupsize) != G) {
      stop(sprintf("groupsize must have length %d (one N_g per group, in the same order as group_vars).", G))
    }
    if (any(!is.finite(groupsize)) || any(groupsize <= 0)) {
      stop("All elements of groupsize must be finite and > 0.")
    }
    N_g <- as.numeric(groupsize)
        names(N_g) <- group_vars

    N_T <- sum(N_g)
    if (N_T <= 0) stop("Sum of groupsize must be > 0.")
    if (N_T > N) {
      warning("Sum(groupsize) > N in the estimation sample. Ensure groupsize corresponds to the model's estimation sample.")
    }

  } else {
    # Validate group indicators are 0/1 in the estimation sample
      for (gv in group_vars) {
        vals <- unique(group_mat[[gv]])
        vals <- vals[!is.na(vals)]
        if (!all(vals %in% c(0, 1))) {
          stop(sprintf("Group indicators must be 0/1. Variable '%s' has values: %s",
                      gv, paste(sort(vals), collapse = ", ")))
        }
      }

      s_i <- rowSums(as.matrix(group_mat), na.rm = TRUE)
      if (any(s_i > 1)) {
        stop("Group indicators must be mutually exclusive within the target sample (found observations with multiple groups).")
      }  
      N_g <- vapply(group_vars, function(gv) sum(!is.na(group_mat[[gv]]) & group_mat[[gv]] == 1), 0L)  
      if (any(N_g == 0)) {
        bad <- group_vars[N_g == 0]
        stop("Some groups have zero treated observations so w=0: ", paste(bad, collapse = ", "),
            ". Drop these groups or redefine group_vars.")
      }
          names(N_g) <- group_vars
    N_T <- sum(N_g)
    if (N_T == 0) stop("No treated observations (all treatment dummies are 0).")  
  }

    p_T <- N_T / N  
  w <- as.numeric(N_g / N_T)
  names(w) <- group_vars  
  
  if (isTRUE(truew)) {
    Sigma_w <- matrix(0, nrow = G, ncol = G)
  } else {
    diag_w <- diag(w, nrow = G)
    ww     <- w %*% t(w)
    Sigma_w <- (1 / N) * (1 / p_T) * (diag_w - ww)
  }
  rownames(Sigma_w) <- colnames(Sigma_w) <- paste0("w_", group_vars)

# Calculate basic statistics 
calc_basic <- function(tau, w, Sigma_w, Sigma_tau) {
    G <- length(tau)     
    # eta0 = log w + tau; Sigma_eta = Var(eta0) = diag(1/w) Sigma_w diag(1/w) + Sigma_tau
    eta0 <- log(w) + tau
    Sigma_eta <- diag(1 / w, nrow = G) %*% Sigma_w %*% diag(1 / w, nrow = G) + Sigma_tau
    
    taubar <- sum(w * tau)
    Var_taubar <- drop(t(w) %*% Sigma_tau %*% w + t(tau) %*% Sigma_w %*% tau)
    
    rho_a <- exp(taubar) - 1
    Var_rho_a <- (exp(taubar))^2 * Var_taubar

    rho_b <- sum(w * exp(tau)) - 1
    exp_eta0 <- exp(eta0)   # = exp(log w + tau) = w * exp(tau)
    Var_b <- drop(t(exp_eta0) %*% Sigma_eta %*% exp_eta0)
    
    # Assemble coefficients
    b <- c(
      taubar = taubar,
      rho_a  = rho_a,
      rho_b  = rho_b
    )
    V <- diag(c(Var_taubar, Var_rho_a, Var_b))
    
    
    rownames(V) <- colnames(V) <- names(b)

    # Guard against small negative variances (numerical / non-PSD vcov):
    vdiag <- diag(V)
    if (any(vdiag < -1e-12)) {
      warning("Some delta-method variances are negative. This suggests a non-PSD vcov matrix (or numerical issues). SEs are computed after truncating variances at 0.")
    }
    vdiag <- pmax(vdiag, 0)
    diag(V) <- vdiag
    se <- sqrt(vdiag)
    names(se) <- c("sd_taubar", "sd_rho_a", "sd_rho_b")

    # Keep named delta pieces without mutating inputs
    w_named   <- w;   names(w_named)   <- paste0("w_",   names(w_named))
    tau_named <- tau; names(tau_named) <- paste0("tau_", names(tau_named))

    Sigma_delta <- as.matrix(bdiag(Sigma_tau, Sigma_w))
    colnames(Sigma_delta) <- rownames(Sigma_delta) <- c(names(tau_named), names(w_named))
    
    list(
      coefficients = b,
      vcov = V,
      resvec = c(b, se),
      delta = c(tau_named, w_named),
      Sigma_delta = Sigma_delta,
      G = G,
      w = w
    )    
  }
  

  result <- calc_basic(tau, w, Sigma_w, Sigma_tau)
  

  # Extra returned components
  result$N   <- N
  result$N_T <- N_T
  result$p_T <- p_T
  result$tau <- tau
  result$Sigma_tau <- Sigma_tau
  result$Sigma_w   <- Sigma_w

  class(result) <- "ate_pct"
  result
}

#' @export
print.ate_pct <- function(x, ...) {
  cat("Average Treatment Effect in Percentage Points\n\n")

  coef <- x$coefficients
  se <- sqrt(diag(x$vcov))
  z <- coef / se
  p <- 2 * stats::pnorm(-abs(z))

  tab <- data.frame(
    Estimate = coef,
    Std.Error = se,
    z.value = z,
    Pr.z = p,
    check.names = FALSE
  )

  print(tab)
  invisible(x)
}

#' @export
summary.ate_pct <- function(object, conf.level = 0.95, ...) {
  coef <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- coef / se
  p <- 2 * stats::pnorm(-abs(z))

  alpha <- 1 - conf.level
  crit  <- stats::qnorm(1 - alpha / 2)
  ci_l  <- coef - crit * se
  ci_u  <- coef + crit * se

  tbl <- data.frame(
    Estimate   = coef,
    Std.Error  = se,
    z.value    = z,
    Pr.z       = p,
    CI.Lower   = ci_l,
    CI.Upper   = ci_u,
    check.names = FALSE
  )

  ans <- list(
    table = tbl,
    conf.level = conf.level,
    N = object$N,
    N_T = object$N_T,
    p_T = object$p_T,
    w = object$w,
    delta_head = utils::head(object$delta),
    Sigma_delta_dim = dim(object$Sigma_delta)
  )
  class(ans) <- "summary.ate_pct"
  ans
}

#' @export
print.summary.ate_pct <- function(x, digits = 4, ...) {
  cat("Average Treatment Effect Percentage Calculation\n\n")
  cat("Coefficients:\n")
  print(round(x$table, digits))
  cat("\n")

  if (!is.null(x$w)) {
    cat("Group weights w (treated share decomposition):\n")
    print(round(x$w, 6))
    cat("\n")
  }

  cat(sprintf("N = %d, N_T = %d, p_T = %.6f\n", x$N, x$N_T, x$p_T))

  if (!is.null(x$delta_head)) {
    cat("\nDelta vector (head):\n")
    print(round(x$delta_head, 6))
  }
  if (!is.null(x$Sigma_delta_dim)) {
    cat("\nSigma_delta (dim): ", paste(x$Sigma_delta_dim, collapse = " x "), "\n", sep = "")
  }

  invisible(x)
}
