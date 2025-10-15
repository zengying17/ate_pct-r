#' Calculation of Average Treatment Effect in Percentage Points
#' Ying Zeng "Estimation and Inference of Average Treatment Effect in Percentage Points"
#' Calculates various average treatment effect percentage measures after running a regression model
#' where group effects are estimated.
#'
#' @param model A fitted model object (e.g., from lm())
#' @param group_vars Character vector of group indicator variable names in the model
#' @param within Logical: indicating whether to calculate ATE in p.c.t allowing for within-group heterogeneity statistics (default: FALSE)
#' @param corr Correlation parameter for within-group calculation (default: 0)
#' @param truew Logical: whether w should be treated as true (variance=0) (default: FALSE)
#' @param VC_model Optional variance-covariance matrix for the fitted model. If NULL, uses vcov(model).
#' @param control Optional character vector of variable names of for dummies of control groups, when calculating the ATE in p.c.t allowing for within-group heterogeneity statistics (default: NULL).  name for the control group. If NULL, control group is a simgle group which has all treatment dummies = 0.
#' @return A list containing:
#'   \item{coefficients}{Named vector of point estimates}
#'   \item{vcov}{Variance-covariance matrix of the reported estimates (delta-method)}
#'   \item{G}{Number of groups}
#'   \item{w}{Group weights (sum to 1)}
#'   \item{s}{Group-specific standard deviations (if within=TRUE)}
#'   \item{Vartau}{Variance vector of varsigma (if within=TRUE)}
#'
#' @examples
#' # model <- lm(log(y) ~ x1 + x2 + d_1 + d_2 + d_3, data = data) #make sure that y is log transformed. 
#' # result <- ate_pct(model, c("d_1", "d_2", "d_3"))
#' # summary(result)
ate_pct <- function(model, group_vars, truew = FALSE, within = FALSE, corr = 1,
                    VC_model = NULL, control = NULL) {
  # ---- Dependencies
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Please install 'Matrix'.")
  bdiag <- Matrix::bdiag
  
  # ---- Basic checks
  beta <- stats::coef(model)
  if (!all(group_vars %in% names(beta))) {
    missing_vars <- group_vars[!group_vars %in% names(beta)]
    stop("The following group variables are not in the model coefficients: ",
         paste(missing_vars, collapse = ", "))
  }
  
  vcov_matrix <- if (is.null(VC_model)) {
    stats::vcov(model)
  } else {
    if (!all(group_vars %in% colnames(VC_model))) {
      missing_vars <- group_vars[!group_vars %in% colnames(VC_model)]
      stop("The following group variables are not in 'VC_model': ",
           paste(missing_vars, collapse = ", "))
    }
    VC_model
  }
  
  # ---- Use the exact estimation sample
  mf <- stats::model.frame(model)                # exact rows used in estimation
  N  <- nrow(mf)
  if (N == 0) stop("Empty model frame after NA handling.")
  G <- length(group_vars)
  
  ## safely fetch a full data.frame to source extra columns (for control), if possible
  original_data <- NULL
  if (!is.null(model$call$data) && is.name(model$call$data)) {
    nm <- as.character(model$call$data)
    if (exists(nm, envir = parent.frame())) {
      original_data <- get(nm, envir = parent.frame())
    } else if (exists(nm, envir = .GlobalEnv)) {
      original_data <- get(nm, envir = .GlobalEnv)
    }
  }
  
  ## Align any extra columns to the estimation sample by row names when possible
  align_to_mf_rows <- function(x) {
    if (is.null(x)) return(NULL)
    if (!is.null(rownames(mf)) && !is.null(rownames(x)) &&
        all(rownames(mf) %in% rownames(x))) {
      x[rownames(mf), , drop = FALSE]
    } else {
      x
    }
  }
  
  ## Build a working data.frame that contains the group dummies used by the model
  # Prefer columns from mf; if missing there (shouldn't be), try to pull from original_data
  get_col <- function(nm) {
    if (nm %in% names(mf)) return(mf[[nm]])
    if (!is.null(original_data) && nm %in% names(original_data)) {
      od <- align_to_mf_rows(as.data.frame(original_data))
      return(od[[nm]])
    }
    stop(sprintf("Column '%s' not found in the model frame or original data.", nm))
  }
  
  group_mat <- as.data.frame(lapply(group_vars, get_col))
  colnames(group_mat) <- group_vars
  tau <- as.numeric(beta[group_vars])
  names(tau) <- group_vars
  Sigma_tau <- as.matrix(vcov_matrix[group_vars, group_vars, drop = FALSE]) 
  
  ## ---- Counts and weights (within treated share)
  s_i <- rowSums(as.matrix(group_mat), na.rm = TRUE)
  if (any(s_i > 1)) warning("Some rows have multiple treatment dummies == 1.")
  
  N_g <- vapply(group_vars, function(gv) sum(!is.na(group_mat[[gv]]) & group_mat[[gv]] == 1), 0L)  
  N_s <- sum(s_i == 1)
  if (N_s == 0) stop("No treated observations (all treatment dummies are 0).")
  p_s <- N_s / N
  
  w <- as.numeric(N_g / N_s)
  names(w) <- group_vars  
  
  # ---- VCV of weights (design-based, multinomial approx) 
  if (isTRUE(truew)) {
    Sigma_w <- matrix(0, nrow = G, ncol = G)
  } else {
    diag_w <- diag(w, nrow = G)
    ww     <- w %*% t(w)
    # Guard against division by zero if p_s == 0 (already checked)
    Sigma_w <- (1 / N) * (1 / p_s) * (diag_w - ww)
  }
  rownames(Sigma_w) <- colnames(Sigma_w) <- paste0("w_", group_vars)
  
  # ---- Residual d.f. for small-sample corrections used in rho_d
  
  # Calculate basic statistics - define calc_basic function inside
  calc_basic <- function(tau, w, Sigma_w, Sigma_tau) {
    G <- length(tau)
    # taubar and its variance
    
    # eta0 = log w + tau; Sigma_eta = Var(eta0) = diag(1/w) Sigma_w diag(1/w) + Sigma_tau
    eta0 <- log(w) + tau
    Sigma_eta <- diag(1 / w, nrow = G) %*% Sigma_w %*% diag(1 / w, nrow = G) + Sigma_tau
    
    taubar <- sum(w * tau)
    Var_taubar <- drop(t(w) %*% Sigma_tau %*% w + t(tau) %*% Sigma_w %*% tau)
    
    # rho_a = exp(taubar) - 1
    rho_a <- exp(taubar) - 1
    Var_rho_a <- (exp(taubar))^2 * Var_taubar
    
    
    # rho_b = sum w * exp(tau) - 1, with delta-method variance
    rho_b <- sum(w * exp(tau)) - 1
    exp_eta0 <- exp(eta0)   # = exp(log w + tau) = w * exp(tau)
    Var_b <- drop(t(exp_eta0) %*% Sigma_eta %*% exp_eta0)
    
    # rho_c = E-bias-corrected via lognormal approximation using diag(Sigma_tau)
    rho_c <- sum(w * exp(tau - 0.5 * diag(Sigma_tau))) - 1   
    
    
    # Assemble coefficients
    b <- c(
      taubar = taubar,
      rho_a  = rho_a,
      rho_b  = rho_b,
      rho_c  = rho_c
    )
    
    se <- c(
      sd_taubar = sqrt(max(0, Var_taubar)),
      sd_rho_a  = sqrt(max(0, Var_rho_a)),
      sd_rho_b  = sqrt(max(0, Var_b)),
      sd_rho_c  = sqrt(max(0, Var_b))
    )
    # Report diagonal VCV for the displayed stats (as in your original)
    V <- diag(se^2)
    rownames(V) <- colnames(V) <- names(b)
    
    # Delta vector and block VCV used
    names(w) <- paste0("w_", names(w))
    names(tau) <- paste0("tau_", names(tau)) 
    
    Sigma_delta <- bdiag(Sigma_tau, Sigma_w)
    Sigma_delta <- as.matrix(Sigma_delta)
    colnames(Sigma_delta) <- rownames(Sigma_delta) <- c(names(tau), names(w))
    
    list(
      coefficients = b,
      vcov = V,
      resvec = c(b, se),
      delta = c(tau, w),
      Sigma_delta = as.matrix(Sigma_delta),
      G = G,
      w = unname(w),
      p_s = p_s
    )    
  }
  
  # ---------- WITHIN-GROUP HETEROGENEITY ----------
calc_within <- function(tau, w, Sigma_w, Sigma_tau, s, p_s, N, corr, control_mode = c("single","matched"), N0 = NULL) {
  control_mode <- match.arg(control_mode)
  G <- length(tau)

  # ----- (unchanged) basic pieces reused in both branches
  eta0      <- log(w) + tau
  Sigma_eta <- diag(1 / w, nrow = G) %*% Sigma_w %*% diag(1 / w, nrow = G) + Sigma_tau
  taubar    <- sum(w * tau)
  Var_taubar <- drop(t(w) %*% Sigma_tau %*% w + t(tau) %*% Sigma_w %*% tau)
  rho_a     <- exp(taubar) - 1
  Var_rho_a <- (exp(taubar))^2 * Var_taubar
  rho_b     <- sum(w * exp(tau)) - 1
  Var_b     <- drop(t(exp(eta0)) %*% Sigma_eta %*% exp(eta0))
  rho_c     <- sum(w * exp(tau - 0.5 * diag(Sigma_tau))) - 1

  # mean-adjusted eta (same as your code)
  eta       <- log(w) + tau - 0.5 * (diag(Sigma_w) / (w^2)) - 0.5 * diag(Sigma_tau)

  if (control_mode == "single") {
    # === Your original single-control path (G + 1)
    s_g <- s[1:G]
    s_0 <- s[G + 1]
    Vartau   <- s_g^2 + s_0^2 - 2 * corr * s_g * s_0
    eta_plus <- eta + 0.5 * Vartau

    rho_b_plus <- sum(w * exp(tau + 0.5 * Vartau)) - 1
    rho_c_plus <- sum(w * exp(tau + 0.5 * Vartau - 0.5 * diag(Sigma_tau))) - 1

    # Sigma_s diagonal, expressed via p_s and N (same numerics as your N_g/N_0 code)
    Sigma_s <- matrix(0, G + 1, G + 1)
    for (i in 1:G) Sigma_s[i, i] <- 2 * s_g[i]^4 / (p_s * w[i] * N)
    Sigma_s[G + 1, G + 1] <- 2 * s_0^4 / (if (!is.null(N0)) N0 else (1 - p_s) * N)

    # block VCV: (tau, w, s_g^2[1..G], s0^2_scalar)
    Sigma_delta <- matrix(0, 3*G + 1, 3*G + 1)
    Sigma_delta[1:G, 1:G] <- Sigma_tau
    Sigma_delta[(G+1):(2*G), (G+1):(2*G)] <- Sigma_w
    Sigma_delta[(2*G+1):(3*G+1), (2*G+1):(3*G+1)] <- Sigma_s

    # gradients
    nabla_s0_Vartau <- 1 - (corr * s_g) / s_0
    nabla_sg_Vartau <- 1 - (corr * s_0) / s_g
    exp_eta_plus    <- exp(eta_plus)

    nabla_delta_rho_plus <- c(
      exp_eta_plus,                  # d/d tau
      exp_eta_plus / w,              # d/d w
      0.5 * nabla_sg_Vartau * exp_eta_plus,             # d/d s_g^2 (length G)
      0.5 * sum(nabla_s0_Vartau * exp_eta_plus)         # d/d s0^2 (scalar)
    )
    Var_b_plus <- drop(t(nabla_delta_rho_plus) %*% Sigma_delta %*% nabla_delta_rho_plus)

    b  <- c(taubar = taubar, rho_a = rho_a, rho_b = rho_b, rho_c = rho_c,
            rho_b_plus = rho_b_plus, rho_c_plus = rho_c_plus)
    se <- sqrt(pmax(0, c(Var_taubar, Var_rho_a, Var_b, Var_b, Var_b_plus, Var_b_plus)))
    V  <- diag(se^2); colnames(V) <- rownames(V) <- names(b)

    names(w)  <- paste0("w_", names(w))
    names(tau) <- paste0("tau_", names(tau))
    names(se) <- paste0("sd_",names(b))
    s2 <- c(s_g^2, s_0^2)
    names(s2) <- c(paste0("s2_", names(tau)), "s2_control")
    colnames(Sigma_delta) <- rownames(Sigma_delta) <- c(names(tau), names(w), names(s2))

    return(list(
      coefficients = b, vcov = V,
      resvec = c(b, se),
      # CHANGE: keep only the named delta (remove duplicate)
      delta = c(setNames(tau, names(tau)),
                setNames(w,   names(w)),
                setNames(s2,  names(s2))),
      Sigma_delta = Sigma_delta,
      G = G, w = unname(w)

    ))
  } else {
    # === Matched controls (2G): s = (s_g[1..G], s0_g[1..G])
    stopifnot(length(s) == 2*G)
    s_g <- s[1:G]
    s_0 <- s[(G+1):(2*G)]
    Vartau   <- s_g^2 + s_0^2 - 2 * corr * s_g * s_0
    eta_plus <- eta + 0.5 * Vartau

    rho_b_plus <- sum(w * exp(tau + 0.5 * Vartau)) - 1
    rho_c_plus <- sum(w * exp(tau + 0.5 * Vartau - 0.5 * diag(Sigma_tau))) - 1

    # diagonal Sigma_s using p_s and N; per-control bins use (1-p_s)
    Sigma_s <- matrix(0, 2*G, 2*G)
    for (i in 1:G) {
      Sigma_s[i, i]       <- 2 * s_g[i]^4 / (p_s * w[i] * N)
      Sigma_s[G+i, G+i]   <- 2 * s_0[i]^4 / ((1 - p_s) * N)
    }

    # block VCV: (tau, w, s_g^2[1..G], s0_g^2[1..G])
    Sigma_delta <- matrix(0, 4*G, 4*G)
    Sigma_delta[1:G, 1:G] <- Sigma_tau
    Sigma_delta[(G+1):(2*G), (G+1):(2*G)] <- Sigma_w
    Sigma_delta[(2*G+1):(3*G), (2*G+1):(3*G)] <- Sigma_s[1:G, 1:G]
    Sigma_delta[(3*G+1):(4*G), (3*G+1):(4*G)] <- Sigma_s[(G+1):(2*G), (G+1):(2*G)]

    # gradients (no scalar sum in the last block)
    nabla_sg_Vartau <- 1 - (corr * s_0) / s_g
    nabla_s0_Vartau <- 1 - (corr * s_g) / s_0
    exp_eta_plus    <- exp(eta_plus)
    nabla_delta_rho_plus <- c(
      exp_eta_plus,                  # d/d tau
      exp_eta_plus / w,              # d/d w
      0.5 * nabla_sg_Vartau * exp_eta_plus,  # d/d s_g^2 (G)
      0.5 * nabla_s0_Vartau * exp_eta_plus   # d/d s0_g^2 (G)
    )
    Var_b_plus <- drop(t(nabla_delta_rho_plus) %*% Sigma_delta %*% nabla_delta_rho_plus)

    b  <- c(taubar = taubar, rho_a = rho_a, rho_b = rho_b, rho_c = rho_c,
            rho_b_plus = rho_b_plus, rho_c_plus = rho_c_plus)
    se <- sqrt(pmax(0, c(Var_taubar, Var_rho_a, Var_b, Var_b, Var_b_plus, Var_b_plus)))
    V  <- diag(se^2); colnames(V) <- rownames(V) <- names(b)

    names(w)  <- paste0("w_", names(w))
    names(tau) <- paste0("tau_", names(tau))
    names(se) <- paste0("sd_", names(b))
    s2 <- c(s_g^2, s_0^2)
    names(s2) <- c(paste0("s2_", names(tau)), paste0("s2c_", names(tau)))
    colnames(Sigma_delta) <- rownames(Sigma_delta) <- c(names(tau), names(w), names(s2))

    return(list(
      coefficients = b, vcov = V,
      resvec = c(b,se),
      delta = c(setNames(tau, names(tau)),
                setNames(w,   names(w)),
                setNames(s2,  names(s2))),
      Sigma_delta = Sigma_delta,
      G = G, w = unname(w),p_s = p_s
    ))
  }
}

  
  # ---------- Compute outputs ----------
  if (isTRUE(within)) {
  resid_sq <- stats::residuals(model)^2

  # s_g per treatment group (unchanged logic)
  s <- numeric(G + 1)  # keep length for single by default
  for (g in 1:G) {
    idx <- !is.na(group_mat[[g]]) & (group_mat[[g]] == 1)
    if (!any(idx)) stop(sprintf("Group '%s' has zero observations in the estimation sample.", group_vars[g]))
    s[g] <- sqrt(mean(resid_sq[idx]))
  }

  # -- NEW (minimal): decide control mode from `control`
  control_mode <- "single"
  ctrl_df <- NULL
  if (is.null(control)) {
    control_mode <- "single"
  } else if (is.character(control)) {
    if (length(control) == 1L) {
      control_mode <- "single"
      ctrl_df <- as.data.frame(lapply(control, get_col))  # 1 column
    } else if (length(control) == G) {
      control_mode <- "matched"
      ctrl_df <- as.data.frame(lapply(control, get_col))  # G columns
    } else {
      stop("`control` must be NULL, one name, or exactly G names (matched controls).")
    }
  } else {
    stop("`control` must be NULL or a character vector of names.")
  }

  if (control_mode == "single") {
    # classic control: all dummies == 0 if ctrl_df is NULL, else use that single control column
    control_indices <-
      if (is.null(ctrl_df)) {
        (rowSums(as.matrix(group_mat), na.rm = TRUE) == 0)
      } else {
        (as.numeric(ctrl_df[[1L]]) != 0)  # CHANGE: robust to 1/TRUE/factors
      }
    if (!any(control_indices)) stop("No observations found for control group.")
    s[G + 1] <- sqrt(mean(resid_sq[control_indices]))

N0 <- sum(control_indices)
result <- calc_within(tau, w, Sigma_w, Sigma_tau, s, p_s, N, corr,
                      control_mode = "single", N0 = N0)
  } else {
    # matched controls: one per group
    if (ncol(ctrl_df) != G) stop("Internal: matched control must have exactly G columns.")
    s0 <- numeric(G)
    for (g in 1:G) {
      idx0 <- (as.numeric(ctrl_df[[g]]) != 0)  # CHANGE: robust
      if (!any(idx0)) stop(sprintf("Matched control '%s' has zero observations.", colnames(ctrl_df)[g]))
      s0[g] <- sqrt(mean(resid_sq[idx0]))
    }
    s <- c(s[1:G], s0)  # length 2G
    N0_g <- vapply(seq_len(G), function(g) sum(as.numeric(ctrl_df[[g]]) != 0), 0L)
result <- calc_within(tau, w, Sigma_w, Sigma_tau, s, p_s, N, corr,
                      control_mode = "matched", N0 = N0_g)
  }
} else {
  result <- calc_basic(tau, w, Sigma_w, Sigma_tau)
}

# NEW: return p_s like Stata e(p_s)
  result$p_s <- p_s
  
  class(result) <- "ate_pct"
  result
}

#' Print method for ate_pct objects
print.ate_pct <- function(x, ...) {
  cat("Average Treatment Effect in Percentage Points \n\n")
  
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

#' Summary method for ate_pct objects
summary.ate_pct <- function(object, conf.level = 0.95, ...) {
  cat("Average Treatment Effect Percentage Calculation\n\n")
  
  # Calculate statistics
  coef <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- coef / se
  p <-  2 * stats::pnorm(-abs(z))
  
  # Calculate confidence intervals
  alpha <- 1 - conf.level
  crit <- qnorm(1 - alpha/2)
  crit  <- stats::qnorm(1 - alpha / 2)
  ci_l  <- coef - crit * se
  ci_u  <- coef + crit * se
  
  # Create formatted table
  result_table <- data.frame(
    Estimate = round(coef, 4),
    `Std. Error` = round(se, 4),
    `z value` = round(z, 3),
    `Pr(>|z|)` = round(p, 4),
    `Lower CI` = round(ci_l, 4),
    `Upper CI` = round(ci_u, 4),
    check.names = FALSE,
    row.names = names(coef)
  )
  
  cat("Coefficients:\n")
  print(result_table, digits = 4)  
  cat("\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  if (!is.null(object$w)) {
    cat("\nGroup weights w (treated share decomposition):\n")
    print(round(object$w, 6))
  }
  if (!is.null(object$s)) {
    cat("\nResidual SDs by group (and control):\n")
    print(round(object$s, 6))
  }
  if (!is.null(object$Vartau)) {
    cat("\nWithin-group Vartau (per group):\n")
    print(round(object$Vartau, 6))
  }
  
  cat("\nDelta vector (subset):\n")
  print(utils::head(round(object$delta, 6)))
  
  cat("\nSigma_delta (shape): ", paste(dim(object$Sigma_delta), collapse = " x "), "\n", sep = "")
  invisible(object)
}