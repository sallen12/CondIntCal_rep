
# wrapper to get prediction intervals from the 6 forecast distributions
get_ints <- function(mu, tau, alpha) {
  n <- length(mu)

  # Ideal forecaster: F = N(mu, 1)
  L_id <- qnorm(alpha/2, mu)
  U_id <- qnorm(1 - alpha/2, mu)
  int_id <<- cbind(L_id, U_id)

  # Unconditional forecaster: F = N(0, 2)
  L_cl <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
  U_cl <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)
  int_cl <<- cbind(L_cl, U_cl)

  # Unfocused forecaster: F = 0.5 N(mu, 1) + 0.5 N(mu + tau, 1), tau = -1/+1 w.e.p
  F_un <- function(q, m, t) 0.5*pnorm(q, m) + 0.5*pnorm(q, m + t)
  Q_un <- function(p, m, t) uniroot(function(x) F_un(x, m, t) - p, interval = c(-10, 10))$root
  L_un <- sapply(1:n, function(i) Q_un(alpha/2, mu[i], tau[i]))
  U_un <- sapply(1:n, function(i) Q_un(1 - alpha/2, mu[i], tau[i]))
  int_un <<- cbind(L_un, U_un)

  # mean biased forecaster: F = N(mu + tau, 1), tau = -1/+1 w.e.p
  F_me <- function(q, m, t) pnorm(q, m + t)
  L_me <- qnorm(alpha/2, mu + tau)
  U_me <- qnorm(1 - alpha/2, mu + tau)
  int_me <<- cbind(L_me, U_me)

  # sign biased forecaster: F = N(-mu, 1)
  L_si <- qnorm(alpha/2, mean = -mu, sd = 1)
  U_si <- qnorm(1 - alpha/2, mean = -mu, sd = 1)
  int_si <<- cbind(L_si, U_si)

  # mixed (unconditional/signed) forecaster: F = N(0, 2)/N(-mu,1) w.e.p
  L_1 <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
  U_1 <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)
  L_2 <- qnorm(alpha/2, mean = -mu, sd = 1)
  U_2 <- qnorm(1 - alpha/2, mean = -mu, sd = 1)
  r <- sample(c(0, 1), n, replace = T)
  L_mi <- r*L_1 + (1 - r)*L_2
  U_mi <- r*U_1 + (1 - r)*U_2
  int_mi <<- cbind(L_mi, U_mi)
}

# function to get marginal coverage
coverage <- function(y, int, closed = T) {
  n <- length(y)
  if (n == 1) {
    if (!is.vector(int)) {
      y <- rep(y, nrow(int))
    } else {
      int <- int |> as.matrix() |> t()
    }
  } else {
    if (is.vector(int)) {
      int <- matrix(int, nrow = n, ncol = 2, byrow = T)
    }
  }
  if (closed) {
    mean(y >= int[, 1] & y <= int[, 2])
  } else {
    mean(y > int[, 1] & y < int[, 2])
  }
}

# function to get average length
ilength <- function(int) {
  if (is.vector(int)) {
    int[2] - int[1]
  } else {
    mean(int[, 2] - int[, 1])
  }
}

# function to perform isotonic interval score decomposition
is_decomp <- function(y, int, alpha, return_fit = T) {
  int_df <- as.data.frame(int)
  fit <- idr(y = y, X = int_df) # fit IDR
  out <- predict(fit, int_df) # get predicted distributions

  int_rc <- cbind(qpred(out, alpha/2), qpred(out, 1 - alpha/2)) # get recalibrated interval forecasts
  int_mg <- c(quantile(y, alpha/2, type = 1), quantile(y, 1 - alpha/2, type = 1)) # get unconditional interval forecasts

  IS <- ints_quantiles(y, int[, 1], int[, 2], 1 - alpha) |> mean() # interval score of original forecast
  IS_mg <- ints_quantiles(y, int_mg[1], int_mg[2], 1 - alpha) |> mean() # interval score of unconditional forecast
  IS_rc <- ints_quantiles(y, int_rc[,1], int_rc[, 2], 1 - alpha) |> mean() # interval score of recalibrated forecast

  MCB <- IS - IS_rc # miscalibration
  DSC <- IS_mg - IS_rc # discrimination

  if (return_fit) {
    return(list(decomp = c("IS" = IS, "UNC" = IS_mg, "DSC" = DSC, "MCB" = MCB), int_rc = int_rc))
  } else {
    return(c("IS" = IS, "UNC" = IS_mg, "DSC" = DSC, "MCB" = MCB))
  }
}

# function to create MCB-DSC plot (adapted from triptych package)
mcbdsc <- function(df_e, n_isolines = 10, colour_values = "black", colour_unc = "#00BF7D", MCBDSC_repel = FALSE, MCB_lim = NA, DSC_lim = NA) {

  # Limits and out-of-bound identification
  default_lims <- function(x) c(0, 1.1 * max(x[is.finite(x)]))
  is_in_range <- function(x, xrange) x >= xrange[1] & x <= xrange[2]
  split_by_out_of_bounds <- function(df, MCB_lim, DSC_lim) {
    oob_state <- factor(
      x = with(df, dplyr::case_when(
        is_in_range(DSC, DSC_lim) & is_in_range(MCB, MCB_lim) ~ "within",
        is_in_range(DSC, DSC_lim) & !is.finite(MCB)           ~ "infty",
        TRUE                                                  ~ "oob"
      )),
      levels = c("within", "infty", "oob")
    )
    res <- split(df, oob_state)
    if (nrow(res$infty)) {
      res$infty$x_geom_text <- MCB_lim[2]
    }
    res
  }
  if (anyNA(MCB_lim)) {
    MCB_lim <- default_lims(df_e$MCB)
  }
  if (anyNA(DSC_lim)) {
    DSC_lim <- default_lims(df_e$DSC)
  }
  df_e_by_state <- split_by_out_of_bounds(df_e, MCB_lim, DSC_lim)
  # Check that the plot is not empty of points!
  if (!nrow(df_e_by_state$within)) {
    warning(paste(
      "The given limits for the MCB-DSC plot exclude all forecasts.",
      "The default choices are used instead."
    ))
    MCB_lim <- default_lims(df_e$MCB)
    DSC_lim <- default_lims(df_e$DSC)
    df_e_by_state <- split_by_out_of_bounds(df_e, MCB_lim, DSC_lim)
  }
  if (nrow(df_e_by_state$oob)) {
    message(paste(
      "The following forecasts are not included in the MCB-DSC plot as their",
      "miscalibration measure is outside the plot limits:",
      paste(df_e_by_state$oob$forecast, collapse = ", ")
    ))
  }

  # Reasonable score values for isolines
  choose_isolines <- function(unc, MCB_lim, DSC_lim) {
    scores <- pretty(x = unc - c(-1.1 * MCB_lim[2], DSC_lim[2]), n = n_isolines)
    tibble::tibble(
      slope = 1,
      intercept = unc - scores,
      label = scores
    ) |>
      # Remove a line if its intercept is too close to zero (less than 1/5 times the line distance).
      dplyr::filter(abs(.data$intercept) > abs(diff(.data$intercept)[1]) / 5)
  }
  df_iso_abline <- choose_isolines(df_e$UNC[1], MCB_lim, DSC_lim)

  # replicate colour_values if it has length 1
  if (length(colour_values) == 1L & nrow(df_e) > 1L) {
    colour_values <- rep(colour_values, nrow(df_e))
  }

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(x = 0, y = 0, xend = .data$max_val, yend = .data$max_val),
      data = tibble::tibble(max_val = 2 * max(MCB_lim, DSC_lim)),
      colour = colour_unc,
      linewidth = 1
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(.data$x, .data$y),
      data = data.frame(x = 0, y = 0),
      colour = colour_unc,
      fill = colour_unc,
      size = 2,
      shape = 15
    ) +
    geomtextpath::geom_labelabline(
      mapping = ggplot2::aes(
        intercept = 0,
        slope = 1,
        label = paste("UNC:", prettyNum(.data$UNC, digits = 3))
      ),
      data = df_e[1, ],
      colour = colour_unc,
      hjust = 0.85,
      size = 7 * 0.36,
      text_only = TRUE,
      boxcolour = NA,
      straight = TRUE
    ) +
    ggplot2::geom_abline(
      data = df_iso_abline,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope),
      colour = "gray"
    ) +
    geomtextpath::geom_labelabline(
      data = df_iso_abline,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope, label = .data$label),
      colour = "gray",
      hjust = 0.85,
      size = 7 * 0.36,
      text_only = TRUE,
      boxcolour = NA,
      straight = TRUE
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, colour = .data$forecast),
      data = df_e_by_state$within
    ) +
    {
      if (!isTRUE(MCBDSC_repel)) {
        ggplot2::geom_text(
          data = df_e_by_state$within,
          mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
          size = 3,
          vjust = 0,
          hjust = 0,
          check_overlap = TRUE,
          position = ggplot2::position_nudge(
            x = diff(MCB_lim) / 80,
            y = -diff(DSC_lim) / 40
          )
        )
      } else if (isTRUE(MCBDSC_repel)) {
        ggrepel::geom_text_repel(
          data = df_e_by_state$within,
          mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
          size = 3
        )
      }
    } +
    {
      if (nrow(df_e_by_state$infty)) {
        list(
          ggplot2::geom_rug(
            data = df_e_by_state$infty,
            mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, colour = .data$forecast),
            sides = "r",
            linewidth = 2
          ),
          ggplot2::geom_text(
            data = df_e_by_state$infty,
            mapping = ggplot2::aes(x = .data$x_geom_text, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
            size = 3,
            hjust = 1,
            check_overlap = TRUE
          )
        )
      }
    } +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::scale_x_continuous(oob = scales::oob_squish_infinite) +
    ggplot2::coord_cartesian(xlim = MCB_lim, ylim = DSC_lim) +
    ggplot2::labs(x = expression(MCB[I]), y = expression(DSC[I])) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      )
    )
}


