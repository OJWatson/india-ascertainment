# additional_functions for rolling
roll_func <- function(x, det) {
  l <- length(det)
  ret <- rep(0, length(x))
  for(i in seq_along(ret)) {
    to_sum <- tail(x[seq_len(i)], length(det))
    ret[i] <- sum(rev(to_sum)*head(det, length(to_sum)),na.rm=TRUE)
  }
  return(ret)
}

#' Create a sero detection probability density function
#' @param ret What antibody are we detecting. Must be one of "igg", "igm" or
#'   "iggm". Default = "igg
#' @param igg_conv Mean days for IgG seroconversion. Default = 13.3
#' @param igm_conv Mean days for IgM seroconversion. Default = 12.3
#' @param igg_sens Max sensitivity for IgG. Default = 0.914
#' @param igm_sens Max sensitivity for IgM. Default = 0.894
#' @param igg_scale Weibull seroreversion decay scale for IgG. Default = 139
#' @param igm_scale Weibull seroreversion decay scale for IgM. Default = 50
#' @param pdfl Length of probability density function. Default = 600 days
#'
#' @importFrom stats dexp dweibull pexp pweibull
sero_det <- function(ret = "igg",
                     igg_conv = 13.3, igm_conv = 12.3,
                     igg_sens = 0.914, igm_sens = 0.894,
                     igg_scale = 139, igm_scale = 50,
                     pdfl = 600) {

  # create igg curve if needed
  if(ret %in% c("igg", "iggm")) {

    # seroconversion and scale by max sens
    sero_conv_igg <- pexp(seq_len(pdfl), rate = 1/igg_conv)
    sero_conv_igg <- (sero_conv_igg/max(sero_conv_igg)) * igg_sens

    # seroreversion and check for below 0s due to insufficiently long pdf
    sero_det_igg <- sero_conv_igg-pweibull(seq_len(pdfl), 3.7, scale = igg_scale)
    sero_det_igg[sero_det_igg < 0] <- 0

  }

  # create igm curve if needed
  if(ret %in% c("igm", "iggm")) {

    # seroconversion and scale by max sens
    sero_conv_igm <- pexp(seq_len(pdfl), rate = 1/igm_conv)
    sero_conv_igm <- (sero_conv_igg/max(sero_conv_igm)) * igm_sens

    # seroreversion and check for below 0s due to insufficiently long pdf
    sero_det_igm <- sero_conv_igm-pweibull(seq_len(pdfl), 3.7, scale = igm_scale)
    sero_det_igm[sero_det_igm < 0] <- 0

  }

  # create iggm curve if needed
  if(ret == "iggm") {

    sero_det_iggm <- 1-((1-sero_det_igg) * (1-sero_det_igm))

  }

  # return requested sero_det curve
  if(ret == "igg"){
    return(sero_det_igg)
  } else if (ret == "igm") {
    return(sero_det_igm)
  } else if (ret == "iggm") {
    return(sero_det_iggm)
  } else {
    stop("Incorrect ret value. Must be one of igg, igm or iggm")
  }

}

#' Create a seropositivity data frame from squire run output
#' @param res squire model run
#' @param sero_det sero detection curve from [[sero_det]]
#'
#' @importFrom dplyr group_by mutate select left_join ungroup
seroprev_df_det <- function(res, sero_det) {

  # get symptom onset data
  date_0 <- max(res$pmcmc_results$inputs$data$date)
  inf <- squire::format_output(res, c("S"), date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
    mutate(S = as.integer(.data$y)) %>%
    group_by(replicate) %>%
    mutate(infections = c(0, diff(max(.data$S,na.rm=TRUE)-.data$S))) %>%
    select(replicate, t, date, .data$S, .data$infections)

  # correctly format
  inf <- left_join(inf,
                   squire::format_output(
                     res, c("infections"),
                     date_0 = max(res$pmcmc_results$inputs$data$date)
                   ) %>%
                     mutate(symptoms = as.integer(.data$y)) %>%
                     select(replicate, t, .data$date, .data$symptoms),
                   by = c("replicate", "t", "date"))

  inf <- inf %>%
    group_by(replicate) %>%
    mutate(sero_positive = roll_func(.data$symptoms, sero_det),
           sero_perc = .data$sero_positive/max(.data$S,na.rm = TRUE)) %>%
    ungroup

  inf$reporting_fraction <- res$pmcmc_results$inputs$pars_obs$phi_death
  return(inf)

}

#' Format output of seroprev output from
#' @param out Seroprev output from [[seroprev_df_det]]
#'
#' @importFrom stats quantile
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr across
format_sero_df <- function(out) {

  out %>% select(replicate, date, sero_perc) %>%
    group_by(date) %>%
    summarise(across(.data$sero_perc, .fns = list(
      min = ~quantile(.x, 0.025, na.rm=TRUE),
      med = ~quantile(.x, 0.5, na.rm=TRUE),
      max = ~quantile(.x, 0.975, na.rm=TRUE)
    ))) %>%
    pivot_longer(-date, names_to = c("stat"), names_pattern = ".*_(.*)") %>%
    pivot_wider(names_from = .data$stat, values_from = .data$value)

}

#' Plot seroprev curve
#' @param out Formatted seroprev output from [[format_sero_df]]
#' @param color Name of column that determines color of curves. Default = NULL
#' @param ci Logical for whether CI is plotted. Default = FALSE
#'
#' @importFrom ggplot2 aes geom_vline geom_line theme_bw xlab theme ylab scale_y_continuous
#' @importFrom ggplot2 scale_color_discrete theme ggplot geom_ribbon
plot_sero_df <- function(out, color = NULL, ci = TRUE) {

  if (is.null(color)) {
    out$color <- 1
  } else {
    out$color <- out[[color]]
  }

  gg <- out %>%
    ggplot(aes(x = date, y = med, ymin = min, ymax = max, color = color)) +
    geom_vline(xintercept = as.Date("2020-08-20")) +
    geom_line() +
    theme_bw() +
    xlab("") +
    ylab("Seroprevalence") + scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "top")

  if (is.null(color)) {
    gg <- gg + theme(legend.position = "none")
  } else {
    gg <- gg + scale_color_discrete(name = color)
  }

  if (ci) {
    gg <- gg + geom_ribbon(alpha = 0.2)
  }

  gg
}


rt_df <- function(out) {

# create the Rt data frame
rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {

  tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change,
                                             change = out$interventions$R0_change,
                                             start_date = out$replicate_parameters$start_date[y],
                                             steps_per_day = 1/out$parameters$dt)

    Rt <- squire:::evaluate_Rt_pmcmc(
      R0_change = tt$change,
      date_R0_change = tt$dates,
      R0 = out$replicate_parameters$R0[y],
      pars = as.list(out$replicate_parameters[y,]),
      Rt_args = out$pmcmc_results$inputs$Rt_args)

  Rt_full <- approx(x = tt$tt, Rt, method = "constant", xout = seq(tt$tt[1], tt$tt[length(tt$tt)]))$y

  df <- data.frame(
    "Rt" = Rt_full,
    "date" = seq.Date(tt$dates[1], tt$dates[length(tt$dates)], 1),
    rep = y,
    stringsAsFactors = FALSE)
  df$pos <- seq_len(nrow(df))
  return(df)
} )

rt <- do.call(rbind, rts)
rt$date <- as.Date(rt$date)

rt <- rt[,c(3,2,1)]

new_rt_all <- rt %>%
  group_by(rep) %>%
  arrange(date) %>%
  complete(date = seq.Date(min(rt$date), date_0, by = "days"))

column_names <- colnames(new_rt_all)[-c(1,2)]
new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))

suppressMessages(sum_rt <- group_by(new_rt_all, date) %>%
                   summarise(Rt_min = quantile(Rt, 0.025),
                             Rt_q25 = quantile(Rt, 0.25),
                             Rt_q75 = quantile(Rt, 0.75),
                             Rt_max = quantile(Rt, 0.975),
                             Rt_median = median(Rt),
                             Rt = mean(Rt)))

return(sum_rt)

}
