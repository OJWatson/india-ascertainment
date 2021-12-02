
get_best <- function(out) {

  if("pmcmc_results" %in% names(out)) {
    if("chains" %in% names(out$pmcmc_results)) {
      mc <- do.call(rbind, lapply(out$pmcmc_results$chains, "[[", "results"))
    } else {
      mc <- out$pmcmc_results$results
    }
    best <- mc[which.max(mc$log_posterior),]
    best <- best[,seq_len(ncol(best)-3)]
    rownames(best) <- NULL
    best$start_date <- as.character(squire:::offset_to_start_date(out$pmcmc_results$inputs$data$date[1], round(best$start_date)))
    best$date_Meff_change <- out$pmcmc_results$inputs$Rt_args$date_Meff_change
    best$Rt_shift_duration <- out$pmcmc_results$inputs$Rt_args$Rt_shift_duration
    best$Rt_rw_duration <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration

    if("chains" %in% names(out$pmcmc_results)) {
      best$covariance_matrix <- out$pmcmc_results$chains$chain1$covariance_matrix[1]
      best$scaling_factor <- mean(
        c(tail(na.omit(out$pmcmc_results$chains$chain1$scaling_factor),1),
          tail(na.omit(out$pmcmc_results$chains$chain2$scaling_factor),1),
          tail(na.omit(out$pmcmc_results$chains$chain3$scaling_factor),1))
      )
    } else {
      best$covariance_matrix <- out$pmcmc_results$covariance_matrix[1]
      best$scaling_factor <- tail(na.omit(out$pmcmc_results$scaling_factor),1)
    }

    best$state <- out$parameters$state

    return(best)
  }
}

stem <- "/home/oj/net/lmic_new/datadrive/lmic/global-lmic-reports-orderly/archive/india_sub_national/"
# if specific wants
reports_all <- function(task = "india_sub_national") {

  wd <- "/home/oj/net/lmic_new/datadrive/lmic/global-lmic-reports-orderly/"
  wdold <- getwd()
  setwd(wd)
  db <- orderly::orderly_db("destination")
  if (is.null(date)) {
    date <- as.character(Sys.Date())
  }

  # get parameters
  pars <- DBI::dbGetQuery(db, 'select * from parameters')
  reports <- DBI::dbGetQuery(db, 'select * from report_version')
  reports <- reports %>% filter(report == task)
  pars <- pars %>% filter(report_version %in% reports$id)
  pars <- pars %>% select(-c("type","id")) %>%
    pivot_wider(names_from = "name", values_from = "value")

  return(pars)
}

pars_init <- readRDS(file.path(here::here(), "src/india_sub_national/pars_init.rds"))

# central pars
reports <- reports_all()
reports <- reports %>% filter(dur_R == 115) %>% group_by(state) %>% filter(report_version == max(report_version[n_mcmc >= 5000]))
lf <- file.path(stem, reports$report_version, "res.rds")
pars_saved <- lapply(lf, function(x){get_best(readRDS(x))})
names(pars_saved) <- unlist(lapply(pars_saved, "[[", "state"))
pars_init$central[match(names(pars_saved), names(pars_init$central))] <- pars_saved

# worst pars
reports <- reports_all()
reports_worst <- reports %>% filter(dur_R == 80) %>% group_by(state) %>% filter(report_version == max(report_version[n_mcmc >= 5000]))
lf <- file.path(stem, reports$report_version, "res.rds")
pars_saved <- lapply(lf, function(x){get_best(readRDS(x))})
names(pars_saved) <- unlist(lapply(pars_saved, "[[", "state"))
pars_init$worst[match(names(pars_saved), names(pars_init$worst))] <- pars_saved


saveRDS(pars_init, file.path(here::here(), "src/india_sub_national/pars_init.rds"))
saveRDS(pars_init, "/home/oj/net/lmic_new/datadrive/lmic/global-lmic-reports-orderly/src/india_sub_national//pars_init.rds")

