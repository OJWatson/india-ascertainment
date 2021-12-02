
gen_summary <- function(res, date_0 = as.Date("2020-02-05"), rf = 0.2, rf2 = 0.1) {
  
  # data
  owid <- read.csv("https://covid.ourworldindata.org/data/owid-covid-data.csv")
  
  # serosurveys
  # May–June, 2020 0·73% (95% CI 0·34–1·13) #  https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30544-1/fulltext
  # Aug 18 and Sept 20, 2020, 7·1% (6·2–8·2) # https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30544-1/fulltext
  # December 17 2020 and January 8 2021 21.5
  sero_df_real <- data.frame("date" = as.Date(c("2020-05-15", "2020-09-01", "2020-12-25")),
                        "prev" = c(0.73*0.01,7.1*0.01, 21.5*0.01),
                        "prev_low" = c(0.34*0.01, 6.2*0.01, 21.5*0.01),
                        "prev_high" = c(1.13*0.01, 8.2*0.01, 21.5*0.01))
  
  out <- squire::format_output(res, c("deaths", "S"), date_0 = date_0)
  out$y[out$compartment == "S"] <- (sum(squire::get_population("India")$n) - out$y[out$compartment == "S"])/(sum(squire::get_population("India")$n))
  
  g1 <- ggplot(out %>% filter(compartment == "deaths"), aes(date, y)) +
    geom_point(aes(as.Date(date), new_deaths), data = owid %>% filter(iso_code == "IND"), inherit.aes = TRUE, alpha = 0.5) +
    #geom_line() +
    geom_line(aes(date, 
                  y*c(rep(rf, as.Date("2021-04-01") - date_0), rep(rf2, nrow(out[out$compartment=="deaths",]) - as.integer(as.Date("2021-04-01") - date_0)))), 
                  color = "red") + ylab("Deaths") + xlab("") +
    geom_vline(xintercept = as.Date("2021-05-8"), linetype = "dashed") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    theme_bw() +
    ggtitle(paste0("Deaths = True Deaths, Red = Reported Deaths at ", rf*100,"% Flat Reporting"))
  
  # quick overview
  ct_df <- readRDS(file.path(here::here(), "analysis/khartoum/ct_df.rds"))
  pcr_days <- length(ct_df$PropDet)
  pcr_det <- ct_df$PropDet
  
  sero_df <- readRDS(file.path(here::here(), "analysis/khartoum/weibull_params.RDS"))
  prob_conversion <-  cumsum(dgamma(0:300,shape = 5, rate = 1/2))/max(cumsum(dgamma(0:300,shape = 5, rate = 1/2)))*0.95
  sero_det <- cumsum(dweibull(0:300, sero_df$wshape, scale = sero_df$wscale))
  sero_det <- cumsum(prob_conversion-sero_det)
  sero_det[sero_det < 0] <- 0
  sero_det <- sero_det/max(sero_det)
  roll_func <- function(x, det) {
    l <- length(det)
    c(NA, 
      zoo::rollapply(x, 
                     list(seq(-l, -1)),
                     function(i) {
                       sum(i*tail(det, length(i)), na.rm = TRUE)
                     },
                     partial = 1
      ))
  }
  
  inf <- squire::format_output(res, c("S"), date_0 = date_0) %>% 
    mutate(S = as.integer(y)) %>% 
    group_by(replicate) %>%  
    mutate(infections = lag(S, 1)-S) %>% 
    select(replicate, t, date, S, infections)
  
  inf <- left_join(inf, 
                   squire::format_output(res, c("infections"), date_0 = date_0) %>% 
                     mutate(symptoms = as.integer(y)) %>% 
                     select(replicate, t, date, symptoms), 
                   by = c("replicate", "t", "date"))
  
  inf <- inf %>% group_by(replicate) %>%
    mutate(pcr_positive = roll_func(infections, pcr_det),
           sero_positive = roll_func(symptoms, sero_det),
           ps_ratio = pcr_positive/sero_positive, 
           sero_perc = sero_positive/max(S,na.rm = TRUE),
           pcr_perc = pcr_positive/max(S,na.rm = TRUE))
  
  sero_df <- data.frame("date" = as.Date(c("2020-05-15", "2020-09-01", "2020-12-25")),
                        "prev" = c(0.73*0.01,7.1*0.01, 21.5*0.01))
  
  g2 <- ggplot(out %>% filter(compartment == "S"), aes(date+21, y)) + geom_line() +
    geom_vline(xintercept = as.Date("2021-05-18"), linetype = "dashed") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ylab("Attack Rate") + xlab("") + theme_bw()  +
    ggtitle("Total Attack Rate Under Assumption of Negligable Waning Immunity")
  
  g3 <- ggplot(inf, aes(date,sero_perc)) + geom_line() + 
    geom_vline(xintercept = as.Date("2021-05-18"), linetype = "dashed") +
    geom_point(aes(as.Date(date), prev), sero_df_real, inherit.aes = FALSE) +
    geom_errorbar(aes(as.Date(date), ymin = prev_low, ymax = prev_high), sero_df_real, inherit.aes = FALSE, width = 0) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw()+ ylab("Seroprevalence") + xlab("") +
    ggtitle("Points Represent Mean of Nationally Representative Seroprevalence Surveys")
    
  
  cowplot::plot_grid( 
    g1, g2, g3, ncol = 1
    )
  
}

r <- squire::run_deterministic_SEIR_model("India", R0 = c(2.65, 1.15, 2.8), tt_R0 = c(0, 70, 410),
                                          prob_severe = squire:::default_probs()$prob_severe*1.2,
                                          day_return = TRUE, time_period = 680, ICU_bed_capacity = 1e10)

gen <- gen_summary(r, date_0 = as.Date("2020-01-21"), rf = 0.2, rf2 = 0.1)
gen 

