RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

system(paste0("echo India Subnational for  ",state))

version_min <- "0.6.7"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least v", version_min)
}

version_min <- "0.1.22"
if(packageVersion("nimue") < version_min) {
  stop("nimue needs to be updated to at least v", version_min)
}

## -----------------------------------------------------------------------------
## 0. Checks and Function Definitions
## -----------------------------------------------------------------------------

# check on the state name
if (!(state %in% c(
  "Andaman and Nicobar Islands","Andhra Pradesh","Arunachal Pradesh","Assam",
  "Bihar","Chandigarh","Chhattisgarh","Dadra and Nagar Haveli and Daman and Diu",
  "Delhi","Goa","Gujarat","Haryana","Himachal Pradesh","Jammu and Kashmir",
  "Jharkhand","Karnataka","Kerala","Ladakh","Lakshadweep","Madhya Pradesh",
  "Maharashtra","Manipur","Meghalaya","Mizoram","Nagaland","Odisha","Puducherry",
  "Punjab","Rajasthan","Sikkim","Tamil Nadu","Telangana","Tripura","Uttar Pradesh",
  "Uttarakhand","West Bengal"))) {
  stop("State is not correct")
}

## -----------------------------------------------------------------------------
## 1. GET INPUT DATA
## -----------------------------------------------------------------------------

## a. Get from local files
## -----------------------------------------------------------------------------

# get pop data from file
demog <- readRDS("demog.rds")
pop <- demog$n[demog$state == state]

# get icu beds from file
icu_beds <- readRDS("icu_beds.rds")
icu_beds <- icu_beds$icu_beds[icu_beds$state == state]

# get hosp beds from file
hosp_beds <- readRDS("hosp_beds.rds")
hosp_beds <- hosp_beds$hosp_beds[hosp_beds$state == state]

# get seroprevalence data points
sero_df <- readRDS("sero.rds")
sero_df <- sero_df[sero_df$state == state,]

# minimum rf allowed
min_rf_df <- readRDS("min_rfs.rds")
min_rf <- max(min_rf_df$min_rf[min_rf_df$state == state], 0.01)
if(nrow(sero_df) > 0) {
max_rf <- min(min_rf_df$max_rf[min_rf_df$state == state]*1.5, 1.00)
} else {
max_rf <- min(min_rf_df$max_rf[min_rf_df$state == state], 1.00)
}


# seroconversion data from brazeau report 34 addjusted in light of more longer term studies
prob_conversion <-  cumsum(dgamma(0:600,shape = 5, rate = 1/2))/max(cumsum(dgamma(0:300,shape = 5, rate = 1/2)))
sero_det <- cumsum(dweibull(0:600, 3.669807, scale = 210.7046))
sero_det <- prob_conversion-sero_det
sero_det[sero_det < 0] <- 0
sero_det <- sero_det/max(sero_det)

## b. Get from remote changing sources
## -----------------------------------------------------------------------------

# get death data
subnat_df <- read.csv("https://data.incovid19.org/csv/latest/states.csv") %>%
  filter(!(State %in% c("India", "State Unassigned"))) %>%
  mutate(Date = as.Date(Date)) %>%
  group_by(State) %>%
  complete(Date = seq.Date(min(as.Date(Date)), max(as.Date(Date)), 1)) %>%
  mutate(State = replace_na(State, unique(!is.na(State)))) %>%
  mutate(cases = replace_na(Confirmed, 0),
         deaths = replace_na(Deceased, 0)) %>%
  select(-c("Tested", "Recovered", "Other", "Confirmed", "Deceased")) %>%
  rename(date = Date, state = State)
df <- subnat_df[subnat_df$state == state, ] %>%
  ungroup %>%
  select(date, deaths, cases) %>%
  arrange(date)

if(state == "Maharashtra") {

  df$deaths[df$date == "2021-11-02"] <- 140274
  df$cases[df$date == "2021-11-02"] <- 6612965
  df <- df[-603,]

}

df$deaths <- c(df$deaths[1], diff(df$deaths))
df$deaths[df$deaths < 0] <- 0
df$cases <- c(df$cases[1], diff(df$cases))
df$cases[df$cases < 0] <- 0

# filter to the current date
df <- df[as.Date(df$date) <= as.Date(date),]

# get vaccination data
subnat_vacc <- read.csv("https://data.incovid19.org/csv/latest/cowin_vaccine_data_statewise.csv")
subnat_vacc <- subnat_vacc %>% rename(region = State, date = Updated.On) %>%
  mutate(date = as.Date(date, "%d/%m/%Y")) %>%
  group_by(region, date) %>%
  summarise(total_vaccinations = sum(Total.Doses.Administered),
            people_vaccinated = sum(First.Dose.Administered),
            people_fully_vaccinated = sum(Second.Dose.Administered)) %>%
  rename(state = region)
subnat_vacc <- subnat_vacc[subnat_vacc$state == state, ]
vacc_inputs <- get_vaccine_inputs(max(df$date), subnat_vacc)

## c. Sort out data issues
## -----------------------------------------------------------------------------

redist_deaths <- function(df, dc, past_days = NULL) {

  dd <- df$deaths[df$date == dc]

  df$deaths[df$date == dc] <- NA
  dd_wk <- df$deaths[df$date %in% seq.Date(dc-3, dc+3, 1)]
  df$deaths[df$date == dc] <- as.integer(mean(dd_wk, na.rm = TRUE))
  dd_new <- df$deaths[df$date == dc]
  to_distribute <- dd - dd_new

  to_add_pos <- which(df$date < dc)
  if(!is.null(past_days)) {
    to_add_pos <- tail(to_add_pos, past_days)
  }
  orig <- df$deaths[to_add_pos]
  to_add <- floor(orig/sum(orig) * to_distribute)
  remaining <- to_distribute - sum(to_add)
  top_up <- which(order(orig)>(length(orig)-remaining))
  to_add[top_up] <- to_add[top_up] + 1

  df$deaths[to_add_pos] <- df$deaths[to_add_pos] + to_add

  return(df)

}

# large death spike fixes for multiple regions
if (state == "Maharashtra") {
  df <- redist_deaths(df, as.Date("2020-06-16"))
  df <- redist_deaths(df, as.Date("2021-07-20"))
}
if (state == "Delhi") {
  df <- redist_deaths(df, as.Date("2020-06-16"))
}
if (state == "Chhattisgarh") {
  df <- redist_deaths(df, as.Date("2020-09-09"))
}
if (state == "Bihar") {
  df <- redist_deaths(df, as.Date("2021-06-09"), 60)
}
if (state == "Madhya Pradesh") {
  df <- redist_deaths(df, as.Date("2021-07-12"), 90)
}
if (state == "Tamil Nadu") {
  df <- redist_deaths(df, as.Date("2020-07-22"))
}
if (state == "West Bengal") {
  df <- redist_deaths(df, as.Date("2020-05-03"))
}
if (state == "Uttarakhand") {
  df <- redist_deaths(df, as.Date("2020-10-17"))
}
if (state == "Goa") {
  df <- redist_deaths(df, as.Date("2021-09-15"), 60)
}
if (state == "Haryana") {
  df <- redist_deaths(df, as.Date("2021-09-13"), 60)
  df <- redist_deaths(df, as.Date("2021-09-28"), 60)
  df <- redist_deaths(df, as.Date("2021-10-14"), 60)
}
if (state == "Nagaland") {
  df <- redist_deaths(df, as.Date("2021-08-29"), 60)
}
## -----------------------------------------------------------------------------
## 2. Fit Model
## -----------------------------------------------------------------------------

# fit model
res <- fit_spline_rt(data = df,
                     country = as.character("India"),
                     pop = pop,
                     min_rf = as.numeric(min_rf),
                     max_rf = as.numeric(max_rf),
                     vacc_inputs = vacc_inputs,
                     n_mcmc = as.numeric(n_mcmc),
                     replicates = as.numeric(replicates),
                     hosp_beds = as.numeric(hosp_beds),
                     icu_beds = as.numeric(icu_beds),
                     sero_df = sero_df,
                     sero_det = sero_det,
                     model = model,
                     pars_obs_dur_R = as.numeric(dur_R),
                     pars_obs_prob_hosp_multiplier = as.numeric(prob_hosp_multiplier),
                     )


# add state for ease and remove the output for memory
res$parameters$state <- state
output <- res$output
res$output <- NULL

# save output without output for memory
saveRDS(res, "res.rds")
res$output <- output

# make a quick plot so we can check fits easily
if (model == "SQUIRE") {
rtp <- rt_plot_immunity(res)
} else {
  rtp <- rt_plot_immunity(res)
}
dp <- dp_plot(res)
cdp <- cdp_plot(res)
sero <- sero_plot(res, sero_df)
ar <- ar_plot(res)

rf_over <- paste0(round(quantile(res$replicate_parameters$rf)[c(2,4)], digits = 2)*100, "%", collapse = " - ")
ggsave("fitting.pdf",width=12, height=15,
       cowplot::plot_grid(
         rtp$plot +
           ggtitle(paste0(state, ". Death Reporting at ", rf_over)) +
           scale_x_date(date_labels = "%b %Y", date_breaks = "3 months"),
         dp, cdp, sero, ar, ncol = 1, align = "v"))

projs <- get_projections(res)
saveRDS(projs, "proj.rds")

