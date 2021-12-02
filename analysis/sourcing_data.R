## --------------------------------
# fetch and format/filter data
## --------------------------------

subnat_vacc <- read.csv("https://raw.githubusercontent.com/sociepy/covid19-vaccination-subnational/main/data/countries/India.csv")
subnat_vacc <- subnat_vacc %>%
  mutate(region = replace(region, which(region_iso %in% c("IN-DN", "IN-DD")), "Dadra and Nagar Haveli and Daman and Diu"),
         region_iso = replace(region_iso, which(region_iso %in% c("IN-DN", "IN-DD")), "IN-DD")) %>%
  group_by(region, date, region_iso) %>%
  summarise(total_vaccinations = sum(total_vaccinations),
            people_vaccinated = sum(people_vaccinated),
            people_fully_vaccinated = sum(people_fully_vaccinated)) %>%
  rename(state = region)

subnat_pop <- read.csv("https://raw.githubusercontent.com/sociepy/covid19-vaccination-subnational/main/data/population.csv")
subnat_pop <- subnat_pop %>%
  filter(grepl("^IN", region_iso)) %>%
  mutate(region_iso = replace(region_iso, which(region_iso %in% c("IN-DN", "IN-DD")), "IN-DD")) %>%
  group_by(region_iso) %>%
  summarise(population = max(population)) %>%
  left_join(unique(subnat_vacc[,c("region_iso", "state")]))

subnat_df <- read.csv("https://api.covid19india.org/csv/latest/states.csv") %>%
filter(!(State %in% c("India", "State Unassigned"))) %>%
  mutate(Date = as.Date(Date)) %>%
  group_by(State) %>%
  complete(Date = seq.Date(min(as.Date(Date)), max(as.Date(Date)), 1)) %>%
  mutate(State = replace_na(State, unique(!is.na(State)))) %>%
  mutate(cases = replace_na(Confirmed, 0),
         deaths = replace_na(Deceased, 0)) %>%
  select(-c("Tested", "Recovered", "Other", "Confirmed", "Deceased")) %>%
  rename(date = Date, state = State)

## --------------------------------
# hospital bed source
# https://www.medrxiv.org/content/10.1101/2020.06.16.20132787v1.full-text
# ------------------------------------------------------------------------------

# -> icu_beds.csv
# -> hosp_beds.csv
icu_beds <- read.csv("analysis/india/icu_beds.csv")
icu_beds <- icu_beds %>%
  mutate(State = gsub("_", " ", State)) %>%
  mutate(State = replace(State, which(State %in% c("Dadra & N Haveli", "Daman & Diu")), "Dadra and Nagar Haveli and Daman and Diu")) %>%
  mutate(State = replace(State, which(State %in% c("Andaman & N. Islands")), "Andaman and Nicobar Islands")) %>%
  mutate(State = gsub("&", "and", State)) %>%
  group_by(State) %>%
  summarise(icu_beds = sum(Total_ICU_Beds_Est)) %>%
  rename(state = State)

# and fill in Ladakh using model based estimate for nation
icu_beds$icu_beds[icu_beds$state == "Ladakh"] <-
  round(squire:::get_healthcare_capacity("India")$ICU_beds*subnat_pop$population[subnat_pop$state=="Ladakh"]/1000)
saveRDS(icu_beds, cp_path("src/india_sub_national/icu_beds.rds"))
saveRDS(icu_beds, cp_path("analysis/data/derived/icu_beds.rds"))

hosp_beds <- read.csv("analysis/india/hosp_beds.csv")
hosp_beds <- hosp_beds %>%
  mutate(State = gsub("_", " ", State)) %>%
  mutate(State = replace(State, which(State %in% c("Dadra & N Haveli", "Daman & Diu")), "Dadra and Nagar Haveli and Daman and Diu")) %>%
  mutate(State = replace(State, which(State %in% c("Andaman & N. Islands")), "Andaman and Nicobar Islands")) %>%
  mutate(State = gsub("&", "and", State)) %>%
  group_by(State) %>%
  summarise(hosp_beds = sum(Total_Hospital_Beds_Est))%>%
  rename(state = State)

# and fill in Ladakh using model based estimate for nation
hosp_beds$hosp_beds[hosp_beds$state == "Ladakh"] <-
  round(squire:::get_healthcare_capacity("India")$hosp_beds*subnat_pop$population[subnat_pop$state=="Ladakh"]/1000)
saveRDS(hosp_beds, cp_path("src/india_sub_national/hosp_beds.rds"))
saveRDS(hosp_beds, cp_path("analysis/data/derived/hosp_beds.rds"))

## --------------------------------
# demography source
# https://censusindia.gov.in/2011census/C-series/C-14.html
#
# DEPRECATED FOR MORE RECENT DEMOGRAPHY BELOW
# ------------------------------------------------------------------------------

# -> indian_demography.csv
# simpleCap <- function(x) {
#   s <- tolower(x)
#   s <- strsplit(s, " ")[[1]]
#   paste(toupper(substring(s, 1,1)), substring(s, 2),
#         sep="", collapse=" ")
# }
#
# demog <- read.csv("analysis/india/indian_demography.csv")
# demog <- filter(demog, !(Ages %in% c("All ages", "Age not stated"))) %>%
#   mutate(State = gsub("State - ", "", State)) %>%
#   filter(State != "India") %>%
#   mutate(State = gsub(" \\(\\d\\d\\)", "", State)) %>%
#   mutate(State = replace(State, which(State == "NCT OF DELHI"), "Delhi")) %>%
#   mutate(State = vapply(State, simpleCap, character(1))) %>%
#   mutate(State = replace(State, which(State %in% c("Dadra & Nagar Haveli", "Daman & Diu")), "Dadra and Nagar Haveli and Daman and Diu")) %>%
#   mutate(State = gsub("&", "and", State)) %>%
#   group_by(State, Ages) %>%
#   summarise(n = sum(n))
# demog$Ages <- factor(demog$Ages, levels = demog$Ages[c(1,10,2:9,11:17)])
# demog <- arrange(demog, State, Ages)
#
# # add in Ladakh and Telangana using national demography
# ind_pop <- squire::get_population("India")$n
# demog <- rbind(demog,
#                data.frame("State" = "Ladakh",
#                           "Ages" = as.character(unique(squire::population$age_group)),
#                           "n" = round((subnat_pop$population[subnat_pop$state == "Ladakh"]/sum(ind_pop))*ind_pop)),
#                data.frame("State" = "Telangana",
#                           "Ages" = as.character(unique(squire::population$age_group)),
#                           "n" = round((subnat_pop$population[subnat_pop$state == "Telangana"]/sum(ind_pop))*ind_pop))
#                ) %>%
# rename(state = State)
# demog <- left_join(demog, subnat_pop, by = "state") %>%
#   group_by(state) %>%
#   mutate(n = as.integer((population[1]/sum(n)) * n)) %>%
#   select(state, Ages, n)
# saveRDS(demog, "src/india_sub_national/demog.rds")

## --------------------------------
## New sourcing demog data
## --------------------------------

location <- "https://nhm.gov.in/New_Updates_2018/Report_Population_Projection_2019.pdf"
tf <- tempfile()
download.file(location, tf)

# Extract the table
out <- tabulizer::extract_tables(tf,pages = 263)

# Total population sizes per state
tot_cell <- out[[1]]
pops <- as.numeric(gsub(",", "", tail(tot_cell[,which(tot_cell[1,]=="2021") - 1], 37)))*1000
states <- tail(tot_cell[,grep("State", tot_cell[2,])], 37) %>%
  gsub("&", "and", .) %>%
  gsub(" \\*", "", .) %>%
  gsub("Dadra and Nagar Haveli", "Dadra and Nagar Haveli and Daman and Diu", .) %>%
  gsub("Daman Diu", "Dadra and Nagar Haveli and Daman and Diu", .) %>%
  gsub("NCT of Delhi", "Delhi", .) %>%
  gsub(" \\(UT\\)", "", .) %>%
  gsub("Uttaranchal", "Uttarakhand", .)
subnat_pop <- data.frame("state" = states, population = pops) %>%
  group_by(state) %>% summarise(population = sum(population))

# now for the breakdown by ages (easier to manually grab from the report)

# get india in total to begin and then get for all state reported and work out for those missing.
ind_age <- c(8.3,8.5,8.7,9.1,9.4,8.8,8.1,7.3,6.5,5.9,5.1,4.2,3.3,2.5,1.9,1.3,1.1)

df_list <- list()

df_list[[1]] <- data.frame(ages = c(6.2,6.7,10.9,10.4,9.9,9,8.3,7.7,6.8,5.9,5.1,4.1,3.1,2.3,1.6,1.1,1),
                           state = "Jammu and Kashmir")

df_list[[2]] <- data.frame(ages = c(6.3,6.8,7.5,8,8.5,8.6,8.4,7.9,7.3,6.7,5.9,5,4.1,3.2,2.4,1.6,1.8),
                           state = "Himachal Pradesh")

df_list[[3]] <- data.frame(ages = c(6.4,6.8,7.1,8,8.9,9.3,9,8.1,7.1,6.3,5.7,4.8,3.8,3,2.4,1.7,1.7),
                           state = "Punjab")

df_list[[4]] <- data.frame(ages = c(7.6,7.6,8.5,9.3,10.1,9.8,8.6,7.2,6.2,5.6,4.8,4,3.2,2.6,2.1,1.4,1.4),
                           state = "Uttarakhand")

df_list[[5]] <- data.frame(ages = c(8.3,8.4,8.3,8.8,9.5,9.5,8.7,7.6,6.6,5.7,4.9,3.9,3,2.4,2,1.3,1.2),
                           state = "Haryana")

df_list[[6]] <- data.frame(ages = c(7,7.3,7.5,8.5,9.8,9.9,9.1,8.4,7.3,6.4,5.3,4.2,3.2,2.3,1.7,1.1,0.9),
                           state = "Delhi")

df_list[[7]] <- data.frame(ages = c(9.7,10,9.5,9.8,9.9,9,7.8,6.7,5.8,5.2,4.4,3.7,2.8,2.1,1.6,1.1,0.9),
                           state = "Rajasthan")

df_list[[8]] <- data.frame(ages = c(10.3,10.5,9.4,10.2,10.4,9.3,7.6,6.2,5.5,5,4.2,3.3,2.6,2,1.6,1.1,0.8),
                           state = "Uttar Pradesh")

df_list[[9]] <- data.frame(ages = c(11,11.2,11.4,11.1,10,7.7,6.4,5.9,5.5,4.9,4.1,3.2,2.5,2,1.6,1,0.6),
                           state = "Bihar")

df_list[[10]] <- data.frame(ages = c(8.6,8.6,9.5,9.7,9.4,8.7,8.1,7.5,6.7,5.9,5,4,3,2.1,1.5,0.9,0.7),
                           state = "Assam")

df_list[[11]] <- data.frame(ages = c(6.5,6.9,7.7,8.3,8.9,9,8.7,7.9,7.1,6.7,6,5,3.9,2.9,2,1.3,1.2),
                            state = "West Bengal")

df_list[[12]] <- data.frame(ages = c(9.5,9.4,10.2,10.4,10.1,8.6,7.4,6.7,6,5.3,4.5,3.6,2.8,2.2,1.7,1.1,0.6),
                            state = "Jharkhand")

df_list[[13]] <- data.frame(ages = c(7.6,7.7,8.3,8.6,8.6,8.2,7.9,7.5,6.9,6.4,5.7,4.8,3.7,3,2.3,1.5,1.2),
                            state = "Odisha")

df_list[[14]] <- data.frame(ages = c(9.3,9.5,9,9.3,9.5,8.9,7.9,6.9,6.2,5.7,5,4.1,3,2.2,1.7,1,0.8),
                            state = "Chhattisgarh")

df_list[[15]] <- data.frame(ages = c(10,10.2,9.2,9.5,9.5,8.7,7.7,6.8,6.1,5.5,4.7,3.8,2.9,2.1,1.6,1,0.9),
                            state = "Madhya Pradesh")

df_list[[16]] <- data.frame(ages = c(8.3,8.3,8.1,8.6,9,8.8,8.3,7.6,6.8,6.1,5.3,4.4,3.5,2.6,1.8,1.2,1.1),
                            state = "Gujarat")

df_list[[17]] <- data.frame(ages = c(6.8,7.1,7.8,8.2,8.9,9,8.8,8.1,7.2,6.4,5.6,4.6,3.6,2.8,2.2,1.6,1.4),
                            state = "Maharashtra")

df_list[[18]] <- data.frame(ages = c(6.5,6.9,7.1,7.9,8.5,8.7,8.7,8.2,7.4,6.8,6,4.9,3.8,3,2.5,1.7,1.3),
                            state = "Andhra Pradesh")

df_list[[19]] <- data.frame(ages = c(7,7.5,7.6,8,8.5,8.8,8.8,8.2,7.3,6.5,5.6,4.7,3.7,2.8,2.2,1.5,1.3),
                            state = "Karnataka")

df_list[[20]] <- data.frame(ages = c(6.5,6.8,7.1,7.2,7.5,7.4,7.3,7.2,7,6.9,6.6,6,5.1,4.2,3.2,2.1,1.9),
                            state = "Kerala")

df_list[[21]] <- data.frame(ages = c(6.2,6.8,6.9,7.3,7.8,8.1,8.3,8.2,7.7,7.1,6.4,5.5,4.4,3.5,2.6,1.7,1.5),
                            state = "Tamil Nadu")

df_list[[22]] <- data.frame(ages = c(6.9,7.1,7.5,8.3,9,9.1,8.9,8.3,7.4,6.6,5.5,4.4,3.4,2.7,2.2,1.5,1.2),
                            state = "Telangana")

# NES States

df_list[[23]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Arunachal Pradesh")

df_list[[24]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Manipur")

df_list[[25]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Meghalaya")

df_list[[26]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Mizoram")

df_list[[27]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Nagaland")

df_list[[28]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Tripura")

df_list[[29]] <- data.frame(ages = c(7.2,6.9,9.5,9.8,9.9,9.5,8.8,7.9,6.8,5.9,5,4,3.1,2.2,1.5,1,1),
                            state = "Sikkim")


# take India for rest

df_list[[30]] <- data.frame(ages = c(ind_age),
                            state = "Chandigarh")

df_list[[31]] <- data.frame(ages = c(ind_age),
                            state = "Dadra and Nagar Haveli and Daman and Diu")

df_list[[32]] <- data.frame(ages = c(ind_age),
                            state = "Goa")

df_list[[33]] <- data.frame(ages = c(ind_age),
                            state = "Lakshadweep")

df_list[[34]] <- data.frame(ages = c(ind_age),
                            state = "Puducherry")

df_list[[35]] <- data.frame(ages = c(ind_age),
                            state = "Andaman and Nicobar Islands")

df_list[[36]] <- data.frame(ages = c(ind_age),
                            state = "Ladakh")


demog <- do.call(rbind, df_list)
demog$Ages <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                "30-34", "35-39", "40-44", "45-49", "50-54",
                "55-59", "60-64", "65-69", "70-74", "75-79", "80+")

demog <- demog %>% select(state, Ages, ages) %>%
  left_join(subnat_pop) %>%
  mutate(n = round(as.integer(ages/100*population))) %>%
  select(state, Ages, n)
saveRDS(demog, cp_path("src/india_sub_national/demog.rds"))
saveRDS(demog, cp_path("analysis/data/derived/demog.rds"))
write.csv(demog, cp_path("analysis/data/derived/indian_demography_2021.csv"), row.names = FALSE)

## --------------------------------
## sero sourcing
## --------------------------------

sero <- read.csv(cp_path("analysis/data/raw/india_serology.csv"))

# point estimate sero
sero$sero_pos <- sero$ns_pos
sero$sero_pos[which(is.na(sero$sero_pos))] <- sero$n_pos[which(is.na(sero$sero_pos))]
sero$sero <- sero$sero_pos/sero$samples

# downsample Delhi as the binomial ll is too strong
sero$samples[sero$state == "Delhi"] <- as.integer(sero$samples[sero$state == "Delhi"] * 0.05)
sero$n_pos[sero$state == "Delhi"] <- as.integer(sero$n_pos[sero$state == "Delhi"] * 0.05)
sero$sero_pos[sero$state == "Delhi"] <- as.integer(sero$sero_pos[sero$state == "Delhi"] * 0.05)

sero$sero_min <- Hmisc::binconf(sero$sero_pos, sero$samples)[,2]
sero$sero_max <- Hmisc::binconf(sero$sero_pos, sero$samples)[,3]
sero$date_end <- as.Date(sero$date_end)
sero$date_start <- as.Date(sero$date_start)

# remove non representative
sero <- sero[sero$notes != "slums",]

# remove selective cities that are not state-wide representative
sero <- sero[-which(sero$district == "Pune" & sero$notes == "Urban"),]
sero <- sero[-which(sero$district == "Hyderabad" & sero$notes == "Urban"),]
sero <- sero[-which(sero$district == "Ahmedabad (Round 1)" & sero$notes == "Urban"),]
sero <- sero[-which(sero$district == "Ahmedabad (Round 2)" & sero$notes == "Urban"),]
sero <- sero[-which(sero$district == "Chennai (Round 1)3" & sero$notes == "Urban"),]
sero <- sero[-which(sero$district == "Mumbai (Round-1)2"),]
sero <- sero[-which(sero$district == "Chennai (Round 2)" & sero$notes == "Urban"),]

# remove volunteer enrollment
sero <- sero[-which(sero$district == "Karnataka (5 regions)5" & sero$state == "Karnataka"),]
sero <- sero[-which(sero$district == "Ganjam" & sero$state == "Odisha" & sero$round == 2),]

# and get the correct round 2 and round 3 dates
sero$date_start[which(sero$round == 2)] <- "2020-08-18"
sero$date_start[which(sero$round == 3)] <- "2020-12-18"
sero$date_end[which(sero$round == 2)] <- "2020-09-20"
sero$date_end[which(sero$round == 3)] <- "2021-01-06"

saveRDS(sero, cp_path("src/india_sub_national/sero.rds"))
saveRDS(sero, cp_path("analysis/data/derived/sero.rds"))

## --------------------------------
## Mortality Reporting Data
## --------------------------------

# https://www.indiatoday.in/magazine/nation/story/20210719-why-official-covid-fatality-figures-grossly-understate-the-reality-1825897-2021-07-09
death_reporting <- data.frame(
  state = c("Jammu and Kashmir","Himachal Pradesh","Punjab","Uttarakhand","Haryana","Delhi","Rajasthan","Uttar Pradesh","Bihar","Assam","West Bengal","Jharkhand","Odisha","Chhattisgarh","Madhya Pradesh","Gujarat","Maharashtra","Andhra Pradesh","Karnataka","Kerala","Tamil Nadu","Telangana","Arunachal Pradesh","Manipur","Meghalaya","Mizoram","Nagaland","Tripura","Sikkim","Chandigarh","Dadra and Nagar Haveli and Daman and Diu","Goa","Lakshadweep","Puducherry","Andaman and Nicobar Islands","Ladakh"
))

death_reporting$death_registration <- c(66.7, 86.4, 100, 95.6, 100, 100, 98.6, 63.3, 51.6, 74, 100, 58.8, 100,
                                        81.5, 89.1, 100, 100, 100, 100, 100, 100, 97.2, 38.6, 21.4, 97.6, 100,
                                        30, 100, 100, 100, 82.1, 100, 88.2, 100, 100, 66.7)
death_reporting$death_certification <- c(NA, 13, 17.5, 8.9, 19.4, 100, 13.9, 6.5, 5.1, 17.2, 14.1, 5.8, 12.6,
                                         21.4, 9.1, 21.3, 38.2, 12.9, 30.4, 11.6, 44, 27.7, 33.4, 67.3, 32.9, 51.9,
                                         12, 33.7, 45.5, 74.4, 50.2, 100, 95.8, 100, 60.1, NA)

death_reporting$death_certification[is.na(death_reporting$death_certification)] <-
  lm(death_certification ~ exp(death_registration),
     death_reporting %>% filter(death_registration>40)) %>%
  predict(newdata = data.frame("death_registration"=66.7))

saveRDS(death_reporting, cp_path("analysis/data/derived/death_reporting.rds"))
