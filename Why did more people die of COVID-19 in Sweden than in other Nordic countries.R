library(tidyverse)
library(lubridate)
library(furrr)
library(progressr)
library(cowplot)
library(zoo)

# see https://arxiv.org/pdf/2005.02090.pdf for the infection-to-death distribution I use by default
infer_infections_draw <- function(
  deaths,
  ifr = 0.01,
  mu_i2d = 3.19,
  sigma_i2d = 0.44
) {
  data_length <- length(deaths)
  
  infections <- rep(0, data_length)
  
  for (i in 1:data_length) {
    if (deaths[i] > 0) {
      infection_lags <- rlnorm(deaths[i], mu_i2d, sigma_i2d)
      infections_per_death <- rgeom(deaths[i], ifr)
      for (j in 1:deaths[i]) {
        infection_date <- i - floor(infection_lags[j])
        
        # make sure we don't impute infections outside of the date range in the data
        if (infection_date > 0) {
          infections[i - infection_lags[j]] <- infections[i - infection_lags[j]] + infections_per_death[j]
        }
      }
    }
  }
  
  infections
}

infer_infections <- function(
  data,
  countries,
  n_sims,
  ifr,
  mu_i2d,
  sigma_i2d,
  rng_seed
) {
  draws <- list()
  
  for (i in 1:length(countries)) {
    c <- countries[i]
    
    country_data <- data %>%
      filter(country == c) %>%
      arrange(date)
    
    print(paste0("Estimation of infections for ", c, " starts"))
    
    args <- list(
      deaths = rep(list(country_data$deaths), n_sims),
      ifr = rep(list(ifr[i]), n_sims),
      mu_i2d = rep(list(mu_i2d[i]), n_sims),
      sigma_i2d = rep(list(sigma_i2d[i]), n_sims)
      )
    
    with_progress({
      p <- progressor(steps = n_sims)
      
      draws[[c]] <- future_pmap_dfr(
        args,
        function(deaths, ifr, mu_i2d, sigma_i2d) {
          p()
          tibble(
            date = country_data$date,
            country = c,
            infections = infer_infections_draw(deaths, ifr, mu_i2d, sigma_i2d)
          )
        },
        .options = furrr_options(seed = rng_seed)
      )
    })
    
    print(paste0("Estimation of infections for ", c, " completed"))
  }
  
  results <- bind_rows(draws) %>%
    group_by(date, country) %>%
    summarize(
      mean = mean(infections),
      median = quantile(infections, probs = 0.5),
      lower1 = quantile(infections, probs = 0.025),
      upper1 = quantile(infections, probs = 0.975),
      lower2 = quantile(infections, probs = 0.25),
      upper2 = quantile(infections, probs = 0.75)
    ) %>%
    ungroup()
}

# see https://www.sciencedirect.com/science/article/pii/S2468042720300634 for the generation time distribution I use by default
simulate_counterfactual_epidemic <- function(
  infections_before_lockdown,
  R_after_lockdown,
  simulation_length,
  mean_gt = 3.99,
  sd_gt = 2.96,
  ifr = 0.01,
  mu_i2d = 3.19,
  sigma_i2d = 0.44
) {
  # parameters of the gamma distribution I use to model the generation time distribution
  gt_alpha <- mean_gt^2 / sd_gt^2
  gt_beta <- mean_gt / sd_gt^2
  
  simulation_start <- length(infections_before_lockdown) + 1
  total_length <- simulation_start + simulation_length - 1
  
  # simulate infections
  infections <- c(infections_before_lockdown, rep(0, simulation_length))
  for (i in 1:total_length) {
    if (infections[i] >= 1) {
      n <- round(infections[i])
      generation_times <- round(rgamma(n, shape = gt_alpha, rate = gt_beta))
      for (j in 1:n) {
        if (i + generation_times[j] >= simulation_start & i + generation_times[j] <= total_length) {
          infections[i + generation_times[j]] <- infections[i + generation_times[j]] + R_after_lockdown
        }
      }
    }
  }
  
  infections <- round(infections)
  
  # simulate deaths
  deaths <- rep(0, total_length)
  for (i in 1:total_length) {
    if (infections[i] >= 1) {
      fatal <- rbinom(infections[i], 1, ifr)
      for (j in 1:infections[i]) {
        if (fatal[j]) {
          death_lag <- round(rlnorm(1, mu_i2d, sigma_i2d))
          if (i + death_lag <= total_length) {
            deaths[i + death_lag] <- deaths[i + death_lag] + 1
          }
        }
      }
    }
  }
  
  tibble(
    infections = infections,
    deaths = deaths
  )
}

simulate_infections_with_counterfactual <- function(
  data,
  countries,
  n_sims,
  lockdowns,
  R_after_lockdown,
  simulation_length,
  mean_gt,
  sd_gt,
  ifr,
  mu_i2d,
  sigma_i2d,
  rng_seed
) {
  draws <- list()
  
  for (i in 1:length(countries)) {
    c <- countries[i]
    
    country_data <- data %>%
      filter(country == c) %>%
      arrange(date)
    
    lockdown_start <- as.integer(lockdowns$date[lockdowns$country == c] - country_data$date[1] + 1)
    
    print(paste0("Estimation of infections for ", c, " starts"))
    
    # I only use deaths that occurred up to 60 days after the start of the lockdown in this country since I will only
    # use infections that occurred before that date for the counterfactual and presumably no one or almost no one who
    # died more than 60 days later should have been infected before that, so this will avoid unnecessary computations
    args1 <- list(
      deaths = rep(list(country_data$deaths[1:(lockdown_start + 60)]), n_sims),
      ifr = rep(list(ifr[i]), n_sims),
      mu_i2d = rep(list(mu_i2d[i]), n_sims),
      sigma_i2d = rep(list(sigma_i2d[i]), n_sims)
      )
    
    with_progress({
      p <- progressor(steps = n_sims)
      
      infections_actual <- future_pmap(
        args1,
        function(deaths, ifr, mu_i2d, sigma_i2d) {
          p()
          infer_infections_draw(deaths, ifr, mu_i2d, sigma_i2d)
        },
        .options = furrr_options(seed = rng_seed)
      )
    })
    
    print(paste0("Estimation of infections for ", c, " completed"))
    
    print(paste0("Estimation of counterfactual for ", c, " starts"))
    
    infections_actual <- map(infections_actual, ~ .x[1:lockdown_start])
    
    args2 <- list(
      infections_before_lockdown = infections_actual,
      R_after_lockdown = rep(list(R_after_lockdown), n_sims),
      simulation_length = rep(list(simulation_length), n_sims),
      mean_gt = rep(list(mean_gt[i]), n_sims),
      sd_gt = rep(list(sd_gt[i]), n_sims),
      ifr = rep(list(ifr[i]), n_sims),
      mu_i2d = rep(list(mu_i2d[i]), n_sims),
      sigma_i2d = rep(list(sigma_i2d[i]), n_sims)
      )
    
    with_progress({
      p <- progressor(steps = n_sims)
      
      draws[[c]] <- future_pmap_dfr(
        args2,
        function(infections_before_lockdown, R_after_lockdown, simulation_length, mean_gt, sd_gt, ifr, mu_i2d, sigma_i2d) {
          p()
          simulate_counterfactual_epidemic(
            infections_before_lockdown,
            R_after_lockdown,
            simulation_length,
            mean_gt,
            sd_gt,
            ifr,
            mu_i2d,
            sigma_i2d
          ) %>%
            add_column(
              date = country_data$date[1:(lockdown_start + simulation_length)],
              country = c
            )
        },
        .options = furrr_options(seed = rng_seed)
      )
    })
    
    draws[[c]] <- bind_rows(draws[[c]])
    
    print(paste0("Estimation of counterfactual for ", c, " completed"))
  }
  
  results <- bind_rows(draws) %>%
    group_by(date, country) %>%
    summarize(
      infections_mean = mean(infections),
      infections_median = quantile(infections, probs = 0.5),
      infections_lower1 = quantile(infections, probs = 0.025),
      infections_upper1 = quantile(infections, probs = 0.975),
      infections_lower2 = quantile(infections, probs = 0.25),
      infections_upper2 = quantile(infections, probs = 0.75),
      deaths_mean = mean(deaths),
      deaths_median = quantile(deaths, probs = 0.5),
      deaths_lower1 = quantile(deaths, probs = 0.025),
      deaths_upper1 = quantile(deaths, probs = 0.975),
      deaths_lower2 = quantile(deaths, probs = 0.25),
      deaths_upper2 = quantile(deaths, probs = 0.75)
    ) %>%
    ungroup()
}

# allow parallel processing
plan(multisession)

# seed for reproducibility
# (see https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/ and https://davisvaughan.github.io/furrr/reference/furrr_options.html#reproducible-random-number-generation-rng-)
rng_seed <- 127

owid_data_url <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"

owid_data <- read_csv(owid_data_url) %>%
  rename(
    country = location
  ) %>%
  mutate(
    deaths = replace_na(new_deaths, 0),
    rolling_average_deaths = round(rollmean(deaths, 7, fill = c(0, 0, 0), align = "right"))
  ) %>%
  select(
    date,
    country,
    deaths,
    rolling_average_deaths
  )

# sometimes the number of deaths for a day is negative because there was a data correction, but it's almost never the
# case for Nordic countries and, even when it happens, the corrections are tiny, so I just set the number of deaths on
# those days to zero instead of trying to deal with that in a fancier way because it wouldn't make any difference anyway
owid_data <- owid_data %>%
  mutate(deaths = ifelse(deaths > 0, deaths, 0))

# make sure the data start on 31 December 2019 for every country
for (c in unique(owid_data$country)) {
  owid_data <- owid_data %>%
    bind_rows(
      tibble(
        date = seq(ymd("2019-12-31"), min(owid_data$date[owid_data$country == c]) - 1, by = "1 day"),
        country = c,
        deaths = rep(0, min(owid_data$date[owid_data$country == c]) - ymd("2019-12-31"))
      )
    )
}

owid_data <- owid_data %>%
  arrange(date, country)

nordic_countries <- c(
  "Denmark",
  "Finland",
  "Norway",
  "Sweden"
)

estimates_nordic <- infer_infections(
  owid_data,
  nordic_countries,
  1e5,
  ifr = c(0.01, 0.01, 0.01, 0.01),
  mu_i2d = c(3.19, 3.19, 3.19, 3.19),
  sigma_i2d = c(0.44, 0.44, 0.44, 0.44),
  rng_seed = rng_seed
)

data_for_plot_nordic <- estimates_nordic %>%
  filter(date >= ymd("2020-02-01") & date <= ymd("2020-06-01"))

# https://en.wikipedia.org/wiki/COVID-19_pandemic_in_Norway#Prevention_measures_and_response
# https://en.wikipedia.org/wiki/COVID-19_pandemic_in_Denmark#Lockdown
# https://en.wikipedia.org/wiki/COVID-19_pandemic_in_Finland#Government_response
lockdowns <- tribble(
  ~country, ~date,
  "Norway", ymd("2020-03-12"),
  "Denmark", ymd("2020-03-18"),
  "Finland", ymd("2020-03-16")
)

ggplot(data_for_plot_nordic, aes(x = date, y = median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower1, ymax = upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = lower2, ymax = upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_vline(data = lockdowns, aes(xintercept = date, color = country), linetype = "dashed", size = 1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 infections in Nordic countries\n(inferred from deaths)") +
  xlab("Date") +
  ylab("Daily number of COVID-19 infections") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  ggsave(
    "Daily number of COVID-19 infections in Nordic countries.png",
    width = 12,
    height = 6
  )

counterfactual_lockdowns <- tribble(
  ~country, ~date,
  "Denmark", ymd("2020-03-15"),
  "Finland", ymd("2020-03-15"),
  "Norway", ymd("2020-03-15"),
  "Sweden", ymd("2020-03-15")
)

counterfactual_estimates <- simulate_infections_with_counterfactual(
  owid_data,
  nordic_countries,
  1e5,
  counterfactual_lockdowns,
  0.6,
  120,
  mean_gt = c(3.99, 3.99, 3.99, 3.99),
  sd_gt = c(2.96, 2.96, 2.96, 2.96),
  ifr = c(0.01, 0.01, 0.01, 0.01),
  mu_i2d = c(3.19, 3.19, 3.19, 3.19),
  sigma_i2d = c(0.44, 0.44, 0.44, 0.44),
  rng_seed = rng_seed
)

data_for_plot_counterfactual <- counterfactual_estimates %>%
  filter(date >= ymd("2020-02-01") & date <= ymd("2020-07-01"))

counterfactual_infections_plot <- ggplot(data_for_plot_counterfactual, aes(x = date, y = infections_median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = infections_lower1, ymax = infections_upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = infections_lower2, ymax = infections_upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 infections in Nordic countries in counterfactual where\nR was suddenly reduced to 0.6 everywhere on March 15") +
  xlab("Date") +
  ylab("Daily number of COVID-19 infections") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

counterfactual_deaths_plot <- ggplot(data_for_plot_counterfactual, aes(x = date, y = deaths_median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = deaths_lower1, ymax = deaths_upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = deaths_lower2, ymax = deaths_upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 deaths in Nordic countries in counterfactual where\nR was suddenly reduced to 0.6 everywhere on March 15") +
  xlab("Date") +
  ylab("Daily number of COVID-19 deaths") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

plot_grid(
  counterfactual_infections_plot,
  counterfactual_deaths_plot,
  labels = c("", ""),
  ncol = 1
) +
  ggsave(
    "Daily number of COVID-19 infections and deaths in Nordic countries in counterfactual where R was suddenly reduced to 0.6 everywhere on March 15.png",
    width = 12,
    height = 12
  )

raw_data_for_plot <- owid_data %>%
  filter(country %in% nordic_countries & date >= ymd("2020-03-01") & date <= ymd("2020-07-01"))

ggplot(raw_data_for_plot, aes(x = date, y = rolling_average_deaths, group = country, color = country)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Seven-day rolling average of daily number of deaths in Nordic countries") +
  xlab("Date") +
  ylab("Daily number of deaths") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "Source : Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  ggsave("Seven-day rolling average of daily number of deaths in Nordic countries.png", width = 12, height = 6)

counterfactual_estimates_sensitivity <- simulate_infections_with_counterfactual(
  owid_data,
  nordic_countries,
  1e5,
  counterfactual_lockdowns,
  0.6,
  120,
  mean_gt = c(3.99, 3.99, 3.99, 3.99),
  sd_gt = c(2.96, 2.96, 2.96, 2.96),
  ifr = c(0.005, 0.005, 0.005, 0.01),
  mu_i2d = c(3.19, 3.19, 3.19, 2.95),
  sigma_i2d = c(0.44, 0.44, 0.44, 0.44),
  rng_seed = rng_seed
)

data_for_plot_counterfactual_sensitivity <- counterfactual_estimates_sensitivity %>%
  filter(date >= ymd("2020-02-01") & date <= ymd("2020-07-01"))

counterfactual_infections_plot_sensitivity <- ggplot(data_for_plot_counterfactual_sensitivity, aes(x = date, y = infections_median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = infections_lower1, ymax = infections_upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = infections_lower2, ymax = infections_upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 infections in Nordic countries in counterfactual where\nR was suddenly reduced to 0.6 everywhere on March 15 (IFR = 0.5% in Denmark/Finland/Norway and 1% in Sweden)") +
  xlab("Date") +
  ylab("Daily number of COVID-19 infections") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

counterfactual_deaths_plot_sensitivity <- ggplot(data_for_plot_counterfactual_sensitivity, aes(x = date, y = deaths_median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = deaths_lower1, ymax = deaths_upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = deaths_lower2, ymax = deaths_upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 deaths in Nordic countries in counterfactual where\nR was suddenly reduced to 0.6 everywhere on March 15 (IFR of 0.5% in Denmark/Finland/Norway and 1% in Sweden)") +
  xlab("Date") +
  ylab("Daily number of COVID-19 deaths") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

plot_grid(
  counterfactual_infections_plot_sensitivity,
  counterfactual_deaths_plot_sensitivity,
  labels = c("", ""),
  ncol = 1
) +
  ggsave(
    "Daily number of COVID-19 infections and deaths in Nordic countries in counterfactual where R was suddenly reduced to 0.6 everywhere on March 15 (first sensitivity analysis).png",
    width = 12,
    height = 12
  )

counterfactual_estimates_sensitivity2 <- simulate_infections_with_counterfactual(
  owid_data,
  nordic_countries,
  1e5,
  counterfactual_lockdowns,
  0.6,
  120,
  mean_gt = c(3.99, 3.99, 3.99, 3.99),
  sd_gt = c(2.96, 2.96, 2.96, 2.96),
  ifr = c(0.005, 0.002, 0.002, 0.01),
  mu_i2d = c(3.19, 3.19, 3.19, 2.95),
  sigma_i2d = c(0.44, 0.44, 0.44, 0.44),
  rng_seed = rng_seed
)

data_for_plot_counterfactual_sensitivity2 <- counterfactual_estimates_sensitivity2 %>%
  filter(date >= ymd("2020-02-01") & date <= ymd("2020-07-01"))

counterfactual_infections_plot_sensitivity2 <- ggplot(data_for_plot_counterfactual_sensitivity2, aes(x = date, y = infections_median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = infections_lower1, ymax = infections_upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = infections_lower2, ymax = infections_upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 infections in Nordic countries in counterfactual where\nR was suddenly reduced to 0.6 everywhere on March 15 (IFR of 0.5% in Denmark, 0.2% in Finland/Norway and 1% in Sweden)") +
  xlab("Date") +
  ylab("Daily number of COVID-19 infections") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

counterfactual_deaths_plot_sensitivity2 <- ggplot(data_for_plot_counterfactual_sensitivity2, aes(x = date, y = deaths_median, group = country, color = country)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = deaths_lower1, ymax = deaths_upper1), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = deaths_lower2, ymax = deaths_upper2), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle("Daily number of COVID-19 deaths in Nordic countries in counterfactual where\nR was suddenly reduced to 0.6 everywhere on March 15 (IFR of 0.5% in Denmark, 0.2% in Finland/Norway and 1% in Sweden)") +
  xlab("Date") +
  ylab("Daily number of COVID-19 deaths") +
  scale_color_discrete(name = "Country") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source of data: Our World in Data - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

plot_grid(
  counterfactual_infections_plot_sensitivity2,
  counterfactual_deaths_plot_sensitivity2,
  labels = c("", ""),
  ncol = 1
) +
  ggsave(
    "Daily number of COVID-19 infections and deaths in Nordic countries in counterfactual where R was suddenly reduced to 0.6 everywhere on March 15 (second sensitivity analysis).png",
    width = 12,
    height = 12
  )

save.image("results.Rdata")
