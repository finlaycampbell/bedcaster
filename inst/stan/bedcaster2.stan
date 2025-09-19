functions {

  // reverse function used in convolution
  vector stan_rev(vector a) {
    int a_size = num_elements(a);
    vector [a_size] a_rev;
    
    for(i in 1:a_size)  
      a_rev[i] = a[a_size - i + 1];
   
    return a_rev;
    
  }

  // convolution function from https://discourse.mc-stan.org/t/
  // dot-products-of-vectors-to-perform-1d-convolution/9053/9 
  vector stan_convolve(int nPoints, vector a,  vector b) {
  
    vector [nPoints] out;
    vector [nPoints] b_rev;
  
    b_rev = stan_rev(b);
  
    for(i in 1:nPoints)  
      out[i] = dot_product(head(a, i), tail(b_rev, i));
 
    return(out);
    
  }


}

data {

  // Number of days with observations
  int<lower = 1> n_obs;

  // Number of days of projection
  int<lower = 1> n_proj;

  // Max of delay distribution
  int<lower = 1> max_delay;

  // Number of reported cases per day
  int cases_reported[n_obs];

  // Number of reported dethas per day
  int deaths_reported[n_obs];

  // ETU occupancy by day
  int etu_reported[n_obs];

  // Number of days with etu occupancy data
  int etu_n;

  // Indices of days with etu occupancy data
  int etu_ind[etu_n];

  // New community alerts by day
  int alerts_reported[n_obs];

  // Number of days with alerts data
  int alerts_n;

  // Indices of days with alerts data
  int alerts_ind[alerts_n];

  // Isolation bed occupancy by day
  int iso_reported[n_obs];

  // Number of days with iso occupancy data
  int iso_n;

  // Indices of days with iso occupancy data
  int iso_ind[iso_n];

  // mean and sd of the normal prior on the logmean of the delay
  // from onset to hospitalisation
  vector[4] prior_onset_to_etu;

  // mean and sd of the normal prior on the logmean of the delay
  // from hospitalisation to survival
  vector[4] prior_etu_to_survival;

  // mean and sd of the normal prior on the logmean of the delay
  // from hospitalisation to death
  vector[4] prior_etu_to_death;

  // mean and sd of the normal prior on the logmean of the delay
  // from case-onset to contact-isolation
  vector[4] prior_onset_to_iso;

  // mean and sd of the normal prior on the logmean of the delay
  // from in lab turnaround
  vector[4] prior_iso_to_release;

  // mean and sd of the normal prior on the plogis transformed CFR
  vector[2] prior_cfr;
  
  // mean and sd of the normal prior on the plogis transformed
  // proportion of alerts isolated
  vector[2] prior_prop_iso;
  
  // mean and sd of the normal prior on the exponentially transformed
  // number of background alerts
  vector[2] prior_alerts_background;

  // mean and sd of the normal prior on the exponentially transformed
  // number of alerts per case
  vector[2] prior_alerts_per_case;

  // delay distribution from onset to case confirmation
  // and entry into database
  vector[n_obs] log_prop_cases_reported;

  // delay distribution from onset to death confirmation
  // and entry into database
  vector[n_obs] prop_deaths_reported;

  // number of alerts_background bins
  int n_alerts_background;
  
  // alerts_background indices
  int alerts_background_ind[n_obs];

  // Number of spline parameters
  int n_spline_param;

  // spline matrix for growth rate fitting
  matrix[n_spline_param, n_obs] spline;

  // asymptote for growth rate  
  vector[n_proj] growthrate_asymptote_weight;

}

transformed data {

  // Build projection index
  vector[n_proj] projection_sq;
  for (i in 1:n_proj) projection_sq[i] = i;

   // Total numbers of days in 
   int n_tot;
   n_tot = n_obs + n_proj;

}

parameters {

  // in-hospital CFR
  real cfr_logit;

  // proportion of alerts isolated
  real prop_iso_logit;

  // logmean delay from onset to hospitalisation
  real onset_to_etu_logmean;
  real<lower = 0> onset_to_etu_sd;  

  // logmean delay from hospitalisation to survival
  real etu_to_survival_logmean;
  real<lower = 0> etu_to_survival_sd;  

  // logmean delay from hospitalisation to death
  real etu_to_death_logmean;
  real<lower = 0> etu_to_death_sd;  

  // logmean delay from case-onset to contact-isolation
  real onset_to_iso_logmean;
  real<lower = 0> onset_to_iso_sd;  

  // logmean delay in lab turnaround
  real iso_to_release_logmean;
  real<lower = 0> iso_to_release_sd;  

  // overdispersion parameter for nbinom sampling of cases, deaths,
  // etu, occupancy, alerts and isolation bed occupancy
  real cases_overdisp_log;
  real deaths_overdisp_log;  
  real etu_overdisp_log;
  real alerts_overdisp_log;
  real iso_overdisp_log;

  // case number at day 0
  real log_cases_intercept;

  // mean number of background of alerts per day (exponentially transformed)
  // real alerts_background;
  vector[n_alerts_background] alerts_background;  

  // mean number of alerts per case
  real alerts_per_case;

  // growth rates
  // vector[n_growth_rate] growth_rate;

  // untransformed spline parameters
  row_vector[n_spline_param] spline_param;

}

transformed parameters {

  // GROWTHRATE CALCULATIONS

  // GLOBAL VARIABLES
  
  // vector of estimated growth rates
  vector[n_obs] growthrate_reported;

  // vector of projected growth rates
  vector[n_proj] growthrate_projected;

  {
  
  // params for growthrate extrapolation
  real growthrate_slope;
  real growthrate_slope_weight;
  real growthrate_slope_weight_min = -0.2;
  real growthrate_slope_weight_max = 0.2;
  real growthrate_asymptote = -0.1;

  // calculate growth rate from spline params
  growthrate_reported = to_vector(spline_param*spline)/10;

  // calculate slope of growth rate at most recent time
  growthrate_slope = growthrate_reported[n_obs] - growthrate_reported[n_obs-1];

  // growth rate slopes are weighted more when the absolute growth rate is high
  growthrate_slope_weight =
    (growthrate_reported[n_obs] - growthrate_slope_weight_min)/
    (growthrate_slope_weight_max - growthrate_slope_weight_min);
  if(growthrate_slope_weight > 1) growthrate_slope_weight = 1;
  if(growthrate_slope_weight < 0) growthrate_slope_weight = 0;

  // minimum slope weighting of 0.2, maximum slope weighting of 1.0
  growthrate_slope_weight = 0.2 + 0.8*growthrate_slope_weight;

  // extrapolate growth rate for projection using slope and slope weighting
  growthrate_projected = growthrate_reported[n_obs] +
    projection_sq * growthrate_slope * growthrate_slope_weight;

  // add growth rate asymptote (weighted mean of extrapolation and asymptote)
  growthrate_projected = growthrate_asymptote * growthrate_asymptote_weight +
    (1 - growthrate_asymptote_weight).*growthrate_projected;	  

  }


  // CASE CALCULATIONS

  // GLOBAL VARIABLES

  // vector of fitted, nowcast and projected cases
  vector[n_obs] cases_truncated_mu;
  vector[n_obs] cases_nowcast_mu;  
  vector[n_proj] cases_projected_mu;

  {

  // calculate in log space
  vector[n_obs] log_cases_truncated_mu;
  vector[n_obs] log_cases_nowcast_mu;  
  vector[n_proj] log_cases_projected_mu;

  // model case observations without reporting delay for first day
  log_cases_nowcast_mu[1] = log_cases_intercept + growthrate_reported[1];

  // model case observations with reporting delay
  log_cases_truncated_mu[1] = log_cases_nowcast_mu[1] + log_prop_cases_reported[1];  
  for(i in 2:n_obs) {
    log_cases_nowcast_mu[i] = log_cases_nowcast_mu[i-1] + growthrate_reported[i];
    log_cases_truncated_mu[i] = log_cases_nowcast_mu[i] + log_prop_cases_reported[i];
  }

  // cap values
  for(i in 1:n_obs) {
    if(log_cases_truncated_mu[i] > 10) log_cases_truncated_mu[i] = 10;
    if(log_cases_nowcast_mu[i] > 10) log_cases_nowcast_mu[i] = 10;    
  }

  // projected cases forward using nowcast values
  log_cases_projected_mu[1] = log_cases_nowcast_mu[n_obs] + growthrate_projected[1];
  for(i in 2:n_proj) {
    log_cases_projected_mu[i] = log_cases_projected_mu[i-1] + growthrate_projected[i];  
    if(log_cases_projected_mu[i] > 10) log_cases_projected_mu[i] = 10;
    if(log_cases_projected_mu[i] < -10) log_cases_projected_mu[i] = -10;
  }

  cases_truncated_mu = exp(log_cases_truncated_mu);
  cases_nowcast_mu = exp(log_cases_nowcast_mu);
  cases_projected_mu = exp(log_cases_projected_mu);  

  }


  // DELAY CALCULATIONS

  // GLOBAL VARIABLES

  // raw delay densities for case-to-etu model
  vector[max_delay] onset_to_etu;
  vector[max_delay] etu_to_survival;
  vector[max_delay] etu_to_death;

  // convolved densities
  vector[max_delay] onset_to_survival;
  vector[max_delay] onset_to_death;

  // raw delay densities for case-to-isolation
  vector[max_delay] onset_to_iso;  
  vector[max_delay] iso_to_release;

  // convolved density
  vector[max_delay] onset_to_release;

  // define cfr and log cfr
  real cfr;
  cfr = inv_logit(cfr_logit);

  {
  
  // calculate delay pmf from lognormal and normalise
  for(i in 1:max_delay) {
    onset_to_etu[i] = exp(
      lognormal_lpdf(i | onset_to_etu_logmean, onset_to_etu_sd)
    );
    etu_to_survival[i] = exp(
      lognormal_lpdf(i | etu_to_survival_logmean, etu_to_survival_sd)
    );
    etu_to_death[i] = exp(
      lognormal_lpdf(i | etu_to_death_logmean, etu_to_death_sd)
    );
    onset_to_iso[i] = exp(
      lognormal_lpdf(i | onset_to_iso_logmean, onset_to_iso_sd)
    );
    iso_to_release[i] = exp(
      lognormal_lpdf(i | iso_to_release_logmean, iso_to_release_sd)
    );    
  }

  // normalise pdf
  onset_to_etu = onset_to_etu ./ sum(onset_to_etu);
  etu_to_survival = etu_to_survival ./ sum(etu_to_survival);
  etu_to_death = etu_to_death ./ sum(etu_to_death);
  onset_to_iso = onset_to_iso ./ sum(onset_to_iso);
  iso_to_release = iso_to_release ./ sum(iso_to_release);

  // convolve onset to survival and include CFR
  onset_to_survival = stan_convolve(max_delay, onset_to_etu, etu_to_survival);
  onset_to_survival = onset_to_survival ./ sum(onset_to_survival) * (1 - cfr);

  // convolve onset to death and include CFR
  onset_to_death = stan_convolve(max_delay, onset_to_etu, etu_to_death);
  onset_to_death = onset_to_death ./ sum(onset_to_death) * cfr;

  // convolve onset to release-from-isolation
  onset_to_release = stan_convolve(max_delay, onset_to_iso, iso_to_release);
  onset_to_release = onset_to_release ./ sum(onset_to_release);

  }
  

  // DEATH CALCULATIONS

  // GLOBAL VARIABLES

  // vector of fitted and projected deaths
  vector[n_obs] deaths_truncated_mu;
  vector[n_obs] deaths_nowcast_mu;  
  vector[n_proj] deaths_projected_mu;

  {

  deaths_truncated_mu = rep_vector(exp(-10), n_obs);
  deaths_nowcast_mu = rep_vector(exp(-10), n_obs);
  deaths_projected_mu = rep_vector(exp(-10), n_proj);

  real density;
  for(i in 1:n_tot) {
    for(j in 1:max_delay) {
      if(i <= n_obs)
        density = cases_nowcast_mu[i] * onset_to_death[j];
      else
        density = cases_projected_mu[i-n_obs] * onset_to_death[j];
      if((i+j) <= n_obs)
        deaths_nowcast_mu[i+j] += density;
      else if((i+j) <= n_tot)
        deaths_projected_mu[i+j-n_obs] += density;
    }
  }

  // reported deaths are nowcast deaths scaled by proportion observed
  for(i in 1:n_obs) deaths_truncated_mu[i] = deaths_nowcast_mu[i] *
    prop_deaths_reported[i];
    
  }
  

  // ALERTS AND HOSPITAL CALCULATIONS

  // GLOBAL VARIABLES

  // numbers of alerts, etu occupancy and iso occupancy
  vector[n_tot] alerts_truncated_mu;
  vector[n_tot] etu_truncated_mu;  
  vector[n_tot] iso_truncated_mu;

  // plogis transform of proportion of alerts isolated
  real prop_iso;
  prop_iso = inv_logit(prop_iso_logit);

  {

  // declare local variables (admission and discharge densities)
  vector[n_tot] etu_admission;
  vector[n_tot] etu_discharge;
  vector[n_tot] iso_admission;
  vector[n_tot] iso_discharge;

  // initialise
  alerts_truncated_mu = rep_vector(0, n_tot);
  etu_admission = rep_vector(0, n_tot);
  etu_discharge = rep_vector(0, n_tot);
  iso_admission = rep_vector(0, n_tot);
  iso_discharge = rep_vector(0, n_tot);

  // calculate admission and discharge
  for (i in 1:(n_tot - 1)) {
    int j_max = (i + max_delay > n_tot) ? (n_tot - i) : max_delay;
  
    // pick the case source (reported vs projected)
    real cases_val;
    real alerts_bg;
    if (i <= n_obs) {
      cases_val = cases_reported[i];
      alerts_bg = exp(alerts_background[alerts_background_ind[i]]);
    } else {
      cases_val = cases_projected_mu[i - n_obs];
      alerts_bg = exp(alerts_background[alerts_background_ind[n_alerts_background]]);
    }

    // precompute alerts term
    real alerts_val = cases_val * exp(alerts_per_case) + alerts_bg;

    for (j in 1:j_max) {
      int t = i + j;

      etu_admission[t]   += cases_val * onset_to_etu[j];
      etu_discharge[t]   += cases_val * (onset_to_survival[j] + onset_to_death[j]);

      alerts_truncated_mu[t] += alerts_val * onset_to_iso[j];
      iso_admission[t]   += alerts_val * onset_to_iso[j] * prop_iso;
      iso_discharge[t]   += alerts_val * onset_to_release[j] * prop_iso;
    }
  }

  // calculate etu and iso occupancy
  etu_truncated_mu = cumulative_sum(etu_admission - etu_discharge);
  iso_truncated_mu = cumulative_sum(iso_admission - iso_discharge);

  for(i in 1:(n_tot)) if(etu_truncated_mu[i] == 0) etu_truncated_mu[i] = 0.00001;
  for(i in 1:(n_tot)) if(alerts_truncated_mu[i] == 0) alerts_truncated_mu[i] = 0.00001;
  for(i in 1:(n_tot)) if(iso_truncated_mu[i] <= 0) iso_truncated_mu[i] = 0.00001;

  }


  // OVERDISPERSION PARAMETERS

  // GLOBAL VARIBLES
  
  // exponential transform of overdispersion parameters
  real cases_overdisp;
  real deaths_overdisp;  
  real etu_overdisp;
  real alerts_overdisp;
  real iso_overdisp;

  cases_overdisp = exp(cases_overdisp_log);
  deaths_overdisp = exp(deaths_overdisp_log);  
  etu_overdisp = exp(etu_overdisp_log);
  alerts_overdisp = exp(alerts_overdisp_log);
  iso_overdisp = exp(iso_overdisp_log);

}

model {

  // prior on delays
  onset_to_etu_logmean ~ normal(prior_onset_to_etu[1], prior_onset_to_etu[2]);
  onset_to_etu_sd ~ normal(prior_onset_to_etu[3], prior_onset_to_etu[4]);
  
  etu_to_survival_logmean ~ normal(prior_etu_to_survival[1], prior_etu_to_survival[2]);
  etu_to_survival_sd ~ normal(prior_etu_to_survival[3], prior_etu_to_survival[4]);
  
  etu_to_death_logmean ~ normal(prior_etu_to_death[1], prior_etu_to_death[2]);
  etu_to_death_sd ~ normal(prior_etu_to_death[3], prior_etu_to_death[4]);
  
  onset_to_iso_logmean ~ normal(prior_onset_to_iso[1], prior_onset_to_iso[2]);
  onset_to_iso_sd ~ normal(prior_onset_to_iso[3], prior_onset_to_iso[4]);

  iso_to_release_logmean ~ normal(prior_iso_to_release[1], prior_iso_to_release[2]);
  iso_to_release_sd ~ normal(prior_iso_to_release[3], prior_iso_to_release[4]);  

  // prior on cfr
  cfr ~ normal(prior_cfr[1], prior_cfr[2]);

  // prior on prop isolated
  prop_iso_logit ~ normal(prior_prop_iso[1], prior_prop_iso[2]);

  // prior on number of background alerts
  alerts_background ~ normal(prior_alerts_background[1], prior_alerts_background[2]);

  // prior on number of alerts per case
  alerts_per_case ~ normal(prior_alerts_per_case[1], prior_alerts_per_case[2]);

  // prior on overdispersion parameter of nbinom sampling distribution
  cases_overdisp_log ~ normal(2, 1);
  etu_overdisp_log ~ normal(0, 1);
  alerts_overdisp_log ~ normal(0, 1);
  iso_overdisp_log ~ normal(0, 1);

  // prior on spline parameters
  spline_param ~ normal(0, 1);

  // likelihood of reported cases
  cases_reported ~ neg_binomial_2(cases_truncated_mu, cases_overdisp);

  // likelihood of reported deaths
  deaths_reported ~ neg_binomial_2(deaths_truncated_mu, deaths_overdisp);

  // likelihood of reported etu
  for(i in etu_ind)
    target += neg_binomial_2_lpmf(
      etu_reported[i] | etu_truncated_mu[i], etu_overdisp
    );

  // likelihood of reported number of alerts
  for(i in alerts_ind)
    target += neg_binomial_2_lpmf(
      alerts_reported[i] | alerts_truncated_mu[i], alerts_overdisp
    );

  // likelihood of reported isolation occupancy
  for(i in iso_ind)
    target += neg_binomial_2_lpmf(
      iso_reported[i] | iso_truncated_mu[i], iso_overdisp
    );

}

generated quantities {

  // SIMULATE CASES AND DEATHS

  // GLOBAL VARIABLES

  int cases_projected_sim[n_proj];
  int cases_truncated_sim[n_obs];
  int cases_nowcast_sim[n_obs];
  
  int deaths_projected_sim[n_proj];
  int deaths_truncated_sim[n_obs];
  int deaths_nowcast_sim[n_obs];

  // generate sample of fitted cases
  cases_truncated_sim = neg_binomial_2_rng(cases_truncated_mu, cases_overdisp);

  // generate sample of fitted cases
  cases_nowcast_sim = neg_binomial_2_rng(cases_nowcast_mu, cases_overdisp);

  // generate sample of projected cases
  cases_projected_sim = neg_binomial_2_rng(cases_projected_mu, cases_overdisp);

  // generate sample of fitted deaths
  deaths_truncated_sim = neg_binomial_2_rng(deaths_truncated_mu, deaths_overdisp);

  // generate sample of fitted deaths
  deaths_nowcast_sim = neg_binomial_2_rng(deaths_nowcast_mu, deaths_overdisp);

  // generate sample of projected deaths
  deaths_projected_sim = neg_binomial_2_rng(deaths_projected_mu, deaths_overdisp);


  // NOWCAST ALERTS AND HOSPITAL OCCUPANCY

  // GLOBAL VARIABLES

  // modelled means
  vector[n_tot] iso_nowcast_mu;
  vector[n_tot] alerts_nowcast_mu;
  vector[n_tot] etu_nowcast_mu;

  {

  // densities for etu admission and discharge
  vector[n_tot] etu_admission_nowcast_mu;
  vector[n_tot] etu_discharge_nowcast_mu;

  // densities for isolation admission and discharge
  vector[n_tot] iso_admission_nowcast_mu;
  vector[n_tot] iso_discharge_nowcast_mu;

  // define empty vectors for nowcast etu and iso densities
  alerts_nowcast_mu = rep_vector(0, n_tot);  
  etu_admission_nowcast_mu = rep_vector(0, n_tot);
  etu_discharge_nowcast_mu = rep_vector(0, n_tot);
  iso_admission_nowcast_mu = rep_vector(0, n_tot);
  iso_discharge_nowcast_mu = rep_vector(0, n_tot);

  // calculate nowcast etu and alerts
  for (i in 1:(n_tot - 1)) {
    int j_max = (i + max_delay > n_tot) ? (n_tot - i) : max_delay;
  
    // pick the case source (observed nowcast vs projected)
    real cases_val;
    real alerts_bg;
    if (i <= n_obs) {
      cases_val = cases_nowcast_mu[i];
      alerts_bg = exp(alerts_background[alerts_background_ind[i]]);
    } else {
      cases_val = cases_projected_mu[i - n_obs];
      alerts_bg = exp(alerts_background[alerts_background_ind[n_alerts_background]]);
    }

    // precompute alerts term
    real alerts_val = cases_val * exp(alerts_per_case) + alerts_bg;

    for (j in 1:j_max) {
      int t = i + j;

      etu_admission_nowcast_mu[t] += cases_val * onset_to_etu[j];
      // discharge comes from survival and deaths
      etu_discharge_nowcast_mu[t] += cases_val *
        (onset_to_survival[j] + onset_to_death[j]);

      alerts_nowcast_mu[t] += alerts_val * onset_to_iso[j];
      iso_admission_nowcast_mu[t] += alerts_val * onset_to_iso[j] * prop_iso;
      iso_discharge_nowcast_mu[t] += alerts_val * onset_to_release[j] * prop_iso;
      
    }
  }

  // calculate etu and iso occupancy
  etu_nowcast_mu = cumulative_sum(
    etu_admission_nowcast_mu - etu_discharge_nowcast_mu
  );
  iso_nowcast_mu = cumulative_sum(
    iso_admission_nowcast_mu - iso_discharge_nowcast_mu
  );

  // correct small values
  for(i in 1:(n_tot))
    if(etu_nowcast_mu[i] == 0) etu_nowcast_mu[i] = 0.00001;
  for(i in 1:(n_tot)) if(alerts_nowcast_mu[i] == 0)
    alerts_nowcast_mu[i] = 0.00001;
  for(i in 1:(n_tot)) if(iso_nowcast_mu[i] <= 0)
    iso_nowcast_mu[i] = 0.00001;

  }


  // SIMULATE ALERTS AND HOSPITAL OCCUPANCY

  // GLOBAL VARIABLES

  // simulate truncated
  int etu_truncated_sim[n_tot];
  int alerts_truncated_sim[n_tot];      
  int iso_truncated_sim[n_tot];

  // simulate nowcast
  int etu_nowcast_sim[n_tot];
  int alerts_nowcast_sim[n_tot];      
  int iso_nowcast_sim[n_tot];

  etu_truncated_sim = neg_binomial_2_rng(etu_truncated_mu, etu_overdisp);
  alerts_truncated_sim = neg_binomial_2_rng(alerts_truncated_mu, alerts_overdisp);
  iso_truncated_sim = neg_binomial_2_rng(iso_truncated_mu, iso_overdisp);

  etu_nowcast_sim = neg_binomial_2_rng(etu_nowcast_mu, etu_overdisp);
  alerts_nowcast_sim = neg_binomial_2_rng(alerts_nowcast_mu, alerts_overdisp);
  iso_nowcast_sim = neg_binomial_2_rng(iso_nowcast_mu, iso_overdisp);


}
