functions {

  // reverse function used in convolution
  vector stan_rev(vector a) {
    int a_size = num_elements(a);
    vector [a_size] a_rev;
    
    for(i in 1:a_size)  
      a_rev[i] = a[a_size - i + 1];
   
    return a_rev;
    
  }

  // convolution function from https://discourse.mc-stan.org/t/dot-products-of-vectors-to-perform-1d-convolution/9053/9 
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

  // Number of days with data
  int<lower = 1> n_days;

  // Max of delay distribution
  int<lower = 1> max_delay;

  // Number of days of projection
  int<lower = 1> days_ahead;

  // Number of reported cases per day
  int cases_observed[n_days];

  // ETU occupancy by day
  int etu_observed[n_days];

  // Number of days with etu occupancy data
  int etu_n;

  // Indices of days with etu occupancy data
  int etu_ind[etu_n];

  // New community alerts by day
  int alerts_observed[n_days];

  // Number of days with alerts data
  int alerts_n;

  // Indices of days with alerts data
  int alerts_ind[alerts_n];

  // Isolation bed occupancy by day
  int iso_observed[n_days];

  // Number of days with iso occupancy data
  int iso_n;

  // Indices of days with iso occupancy data
  int iso_ind[iso_n];

  // mean and sd of the normal prior on the logmean of the delay from onset to hospitalisation
  vector[4] prior_onset_to_etu;

  // mean and sd of the normal prior on the logmean of the delay from hospitalisation to survival
  vector[4] prior_etu_to_survival;

  // mean and sd of the normal prior on the logmean of the delay from hospitalisation to death
  vector[4] prior_etu_to_death;

  // mean and sd of the normal prior on the logmean of the delay from case-onset to contact-isolation
  vector[4] prior_onset_to_iso;

  // mean and sd of the normal prior on the logmean of the delay from in lab turnaround
  vector[4] prior_iso_to_release;

  // mean and sd of the normal prior on the plogis transformed CFR
  vector[2] prior_cfr;
  
  // mean and sd of the normal prior on the plogis transformed proportion of alerts isolated
  vector[2] prior_prop_iso;
  
  // mean and sd of the normal prior on the exponentially transformed number of background alerts
  vector[2] prior_alerts_background;

  // mean and sd of the normal prior on the exponentially transformed number of alerts per case
  vector[2] prior_alerts_per_case;

  // delay distribution from onset to confirmation and entry into database
  vector[n_days] log_prop_cases_reported;

  // number of alerts_background bins
  int n_alerts_background;
  
  // alerts_background indices
  int alerts_background_ind[n_days];

  // Number of spline parameters
  int n_spline_param;

  // spline matrix for growth rate fitting
  matrix[n_spline_param, n_days] spline;

}

parameters {

  // in-hospital CFR
  real cfr;

  // proportion of alerts isolated
  real prop_iso;

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

  // overdispersion parameter for nbinom sampling of cases, etu occupancy, alerts and isolation bed occupancy
  real cases_overdisp;  
  real etu_overdisp;
  real alerts_overdisp;  
  real iso_overdisp;  

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

  // // cases per day that are missed due to delayed reporting
  // vector[n_days] log_cases_missed;

}

transformed parameters {

  // vector of fitted and projected cases
  vector[n_days] log_cases_fitted;
  vector[n_days] log_cases_inflated;  
  vector[days_ahead] log_cases_projected;

  // raw delay densities for case-to-etu model
  vector[max_delay] onset_to_etu;
  vector[max_delay] etu_to_survival;
  vector[max_delay] etu_to_death;

  // convolved densities
  vector[max_delay] onset_to_survival;
  vector[max_delay] onset_to_death;

  // densities for etu admission, discharge, change in occupancy
  vector[n_days + days_ahead] etu_admission;
  vector[n_days + days_ahead] etu_discharge;
  vector[n_days + days_ahead] etu_modelled;
  vector[n_days] etu_modelled_for_fit;

  // raw delay densities for case-to-isolation
  vector[max_delay] onset_to_iso;  
  vector[max_delay] iso_to_release;

  // convolved density
  vector[max_delay] onset_to_release;

  // densities for alerts, isolation admission, discharge, change in occupancy
  vector[n_days + days_ahead] alerts_modelled;
  vector[n_days + days_ahead] iso_admission;
  vector[n_days + days_ahead] iso_discharge;
  vector[n_days + days_ahead] iso_modelled;
  vector[n_days] alerts_modelled_for_fit;  
  vector[n_days] iso_modelled_for_fit;

  // plogis transform of proportion of alerts isolated
  real prop_iso_trans;

  // vector of estimated growth rates
  vector[n_days] growth_rate_vec;

  // plogis transform of proportion of alerts isolated
  prop_iso_trans = logistic_cdf(prop_iso, 0, 1);

  // calculate growth rate from spline params
  growth_rate_vec = to_vector(spline_param*spline)/10;

  // model case observations with reporting delay for first day
  log_cases_fitted[1] = log_cases_intercept + growth_rate_vec[1] + log_prop_cases_reported[1];
  
  // model case observations without reporting delay for first day
  log_cases_inflated[1] = log_cases_intercept + growth_rate_vec[1];

  // model case observations for all remaining observed days
  if(log_cases_fitted[1] > 10) log_cases_fitted[1] = 10;
  if(log_cases_inflated[1] > 10) log_cases_inflated[1] = 10;  
  for(i in 2:n_days) {
    log_cases_fitted[i] = log_cases_inflated[i-1] + growth_rate_vec[i] + log_prop_cases_reported[i];
    log_cases_inflated[i] = log_cases_inflated[i-1] + growth_rate_vec[i];
    if(log_cases_fitted[i] > 10) log_cases_fitted[i] = 10;
    if(log_cases_inflated[i] > 10) log_cases_inflated[i] = 10;    
  }

  // projected cases forward using inflated values
  for(i in 1:days_ahead) {
    log_cases_projected[i] = log_cases_inflated[n_days] + growth_rate_vec[n_days]*i;
    if(log_cases_projected[i] > 10) log_cases_projected[i] = 10;
    if(log_cases_projected[i] < -10) log_cases_projected[i] = -10;
  }

  alerts_modelled = rep_vector(0, n_days + days_ahead);
  etu_admission = rep_vector(0, n_days + days_ahead);
  etu_discharge = rep_vector(0, n_days + days_ahead);
  iso_admission = rep_vector(0, n_days + days_ahead);
  iso_discharge = rep_vector(0, n_days + days_ahead);

  // calculate delay pmf matrix and normalise
  for(i in 1:max_delay) {
    // calculate pmf from lognormal 
    onset_to_etu[i] = exp(lognormal_lpdf(i | onset_to_etu_logmean, onset_to_etu_sd));
    etu_to_survival[i] = exp(lognormal_lpdf(i | etu_to_survival_logmean, etu_to_survival_sd));
    etu_to_death[i] = exp(lognormal_lpdf(i | etu_to_death_logmean, etu_to_death_sd));
    onset_to_iso[i] = exp(lognormal_lpdf(i | onset_to_iso_logmean, onset_to_iso_sd));
    iso_to_release[i] = exp(lognormal_lpdf(i | iso_to_release_logmean, iso_to_release_sd));    
  }

  // normalise pdf
  onset_to_etu = onset_to_etu ./ sum(onset_to_etu);
  etu_to_survival = etu_to_survival ./ sum(etu_to_survival);
  etu_to_death = etu_to_death ./ sum(etu_to_death);
  onset_to_iso = onset_to_iso ./ sum(onset_to_iso);
  iso_to_release = iso_to_release ./ sum(iso_to_release);

  // convolve onset to survival and include CFR
  onset_to_survival = stan_convolve(max_delay, onset_to_etu, etu_to_survival);
  onset_to_survival = onset_to_survival ./ sum(onset_to_survival) * (1 - logistic_cdf(cfr, 0, 1));

  // convolve onset to death and include CFR
  onset_to_death = stan_convolve(max_delay, onset_to_etu, etu_to_death);
  onset_to_death = onset_to_death ./ sum(onset_to_death) * logistic_cdf(cfr, 0, 1);

  // convolve onset to release-from-isolation
  onset_to_release = stan_convolve(max_delay, onset_to_iso, iso_to_release);
  onset_to_release = onset_to_release ./ sum(onset_to_release);

  // calculate admission and discharge
  for(i in 1:(n_days + days_ahead - 1)) {
    if((i + max_delay) > (n_days + days_ahead)) {
      // in case of right truncation move until end of analysis period
      for(j in 1:(n_days + days_ahead - i)) {
        // pull from observed cases
        if(i <= n_days) {
          etu_admission[i + j] +=
	    cases_observed[i]*onset_to_etu[j];
          etu_discharge[i + j] +=
	    cases_observed[i]*(onset_to_survival[j] + onset_to_death[j]);
	  // here we are adding alerts from cases and alerts from background together
          alerts_modelled[i + j] +=
	    (cases_observed[i]*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j];
	  // iso admission density is alert density multiplied by prop isolated
          iso_admission[i + j] +=
	    (cases_observed[i]*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j]*prop_iso_trans;
	  // iso discharge density is discharge density multiplied by prop isolated	  
          iso_discharge[i + j] +=
	    (cases_observed[i]*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_release[j]*prop_iso_trans;
	// pull from projected cases
	} else {
          etu_admission[i + j] +=
	    exp(log_cases_projected[i - n_days])*onset_to_etu[j];
          etu_discharge[i + j] +=
	    exp(log_cases_projected[i - n_days])*(onset_to_survival[j] + onset_to_death[j]);
          alerts_modelled[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j];
          iso_admission[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_release[j]*prop_iso_trans;
        }
      }
    } else {
      // if no right truncation go to end of delay distribution
      for(j in 1:max_delay) {
        // pull from observed cases
        if(i <= n_days) {
          etu_admission[i + j] +=
	    cases_observed[i]*onset_to_etu[j];
          etu_discharge[i + j] +=
	    cases_observed[i]*(onset_to_survival[j] + onset_to_death[j]);
          alerts_modelled[i + j] +=
	    (cases_observed[i]*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j];
          iso_admission[i + j] +=
	    (cases_observed[i]*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge[i + j] +=
	    (cases_observed[i]*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_release[j]*prop_iso_trans;
	// pull from projected cases
	} else {
          etu_admission[i + j] +=
	    exp(log_cases_projected[i - n_days])*onset_to_etu[j];
          etu_discharge[i + j] +=
	    exp(log_cases_projected[i - n_days])*(onset_to_survival[j] + onset_to_death[j]);
          alerts_modelled[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j];
          iso_admission[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_release[j]*prop_iso_trans;
	}
      }
    }
  }

  // calculate etu and iso occupancy
  etu_modelled = cumulative_sum(etu_admission - etu_discharge);
  iso_modelled = cumulative_sum(iso_admission - iso_discharge);

  for(i in 1:(n_days + days_ahead)) if(etu_modelled[i] == 0) etu_modelled[i] = 0.00001;
  for(i in 1:n_days) etu_modelled_for_fit[i] = etu_modelled[i];

  for(i in 1:(n_days + days_ahead)) if(alerts_modelled[i] == 0) alerts_modelled[i] = 0.00001;
  for(i in 1:n_days) alerts_modelled_for_fit[i] = alerts_modelled[i];

  for(i in 1:(n_days + days_ahead)) if(iso_modelled[i] <= 0) iso_modelled[i] = 0.00001;
  for(i in 1:n_days) iso_modelled_for_fit[i] = iso_modelled[i];

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
  prop_iso ~ normal(prior_prop_iso[1], prior_prop_iso[2]);

  // prior on number of background alerts
  alerts_background ~ normal(prior_alerts_background[1], prior_alerts_background[2]);

  // prior on number of alerts per case
  alerts_per_case ~ normal(prior_alerts_per_case[1], prior_alerts_per_case[2]);

  // prior on overdispersion parameter of nbinom sampling distribution
  cases_overdisp ~ normal(2, 1);
  etu_overdisp ~ normal(0, 1);
  alerts_overdisp ~ normal(0, 1);
  iso_overdisp ~ normal(0, 1);

  // prior on spline parameters
  spline_param ~ normal(0, 1);

  // likelihood of observed cases
  cases_observed ~ neg_binomial_2(exp(log_cases_fitted), exp(cases_overdisp));

  // likelihood of observed etu
  for(i in etu_ind)
    target += neg_binomial_2_lpmf(etu_observed[i] | etu_modelled_for_fit[i], exp(etu_overdisp));

  // likelihood of observed number of alerts
  for(i in alerts_ind)
    target += neg_binomial_2_lpmf(alerts_observed[i] | alerts_modelled_for_fit[i], exp(alerts_overdisp));

  // likelihood of observed isolation occupancy
  for(i in iso_ind)
    target += neg_binomial_2_lpmf(iso_observed[i] | iso_modelled_for_fit[i], exp(iso_overdisp));

}

generated quantities {

  int cases_projected_dist[days_ahead];
  int cases_fitted_dist[n_days];
  int cases_inflated_dist[n_days];    
  int etu_modelled_dist[n_days + days_ahead];
  int alerts_modelled_dist[n_days + days_ahead];      
  int iso_modelled_dist[n_days + days_ahead];

  // densities for etu admission, discharge, change in occupancy
  vector[n_days + days_ahead] etu_admission_inflated;
  vector[n_days + days_ahead] etu_discharge_inflated;
  vector[n_days + days_ahead] etu_modelled_inflated;

  // densities for isolation admission, discharge, change in occupancy
  vector[n_days + days_ahead] alerts_modelled_inflated;  
  vector[n_days + days_ahead] iso_admission_inflated;
  vector[n_days + days_ahead] iso_discharge_inflated;
  vector[n_days + days_ahead] iso_modelled_inflated;

  // generate sample of fitted cases
  cases_fitted_dist = neg_binomial_2_rng(exp(log_cases_fitted), exp(cases_overdisp));

  // generate sample of fitted cases
  cases_inflated_dist = neg_binomial_2_rng(exp(log_cases_inflated), exp(cases_overdisp));

  // generate sample of projected cases
  cases_projected_dist = neg_binomial_2_rng(exp(log_cases_projected), exp(cases_overdisp));

  // generate sample of etu occupancy
  etu_modelled_dist = neg_binomial_2_rng(etu_modelled, exp(etu_overdisp));

  // generate sample of alerts
  alerts_modelled_dist = neg_binomial_2_rng(alerts_modelled, exp(alerts_overdisp));
  
  // generate sample of isolation occupancy
  iso_modelled_dist = neg_binomial_2_rng(iso_modelled, exp(iso_overdisp));

  // define empty vectors for inflated etu and iso densities
  alerts_modelled_inflated = rep_vector(0, n_days + days_ahead);  
  etu_admission_inflated = rep_vector(0, n_days + days_ahead);
  etu_discharge_inflated = rep_vector(0, n_days + days_ahead);
  iso_admission_inflated = rep_vector(0, n_days + days_ahead);
  iso_discharge_inflated = rep_vector(0, n_days + days_ahead);

  // calculate etu_admission_inflated and etu_discharge_inflated density
  for(i in 1:(n_days + days_ahead - 1)) {
    if((i + max_delay) > (n_days + days_ahead)) {
      // in case of right truncation move until end of analysis period
      for(j in 1:(n_days + days_ahead - i)) {
        // pull from observed cases
        if(i <= n_days) {
          etu_admission_inflated[i + j] +=
	    exp(log_cases_inflated[i])*onset_to_etu[j];
          etu_discharge_inflated[i + j] +=
	    exp(log_cases_inflated[i])*(onset_to_survival[j] + onset_to_death[j]);
	  // here we are adding alerts from cases and alerts from background together
          alerts_modelled_inflated[i + j] +=
	    (exp(log_cases_inflated[i])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j];
          iso_admission_inflated[i + j] +=
	    (exp(log_cases_inflated[i])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge_inflated[i + j] +=
	    (exp(log_cases_inflated[i])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_release[j]*prop_iso_trans;
	// pull from projected cases
	} else {
          etu_admission_inflated[i + j] +=
	    exp(log_cases_projected[i - n_days])*onset_to_etu[j];
          etu_discharge_inflated[i + j] +=
	    exp(log_cases_projected[i - n_days])*(onset_to_survival[j] + onset_to_death[j]);
          alerts_modelled_inflated[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j];
          iso_admission_inflated[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge_inflated[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_release[j]*prop_iso_trans;
        }
      }
    } else {
      // if no right truncation go to end of delay distribution
      for(j in 1:max_delay) {
        // pull from observed cases
        if(i <= n_days) {
          etu_admission_inflated[i + j] +=
	    exp(log_cases_inflated[i])*onset_to_etu[j];
          etu_discharge_inflated[i + j] +=
	    exp(log_cases_inflated[i])*(onset_to_survival[j] + onset_to_death[j]);
          alerts_modelled_inflated[i + j] +=
	    (exp(log_cases_inflated[i])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j];
          iso_admission_inflated[i + j] +=
	    (exp(log_cases_inflated[i])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge_inflated[i + j] +=
	    (exp(log_cases_inflated[i])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[i]]))*onset_to_release[j]*prop_iso_trans;
	// pull from projected cases
	} else {
          etu_admission_inflated[i + j] +=
	    exp(log_cases_projected[i - n_days])*onset_to_etu[j];
          etu_discharge_inflated[i + j] +=
	    exp(log_cases_projected[i - n_days])*(onset_to_survival[j] + onset_to_death[j]);
          alerts_modelled_inflated[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j];
          iso_admission_inflated[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_iso[j]*prop_iso_trans;
          iso_discharge_inflated[i + j] +=
	    (exp(log_cases_projected[i - n_days])*exp(alerts_per_case) +
	    exp(alerts_background[alerts_background_ind[n_alerts_background]]))*onset_to_release[j]*prop_iso_trans;
	}
      }
    }
  }

  // calculate etu and iso occupancy
  etu_modelled_inflated = cumulative_sum(etu_admission_inflated - etu_discharge_inflated);
  iso_modelled_inflated = cumulative_sum(iso_admission_inflated - iso_discharge_inflated);

  for(i in 1:(n_days + days_ahead)) if(etu_modelled_inflated[i] == 0) etu_modelled_inflated[i] = 0.00001;
  for(i in 1:(n_days + days_ahead)) if(alerts_modelled_inflated[i] == 0) alerts_modelled_inflated[i] = 0.00001;
  for(i in 1:(n_days + days_ahead)) if(iso_modelled_inflated[i] <= 0) iso_modelled_inflated[i] = 0.00001;

}
