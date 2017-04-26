library(magrittr)
library(purrr)
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)

# Generate range of parameters for beta-binomial-hypergeometric distribution with N = 50 and n at most 25
sample_params <- data_frame(x = 0:25)
sample_params <- mutate(sample_params, n = 25)
sample_params <- map_df(.x = 0:50, .f = ~ mutate(sample_params, M = .x) %>%
                          mutate(x = ifelse(x > .x, .x, x)))
sample_params <- distinct(sample_params)
sample_params <- mutate(sample_params, N = 50)
sample_params <- map_df(.x = 1:50, .f = ~ mutate(sample_params, alpha = .x))
sample_params <- map_df(.x = 1:50, .f = ~ mutate(sample_params, beta = .x))
sample_params <- distinct(sample_params)




memory.limit(size = 25000)
beta_binomial_density <- function(x, n, alpha, beta) {choose(n, x)*beta(alpha+x,beta+n-x)/beta(alpha,beta)}

params_with_posterior_probs <- sample_params %>%
  mutate(Likelihood_of_x_for_MN = dhyper(x, M, N - M, n),
         Likelihood_of_MN_for_Alpha_Beta = beta_binomial_density(M, N, alpha, beta),
         Likelihood_of_x_for_MN_Alpha_Beta = Likelihood_of_x_for_MN*Likelihood_of_MN_for_Alpha_Beta) %>%
  group_by(N, n, x) %>%
  mutate(PP_Alpha_Beta_MN_given_xnN = Likelihood_of_x_for_MN_Alpha_Beta/sum(Likelihood_of_x_for_MN_Alpha_Beta)) %>%
  ungroup %>%
  group_by(alpha, beta, x, n, N) %>%
  summarize(PP_Alpha_Beta_given_xnN = sum(PP_Alpha_Beta_MN_given_xnN)) %>%
  ungroup %>%
  group_by(N, n, x) %>%
  mutate(PP_Alpha_Beta_given_xnN = PP_Alpha_Beta_given_xnN/sum(PP_Alpha_Beta_given_xnN))

likelihood_of_x_for_alpha_beta <- sample_params %>%
  mutate(Likelihood_of_x_for_MN = dhyper(x, M, N - M, n),
         Likelihood_of_MN_for_Alpha_Beta = beta_binomial_density(M, N, alpha, beta),
         Likelihood_of_x_for_MN_Alpha_Beta = Likelihood_of_x_for_MN*Likelihood_of_MN_for_Alpha_Beta) %>%
  group_by(N, n, x) %>%
  

isoquant_sets <- params_with_posterior_probs %>%
  group_by(N, n, x, PP_Alpha_Beta_given_xnN) %>%
  filter(n() > 1) %>%
  mutate(isoquant = group_indices(1))

pp_plot <- params_with_posterior_probs %>%
  ggplot(aes(x = alpha, y = beta, 
             alpha = PP_Alpha_Beta_given_xnN,
             frame = x)) +
  geom_raster(fill = "dodgerblue4") +
  coord_equal()

gg_animate(pp_plot)

# Create a population composed of equal-sized sub-populations
  population <- data_frame(x = c(rbinom(500, 1, .1), rbinom(500, 1, .95)),
                           subpopulation = map(1:20, ~rep(.x, 50)) %>% unlist)
  population_summaries <- population %>%
    group_by(subpopulation) %>%
    summarize(M = sum(x))

# Sample equally from each subpopulation
  sample <- population %>%
    group_by(subpopulation) %>%
    sample_n(size = 25, replace = FALSE)
# Sumamrize the samples from each subpopulation
  sample_summaries <- sample %>%
    group_by(subpopulation) %>%
    summarize(x = sum(x)) %>%
    mutate(n = 25, N = 50)


# Obtain likelihoods of alpha-and-beta parameters for a beta distribution
# in each sub-sample, combine into a single likelihood for the aggregation of all but one sub-population's sample.
# Then convert this to posterior probabilities of alpha and beta given a uniform prior on their joint distribution.
  alpha_beta_df <- map2_df(.x = map(1:1000, ~data_frame(alpha = 1:1000)), .y = 1:1000,
                           .f = ~ mutate(.x, beta = .y))
  alpha_beta_df <-  map_df(.x = 1:20,
                           .f = ~ mutate(alpha_beta_df, subpopulation = .x))
  
  sample_summaries_w_alpha_beta_densities <- sample_summaries %>%
    group_by(subpopulation) %>%
    left_join(x = .,
              y = alpha_beta_df) %>%
    mutate(alpha_beta_density = beta_binomial_density(x, n, alpha, beta))

# For a given sub-population, use the combined sub-samples from the other populations
# to develop a prior on the combination of alpha-and-beta (i.e. prior to incorporating the sample for that sub-population).
# Use this prior on alpha-and-beta to generate a prior for M, the number of successes in the finite sub-population.

  # Demonstration of generating a prior and posterior for a single sub-population
  alpha_beta_prior_for_subpop_1 <- sample_summaries_w_alpha_beta_densities %>%
    ungroup %>%
    filter(subpopulation != 1) %>%
    group_by(alpha, beta) %>%
    summarize(alpha_beta_density = prod(alpha_beta_density)) %>%
    ungroup %>%
    mutate(alpha_beta_pp = 1000*alpha_beta_density/sum(1000*alpha_beta_density, na.rm = TRUE),
           alpha_beta_pp = alpha_beta_pp/sum(alpha_beta_pp, na.rm = TRUE))
  


  M_prior_for_subpop_1 <- map_df(.x = 1:50,
                                 .f = ~ mutate(alpha_beta_prior_for_subpop_1,
                                               M = .x)) %>%
    mutate(M_Likelihood = beta_binomial_density(x = M,
                                                n = 50,
                                                alpha = alpha, beta = beta),
           Probability_of_M = alpha_beta_pp*M_Likelihood) %>%
    group_by(M) %>%
    summarize(Probability_of_M = sum(Probability_of_M, na.rm = TRUE)) %>%
    mutate(Probability_of_M = Probability_of_M/sum(Probability_of_M, na.rm = TRUE)) %>%
    rename(Prior_Probability_of_M = Probability_of_M)

  subpop_1_posterior <- M_prior_for_subpop_1 %>%
    mutate(Sample_Summary = map(.x = M,
                                .f = ~ filter(sample_summaries,
                                              subpopulation == 1))) %>%
    unnest %>%
    mutate(Likelihood_of_M = dhyper(x = x, m = M, n = N - M, k = n),
           Posterior_of_M = Prior_Probability_of_M*Likelihood_of_M,
           Posterior_of_M = Posterior_of_M/sum(Posterior_of_M))
  
  subpop_1_posterior %>%
    gather(key = Component, value = Value, -M) %>%
    filter(Component %in% c("Prior_Probability_of_M",
                            "Likelihood_of_M",
                            "Posterior_of_M")) %>%
    mutate(Component = factor(Component, levels = c("Prior_Probability_of_M",
                                                    "Likelihood_of_M",
                                                    "Posterior_of_M"), ordered = TRUE)) %>%
    ggplot(aes(x = M, y = Value, fill = Component)) +
    geom_col() +
    facet_wrap(~ Component)
  
  # Generate posteriors for each subpopulation
  
    # Define function to produce posterior for a given sub-population
    get_subpop_posterior <- function(subpop_id) {
      alpha_beta_prior_for_subpop <- sample_summaries_w_alpha_beta_densities %>%
        ungroup %>%
        filter(subpopulation != subpop_id) %>%
        group_by(alpha, beta) %>%
        summarize(alpha_beta_density = prod(alpha_beta_density)) %>%
        ungroup %>%
        mutate(alpha_beta_pp = alpha_beta_density/sum(alpha_beta_density, na.rm = TRUE))
      
      M_prior_for_subpop <- map_df(.x = 1:50,
                                     .f = ~ mutate(alpha_beta_prior_for_subpop,
                                                   M = .x)) %>%
        mutate(M_Likelihood = beta_binomial_density(x = M,
                                                    n = 50,
                                                    alpha = alpha, beta = beta),
               Probability_of_M = alpha_beta_pp*M_Likelihood) %>%
        group_by(M) %>%
        summarize(Probability_of_M = sum(Probability_of_M, na.rm = TRUE)) %>%
        mutate(Probability_of_M = Probability_of_M/sum(Probability_of_M, na.rm = TRUE)) %>%
        rename(Prior_Probability_of_M = Probability_of_M)
      
      M_prior_for_subpop %>%
        mutate(Sample_Summary = map(.x = M,
                                    .f = ~ filter(sample_summaries,
                                                  subpopulation == subpop_id))) %>%
        unnest %>%
        mutate(Likelihood_of_M = dhyper(x = x, m = M, n = N - M, k = n),
               Posterior_of_M = Prior_Probability_of_M*Likelihood_of_M,
               Posterior_of_M = Posterior_of_M/sum(Posterior_of_M))
}
    
    # Apply the function to each subpopulation
    subpopulation_posteriors <- vector(mode = "list", length = 20)
    for(subpopulation in 1:20) {
      subpopulation_posteriors[[subpopulation]] <- get_subpop_posterior(subpopulation)
    }
    subpopulation_posteriors <- bind_rows(subpopulation_posteriors)
    

    
    # Obtain MAP and ML estimates for each sub-population
    # place alongside the true population values
    Subpopulation_Estimates <- subpopulation_posteriors %>%
      group_by(subpopulation) %>%
      summarize(ML_Estimate = unique(M[Likelihood_of_M == max(Likelihood_of_M)]),
                MAP_Estimate = unique(M[Posterior_of_M == max(Posterior_of_M)])) %>%
      left_join(x = .,
                y = population_summaries,
                by = "subpopulation")
    
    # Visualize the priors and posteriors for each sub-population estimate
    subpopulation_posteriors %>%
      group_by(subpopulation) %>%
      do({gather(., key = Component, value = Value, -M)}) %>%
      filter(Component %in% c("Prior_Probability_of_M",
                              "Likelihood_of_M",
                              "Posterior_of_M")) %>%
      mutate(Component = factor(Component, levels = c("Prior_Probability_of_M",
                                                      "Likelihood_of_M",
                                                      "Posterior_of_M"), ordered = TRUE)) %>%
      ggplot(aes(x = M, y = Value, color = Component, fill = Component)) +
      geom_area(alpha = 0.5) +
      geom_vline(data = Subpopulation_Estimates,
                 aes(xintercept = ML_Estimate), color = "green") +
      geom_vline(data = Subpopulation_Estimates,
                 aes(xintercept = MAP_Estimate), color = "blue") +
      geom_vline(data = Subpopulation_Estimates,
                 aes(xintercept = M), color = "black") +
      facet_wrap(~ subpopulation)

    subpopulation_posteriors %>%
      group_by(subpopulation) %>%
      summarize(ML_Estimate = mean(N*(mean(x)/mean(n))),
                MAP_Estimate = unique(M[Posterior_of_M == max(Posterior_of_M)])) %>%
      summarize(ML_Average_Error = mean(abs(ML_Estimate - M)),
                MAP_Average_Error = mean(abs(MAP_Estimate - M)))

