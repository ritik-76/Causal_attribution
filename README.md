# Causal_attribution

This repository contains the code and data underlying the article :

 ## Estimating Causal Attribution of Anthropogenic Forcing on High-Temperature Extremes Usinga Latent Gaussian Spatial Model 

### 📂 The description of all attached files is provided below :

#### Likelihood_estimation.R

Generalized Likelihood estimates of the latent variables.

#### Model_components.R

For creating the necessary 

#### Update_MCMC.R

This file contains the function for the updation steps in the Gibbs Sampling Method.

#### MCMC_Out.R

For generating posterior samples corresponding to the parameters in the latent layer and hyperparameters.

#### Table_1_linear_model.R

Generating the entries for the table 1 in the article.

#### Figure_1_corrmat_EDA.R

Plot of the heatmap of the correlation matrix in Figure 1 of the artice.

#### Figure_1_variogram_EDA.R

Plot of variogram corresponding to the latent variables in Figure 1 of the article.

#### Figure_2_trans_xi.R

Plot of the characterisation of transformed shape parameter.

#### Table_2_gamma_posterior.R

Posterior summaries of gamma in Table 2 of the article.

#### Table_3_Posterior_diag.R

Diagnostics of approximate posterior samples of gamma in Table 3 of the article. 

#### Figure_3_causal_effect.R

Plot of grid wise estimated causal effect along with the standard errors in Figure 3 of the article.

#### Figure_4_causal_effect.R

Plot of grid wise estimated causal effect along with the standard errors during the pre-industrial period in Figure 4 of the article.

#### Figure_5_trend_plot.R

Plot of trend and trend difference in Figure 5 of the article.

#### Figure_6_Hotspot.R

Plot of the hotspot regions in Figure 6 of the article.


###  ⚙️ Requirements

#### Install required packages :

install.packages(c("extRemes", "fields", "spam", "Matrix", "LaplacesDemon", "ggplot", "cowplot", "ggmap", "ExceedanceTools", "coda"))

#### Instruction for installing evd package :

- Download the source package.
- Replace the provided bvfit.R with the source code of evd package.
- Install the modified evd package from the assigned local path.
