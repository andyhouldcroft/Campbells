# Campbells

This repository contains code to perform the wildlife density modelling approach outlined in <a href="https://doi.org/10.1111/ecog.07399">Houldcroft et al. (2024)</a>, with the associated survey data available from the Dryad Digital Repository: https://doi.org/10.5061/dryad.xd2547dqx

[H24_SETUP.R](H24_SETUP.R) contains the code to produce spatial covariates, a barrier SPDE mesh and format all data as sf and terra objects compatible with inlabru.

[H24_DENSITIES.R](H24_DENSITIES.R) contains the code to run the point process model with a spatial cluster size distribution, cluster size detection function and barrier SPDE models.

Please ensure that the demonstrated hazard-rate detection function and discretised truncated log-Normal cluster size distribution are appropriate for your specific data and research circumstances. In the event of model convergence issues, consider including additional priors on the precision and mean of model parameters. 

Code development was led by <a href="https://github.com/andyhouldcroft" target="_blank" rel="noopener noreferrer">Andrew Houldcroft</a>, with <a href="https://github.com/finnlindgren" target="_blank" rel="noopener noreferrer">Finn Lindgren</a> contributing code for the spatial cluster size distribution and cluster size effect in the hazard-rate detection function.

Please see the <a href="https://github.com/inlabru-org/inlabru/tree/devel/vignettes/articles" target="_blank" rel="noopener noreferrer">vignettes</a> in the <a href="https://github.com/inlabru-org/inlabru" target="_blank" rel="noopener noreferrer">inlabru</a> repository for further guidance on using inlabru for spatial modelling.
Code in this repository are free to use, but please cite Houldcroft et al. (2024) if you do so.
