# Campbells

This repository contains code to perform the wildlife density modelling approach outlined in Houldcroft et al. (2023). Survey data may be requested but access cannot be guaranteed due to the high vulnerability of Campbell’s monkeys to hunting at present.

<b>H23_SETUP.R</b> contains the code to produce spatial covariates, a barrier SPDE mesh and format all data as sf and terra objects compatible with inlabru.

<b>H23_DENSITIES.R</b> contains the code to run the point process model with a spatial cluster size distribution, cluster size detection function and barrier SPDE models.

Code development was led by <strong><a href="https://github.com/andyhouldcroft" target="_blank" rel="noopener noreferrer">Andrew Houldcroft</a></strong>, with <strong><a href="https://github.com/finnlindgren" target="_blank" rel="noopener noreferrer">Finn Lindgren</a></strong> contributing code for the spatial cluster size distribution and cluster size effect in the hazard-rate detection function.

Please see the <strong><a href="https://github.com/inlabru-org/inlabru/tree/devel/vignettes/articles" target="_blank" rel="noopener noreferrer">vignettes</a></strong> in the <strong><a href="https://github.com/inlabru-org/inlabru" target="_blank" rel="noopener noreferrer">inlabru</a></strong> repository for further guidance on using inlabru for spatial modelling.
Code in this repository are free to use, but please cite Houldcroft et al. (2023) if you do so.
