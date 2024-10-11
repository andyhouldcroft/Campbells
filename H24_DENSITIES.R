################################################################################
# Article:  Joint spatial modeling of cluster size and density for a heavily
#           hunted primate persisting in a heterogeneous landscape
# Contact:  Andrew Houldcroft
#           ah1101@exeter.ac.uk
# Notes:    Please replace "..." with your own values.
################################################################################
library(INLA)
library(INLAspacetime)
library(inlabru)
library(fmesher)
library(sf)
library(terra)

# Set INLA mode
bru_options_set(inla.mode = "experimental")

# Import data
data <- readRDS("data.rds")

# Unwrap covariates
covariate.1 <- unwrap(data$covariates$covariate.1)
covariate.2 <- unwrap(data$covariates$covariate.2)

# Define key parameters
max.dist <- max(data$points$distance)
max.cluster <- max(data$points$cluster)

# Hazard-rate detection model
log_hr <- function(distance, log.sigma, log.gamma) {
  log1p(-exp(-(distance/exp(log.sigma))^-exp(log.gamma)))
}

# Discretized truncated log-Normal cluster size distribution
log_pd <- function(cluster, meanlog, log.sdlog) {
  bounds <- c(0, 1:max.cluster)
  log(plnorm(bounds[cluster+1], meanlog = meanlog, sdlog = exp(log.sdlog)) - 
        plnorm(bounds[cluster], meanlog = meanlog, sdlog = exp(log.sdlog))) -
    plnorm(max(bounds), meanlog = meanlog, sdlog = exp(log.sdlog), log.p = TRUE)
}

# Define barrier SPDE Matérn models
barrier.matern <- barrierModel.define(data$mesh, 
                                      data$barrier,
                                      prior.sigma = c("...", "..."),
                                      prior.range = c("...", "..."))

# Define the model components
components <- ~ 0 + 
  BSPDE(main = geometry, model = barrier.matern) +
  meanlog(main = geometry, model = barrier.matern) +
  beta.cov.1(covariate.1, model = "linear") +
  beta.cov.2(covariate.2, model = "linear") +
  beta.cluster(seq_len(max.cluster), model = "linear") +
  log.sigma(1) +
  log.gamma(1) +
  log.sdlog(1) +
  intercept(1) +
  intercept.2(1)

# Define the model formula
formula <- geometry + distance + cluster ~
  BSPDE +
  log_hr(distance, log.sigma + beta.cluster[cluster], log.gamma) +
  log_pd(cluster, meanlog + intercept.2, log.sdlog) +
  beta.cov.1 +
  beta.cov.2 +
  intercept +
  log(2)

# Run the model
model <- lgcp(
  components = components,
  data = data$points,
  samplers = data$samplers,
  domain = list(
    geometry = data$mesh,
    distance = fm_mesh_1d(seq(0, max.dist, length.out = 30)),
    cluster = seq_len(max.cluster)),
  formula = formula,
  options = list(
    bru_initial = list(beta.cov.1 = "...",
                       beta.cov.2 = "...",
                       beta.cluster = "...", 
                       log.sigma = "...", 
                       log.gamma = "...",
                       log.sdlog = "...",
                       intercept = "...", 
                       intercept.2 = "..."), 
    bru_verbose = 4, 
    bru_method = list(rel_tol = 0.15), 
    control.compute = list(dic = TRUE)))

# View model summary
summary(model)

# Define integration points
pred.points <- fm_int(data$mesh, data$site, format = "sf")

# Predict individual abundance
abundance <- predict(model, pred.points, ~ {
  bounds <- c(0, 1:max.cluster)
  cluster <- 1:max.cluster
  cluster.expectation <- vapply(
    seq_along(meanlog),
    function(k){
      prob.vector <- plnorm(
        bounds[cluster + 1], 
        meanlog = meanlog[k] + intercept.2, 
        sdlog = exp(log.sdlog)) - 
        plnorm(bounds[cluster], 
               meanlog = meanlog[k] + intercept.2, 
               sdlog = exp(log.sdlog))
      prob.vector <- prob.vector/plnorm(max(bounds), 
                                        meanlog = meanlog[k] + intercept.2, 
                                        sdlog = exp(log.sdlog))
      sum((1:max.cluster) * prob.vector)
    }, 0.0)
  sum(weight * exp(BSPDE + beta.cov.1 + beta.cov.2 + intercept) * 
        cluster.expectation)
}, 
n.samples = "...")
abundance

# Predict individual abundance
average.density <- predict(model, pred.points, ~ {
  bounds <- c(0, 1:max.cluster)
  cluster <- 1:max.cluster
  cluster.expectation <- vapply(
    seq_along(meanlog),
    function(k){
      prob.vector <- plnorm(
        bounds[cluster + 1], 
        meanlog = meanlog[k] + intercept.2, 
        sdlog = exp(log.sdlog)) - 
        plnorm(bounds[cluster], 
               meanlog = meanlog[k] + intercept.2, 
               sdlog = exp(log.sdlog))
      prob.vector <- prob.vector/plnorm(max(bounds), 
                                        meanlog = meanlog[k] + intercept.2, 
                                        sdlog = exp(log.sdlog))
      sum((1:max.cluster) * prob.vector)
    }, 0.0)
  sum(weight * exp(BSPDE + beta.cov.1 + beta.cov.2 + intercept) * 
        cluster.expectation) / sum(weight)
}, 
n.samples = "...")
average.density

# Predict cluster abundance
cluster.abundance <- predict(model, pred.points, ~ {
  sum(weight * exp(BSPDE + beta.cov.1 + beta.cov.2 + intercept))
}, 
n.samples = "...")
cluster.abundance

# Predict mean cluster density
average.cluster.density <- predict(model, pred.points, ~ {
  sum(weight * exp(BSPDE + beta.cov.1 + beta.cov.2 + intercept)) / sum(weight)
}, 
n.samples = "...")
average.cluster.density

# Predict cluster expectation
average.cluster.size <- predict(model, pred.points, ~ {
  bounds <- c(0, 1:max.cluster)
  cluster <- 1:max.cluster
  cluster.expectation <- vapply(
    seq_along(meanlog),
    function(k){
      prob.vector <- plnorm(
        bounds[cluster + 1], 
        meanlog = meanlog[k] + intercept.2, 
        sdlog = exp(log.sdlog)) - 
        plnorm(bounds[cluster], 
               meanlog = meanlog[k] + intercept.2, 
               sdlog = exp(log.sdlog))
      prob.vector <- prob.vector/plnorm(max(bounds), 
                                        meanlog = meanlog[k] + intercept.2, 
                                        sdlog = exp(log.sdlog))
      sum((1:max.cluster) * prob.vector)
    }, 0.0)
  sum(weight * cluster.expectation) / sum(weight)
}, 
n.samples = "...")
average.cluster.size

# Define prediction grid
grid <- fm_pixels(data$mesh, mask = TRUE)

# Predict individual density surface
density.surface <- predict(model, grid, ~ {
  bounds <- c(0, 1:max.cluster)
  cluster <- 1:max.cluster
  cluster.expectation <- vapply(
    seq_along(meanlog),
    function(k){
      prob.vector <- plnorm(
        bounds[cluster + 1], 
        meanlog = meanlog[k] + intercept.2, 
        sdlog = exp(log.sdlog)) - 
        plnorm(bounds[cluster], 
               meanlog = meanlog[k] + intercept.2, 
               sdlog = exp(log.sdlog))
      prob.vector <- prob.vector/plnorm(max(bounds), 
                                        meanlog = meanlog[k] + intercept.2, 
                                        sdlog = exp(log.sdlog))
      sum((1:max.cluster) * prob.vector)
    }, 0.0)
  exp(BSPDE + beta.cov.1 + beta.cov.2 + intercept) * cluster.expectation
}, 
n.samples = "...")
ggplot() + gg(density.surface, geom = "tile")

# Predict cluster density surface
cluster.density.surface <- predict(model, grid, ~ {
  exp(BSPDE + beta.cov.1 + beta.cov.2 + intercept)
}, 
n.samples = "...")
ggplot() + gg(cluster.density.surface, geom = "tile")

# Predict spatial cluster distribution
cluster.size.surface <- predict(model, grid, ~ {
  bounds <- c(0, 1:max.cluster)
  cluster <- 1:max.cluster
  cluster.expectation <- vapply(
    seq_along(meanlog),
    function(k){
      prob.vector <- plnorm(
        bounds[cluster + 1], 
        meanlog = meanlog[k] + intercept.2, 
        sdlog = exp(log.sdlog)) - 
        plnorm(bounds[cluster], 
               meanlog = meanlog[k] + intercept.2, 
               sdlog = exp(log.sdlog))
      prob.vector <- prob.vector/plnorm(max(bounds), 
                                        meanlog = meanlog[k] + intercept.2, 
                                        sdlog = exp(log.sdlog))
      sum((1:max.cluster) * prob.vector)
    }, 0.0)
  cluster.expectation
}, 
n.samples = "...")
ggplot() + gg(cluster.size.surface, geom = "tile")

# Predict barrier SPDE spatial effect
BSPDE.spatial.effect <- predict(model, grid, ~ {
  exp(BSPDE)
}, 
n.samples = "...")
ggplot() + gg(BSPDE.spatial.effect, geom = "tile")

# Predict covariate.1 spatial effect
cov.1.spatial.effect <- predict(model, grid, ~ {
  exp(beta.cov.1)
}, 
n.samples = "...")
ggplot() + gg(cov.1.spatial.effect, geom = "tile")

# Predict covariate.2 spatial effect
cov.2.spatial.effect <- predict(model, grid, ~ {
  exp(beta.cov.2)
}, 
n.samples = "...")
ggplot() + gg(cov.2.spatial.effect, geom = "tile")
