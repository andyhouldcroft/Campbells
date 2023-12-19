################################################################################
# Article:  Using inlabru to predict and map wildlife densities in heterogenous 
#           landscapes
# Contact:  Andrew Houldcroft
#           ah1101@exeter.ac.uk
# Notes:    Please replace "..." with your own values.
################################################################################
library(INLA)
library(fmesher)
library(sf)
library(terra)

# Set INLA mode
bru_options_set(inla.mode = "experimental")

# Define project CRS in kilometers
crs.km <- "..."

# Import detection locations, distances and cluster sizes
points <- read.csv("...", header = TRUE)
points <- st_as_sf(points, coords = c("longitude", "latitude"), crs = 4326)
points <- st_transform(points, crs.km)

# Import study site
site <- st_read(dsn = "~/...", layer = "...")
site <- st_transform(site, crs.km)

# Import mesh hull
hull <- st_read(dsn = "~/...", layer = "...")
hull <- st_transform(hull, crs.km)

# Import barrier polygon
barrier <- st_read(dsn = "~/...", layer = "...")
barrier <- st_transform(barrier, crs.km)

# Import line transect repeats
samplers <- st_read(dsn = "~/...", layer = "...")
samplers <- st_transform(samplers, crs.km)

# Import your vector data
vector <- st_read(dsn = "~/...", layer = "...")
vector <- st_transform(vector, crs.km)

# Import your raster data
raster <- rast("~/...")
raster <- project(raster, crs.km)

# Define triangulation mesh
mesh <- inla.mesh.2d(boundary = hull,
                     max.edge = c("...", "..."),
                     crs = crs.km)

# Define barrier triangles
barrier <- unlist(fm_contains(
  x = barrier, 
  y = mesh, 
  type = "centroid"))

# Compute covariate grid with 25-m resolution
grid <- rast(st_as_sf(data.frame(mesh$loc[,1:2]), 
                      coords = c("X1", "X2"), 
                      crs = crs.km),
             resolution = c(0.025, 0.025))

# Example 1: Proximity to vector in kilometers
covariate.1 <- distance(grid, vect(vector), unit = "km")

# Example 2: Mean of raster per 500-m buffer
buffer <- focalMat(raster, 0.5, type = "circle")
covariate.2 <- focal(raster, w = buffer, fun = "mean", na.rm = TRUE)
covariate.2 <- resample(covariate.2, grid, method = "average")

# Export data
saveRDS(list(points = points, 
             samplers = samplers, 
             site = site,
             hull = hull,
             mesh = mesh,
             barrier = barrier,
             covariates = list(covariate.1 = wrap(covariate.1),
                               covariate.2 = wrap(covariate.2))),
        "data.rds")
