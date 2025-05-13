# Load required library
library(terra)
library(ggplot2)
# also rnaturalearth for world geometry 

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
a <- world[world$name %in% c("Romania", "Serbia", "Bulgaria"), ]
b <- sf::st_union(a)
c <- sf::st_make_valid(b)

# Plot
ggplot(data = c) +
    geom_sf(fill = "white") +
    theme_bw() +
    theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
    )

wd <- "~/wallacean" # replace with your working directory
env_dir <- "~/env" # replace with your env directory
env_files <- c(
    "CHELSA_bio_4.tif", # mean annual temp
    "CHELSA_bio_13.tif", # precip of wettest 1/4
    NULL
)
env.crs <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # nolint
env_paths <- file.path(env_dir, env_files)

env <- raster::stack()
for (env_path in env_paths) {
    print(paste0("Loading ", env_path, " ..."))
    env_layer <- raster::raster(env_path)

    # If needed, reproject the env layer
    if (!raster::compareCRS(raster::crs(env_layer), env.crs)) {
        print(paste0(" Reprojecting..."))
        env_layer <- raster::projectRaster(env_layer, env.crs)
    }

    env <- raster::addLayer(env, env_layer)
}

# Temp
tras <- raster::crop(env$CHELSA_bio_4, raster::extent(a))
tras2 <- raster::mask(tras, a)
tras3 <- raster::scale(tras2)

# Precip
pras <- raster::crop(env$CHELSA_bio_13, raster::extent(a))
pras2 <- raster::mask(pras, a)
pras3 <- raster::scale(pras2)

ex <- c(25, 26, 46.5, 47.5)
tras4 <- raster::crop(tras3, ex)

# Create a species distribution
# Define logistic relationship
beta0 <- -1 # intercept
beta1 <- 2.5 # temperature effect
temp_scaled <- terra::values(pras3) # scaled temperature values
occ_prob <- plogis(beta0 + beta1 * temp_scaled)

# Simulate true presence/absence
set.seed(42)
true_presence <- pras3
values(true_presence) <- rbinom(ncell(pras3), size = 1, prob = occ_prob)
true_presence <- terra::rast(true_presence)

# Define the sampling process

# Helper: extract presence at points
get_obs <- function(points, truth_raster) {
    extracted <- extract(truth_raster, points)
    points$obs <- extracted[, 2]
    return(points)
}

# a. Random sampling
set.seed(1)
pras3_terra <- terra::rast(pras3) # Convert to SpatRaster

# Ensure raster has valid values (e.g., temperature or presence)
valid_mask <- !is.na(true_presence) # Mask of non-NA cells

# Mask your raster to just valid cells
pras3_terra <- mask(pras3_terra, valid_mask)

# Now sample only from valid cells

random_pts <- spatSample(pras3_terra, size = 300, method = "random", as.points = TRUE, na.rm = TRUE)

# Extract observed values from the true presence raster
random_pts <- get_obs(random_pts, true_presence)

# Vis
terra::plot(random_pts)
terra::lines(a)

# b. Systematic sampling
set.seed(2)
grid_pts <- spatSample(pras3_terra, size = 300, method = "regular", as.points = TRUE)
grid_pts <- get_obs(grid_pts, true_presence)

# Vis
terra::plot(grid_pts)
terra::lines(a)

# TODO: Discard points outside the study area

# c. Clustered sampling
set.seed(3)
# Sample many random points, keep only those within distance of a center
candidate_pts <- spatSample(pras3_terra, size = 2000, method = "random", as.points = TRUE)

# Choose a bias center in raster space (e.g., lower-left corner)
bias_center <- matrix(c(
    ext(pras3_terra)[1] + 0.2 * (ext(pras3_terra)[2] - ext(pras3_terra)[1]),
    ext(pras3_terra)[3] + 0.2 * (ext(pras3_terra)[4] - ext(pras3_terra)[3])
), ncol = 2)
bias_center_vect <- vect(bias_center, type = "points", crs = crs(pras3_terra))

# Calculate distance from bias center
dists <- distance(candidate_pts, bias_center_vect)
cluster_pts <- candidate_pts[dists < 0.2 * max(dists), ][1:300]

# Remove points that do not intersect with raster pras3_terra
valid_mask <- !is.na(terra::extract(pras3_terra, cluster_pts)[, 2])
cluster_pts <- cluster_pts[valid_mask, ]

cluster_pts <- get_obs(cluster_pts, true_presence)

# Vis
terra::plot(a$geometry)
terra::points(cluster_pts)


# d. Repeated samples from two sites 
set.seed(4)
# Define two fixed coordinates
site_coords <- matrix(
    c(
        ext(pras3_terra)[1] + 0.25 * (ext(pras3_terra)[2] - ext(pras3_terra)[1]),
        ext(pras3_terra)[3] + 0.25 * (ext(pras3_terra)[4] - ext(pras3_terra)[3]),
        ext(pras3_terra)[1] + 0.75 * (ext(pras3_terra)[2] - ext(pras3_terra)[1]),
        ext(pras3_terra)[3] + 0.75 * (ext(pras3_terra)[4] - ext(pras3_terra)[3])
    ),
    ncol = 2, byrow = TRUE
)

# Jitter around these locations
n_repeats <- 150
rep_coords <- do.call(rbind, lapply(1:n_repeats, function(i) {
    site <- site_coords[sample(1:2, 1), ]
    jitter <- rnorm(2, mean = 0, sd = 0.005 * (ext(pras3_terra)[2] - ext(pras3_terra)[1]))
    return(site + jitter)
}))
rep_pts <- vect(rep_coords, type = "points", crs = crs(pras3_terra))
rep_pts <- get_obs(rep_pts, true_presence)

# Vis 
terra::plot(a$geometry)
terra::points(rep_pts)
