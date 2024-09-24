# Load necessary packages


# Function to install and load packages
install_and_load <- function(package) {
  # Check if the package is installed
  if (!require(package, character.only = TRUE)) {
    # Install the package if it's not installed
    install.packages(package)
    # Load the package
    library(package, character.only = TRUE)
  }
}

# List of required packages
packages <- c("sf", "maps", "concaveman", "raster", "terra")

# Loop through each package and install/load it
for (pkg in packages) {
  install_and_load(pkg)
}

##Packages:
#sf:          For spatial data handling, buffering, and transformations.
#maps:        To load geographic boundary data, such as Brazil's map.
#concaveman:  To generate concave hull polygons from geographic points.
#raster:      For manipulating and rasterizing geospatial grids and polygons.
#terra:       For raster data handling and saving in .tif format.

# Function to generate random coordinates within any country
generate_coordinates <- function(n, country_name) {
  
  # Load the shapefile for the specified country and transform CRS
  country <- st_geometry(st_transform(st_as_sf(maps::map("world", country_name, plot = FALSE, fill = TRUE)), 4326))
  
  # Get the bounding box of the country (longitude and latitude boundaries)
  bbox <- st_bbox(country)
  long_min <- bbox["xmin"]
  long_max <- bbox["xmax"]
  lat_min <- bbox["ymin"]
  lat_max <- bbox["ymax"]
  
  # Initialize dataframe to store valid coordinates
  coords <- data.frame(Long = numeric(0), Lat = numeric(0))
  
  # Generate coordinates until the required amount is reached
  while (nrow(coords) < n) {
    # Generate random points within the defined boundaries
    points_sf <- st_as_sf(data.frame(Long = runif(n, long_min, long_max), Lat = runif(n, lat_min, lat_max)), 
                          coords = c("Long", "Lat"), crs = 4326)
    
    # Filter points that are within the specified country
    valid_points <- st_intersection(points_sf, country)
    coords <- rbind(coords, st_coordinates(valid_points))
  }
  
  # Return only the required number of coordinates
  return(coords[1:n, ])
}

###########
##TPEmap ##
###########

TPEmap <- function(coordinates, point_buffer = 100, concavity = 2, length_threshold = 10, 
                   expansion_buffer = 100, simplify_tolerance = 0.01, pixel_size = 0.01) {
  
  # Define the target CRS (WGS84)
  target_crs <- 4326
  
  # Convert the dataframe of coordinates into an sf object with the target CRS
  coords_sf <- st_as_sf(coordinates, coords = c("X", "Y"), crs = target_crs)
  
  # Apply buffer to each point
  coords_buffered <- st_buffer(coords_sf, dist = point_buffer * 1000) # Convert km to meters
  
  # Calculate the concave hull
  concave_polygon <- concaveman::concaveman(coords_buffered, concavity = concavity, 
                                            length_threshold = length_threshold)
  
  # Apply an additional buffer around the concave hull
  final_polygon <- st_buffer(concave_polygon, dist = expansion_buffer * 1000) # Convert km to meters
  
  # Simplify the final polygon if needed
  final_polygon_simplified <- st_simplify(final_polygon, dTolerance = simplify_tolerance * 1000)
  
  # Ensure the final_polygon_simplified uses the correct CRS
  final_polygon_simplified <- st_transform(final_polygon_simplified, crs = target_crs)
  
  # Load the map for all the Americas (South, Central, and North)
  americas <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico", 
                                                      "Guatemala", "Honduras", "El Salvador", 
                                                      "Nicaragua", "Costa Rica", "Panama", 
                                                      "Cuba", "Haiti", "Dominican Republic",
                                                      "Jamaica", "Bahamas", "Brazil", "Argentina", 
                                                      "Bolivia", "Chile", "Colombia", "Ecuador", 
                                                      "Guyana", "Paraguay", "Peru", "Suriname", 
                                                      "Uruguay", "Venezuela"), 
                                 plot = FALSE, fill = TRUE))
  
  # Unite all countries in the Americas into one polygon (representing the landmass)
  americas_union <- st_union(americas)
  
  # Ensure americas_union uses the same CRS
  americas_union <- st_transform(americas_union, crs = target_crs)
  
  # Intersect the polygon with the Americas landmass to exclude ocean areas
  final_polygon_land <- st_intersection(final_polygon_simplified, americas_union)
  
  # Create the raster base based on the land-intersected polygon
  bbox <- st_bbox(final_polygon_land)
  raster_template <- raster(xmn = bbox["xmin"], xmx = bbox["xmax"], 
                            ymn = bbox["ymin"], ymx = bbox["ymax"], 
                            res = pixel_size, crs = st_crs(final_polygon_land)$proj4string)
  
  # Rasterize the polygon
  raster_polygon <- rasterize(final_polygon_land, raster_template, field = 1, background = NA)
  
  # Count the number of pixels in the raster and those within the polygon
  total_pixels <- ncell(raster_polygon)
  rasterized_pixels <- sum(!is.na(values(raster_polygon)))
  
  cat("Total number of pixels in the raster:", total_pixels, "\n")
  cat("Number of rasterized pixels (non-NA):", rasterized_pixels, "\n")
  
  # Plot the polygon and the rasterized version
  plot(raster_polygon, main = "Rasterized TPE (Excluding Ocean in the Americas)")
  
  # Convert the raster object to 'terra' format
  raster_polygon_terra <- rast(raster_polygon)
  
  return(list(polygon = final_polygon_land, raster = raster_polygon))
}

Countries_of_the_Americas <-c("Canada", "USA", "Mexico", "Guatemala", 
                  "Honduras", "El Salvador", "Nicaragua", 
                  "Costa Rica", "Panama", "Cuba", "Haiti", 
                  "Dominican Republic", "Jamaica", "Bahamas", 
                  "Brazil", "Argentina", "Bolivia", "Chile", 
                  "Colombia", "Ecuador", "Guyana", "Paraguay", 
                  "Peru", "Suriname", "Uruguay", "Venezuela")
