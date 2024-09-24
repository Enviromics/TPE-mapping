# TPE-mapping

# Defining Target Population of Environments to Enviromics Studies Using R-based GIS Tools

The “Target Population of Environments” (TPE) is a concept used to
describe specific geographic regions where plant varieties or cultivars
are cultivated. It aids in identifying and characterizing these
environments, guiding the selection of genotypes that are better adapted
and more likely to perform well. Defining the TPE involves integrating
high-resolution environmental and phenotypic data, allowing for the
prediction of genotype performance under varying conditions and
tailoring selection to the unique attributes of each location.

To operationalize this, the objective is to design an R function called
**TPEmap**. This function will serve to clip the TPE region for any
crop, enabling users to clip the TPE boundaries based on a specified
geoprocessing area. By doing so, breeding researchers can focus their
analysis on the relevant environmental zones within the target area,
ensuring more accurate and region-specific selection strategies.

# 1. **TPEmap** Function

The **TPEmap** function assists in creating a TPE for plant breeding
using geospatial data from MET or on-farm trials. Key arguments include
adjustments for buffers, polygon concavity, and pixel size for
rasterization.

The following are the parameters that can be detailed in a TPE (Target
Population of Environments):

• `coordinates`: A dataframe with X (longitude) and Y (latitude)
columns, representing the geographic coordinates of MET or on-farm trial
points.

• `point_buffer`: A numeric value that defines the buffer distance to be
applied around each point (in kilometers). This argument allows the user
to adjust the size of the area of influence around the points.

• `concavity`: A numeric value that defines the degree of polygon
concavity. Lower values create a more detailed polygon, while higher
values (up to infinity) result in a Convex Hull.

• `length_threshold`: A numeric value that defines the edge length
threshold in the concavity algorithm. Edge segments with lengths below
this value are not considered for additional detail. This argument
allows the user to control the level of detail in the polygon.

• `expansion_buffer`: A numeric value that defines the additional buffer
distance to be applied after the concave polygon is generated, expanding
the final TPE (in kilometers).

• `simplify_tolerance`: A numeric value that defines the simplification
tolerance of the final polygon. This option allows the user to smooth
the polygon, removing excessive detail and small irregularities.

• `pixel_size`: A numeric value that defines the pixel size for
rasterizing the TPE. This argument is used in the conversion of the
final polygon to a raster.

## 1.1. Source the functions; and the **ggplot2** pack:

The first step is to execute the **TPEmap.R** file, which includes the
generate\_coordinates and **TPEmap** functions. Additionally, the
**ggplot2** package is used visualizing data, especially creating maps
with layers.

    source("TPEmap.R")

    ## pacote 'proxy' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'e1071' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'wk' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'classInt' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 's2' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'units' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'sf' desempacotado com sucesso e somas MD5 verificadas
    ## 
    ## Os pacotes binários baixados estão em
    ##  C:\Users\Administrador\AppData\Local\Temp\RtmpmAzNLV\downloaded_packages
    ## pacote 'V8' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'concaveman' desempacotado com sucesso e somas MD5 verificadas
    ## 
    ## Os pacotes binários baixados estão em
    ##  C:\Users\Administrador\AppData\Local\Temp\RtmpmAzNLV\downloaded_packages
    ## pacote 'sp' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'terra' desempacotado com sucesso e somas MD5 verificadas
    ## pacote 'raster' desempacotado com sucesso e somas MD5 verificadas
    ## 
    ## Os pacotes binários baixados estão em
    ##  C:\Users\Administrador\AppData\Local\Temp\RtmpmAzNLV\downloaded_packages

    library("ggplot2")

## 1.2. Example usage of the **generate\_coordinates** function:

Here is an example of using this function, which can be adjusted
according to your data.

To view the list of countries in the Americas for which the
**generate\_coordinates** and **TPEmap** functions will work

    Countries_of_the_Americas 

    ##  [1] "Canada"             "USA"                "Mexico"            
    ##  [4] "Guatemala"          "Honduras"           "El Salvador"       
    ##  [7] "Nicaragua"          "Costa Rica"         "Panama"            
    ## [10] "Cuba"               "Haiti"              "Dominican Republic"
    ## [13] "Jamaica"            "Bahamas"            "Brazil"            
    ## [16] "Argentina"          "Bolivia"            "Chile"             
    ## [19] "Colombia"           "Ecuador"            "Guyana"            
    ## [22] "Paraguay"           "Peru"               "Suriname"          
    ## [25] "Uruguay"            "Venezuela"

### 1.3. Generate Coordinates

Here, you can generate the number of desired coordinates based on the
selected country.

    coords <- generate_coordinates(n = 20, country_name = "Brazil")
    head(coords)

    ##           X          Y
    ## 1 -46.47055  -4.497927
    ## 2 -53.94479 -24.832044
    ## 3 -47.57103  -6.451637
    ## 4 -65.75033  -2.362546
    ## 5 -45.67567  -3.876433
    ## 6 -47.89452  -1.515029

**Your own coordinates data points should be inserted at this stage.**
Instead, to illustrate the use of the **TPEmap** function, you can use
the **generate\_coordinates** function to generate simulated data
applicable to **TPEmap**. The `coords` object represents a set of
locations (multiple geographic coordinates) where your actual data,
originating from experiments or on-farm trials, will be inserted for
analysis.

### 1.4. Example usage of the TPEmap function to generate the concave hull and raster

#### 1.4.1. TPEmap default arguments:

    # default arguments. Not run: 
    # TPEmap(coordinates =,  point_buffer = 100, concavity = 2, length_threshold = 10, 
          # expansion_buffer = 100, simplify_tolerance = 0.01, pixel_size = 0.01)

    TPEresult <- TPEmap(coords,simplify_tolerance = 10,pixel_size = 0.1)

    ## Total number of pixels in the raster: 109525 
    ## Number of rasterized pixels (non-NA): 48764

![](TPEmap_applying_files/figure-markdown_strict/unnamed-chunk-6-1.png)

##### Pixel=11.1km

The table below shows the approximate correspondence between variations
in geographic coordinate degrees and the equivalent distance in
kilometers, considering 1º of latitude or longitude to be equivalent to
111 km. This can be used to inform the `pixel_size argument`.

<table>
<thead>
<tr class="header">
<th>Degrees (°)</th>
<th>Distance (km)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>0.01</td>
<td>1.11</td>
</tr>
<tr class="even">
<td>0.10</td>
<td>11.10</td>
</tr>
<tr class="odd">
<td>0.50</td>
<td>55.50</td>
</tr>
<tr class="even">
<td>1.00</td>
<td>111.00</td>
</tr>
<tr class="odd">
<td>5.00</td>
<td>555.00</td>
</tr>
<tr class="even">
<td>10.00</td>
<td>1110.00</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
</tr>
</tbody>
</table>

#### 1.4.2. List of the TPEmap results:‘polygon’ ‘raster’

    names(TPEresult)

    ## [1] "polygon" "raster"

### 1.5. Extract the final polygon from TPEresult

    final_polygon <- TPEresult$polygon

# 2. MAP

## 2.1. Load shapefiles of Brazil and surrounding Latin American countries:

    americas <- st_as_sf(maps::map("world", regions = Countries_of_the_Americas, 
                                plot = FALSE, fill = TRUE))

## 2.2. Transform CRS to WGS84

    americas <- st_transform(americas, crs = 4326) #CRS=Coordinate Reference System

## 2.3. Define zoom limits (around the Country)

    x_range <- range(coords$X)
    y_range <- range(coords$Y)

### 2.3.1. Setting a zoom % for the clipped TPE.

    zoom_perc<- 0.25 

The lower this `zoom_pack` value, the greater the zoom applied.

    zoom_lims <- list(
      x = c(x_range[1] - (zoom_perc * (x_range[2] - x_range[1])),  # Min X with 25% decrease
            x_range[2] + (zoom_perc * (x_range[2] - x_range[1]))), # Max X with 25% increase
      y = c(y_range[1] - (zoom_perc * (y_range[2] - y_range[1])),  # Min Y with 25% decrease
            y_range[2] + (zoom_perc * (y_range[2] - y_range[1])))  # Max Y with 25% increase
    )

## 2.4. Convert coords dataframe into an sf object

    coords_sf <- st_as_sf(coords, coords = c("X", "Y"), crs = 4326)

## 2.5. Create the map with ggplot2

    map <- ggplot() +
      geom_sf(data = americas, fill = "gray90", color = "white") +  
      geom_sf(data = final_polygon, fill = "lightblue", color = "blue", alpha = 0.5) +  
      geom_sf(data = coords_sf, color = "slateblue", size = 2) +  
      coord_sf(xlim = zoom_lims$x, ylim = zoom_lims$y, expand = FALSE) +  
      theme_minimal()

Latin America with light gray

Plot the final polygon

Plot the original points

Apply zoom limits

## 2.6. Show the map

    print(map)

![](TPEmap_applying_files/figure-markdown_strict/unnamed-chunk-16-1.png)

# 3. RASTER BASE

## 3.1. Load the raster file if needed (Optional)

## 3.2. Optionally, you could visualize the raster directly if needed using:

    raster_polygon<-TPEresult$raster

## 3.3. Saving the raster (TPE) as a .tif file

    writeRaster(raster_polygon, filename = "raster_TPE.tif", overwrite = TRUE)
    getwd() #Access your TPE here

By applying the **TPEmap** function, we will generate a base raster file
in **.tif** format that can be used as a mask to retrieve numerous
envirotyping data. This allows plant breeders to associate these data
with phenotypic information and subsequently perform enviromic analyses.
