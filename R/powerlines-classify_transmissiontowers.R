#' Classify the transmission towers
#'
#' Attribute the class 15 to points that correspond to transmission towers using the positionning
#' of the towers given by \link{find_transmissiontowers}. Using the transmission towers types and their
#' orientation it computes the spatial box that emcompasses the towers and classify points within this
#' rectangle as 'transmission tower'.
#'
#' @param las An object of class LAS
#' @param towers SpatialPointsDataFrame returned by \link{find_transmissiontowers}.
#' @param dtm A RasterLayer. The digital terrain model is useful to find the bottom of the towers
#' @param threshold numeric. Height above ground. Points below this elevation are not classified
#' as transmission towers.
#'
#' @return A LAS object with an updated classification.
#'
#' @references
#' Roussel J, Achim A, Auty D. 2021. Classification of high-voltage power line structures in low density
#' ALS data acquired over broad non-urban areas. PeerJ Computer Science 7:e672 https://doi.org/10.7717/peerj-cs.672
#'
#' @examples
#' \donttest{
#' # A simple file with wires already clipped from 4 files + shapefile
#' # of the network
#' LASfile <- system.file("extdata", "wires.laz", package="lidRplugins")
#' wireshp <- system.file("extdata", "wires.shp", package="lidRplugins")
#' dtmtif  <- system.file("extdata", "wire-dtm.tif", package="lidRplugins")
#' las <- readLAS(LASfile, select = "xyzc")
#' network <- sf::st_read(wireshp)
#' dtm <- raster::raster(dtmtif)
#'
#' towers <- find_transmissiontowers(las, network, dtm, "waist-type")
#' las <- classify_transmissiontowers(las, towers, dtm)
#'
#' plot(las, color = "Classification")
#' }
#' @family electrical network
#' @export
classify_transmissiontowers = function(las, towers, dtm, type = NULL, threshold = 2)
{
  UseMethod("classify_transmissiontowers", las)
}

#' @export
classify_transmissiontowers.LAS = function(las, towers, dtm, type = NULL, threshold = 2)
{
  towers <- sf::st_crop(towers, lidR::st_bbox(las))
  towers <- tower.boundingbox(towers, type)
  sf::st_crs(towers) <- lidR::st_crs(las)

  tmp <- lidR::merge_spatial(las, towers, "towers")
  tmp <- lidR::normalize_height(tmp, dtm)
  las@data$Classification[tmp$Z > threshold & tmp$towers == TRUE] <- lidR::LASTRANSMISSIONTOWER
  return(las)
}

tower.boundingbox = function(towers, type = NULL)
{
  if (length(towers) == 0L)
  {
    data = data.frame(id = integer(0), maxZ = numeric(0), minZ = numeric(0))
    out = sf::st_sf(data, geometry = sf::st_sfc())
    return(out)
  }

  lines <- vector("list", nrow(towers))
  for (i in 1:nrow(towers))
  {
    tower <- towers[i,]
    if (is.null(type)){
      tower.spec <- get_tower_spec(tower$type)
    } else {
      tower.spec <- get_tower_spec(type)
    }
    
    width <- tower.spec$width[2]
    hwidth <- width/2
    p1 <- sf::st_coordinates(tower)
    ux <- tower$ux
    uy <- tower$uy

    p2 <- p1
    p2[,1] <- p2[,1] + ux * hwidth
    p2[,2] <- p2[,2] + uy * hwidth
    p1[,1] <- p1[,1] - ux * hwidth
    p1[,2] <- p1[,2] - uy * hwidth

    lines[[i]] <- sf::st_sf(id = as.character(i), geometry = sf::st_sfc(sf::st_linestring(rbind(p1, p2))))
  }

  tower.orientation <- do.call(rbind, lines)
  #plot(tower.orientation, add = T)

  # Compute an extent for the tower by buffering the lines
  towers.extent <- sf::st_buffer(tower.orientation, dist = tower.spec$length[2]/2, endCapStyle = "FLAT")
  towers.extent$maxZ <- towers$Z
  towers.extent$minZ <- towers$dtm

  return(towers.extent)
}
