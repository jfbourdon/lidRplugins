#' Track the wires of the powerlines using the tower positions
#'
#' Assuming the coordinates, the elevation and the type of transmission tower are known the functions
#' tracks the wires of the powerlines. To achieve this task it computes the catenary equation of the
#' wires between two consecutive towers. Strictly speaking the function is not able to find the wires
#' instead it guesses their equation that can be computed deterministically from the tower coordinates,
#' their height and their type.
#'
#' @param towers \code{sf POINT} or \code{sfc POINT} containing the positions of the towers as returned by
#' \link{find_transmissiontowers}
#' @param powerline \code{sf LINESTRING} or \code{sfc LINESTRING} that map the electrical network accurately. The method may
#' be improved later to get rid of this information.
#' @param dtm A \code{SpatRaster}. Digital Terrain Model is useful to find a relevant elevation for virtual
#' towers
#' @param type One of "waist-type" or "double-circuit" according to
#' \href{http://www.hydroquebec.com/learning/transport/types-pylones.html}{Hydro-Quebec}. Can also
#' be a list with custom specifications. See examples.
#' @param debug logical. Plot the different steps of the algorithm so one can try to figure out what
#' is going wrong.
#'
#' @return  A \code{sf POINT} that actually represents 3D lines with several attributes
#' per points. \code{Z} the elevation, \code{virtual} tells if the wire has been found with two
#' consecutive tower and is thus accurate or if only one tower was used and in this case the wire is
#' a pure guess, \code{section} attributes an ID to each wire section i.e. between two towers, \code{ID}
#' attributes and ID to each powerline \code{type} store the transmission tower type.
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
#' dtm <- terra::rast(dtmtif)
#'
#' towers <- find_transmissiontowers(las, network, dtm, "waist-type")
#' wires <- track_wires(towers, network, dtm, "waist-type")
#'
#' col <- c("red", "blue", "forestgreen", "darkorchid", "darkorange")[wires$ID]
#' col[wires$virtual & col == "red"] <- "pink"
#' col[wires$virtual & col == "blue"] <- "lightblue"
#' col[wires$virtual & col == "forestgreen"] <- "lightgreen"
#' col[wires$virtual & col == "darkorchid"] <- "plum"
#' col[wires$virtual & col == "darkorange"] <- "goldenrod1"
#'
#' plot(header(las))
#' plot(towers, add = TRUE, col = towers$deflection + 1)
#' plot(wires, col = col, add = TRUE, cex = 0.1)
#'
#' plot(las, clear_artifacts = FALSE)
#' rgl::points3d(wires@coords[,1], wires@coords[,2], wires$z, col = col, size = 5)
#' }
#' @family electrical network
#' @export
track_wires <- function(towers, powerline, dtm, type = c("waist-type", "double-circuit"), debug = FALSE)
{
  if (!(is(towers, "sf") | is(towers, "sfc"))) stop("towers must be object of class sf of sfc")
  if (!(is(powerline, "sf") | is(powerline, "sfc"))) stop("powerline must be object of class sf of sfc")

  tower.spec <- get_tower_spec(type)

  tlocation <- towers
  proj <- sf::st_crs(tlocation)

  # If 0 tower the question is closed: return nothing
  if (length(tlocation) == 0)
  {
    data   <- data.frame(z = numeric(0), virtual = integer(0), ID = integer(0), type = character())
    output <- sf::st_sf(data, geometry = sf::st_sfc(), crs = sf::st_crs(towers))
    return(output)
  }

  if (debug)
  {
    opar = graphics::par("mfrow")
    graphics::par(mfrow = c(2,3))
    on.exit(graphics::par(mfrow = opar))
  }

  # Crop the lines to the extent of the ROI
  pwll <- sf::st_crop(powerline, terra::ext(dtm))
  if (nrow(pwll) == 0)
  {
    data   <- data.frame(z = numeric(0), virtual = integer(0), ID = integer(0), type = character())
    output <- sf::st_sf(data, geometry = sf::st_sfc(), crs = sf::st_crs(towers))
    return(output)
  }

  if (debug)
  {
    plot(terra::ext(dtm), main = paste0("Raw powerline network"), asp = 1)
    plot(pwll, add = T, col = 1:length(pwll))
  }

  # This starts like the tower detection by fixing the shapefile
  pwll <- pwll |>
    sf::st_union() |>
    sf::st_line_merge() |>
    sf::st_cast('LINESTRING') |>
    sf::st_as_sf() |>
    sf::st_simplify(preserveTopology = FALSE, dTolerance = 40)
  spwll <- gSplitLines(pwll)
  spwll <- sf::st_as_sf(spwll)
  spwlp <- sf::st_buffer(spwll, dist = 125, endCapStyle = 'SQUARE')

  if (debug)
  {
    plot(terra::ext(dtm), main = "Post-processed lines", asp = 1)
    plot(spwll, add = T, col =  1:length(spwll))
    plot(spwlp, add = T, border =  1:length(spwlp), lty = 3)
  }

  # Decompose the wires with different orientations into linear sections
  if (any(towers$deflection))
  {
    angles <- unique(round(tlocation$theta, 2))
    tlocations <- vector("list", length(angles))
    for (i in 1:length(angles)) tlocations[[i]] <- tlocation[round(tlocation$theta,2) == angles[i],]
  } else {
    tlocations <- list(tlocation)
  }

  if (length(tlocations) != nrow(spwll))
    stop("Internal error: different number of sections.", call. = FALSE)

  if (debug)
  {
    plot(terra::ext(dtm), main = "Detection of the powerlines", asp = 1)
  }

  # Loop on each section
  ID <- 1                    # ID for each line
  SECTION <- 0               # ID for each section
  HXY <- vector("list", length(tlocations))
  for (kk in 1:length(tlocations))
  {
    # Get the section
    tlocation <- tlocations[[kk]]
    pwlp <- spwlp[kk,]
    #plot(tiles)
    #plot(tlocation, add = T)
    #plot(pwlp, add = T)

    # Initialize vars
    n <- nrow(tlocation)
    posx <- sf::st_coordinates(tlocation)[,1]
    posy <- sf::st_coordinates(tlocation)[,2]
    wire <- vector("list", n)
    #plot(tiles)
    #plot(tlocation, add = T, col = tlocation$deflection +1)
    #arrows(posx, posy, posx + 100*tlocation$ux, posy + 100*tlocation$uy, length = 0.05)

    # Draw 1000 m lines passing throught each tower
    # Then buffer to create a polygons that encommpass each wire line
    for (i in 1:n)
    {
      l <- 500
      x <- posx[i]
      y <- posy[i]
      ux <- tlocation$ux[i]
      uy <- tlocation$uy[i]
      coords <- matrix(c(x - l*ux, y - l*uy, x + l*ux, y + l*uy), ncol = 2, byrow = T)
      wire[[i]] <- sf::st_sf(ID = as.character(i), geometry = sf::st_sfc(sf::st_linestring(coords)))
    }

    lwires <- do.call(rbind, wire)
    sf::st_crs(lwires) <- proj
    crlwires <- sf::st_crop(lwires, terra::ext(dtm) - 1)
    pwires <- sf::st_buffer(lwires, dist = 0.3*tower.spec$length[2], endCapStyle = "SQUARE")
    pwires_crop <- sf::st_crop(pwires, pwlp)
    pwires_crop <- sf::st_cast(pwires, 'POLYGON')
    if (nrow(pwires_crop) == 0)
    {
      data   <- data.frame(z = numeric(0), virtual = integer(0), ID = integer(0), type = character())
      output <- sf::st_sf(data, geometry = sf::st_sfc(), crs = proj)
      return(output)
    }
    
    if (debug)
    {
      plot(pwires, add = T, col = 1:n+1)
      plot(lwires, add = T)
      plot(tlocation, add = T, col = kk)
    }
    
    # Groups tower series
    sees <- sf::st_intersects(pwires, tlocation, sparse=FALSE)
    groups <- vector(mode = 'list', length = nrow(sees))
    for (i in 1:nrow(sees)){
      groups[i] <- list(which(sees[i,]))
    }
    while (sum(sapply(groups, length)) != n){
      for (i in 1:length(groups)){
        for (j in 1:length(groups)){
          if (length(intersect(groups[[i]], groups[[j]]))>0){
            groups[i] <- list(union(groups[[i]], groups[[j]]))
          }
        }
      }
      groups <- lapply(groups, sort)
      groups <- unique(groups)
    }

    # Generate catenary between two consecutive towers
    nlines <- length(groups)
    Hxy <- vector("list", nlines)
    for (i in 1:nlines)
    {
      # Get the towers for the processing line
      line <- sf::st_union(pwires[groups[[i]],])
      toww <- tlocation[groups[[i]],]
      tow  <- toww[, c("Z", "deflection")]
      tow$virtual = FALSE

      # Generate virtual towers. Virtual towers are non existing
      # towers that help to prolongate the lines when there is not
      # a second towers
      lin <- sf::st_intersection(crlwires, line)
      x1 <- sapply(lin$geometry, function(x) {
        m <- sf::st_coordinates(x)
        m <- m[,1:2]
        i <- which.max(m[,1])
        m[i,]
      })
      x1 <- matrix(x1[,which.max(x1[1,])], ncol = 2)

      x2 <- sapply(lin$geometry, function(x) {
        m <- sf::st_coordinates(x)
        m <- m[,1:2]
        i <- which.min(m[,1])
        m[i,]
      })
      x2 <- matrix(x2[, which.min(x2[1,])], ncol = 2)
      
      vtowers1 <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(x1), crs = proj))
      vtowers2 <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(x2), crs = proj))
      vtowers <- rbind(vtowers1, vtowers2)
      sf::st_crs(vtowers) <- sf::st_crs(tow)
      Z <- terra::extract(dtm, vtowers)[,2] + mean(tlocation$Z - tlocation$dtm)

      if (anyNA(Z)) stop("Impossible to find DTM value at the edge of the raster. The DTM is not large enough.")

      vtowers <- sf::st_sf(geometry = sf::st_sfc(vtowers$geometry), data.frame(Z, deflection = FALSE, virtual = TRUE))
      tow <- rbind(tow, vtowers)
      #plot(tow, add = T, col = tow$virtual + 1, cex = 2)

      # Order the towers by distance to an arbitrary point so they are in good order in the SPDF
      atower <- toww[1,]
      pmin <- matrix(c(sf::st_coordinates(atower)[,1] + atower$ux * 5000, sf::st_coordinates(atower)[,2] + atower$uy * 5000), ncol = 2)
      pmin <- sf::st_sfc(sf::st_point(pmin), crs = proj)
      d <- sf::st_distance(pmin, tow)
      j <- order(d)
      tow <- tow[j,]
      #plot(tow, add = T, col = "blue")
      #plot(tow, add = T, col = tow$deflection + 2)

      # Remove the virtual towers if not needed.
      # VT are not needed if they prolongate a line at a deflection point.
      # Special case if there is only a deflection. In this case we are scrapped
      if (sum(tow$deflection) <= sum(!tow$virtual))
      {
        rm = c(FALSE, tow$deflection[-nrow(tow)] & tow$virtual[-1]) | c(tow$deflection[-1] & tow$virtual[-nrow(tow)], FALSE)
        tow = tow[!rm,]
      }
      #plot(tiles)
      #plot(tow, add = T, col = "blue")

      # If we have more than a tower we can compute the catenary between two consecutive towers
      # else we do not compute anything
      if (nrow(tow) > 1)
      {
        hxy <- vector("list", nrow(tow) - 1)
        for (k in 1:(nrow(tow) - 1))
        {
          p1 <- list(x = sf::st_coordinates(tow)[k,1], y =sf::st_coordinates(tow)[k,2], z = tow$Z[k] - tower.spec$wire.distance.to.top, virtual = tow$virtual[k])
          p2 <- list(x = sf::st_coordinates(tow)[k + 1, 1], y = sf::st_coordinates(tow)[k + 1, 2], z = tow$Z[k + 1] - tower.spec$wire.distance.to.top, virtual = tow$virtual[k + 1])
          hxy[[k]] <- catenary(p1, p2, tower.spec$tension)
          hxy[[k]]$virtual <- p1$virtual | p2$virtual
          hxy[[k]]$section <- SECTION
          SECTION <- SECTION + 1
        }

        hxy <- data.table::rbindlist(hxy)
        hxy$ID <- ID
        ID <- ID + 1
        Hxy[[i]] <- hxy
      } else {
        cc <- list(x = 0, y = 0, z = 0)
        res <- catenary(cc, cc, tower.spec$tension)
        res$virtual = TRUE
        res$section <- 0
        res$ID = 0
        Hxy[[i]] <- res[0,]
      }
    }

    Hxy = data.table::rbindlist(Hxy)
    points <- sf::st_sfc(sf::st_multipoint(as.matrix(Hxy[,c('x', 'y')])))
    points <- sf::st_cast(points, 'POINT')
    if (nrow(Hxy) == 0) {
      points <- points[0]
    }
    Hxy <- sf::st_sf(Hxy[,c('z', 'virtual', 'section', 'ID')], geometry = points)
    # rgl::points3d(Hxy@coords[,1],Hxy@coords[,2], Hxy$z)
    sf::st_crs(Hxy) <- proj
    Hxy <- Hxy[!is.na(sf::st_within(Hxy, pwlp, sparse = FALSE)),]
    HXY[[kk]] <- Hxy
  }

  HXY <- do.call(rbind, HXY)
  HXY$type = type$name
  wires <- HXY

  if (debug)
  {
    col <- c("red", "blue", "forestgreen", "darkorchid", "darkorange", "yellow")[wires$ID]
    col[wires$virtual & col == "red"] <- "pink"
    col[wires$virtual & col == "blue"] <- "lightblue"
    col[wires$virtual & col == "forestgreen"] <- "lightgreen"
    col[wires$virtual & col == "darkorchid"] <- "plum"
    col[wires$virtual & col == "darkorange"] <- "goldenrod1"
    col[wires$virtual & col == "yellow"] <- "white"

    plot(terra::ext(dtm), main = paste0("Final extraction"))
    plot(towers, add = T, col = towers$deflection + 1)
    #plot(textent, add = T,  border = textent$deflection + 1)
    plot(wires, col = col, add = T, cex = 0.1)
  }

  # clean wires below ground
  z0 <- terra::extract(dtm, wires)[,2]
  invalid <- unique(wires$section[which(wires$z - z0 < 0)])
  wires <- wires[!wires$section %in% invalid,]
  return(wires)
}

# Derivation of Equations for Conductor and Sag Curves of an Overhead Line Based on a Given Catenary Constant
# Alen Hatibovic
# Electrical Engineering and Computer Science 58/1 (2014) 23–27 doi: 10.3311/PPee.6993
catenary = function(p1, p2, c = 1500)
{
  x1 = p1$x
  x2 = p2$x
  y1 = p1$y
  y2 = p2$y
  h1 = p1$z
  h2 = p2$z
  dx = x2-x1
  dy = y2-y1
  dz = h2-h1
  S  = sqrt(dx*dx+dy*dy+dz*dz)
  x  = seq(0, S, by = 1)

  term0 = asinh((h2-h1)/(2*c*sinh(S/(2*c))))
  term1 = x - S/2
  term2 = c*term0
  term3 = S/(2*c)
  term4 = term0

  term5 = sinh(1/(2*c)*(term1+term2))^2
  term6 = sinh(0.5*(term3-term4))^2

  hx = 2*c * (term5 - term6) + h1

  x = seq(x1, x2, length.out = length(hx))
  y = seq(y1, y2, length.out = length(hx))
  z = hx
  return(data.frame(x,y,z))
}

