#' Plot an object of class PresenceAbsence
#' 
#' @author Bruno Vilela
#' 
#' @description Plots species richness map from an object of class PresenceAbsence or a
#' particular species' map.
#'
#' @param x An object of class \code{\link{PresenceAbsence}}.
#' @param name A character specifying a species to be ploted instead of the complete species richness map.
#' @param world If \code{TURE} a map of political divisions (countries) is added to the plot.
#' @param col_rich Color function (e.g. \code{\link{rainbow}}, \code{\link{heat.colors}}, \code{\link{colorRampPalette}}) to be used in the richness map.
#' @param col_name The color to use when ploting single species.
#' @param ... Other parameters pass to the plot function.
#' 
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.presab.birds}}
#' 
#' @examples \dontrun{
#' data(PAM)
#' plot(PAM)
#' plot(PAM, xlab = "Longitude", ylab = "Latitude",
#'      main = "Phyllomedusa species richness")
#' plot(PAM, name = "Phyllomedusa atelopoides")
#' plot(PAM, name = "Phyllomedusa azurea")
#' }
#' 
#' @method plot PresenceAbsence
#' @export

plot.PresenceAbsence <- function(x, name = NULL, world = TRUE,
                                 col_rich = NULL, col_name = "red",
                                 ...) {
  x <- .check_pam(x)
  # Richness plot
  if (is.null(name)) {
    
    # Creating the color function
    if (is.null(col_rich)) {      
      # Colour ramp from colorbrewer (T. Lucas suggestion)
      colfunc <- grDevices::colorRampPalette(c("#fff5f0", "#fb6a4a", "#67000d"))
    } else {
      colfunc <- col_rich
    }
    
    # Getting values  
    v <- terra::values(x$Rich)
    c <- max(v, na.rm = TRUE)
    
    # Set zero to NA
    v[(v == 0)] <- NA
    terra::values(x$Rich) <- v
    
    # Plot, add one and remove the first to not be white.
    plot(x$Rich, col = colfunc(c + 1)[-1], ...)  
  }
  
  # Individual plot
  if (!is.null(name)) {
    # Species position in the PAM
    pos <- which(x$Sp == name)
    # Transform the one species in raster
    r <- terra::rasterize(x$Presen[ , 1:2], 
                          x$Rich, 
                          x$Presen[ , (pos + 2)])
    # Plot
    plot(r, col = c("white", col_name), legend = FALSE, ...)
  }
  if (world) {
    wrld_simpl <- get(utils::data("wrld_simpl", package = "letsR"))
    plot(sf::st_geometry(wrld_simpl), add = TRUE)
  }

}

