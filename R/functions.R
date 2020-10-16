create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }
  
  return(labs)
}

floor_dec <- function(x, digits=1) round(x - 5*10^(-digits-1), digits)
ceiling_dec <- function(x, digits=1) round(x + 5*10^(-digits-1), digits)

map <- function(x, shape, palette, labs = NULL, breaks = NULL, n, digits = 0, 
                legend_title, panel_title = NULL, labg = F, ...) {
  
  # Generate breaks for color legend
  if (is.null(breaks)) {
    breaks <- as.numeric(quantile(shape[[x]], probs = seq(0, 1, l = n + 1), na.rm = T))
    breaks[2:n] <- round(breaks[2:n], digits)
    breaks[1] <- floor_dec(breaks[1], digits)
    breaks[n + 1] <- ceiling_dec(breaks[n + 1], digits) 
  }
  
  # Generate labels for color legend
  if (is.null(labs)) labs <- create_labels(breaks, greater = labg)
  
  
  
  # Generate color palette
  if (is.character(palette)) {
    pal <- tmaptools::get_brewer_pal(palette, n = length(labs), contrast = c(0, 1), plot = F)
  } else {
    pal <- palette(length(labs))
  }
  
  # Produce map
  tm_shape(shape) +
    tm_polygons(col = x, palette = pal, breaks = breaks, style = "fixed", 
                legend.show = T, labels = labs, title = legend_title,
                border.col = "black", colorNA = "grey", textNA = "Not analyzed") +
    tm_compass(position = c("right", "top")) +
    tm_scale_bar(position = c("right", "bottom"), text.size = .6) +
    tm_layout(outer.margins = 0, asp = 0, 
              # =legend.bg.color = "white", legend.frame = "black",
              legend.text.size = 0.95, legend.title.size = 1.35,
              panel.show = T, panel.labels = panel_title, 
              panel.label.size = 1.4, 
              ...)
}
