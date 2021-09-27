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

extrac_var_climate <- function(var, data) {
  data_clean <- data %>% 
    select(Date, name, var) %>% 
    filter(name %in% colnames(counts)) %>% 
    tidyr::spread(key = name, value = var) %>% 
    select(-Date) %>% 
    as.matrix()
  rownames(data_clean) <- as.character(unique(data$Date))
  return(data_clean)
}

reshape_df <- function(x, value) {
  y <- as_tibble(t(x))
  y$COUNTRY <- colnames(x)
  y <- y %>% 
    tidyr::gather(key = "time", value = !!value, -COUNTRY) %>% 
    mutate(time = as.Date(time))
  return(y)
}

tidymat <- function(x, type, var) {
  colnames(x) <- colnames(fit$stsObj)
  x <- as_tibble(x)
  x$time <- as.Date(final_dates)
  x$type <- type
  x$var <- var
  tidyr::gather(x, key = "country", value = "contrib", -time, -type, -var)
} 

plot_preds <- function(predictions, plot_opts=NULL){
  
  
  fits <- predictions[predictions$type == "IS", ]
  preds <- predictions[predictions$type == "OUS", ]
  
  if(is.null(plot_opts)) plot_opts = list()
  
  xaxis = list(title="", 
               range = c(min(predictions$time), max(predictions$time) + 1))
  yaxis = list(title="Count")
  
  ## make the first prediction join up with the last data:
  # preds = rbind(fits[nrow(fits),], preds)
  
  fig = plot_ly(fits, x=~time)
  
  if(isFALSE(plot_opts$displayModeBar)) fig = config(fig, displayModeBar=FALSE)
  
  if(is.null(plot_opts$showlegend)) plot_opts$showlegend=TRUE
  
  fig <- fig %>% add_ribbons(ymin=~low95, ymax=~up95,
                             legendgroup="Model",
                             line=list(color="transparent"),
                             fillcolor='rgb(208,227,245)',
                             showlegend=plot_opts$showlegend, name="95% CI") # name doesn't show
  fig <- fig %>% add_ribbons(ymin=~low50, ymax=~up50,
                             legendgroup="Model",
                             line=list(color="transparent"),
                             fillcolor='rgb(171,205,237)',
                             showlegend=plot_opts$showlegend, name="50% CI") # name doesn't show
  fig  <- fig %>% add_trace(y=~mean, type="scatter",mode="lines",
                            legendgroup="Model",
                            line=list(color="black"),
                            showlegend=plot_opts$showlegend,
                            fillcolor='rgba(100,100,80,.2)', name="Mean")
  
  
  fig <- fig %>% add_ribbons(data=preds,
                             ymin=~low95, ymax=~up95,
                             legendgroup="Forecast",
                             line=list(color="transparent"),
                             fillcolor='rgb(255,237,204)',
                             showlegend=plot_opts$showlegend, 
                             name="95% CI") # name doesn't show
  fig <- fig %>% add_ribbons(ymin=~low50, ymax=~up50,
                             legendgroup="Forecast",
                             line=list(color="transparent"),
                             fillcolor='rgb(255,201,102)',
                             showlegend=plot_opts$showlegend, 
                             name="50% CI") # name doesn't show
  
  fig <- fig %>% add_trace(x=~time, y=~mean, data=preds,
                           legendgroup="Forecast",
                           line=list(color="red"),
                           type="scatter", mode="lines",
                           showlegend=plot_opts$showlegend, name="Forecast")
  
  fig <- fig %>% add_markers(x=~time, y=~observed, data=preds,
                             legendgroup="New Data",
                             marker = list(color="red", symbol="cross",
                                           line=list(color="white", width=1)),
                             showlegend=plot_opts$showlegend, 
                             name="New Data")
  
  fig <- fig %>% add_markers(data=fits, y=~observed, x=~time,
                             marker=list(color="grey",  symbol="circle",
                                         line=list(width=1, color="white")
                             ),
                             showlegend=plot_opts$showlegend,
                             name="Cases")
  
  fig <- fig %>% layout(xaxis=xaxis, 
                        yaxis=yaxis,
                        title=preds$country[1])
  return(fig)
  
}


plot_contributions <- function(pf, fitted, units){
  
  units_names <- units
  
  ids <- match(units, colnames(fitted$mean))
  results <- cbind(fitted$endemic[, ids[1]], fitted$epi.own[, ids[1]], fitted$epi.neighbours[, ids[1]])
  
  results <- as.data.frame(results)
  names(results) <- c("Endemic", "Within", "Between")
  results$cname <- rep(units_names, each = length(fitted$inference_dates))
  results$COUNTRY <- rep(units, each = length(fitted$inference_dates))
  results$time <- fitted$inference_dates
  
  cases_to_plot <- pf %>%
    filter(country %in% units) %>%
    rename(COUNTRY = country) %>%
    mutate(time = as.Date(time,format="%Y-%m-%d"))
  
  
  fig =  results %>%
    inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
    select(-COUNTRY) %>% 
    tidyr::gather(key = "Contribution", value = "value", -cname, -time, -observed) %>% 
    mutate(Contribution = factor(Contribution, 
                                 levels = c("Between", "Within", "Endemic"))) %>% 
    ggplot(aes(x = time)) +
    geom_area(aes(y = value, fill = Contribution), alpha = 1) +
    geom_line(aes(y = observed), linetype = 1, size = .3) +
    geom_point(aes(y = observed), size = 1) +
    facet_wrap(~ cname, scales = "free_y", ncol = 1) +
    scale_fill_brewer(type = "q", palette= 7)+
    scale_x_date(expand = c(0, 1)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(x = "", y = "Reported cases", fill = "") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))
  
  
  
  return(fig)
  
}