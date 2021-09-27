

##Loading Packages
library(dplyr)
library(amt)
library(shiny)
library(htmlwidgets)
library(htmltools)
library(leaflet)
library(geojsonio)
library(sf)
library(leafletCN)
library(ggplot2)
library(ggrepel)
library(plotly)
    library(knitr)
    library(kableExtra)
library(viridis)

##Sourcing Functions


source("R/functions.R")

#Check if RE data is available
if (file.exists("output/models/predictions_RE.rds")){
  suffix = "_RE"
} else {
  suffix = ""
}


##Reading in Data

#Shinydata.rds is saved in model_fitting
pf = readRDS("output/models/shinydata.rds")
pf_fitted = pf;


md = pf[[2]]
md$Pop2020_Mi = md$Pop2020/1000000
pf = pf[[1]] %>%
  filter(as.Date(time) %in% seq(as.Date(pf[[3]]),as.Date(pf[[4]]),1))
pf_fitted = pf_fitted[[1]] %>%
  select("time","country","observed") %>%
  filter(as.Date(time) %in% seq(as.Date(pf_fitted[[3]]),as.Date(pf_fitted[[4]]),1))
pf$HDI_2018 = rep(md$HDI_2018, each=length(unique(pf$time)))
pf$Pop2020 = rep(md$Pop2020, each=length(unique(pf$time)))
pf$Pop2020_Mi = rep(md$Pop2020/1000000, each=length(unique(pf$time)))
pf$SSA = rep(md$SSA, each=length(unique(pf$time)))


##Data Processing

#creating data vector
date_vector <- as.Date(unique(pf$time))

#adding final date observed cases to md
md$observed  <- pf$observed[as.Date(pf$time)==date_vector[length(date_vector)]]
md$rain <- pf$rain[as.Date(pf$time)==date_vector[length(date_vector)]]
md$temperature  <- pf$temperature[as.Date(pf$time)==date_vector[length(date_vector)]]
md$humidity  <- pf$humidity[as.Date(pf$time)==date_vector[length(date_vector)]]
md$stringency  <- pf$stringency[as.Date(pf$time)==date_vector[length(date_vector)]]
md$testing  <- pf$testing[as.Date(pf$time)==date_vector[length(date_vector)]]

#Layermap
LayerMapping =  c(
  "Reported Cases" = "observed",
  "SSA" = "SSA",
  "Population (millions)" = "Pop2020_Mi",
  "HDI in 2018" = "HDI_2018",
  "Rainfall (mm)" = "rain",
  "Temperature (\u00B0C)" = "temperature",
  "Humidity (g/kg)" = "humidity",
  "Strigency" = "stringency",
  "Testing" = "testing"
)

#creating afmap
afmap <- md
afmap_geo <- geojson_json(afmap, geometry = "polygon")

#creating fit (THIS DOES NOT EXIST)
all_fitted = readRDS(paste0("output/models/fitted_mean",suffix,".rds"))


all_preds = readRDS(paste0("output/models/predictions",suffix,".rds"))
all_preds$time <- as.Date(all_preds$time,"%Y-%m-%d")
all_preds = all_preds %>%
  rename(country = cname)
if (!any(grepl("type",colnames(all_preds)))){
  all_preds$type = matrix("IS",dim(all_preds)[1],1)
  pred_ind = logical(dim(all_preds)[1])
  for(i in 0:6){
    pred_ind = pred_ind | grepl(as.Date(last_date)-i,all_preds$time)
  }
  all_preds$type[pred_ind] = "OUS"
}

#color palette
cpal = cividis(256)

##Creating Interface

ui <- navbarPage("COVID-19 in Africa",
                tabPanel("Explore Data", 
                          sidebarLayout(             #creates sidebar of inputs and countries
                            sidebarPanel(
                              width = 3,
                              selectInput("layer", "Layer",
                                          LayerMapping,
                                          selected = "observed"),
                              selectInput("countries",
                                          "Countries",
                                          md$name, 
                                          selected = "Algeria",
                                          multiple = TRUE), 
                              h4(HTML(paste("This is an interactive companion to the results presented in ",
                                            tags$a(href="https://www.pnas.org/content/118/28/e2026664118", 
                                                              "Pan-African Evolution of Within and Between Country COVID-19 Dynamics"),
                                                       " (Ssentongo et al., 2021) which has been updated to include the the most recent",
                                                       "90 days of data. The code to reproduce",
                                            "these results is available at the",
                                            tags$a(href="https://github.com/Schiff-Lab/COVID19-HHH4-Africa", 
                                                   "Schiff Lab github.")))) 
                            ),
                            mainPanel(  #offsetting map and time series
                              fluidRow(column(11, offset = 0, leafletOutput("mymap"))),
                              fluidRow(column(11, offset = 0, plotOutput("timeseries")))),
                          )),
                 tabPanel("Model Summary", 
                          sidebarLayout(
                            "",
                            mainPanel(
                              strong("Table with model summary"),
                              tableOutput("table"),
                            )
                          )),
                 tabPanel("Predictions",
                          sidebarLayout(             #creates sidebar of inputs and countries
                            sidebarPanel(
                              width = 3,
                              selectInput("pred_countries",
                                          "Countries",
                                          md$name, 
                                          selected = "Algeria", 
                                          multiple=FALSE),
                            ),
                            mainPanel(  #offsetting map and time series
                              strong("Fit and 7-day prediction"),
                              plotlyOutput("fandpplot"),
                            br(),
                            strong("Within- and between-country contributions"),
                            plotlyOutput("wb"),
                          )
                          ))
)


# Define server logic
server <- function(input, output) {
  

  
  # Generate an HTML table view of the data ----
  output$table <- function(){
    fit <- readRDS(paste0("output/models/fitted_model_LAG7",suffix,".rds"))
    beta_hat <- fit$coefficients[1:16]
    sd_hat <- fit$se[1:16]
    zscores <- beta_hat / sd_hat
    pvalues <- 2 * pnorm(abs(zscores), lower.tail = F)
    pvalues <- as.character(ifelse(pvalues < 0.001, "< 0.001", round(pvalues, 3)))
    pvalues <- ifelse(nchar(pvalues) < 5, paste0(pvalues, 0), pvalues)
    pvalues <- tibble(Params = names(beta_hat), pvalues)
    
    tab <- readr::read_csv(paste0("output/tables/tab_params_LAG7",suffix,".csv"))
    tab <- tab %>%
      left_join(pvalues)
    tab <- tab[c(16, 9, 1:8, 15, 10:14, 17, 18), ]
    tab$pvalues[c(1, 2, 11, 17:18)] <- "-"
    tab[, 2:4] <- apply(tab[,2:4], 2,
                        function(x) ifelse(nchar(x) < 5, paste0(x, 0), x))
    tab$Params <- c("Intercept", "Intercept", "log(population)", "HDI",
"Landlocked", "Stringency(t-7)", "Testing(t-7)",
                    "Rain(t-7)", "Temperature(t-7)", "Humidity(t-7)",
                    "Intercept", "log(population)",
                    "HDI", "Landlocked", "Stringency(t-7)", "Testing(t-7)",
                   "rho", "psi")
    CI <- paste0("(", tab$`2.5 %`, ", ", tab$`97.5 %`, ")")
    tab$`97.5 %` <- NULL
    tab$`2.5 %` <- CI
    names(tab)[3] <- "CI"
    boldID <- which(tab$pvalues == "< 0.001" | as.numeric(tab$pvalues) <= 0.1)
    kbl(tab, col.names = c("Parameter", "Relative Risk", "95% CI", "p-value"),
        align = c("lccc") ) %>%
      kable_material(full_width = T, font_size = 20) %>%
      row_spec(0, bold = T) %>%
      row_spec(boldID, background = "#cdf7d6", color = "Black") %>%
      pack_rows("Endemic", 1, 1, background = "#e8e8e8") %>%
      pack_rows("Within-country", 2, 10, background = "#e8e8e8") %>%
      pack_rows("Between-country", 11, 16, background = "#e8e8e8") %>%
      pack_rows("", 17, 18, background = "#e8e8e8")
  }
  
  
  output$mymap <- renderLeaflet({
    print("in leaflet")
    leaflet() %>%
      addProviderTiles("CartoDB")
  })
  
  
  #This is called within the observer that changes the bounds of the map
  geom_vec <- reactive({
    if (is.null(input$countries)){
      st_coordinates(md$geom) 
    }
    else{
      st_coordinates(md$geom[is.element(md$name,input$countries)]) 
    }
  })
  
  output$wb <-renderPlotly({
    plot_contributions(pf_fitted, all_fitted, input$pred_countries)
  })
  
  
  #This is called when the prediction plot is generated -- calls plot_preds in functions.r
  output$fandpplot <-renderPlotly({
     ap = all_preds %>%
       filter(country == input$pred_countries)
     plot_preds(ap)
  })
  
  
  observeEvent(input$countries,{
    print("in leafbound")
    geom <- geom_vec()
    leafletProxy("mymap") %>%
      fitBounds(
        lat1 = min(geom[, 2]),
        lng1 = min(geom[, 1]),
        lat2 = max(geom[, 2]),
        lng2 = max(geom[, 1])
      )
  })
  
  #Leaflet Map
  observe({
    lidx <- which(colnames(md)==input$layer)
    
    if (input$layer=="observed" | input$layer=="Pop2020_Mi"){ # draw colors on a log scale for and population
      dmn <- log(md[[lidx]]+1)
      binnum = min(c(10,length(unique(dmn))))
      binlab = round(exp(seq(max(dmn),min(dmn),length.out=binnum+1))-1)
      labels = paste0(binlab[2:(binnum+1)],"-", binlab[1:binnum])
    }
    else if (input$layer=="SSA" | input$layer=="testing"){ # show single value for category
      dmn <- md[[lidx]]
      binnum = length(unique(dmn))
      labels = sort(unique(dmn),decreasing=T)
    }
    else if (input$layer=="humidity" | input$layer=="HDI_2018"){ # round to three decimal place ranges
      dmn <- md[[lidx]]
      binnum = min(c(10,length(unique(dmn))))
      binlab = round(seq(max(dmn),min(dmn),length.out=binnum+1),3)
      labels = paste0(binlab[2:(binnum+1)],"-", binlab[1:binnum])
    }
    else{ # round to integers range
      dmn <- md[[lidx]]
      binnum = min(c(10,length(unique(dmn))))
      binlab = round(seq(max(dmn),min(dmn),length.out=binnum+1))
      labels = paste0(binlab[2:(binnum+1)],"-", binlab[1:binnum])
      
    }
    

    pal <- colorBin(
      bins=binnum,  
      palette = cpal,
      domain = dmn, reverse = T, pretty=FALSE)
    
    pal_rev <- colorBin(
      bins=binnum,  
      palette = cpal,
      domain = dmn, reverse = F, pretty=FALSE)
    
    StrokeOp <- matrix(0,nrow(md),1)
    StrokeOp[match(input$countries,md$name)] = 1
    FillOp <- matrix(.2, nrow(md),1)
    FillOp[match(input$countries, md$name)] = .6
    StrokeWgt <- matrix(2,nrow(md),1)
    #StrokeWgt[match(input$countries,md$name)] = 3
    
    leafletProxy("mymap") %>%
      clearShapes() %>%
      removeControl("legend") %>%
      addPolygons(
        data=md$geom,
        weight=StrokeWgt,
        color= pal(dmn),
        opacity=StrokeOp,
        label = sprintf("<strong> %s </strong>
                  : %s",
                  as.character(md$name),
                  md[[lidx]]) %>%
          lapply(htmltools::HTML),
        fillColor=pal(dmn),
        fillOpacity = FillOp
      ) %>%
      
      #       addLegend(pal=pal_rev, values = dmn, opacity = 1,
      #      title = paste(names(LayerMapping[LayerMapping %in% input$layer]),"<br>on",pf$time[length(pf$time)]), 
      #     position="bottomright", layerId = "legend",
      #      labFormat = labelFormat(digits=lf.digits, transform=lf.transform, between=lf.between))
      addLegend(pal=pal_rev, values = dmn, opacity = 1,
                title = paste(names(LayerMapping[LayerMapping %in% input$layer]),"<br>on",pf$time[length(pf$time)]), 
                position="bottomright", layerId = "legend",
                labFormat = function(type,cuts,p) {
                  paste0(labels)
                })
    
    #addLegend(pal=pal(seq(min(dmn),max(dmn),length.out=10)), values = ~ seq(min(md[[lidx]]),max(md[[lidx]]),length.out=10),
    #     opacity = 1, title = "Est. GDP (2010)", position="bottomright", layerId = "legend")
  })
  
  
  #Time series
  
  output$timeseries <- renderPlot({
    #creating palette
    lidx <- which(colnames(md)==input$layer)
    
    if (input$layer=="observed" | input$layer=="Pop2020_Mi"){ # draw colors on a log scale for and population
      dmn <- log(md[[lidx]]+1)
      binnum = min(c(10,length(unique(dmn))))
      binlab = round(exp(seq(max(dmn),min(dmn),length.out=binnum+1))-1)
      labels = paste0(binlab[2:(binnum+1)],"-", binlab[1:binnum])
    }
    else if (input$layer=="SSA" | input$layer=="testing"){ # show single value for category
      dmn <- md[[lidx]]
      binnum = length(unique(dmn))
      labels = sort(unique(dmn),decreasing=T)
      
      
    }
    else if (input$layer=="humidity" | input$layer=="HDI_2018"){ # round to three decimal place ranges
      dmn <- md[[lidx]]
      binnum = min(c(10,length(unique(dmn))))
      binlab = round(seq(max(dmn),min(dmn),length.out=binnum+1),3)
      labels = paste0(binlab[2:(binnum+1)],"-", binlab[1:binnum])
    }
    else{ # round to integers range
      dmn <- md[[lidx]]
      binnum = min(c(10,length(unique(dmn))))
      binlab = round(seq(max(dmn),min(dmn),length.out=binnum+1))
      labels = paste0(binlab[2:(binnum+1)],"-", binlab[1:binnum])
      
    }
    
    pal <- colorBin(
      bins=binnum,  
      palette = cpal,
      domain = dmn, reverse = T, pretty=FALSE)
    
    
    ts_tmp <- pf %>%
      select(time, country, input$layer) %>%
      mutate(time=as.Date(time, format="%Y-%m-%d")) %>%
      filter(pf$country %in% input$countries)
    
    ts_last <- aggregate(ts_tmp, by=list(ts_tmp$country), FUN=last)
    
    
    
    country_label <- ts_tmp %>%
      filter(time %in% last(time))
    ggplot(data = ts_tmp,
           mapping = aes(x = time, group = country, color=country)
    ) +
      aes_string(y = input$layer) +
      ylab(names(LayerMapping[LayerMapping %in% input$layer])) +
      xlab("") +
      scale_x_date(date_labels = ("%b %Y")) +
      geom_line(size=1) +
      scale_color_manual(values=pal(dmn[md$name %in% input$countries]))+
      #thematic elements
      #ggtitle(paste(names(LayerMapping[LayerMapping %in% input$layer]), pf$time[1], "-", pf$time[length(pf$time)]))+
      geom_label_repel(data = country_label,
                       mapping = aes(x = time, group = country, label=country),
                       nudge_x = -10,
                       na.rm = TRUE) +
      theme(legend.position = "none",
            plot.title = element_text(size=20, face="bold", margin=margin(10,0,10,0)))
    #legend.justification=c(0,0), legend.position=c(-.1,0))
    
  })
}


# Run the application 
shinyApp(ui = ui, server = server)


