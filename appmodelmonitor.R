library(shiny)
library(magrittr)
library(tidyverse)
library(lubridate)
library(scales)
library(NMF)
library(plotly)
library(TSclust)
library(dynamicTreeCut)
library(fdasrvf)
library(lmerTest)
library(ecp)


setwd("~/Downloads/apply/shiny/present")
load("data.RData")

setwd("~/Downloads/apply/shiny/present")
source("monitoring.R")

app.name <- "Monthly Energy Consumptions"


building.choice <- c()
for(s in unique(df$Facility_Name)){
  building.choice <- c(building.choice, 
                       unlist(str_split(unlist(str_split(
                         unlist(str_split(s, " "))[1], 
                         "-"))[1], ","))[1] )
}
building.choice <- sort(unique(toupper(building.choice)))
function.type <- unique(df$AccountType)
# ui <- fluidPage(
#   titlePanel(app.name),
#   
#   navlistPanel(
#     "Header",
#     tabPanel("EDA",
#              h3("Missing Values"),
#              uiOutput("Missing")
#              ),
#     tabPanel("Monitoring",
#              h3("Monthly Monitoring"),
#              uiOutput("body")
#              )
#   )
# )

ui <- navbarPage(app.name,
             tabPanel("EDA",
                      h3("Data Description"),
helpText("The data were collected by the University’s Facilities Operations division 
for Utility Operations and Energy Management. Monthly gas consumption was measured by the 
 utility provider in 238 separate accounts located across 115 buildings. The number of 
 accounts varies between buildings. In some cases, one physical building contains multiple 
 accounts, for example, multiple individual units within a single apartment building. The 
energy consumption is measured in units of hundred cubic feet(CCF). Each observation was 
collected form an individual utility invoice for a single account in a single billing period,
and measures the difference between two meter readings."),   
                      h3("Missing Values"),
                      uiOutput("Missing"),
                      h3("Visualization"),
                      uiOutput("body_visual"),
                      column(3, selectInput("bldgeda", "Building Type",
                                            choices=building.choice,
                                            selected = "NORTHWOOD",
                                            width="100%"),
                             verbatimTextOutput("plotlyaccinfo")
                             )
                      ),
             tabPanel("Modeling",
helpText("It has two parts, the clustering on different types of buildings and the modeling 
for doing forecasting and change point detection."),
                      fluidRow(
                      h3("Hierarchical Clustering using Karchar Mean"),
helpText("Since the time series have different lengths, starting at different time point, and 
have different scales (shape>scale). Therefore, we use the Karchar mean.
The Karchar mean is essentially embedding the monthly energy consumptions into the 
functional domain and find the Frechet mean on each functional(account). Each Karchar
mean from the account ends up total 12 value points, which is used as features for clustering."),
                      uiOutput("body_clustering"),
                      # column(3, selectInput("n.clu", "# Cluster",
                      #                       choices=2:8,
                      #                       selected = "6",
                      #                       width="100%")),
                      # actionButton(inputId="cluster.go", label="Refresh Clusters"),
                      column(3, selectInput("clu", "Cluster",
                                            choices=1:6,
                                            selected = "2",
                                            width="100%")),
                      column(6,  tableOutput("clustertext"))
                      ),
                      fluidRow(
                      h3("Forecasting and Change Point"),
helpText("The estimation model is using the mixed ARIMAX model WITHIN each cluster, specifically 
CCF~(L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12++0)*(CCF1+CCF6+CCF12)+avgt+snow+cdd +
(1|id)+(0+CCF1|id)+(0+CCF6|id)+(0+CCF12|id). There are three lags/delay consumptions, chosen by 
the ACF and PACF plots, monthly indicators, and three weather related data.
The forecasting is done for the most recent whole year and the change point detection is done on
the historical years. 

The change point detection is using the E-statistic, which is essentially a goodness-of-fit 
measure of difference between segments. There are two types of change points, one is using the
CCF consumptions, the other one is using the residuals from the estimated model. Specifically,
the red line is using the CCF and blue dotted line is using the residuals."),
                      uiOutput("body_ForcastChange")
                      )
                      ),
             tabPanel("Monitoring",
helpText("Monthly monitoring is for the incoming data. It is the ON-GOING project. The motivation
is that our clients want to know periodically which account is behaving out-of-form, in a sense 
that consuming too much energy. The idea is using a regression on the `in-control` historical data
and use estimates to predict the future. Then use the predictions to decide whether the incoming
data is out-of-form in a statistical hypothesis testing framework. 
But the historical data may not be good enough either and so the estimates will be refined by 
using a Bayesian analysis framework."),
                      h3("Monthly Monitoring"),
                      column(3, selectInput("bldg", "Building Type",
                                            choices=building.choice,
                                            selected = "NORTHWOOD",
                                            width="100%")),
                      uiOutput("body")
                      )
)

server <- function(input, output, session) {
  storage <- reactiveValues()
  observeEvent(c(input$bldg, input$iter),{
    req(input$bldg)
    storage$df <- df[apply(df%>% select(Facility_Name), 1, FUN = function(x){
      unlist(str_split(x, " "))[1] == input$bldg
    }), ]
  })
  # 
  observeEvent(c(input$bldg, input$iter, input$monitor.year),{
    req(storage$df, input$bldg, input$monitor.year, input$iter)
    
    monitor.year <- as.numeric(input$monitor.year)
    iter <- as.numeric(input$iter)
    storage$df.test <- storage$df %>% filter(Year >= monitor.year - iter,
                                             Year < monitor.year)
    storage$df.future <- storage$df %>% filter(Year == monitor.year) %>%
      select(Date, CNGAccount, CCF) %>% spread(CNGAccount, CCF)
    storage$df.future <- storage$df.future[, -1]

  }, ignoreNULL = F)
  
  output$Missing <- renderUI({
    plotOutput("Missing.heapmap", height="600px")
  })
  output$Missing.heapmap <- renderPlot({
    par(mfrow = c(2,1))
    aheatmap(1-is.na(cng), Rowv=NA, Colv=NA,
             color="Greys", main = "Data Pattern")
    aheatmap(cng, Rowv=NA, Colv=NA, color="Reds", main = "Value Pattern")
  })
  
  output$body_visual <- renderUI({
    tabPanel("",
             fluidRow(
                  plotlyOutput("plotly.accvisual", height="600px")
             )
    )
  })
  output$plotlyaccinfo <- renderPrint({paste0("There are total ",
                                             length(unique(storage$df$Meter)),
                                             " accounts in ", 
                                             input$types)})
  
  output$plotly.accvisual <- renderPlotly({
    df1 <- df[apply(df%>% select(Facility_Name), 1, FUN = function(x){
      unlist(str_split(x, " "))[1] == input$bldgeda
    }), ]
    df1$READ.TO <- ymd(df1$READ.TO)
    p <- plot_ly(df1, x = ~ READ.TO, y = ~ CCF, split = df1$Meter) %>%
      add_lines() %>%
      layout(title = paste(input$types, "Acc/Meters"),
             xaxis = list(title = 'Time',
              rangeselector=list(
                buttons = list(
                    list(
                      count = 6,
                      label = "6 mo",
                      step = "month",
                      stepmode = "backward"),
                    list(
                      count = 1,
                      label = "1 yr",
                      step = "year",
                      stepmode = "backward"),
                    list(
                      count = 1,
                      label = "YTD",
                      step = "year",
                      stepmode = "todate"),
                    list(step = "all")
              )),
              rangeslider = list(type = "READ.TO")),
             yaxis = list (title = 'CCF'))
    p$elementId <- NULL
    
      p
  })
  
  
  output$body_clustering <- renderUI({
    tabPanel("",
             fluidRow(
          plotOutput("clustering", height="500px", width = "95%")
             )
    )
  })
  
  # observeEvent(c(input$cluster.go), {
  #   n.clu = as.numeric(input$n.clu)
  #   storage$fmean.memb <- cutree(fmean.hc, k = n.clu)
  # }, ignoreNULL = FALSE)
  
  output$clustertext <- renderTable({
    tb <- data.frame(t(t(table(fmean.memb))))[, -2]
    colnames(tb) <- c("cluster", "count")
    tb
  })
  
  
  output$clustering <- renderPlot({
    clu.use = as.numeric(input$clu)
    if(sum(fmean.memb==clu.use) > 20){
      ggplot(data = data.use %>% filter(Type == "CCF",
                                        clu == clu.use),
             aes(x = ymd(READ.TO), y = val, color = factor(Meter))) +
        geom_line() + xlab("Time")
    }else{
      ggplot(data = data.use %>% filter(Type == "CCF",
                                        clu == clu.use),
             aes(x = ymd(READ.TO), y = val)) +
        geom_line() + xlab("Time")+ facet_wrap(~Meter, scales = "free")
    }
  })
  
  output$body_ForcastChange <- renderUI({
    tabPanel("",
             fluidRow(
        plotOutput("ForcastChange", height="600px")
             )
    )
  })
  output$ForcastChange <- renderPlot({
    # dev.off()
    series.use <- 1:16
    ggplot(data = data.use %>% filter(Meter %in% Meters[series.use]),
           aes(x=ymd(READ.TO), y = val, color = Type)) +
      geom_line() + facet_wrap(~Meter.clu, scales = "free") + xlab("Year") +
      geom_vline(data = change.use%>% filter(Meter %in% Meters[series.use]),
                                             aes(xintercept = change), color = "red") +
      geom_vline(data = change.use%>% filter(Meter %in% Meters[series.use]),
                     aes(xintercept = res.change), color = "blue", linetype = "dotted")
  })
  
  
  output$body <- renderUI({
    req(storage$df)
    # tabsetPanel(
      # tabPanel("Plots",
      #          plotOutput("plot_tag", height="750px")),
      tabPanel("",
               # h2("Model Summary", style="text-align: center; margin-bottom: 30px;"),
               column(3, selectInput("iter", "# Historical Years",
                                      # min=2, max=(length(unique(storage$df$Year))-1),
                                      # step = 1, value = 3
                                      2:(length(unique(storage$df$Year))-1),
                                     selected = 3
                                      )),
               column(3, selectInput("monitor.year", "The year to monitor",
                                     unique(storage$df$Year)[-c(1,2)],
                                     selected = 2018
                                     )
                      ),
                 column(6,  actionButton("evReactiveButton", "Start to monitor"),
                        verbatimTextOutput("model_summary")),
               fluidRow(
                  plotOutput("ggmonitor.building", height="600px")
               )
            # )
        )
  })
  output$ggmonitor.building <- renderPlot({
    monitor.year <- as.numeric(input$monitor.year)
    iter <- as.numeric(input$iter)
    p.null <- is.null(storage$p.alpha)
    p.true <- F
    if(!p.null)
      p.true <- storage$p.alpha
    
    if(input$evReactiveButton[1] & !p.null & p.true ){
      df1 <- storage$df
      df1$lbl <- "in-control"
      df1$lbl[df1$CNGAccount %in% names(storage$p.each)[storage$p.each]] <- "out-control"
      ggplot(data = df1 %>% filter(Year >= monitor.year - iter, Year <= monitor.year),
             aes(x = Year + (Month-1)/12, y = CCF, color = factor(Meter))) +
        geom_line() + xlab("Year") + facet_grid(~lbl)
    }else{
      ggplot(data = storage$df %>% filter(Year >= monitor.year - iter,
                                          Year <= monitor.year), 
             aes(x = Year + (Month-1)/12, y = CCF, color = factor(Meter))) +
        geom_line() + xlab("Year")
    }
  })
  etext <- eventReactive(input$evReactiveButton, {
    req(storage$df.test, storage$df.future, input$iter, input$monitor.year)
    
    monitor.year <- as.numeric(input$monitor.year)
    iter <- as.numeric(input$iter)
    out <- real.monitoring(storage$df.test, storage$df.future, iter, monitor.year)
    storage$coef.sigma <- out[[1]]
    storage$p.alpha <- out[[2]]
    if(storage$p.alpha){
      print("There is at least one acc `out-of-form`")
    }else{
      print("All acc `in-control`")
    }
    
    storage$p.each <- out[[3]]
  })
  output$model_summary <- renderPrint({
    etext()
  })
}

shinyApp(ui, server)
1
