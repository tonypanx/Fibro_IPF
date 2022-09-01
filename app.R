library(shiny)
library(Seurat)
library(tibble)

fibro.markers <- readRDS("/project2/jhuangime/tony/Fibro_CAR/Human_IPF_Atlas_rebuilt/fibro_markers.rds")
ipf.markers <- readRDS("/project2/jhuangime/tony/Fibro_CAR/Human_IPF_Atlas_rebuilt/ipf_markers.rds") 
onlyfib.markers <- readRDS("/project2/jhuangime/tony/Fibro_CAR/Human_IPF_Atlas_rebuilt/onlyfib_markers.rds")
myofib.markers <- readRDS("/project2/jhuangime/tony/Fibro_CAR/Human_IPF_Atlas_rebuilt/myofib_markers.rds")

fib.merge <- readRDS("/project2/jhuangime/tony/Fibro_CAR/Human_IPF_Atlas_rebuilt/fib_merge.rds")
Idents(fib.merge) <- fib.merge$condition

fibro.markers <- add_column(fibro.markers, Gene = rownames(fibro.markers), .before = 1)
ipf.markers <- add_column(ipf.markers, Gene = rownames(ipf.markers), .before = 1)
onlyfib.markers <- add_column(onlyfib.markers, Gene = rownames(onlyfib.markers), .before = 1)
myofib.markers <- add_column(myofib.markers, Gene = rownames(myofib.markers), .before = 1)

ui <- fluidPage(
  titlePanel("IPF Marker Expression"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("var", 
                  label = "Choose a marker list to display",
                  choices = c("Fibro vs Non-fibro" = "1", 
                              "IPF vs Healthy, Fibroblasts Combined" = "2", 
                              "IPF vs Healthy, Regular Fibroblasts" = "3", 
                              "IPF vs Healthy, Myofibroblasts" = "4"),
                  selected = "Fibro vs Non-fibro"),
      
      column(8, textInput("gene", h3("Investigate Gene"), value = "COL1A1")),
      
      column(8, radioButtons('radio', h3("Select plot"), choices = c("IPF vs Healthy, all fibroblasts" = "1", "IPF vs Healthy, regular fibroblasts"="2", "IPF vs Healthy, myofibroblasts"="3"))),
      
      column(3,
             h3("Plot Gene"),
             actionButton("do", "Submit")),
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("DEG lists", dataTableOutput('table')),
        tabPanel("Gene comparison", dataTableOutput('table2')),
        tabPanel("Plot", plotOutput('plt'))
      )
    )

      

    )
)



server <- function(input, output, session) {

  output$table <- renderDataTable({
    if(input$var == "1"){
      fibro.markers
    }
    else if (input$var == "2"){
      ipf.markers
    } else if (input$var == "3"){
      onlyfib.markers
    } else{
      myofib.markers
    }
  }
  )
  
  temp <- data.frame()
  
  output$table2 <- renderDataTable(setNames(cbind(c("Fibro vs NonFibro", "IPF Fib vs Healthy Fib", "IPF Reg Fib vs Healthy Reg Fib", "IPF Myofib vs Healthy Myofib"), 
                                         rbind(fibro.markers[input$gene, c(3:6)], 
                                         ipf.markers[input$gene, c(3:6)], 
                                         onlyfib.markers[input$gene, c(3:6)], 
                                         myofib.markers[input$gene, c(3:6)])), c("Gene List", "Avg_Log2FC", "% Expr in A", "% Expr in B", "Adj. P val")))
  
  
  gene <- reactiveVal()
  
  observeEvent(input$do,{
    gene(input$gene)
    
  })
  
  output$plt <- renderPlot({
    if (input$radio == '1'){
      VlnPlot(fib.merge, features = gene(), pt.size = .1)
    }else if (input$radio == '2'){
      VlnPlot(subset(fib.merge, celltype == "Fibroblast"), features = gene(), pt.size = .1)
    }else{
      VlnPlot(subset(fib.merge, celltype == "Myofibroblast"), features = gene(), pt.size = .1)
    }
  })

}

shinyApp(ui = ui, server = server)



