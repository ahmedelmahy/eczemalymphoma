library(shiny)
shiny_pheatmap <- function(datalist){
    ui <- fluidPage(
        sidebarLayout(
            sidebarPanel(
                # select variable to plot
                sliderInput("number_of_clusters_cytokines", "number of clusters cytokines:", 3, min = 1, max = 20),
                sliderInput("number_of_clusters_patients", "number of clusters persons:", 3, min = 1, max = 20),
                selectInput(inputId = "datasets",label="What data to include ? ",choices = names(datalist), multiple = TRUE)
            ),
            mainPanel(
                # plot ggplot2 out
                plotOutput("plotDisplay",height = "670px"),
                textOutput("dims")

            )))

    server <- function(input, output) {

        df <- reactive({as.data.frame(do.call("rbind",datalist[input$datasets]))})
        output$plotDisplay <- renderPlot({
            annotation_row2 <- data.frame(status = df()$class)
            rownames(annotation_row2) = rownames(df())
            pheatmap(df()[,1:dim(df())[2]-1], cutree_cols = input$number_of_clusters_cytokines,
                     cutree_rows = input$number_of_clusters_patients,  annotation_row  = annotation_row2, fontsize_row = 3)
        }
        )
        output$dims <- renderPrint({
            print(dim(df()))
        })
    }
    shinyApp(ui = ui, server = server)
}
