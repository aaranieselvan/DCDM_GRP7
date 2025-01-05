library(shiny)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#Plot3:Plot 3 shows hierarchical clustering of gene symbols based on log-transformed 
#p-values, with the top 25 most significant genes displayed by default. Users can adjust the 
#number of clusters and select specific genes to explore their clustering behavior, with a color-coded heatmap and a list of genes in each cluster.

# Import cleaned data file
data <- read.csv("~/Desktop/DCDM_GRP7/outputs/clean_final_data.csv")
print(nrow(data))  # Check the total number of rows (should be 140931)

#Filter to  retain statistically significant pvalues (p <= 0.05) 
#log transform significant pvalues to prepare for clustering of gene_symbols 
#and  exclude all infinite or missing log_p_values
filtered_p_value_data <- data %>%
  mutate(pvalue = as.numeric(pvalue)) %>%
  filter(pvalue > 0 & pvalue <= 0.05) %>%
  mutate(log_p_value = -log10(pvalue)) %>%
  filter(!is.na(log_p_value) & is.finite(log_p_value))
print(nrow(filtered_p_value_data))  # Check the filtered number of rows (should be 910)

# Combine multiple identical gene_symbols (who have different p values) into a single value for clustering
clustering_data <- filtered_p_value_data %>%
  group_by(gene_symbol) %>%
  summarise(log_p_value = mean(log_p_value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(log_p_value))  # Sort in descending order by log_p_value
print(nrow(clustering_data))  # Check the number of unique gene_symbols (should be 123)
View(clustering_data)

# Extract top 25 most significant gene_symbols, these will be used as default plot heatmap as they are most significant
top_25_genes <- clustering_data %>%
  top_n(25, log_p_value)

# Define user interface with number of clusters, gene selection & download heatmap button
ui <- fluidPage(
  titlePanel("Hierarchical Clustering of Gene_Symbols Based normalised log p-value"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("clusters", "Number of Clusters:", min = 2, max = 10, value = 3),
      selectizeInput("selected_genes", label = "Select gene_symbols:",
                     choices = clustering_data$gene_symbol,
                     selected = top_25_genes$gene_symbol, 
                     multiple = TRUE,
                     options = list(create = TRUE, maxItems = 100, placeholder = 'Start typing a gene symbol')),
      textOutput("selected_genes_count"),
      actionButton("reset_button", "Reset to Top 25 gene_symbols"),
      downloadButton("download_heatmap", "Download Heatmap"),
      # Show gene symbols for each cluster
      uiOutput("cluster_gene_symbols")
    ),
    mainPanel(
      plotOutput("heatmap", height = "750px", width = "750px")
    )
  )
)

# Server logic to handle reactive inputs and clustering
server <- function(input, output, session) {
  
  # Reactive value to store the selected genes (default to top 25)
  selected_genes_reactive <- reactiveVal(top_25_genes$gene_symbol)
  
  # Reactive to dynamically update label for selectInput
  selected_genes_label <- reactive({
    num_selected <- length(input$selected_genes)
    paste("Select gene_symbols (", num_selected, " selected):", sep = "")
  })
  
  # Update label for gene selection input dynamically
  observe({
    updateSelectInput(session, "selected_genes", label = selected_genes_label())
  })
  
  # Observe the reset button to reset selections to the top 25 genes
  observeEvent(input$reset_button, {
    selected_genes_reactive(top_25_genes$gene_symbol)
    updateSelectizeInput(session, "selected_genes", selected = top_25_genes$gene_symbol)
  })
  
  # Reactive to update heatmap based on selected genes
  reactive_heatmap <- reactive({
    selected_genes <- selected_genes_reactive()  # Get selected genes
    
    # Ensure that only genes from clustering_data are used
    selected_genes_data <- clustering_data %>%
      filter(gene_symbol %in% selected_genes)
    
    # Convert dataframe into a matrix for clustering, set row names to gene_symbols
    clustering_matrix <- as.matrix(selected_genes_data[, -1])
    rownames(clustering_matrix) <- selected_genes_data$gene_symbol
    
    # Normalize log_p_values in matrix into normalised log p values
    clustering_matrix <- scale(clustering_matrix)
    
    # Rename the column to normalised log p-value)
    colnames(clustering_matrix) <- "normalised log p-value"
    
    #calculate pairwise euclidean distances between rows of matrix which is used for
    #hierarchical clustering to construct a dendogram which visually represents 
    #clustering of multiple gene_symbols
    distance_matrix <- dist(clustering_matrix, method = "euclidean")
    hc <- hclust(distance_matrix, method = "complete")
    
    # Dendrogram is cut into clusters "k" which is determined by the input slider
    kmeans_result <- cutree(hc, k = input$clusters)
    
    # Dataframe contains clusters with their assigned gene_symbols
    annotation_row <- data.frame(Cluster = factor(kmeans_result))
    rownames(annotation_row) <- rownames(clustering_matrix)
    
    # Define cluster colours
    cluster_levels <- unique(kmeans_result)
    cluster_colors <- brewer.pal(min(length(cluster_levels), 12), "Set3")
    names(cluster_colors) <- as.character(cluster_levels)
    annotation_colors <- list(Cluster = cluster_colors)
    
    # Display the gene_symbols in each cluster
    cluster_genes <- lapply(cluster_levels, function(cluster_id) {
      genes_in_cluster <- rownames(clustering_matrix)[kmeans_result == cluster_id]
      paste(genes_in_cluster, collapse = ", ")
    })
    names(cluster_genes) <- paste("Cluster", cluster_levels)
    
    # Update UI to display gene_symbols (genes) per cluster
    output$cluster_gene_symbols <- renderUI({
      tagList(
        lapply(names(cluster_genes), function(cluster_name) {
          div(
            strong(cluster_name), 
            ": ", cluster_genes[[cluster_name]]
          )
        })
      )
    })
    
    # Dynamically adjust font size and plot dimensions
    num_genes <- length(selected_genes)
    fontsize_row <- ifelse(num_genes < 20, 15, 10)
    fontsize_number <- ifelse(num_genes < 20, 15, 10)
    plot_height <- ifelse(num_genes < 20, 750, 100 + num_genes * 20)
    plot_width <- ifelse(num_genes < 20, 750, 100 + num_genes * 20)
    
    pheatmap(
      clustering_matrix,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      color = colorRampPalette(c("plum1", "maroon2", "cadetblue1"))(50),
      fontsize_row = fontsize_row,
      fontsize_col = 12,  # Adjust size for better readability
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = fontsize_number,
      treeheight_row = 100,
      treeheight_col = 100,
      main = paste("Clustering of genes by ", input$clusters, " clusters", sep = ""),
      height = plot_height,
      width = plot_width,
      angle_col = 0  # This will force the column label to remain horizontal (default)
    )
  })
  
  # Render the heatmap
  output$heatmap <- renderPlot({
    reactive_heatmap()
  })
  
  # Update the selected genes based on user input
  observe({
    selected_genes_reactive(input$selected_genes)
  })
  
  # Display the count of selected gene_symbol data
  output$selected_genes_count <- renderText({
    num_genes <- length(input$selected_genes)
    paste("Number of selected gene_symbols: ", num_genes)
  })
  
  # Download the heatmap as PNG
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("heatmap.png")
    },
    content = function(file) {
      png(file, width = 750, height = 750)
      print(reactive_heatmap())
      dev.off()
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)