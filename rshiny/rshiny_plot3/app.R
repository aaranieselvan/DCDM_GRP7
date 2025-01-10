# Plot 3.3 heatmap 

library(shiny)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Import cleaned data file
data <- read.csv("~/Desktop/DCDM_GRP7/outputs/clean_final_data.csv")
print(nrow(data))  # Check the total number of rows (should be 140931)
View(data)

# Step 1: Filter to retain statistically significant p-values (p <= 0.05)
filtered_p_value_data1 <- data %>%
  mutate(pvalue = as.numeric(pvalue)) %>%  # Convert the 'pvalue' column to numeric
  filter(pvalue > 0 & pvalue <= 0.05) %>%  # Retain rows where 'pvalue' is greater than 0 and less than or equal to 0.05
  filter(!is.na(pvalue) & is.finite(pvalue)) # Remove rows with NA or non-finite p-values
print(nrow(filtered_p_value_data1))  # Check the filtered number of rows (should be 910) 

# Step 2: Apply FDR correction (Benjamini-Hochberg procedure)
negative_or_zero_fdr <- filtered_p_value_data1 %>%
  mutate(fdr = p.adjust(pvalue, method = "BH"))  # Add the FDR-adjusted p-values

# Step 3: Check for negative or zero FDR values
negative_or_zero_fdr_check <- negative_or_zero_fdr %>%
  filter(fdr <= 0)  # Filter out rows with FDR <= 0

if (nrow(negative_or_zero_fdr_check) > 0) {
  stop("Error: There are negative or zero FDR values. Please investigate the data before proceeding.")
} else {
  print("FDR values are valid (none are negative or zero). Proceeding with log transformation.")
}

#Output = FDR values are valid (none are negative or zero). Proceeding with log transformation.

# Step 4: Log-transform FDR-adjusted p-values (instead of raw p-values)
log_transformed_FDR <- negative_or_zero_fdr %>%
  mutate(log_fdr = -log10(fdr)) %>%  # Log-transform the FDR-adjusted p-values
  filter(!is.na(log_fdr) & is.finite(log_fdr))  # Remove rows where log_fdr is NA or non-finite

print(nrow(log_transformed_FDR))  # Check the number of rows after transformation (should be 910)

# Combine multiple identical gene_symbols (who have different p values) into a single value for clustering
clustering_data <- log_transformed_FDR %>%
  group_by(gene_symbol) %>%  # Group data by unique 'gene_symbol'
  summarise(log_fdr = mean(log_fdr, na.rm = TRUE),  # Calculate the mean log FDR value for each gene
            mouse_strain = first(mouse_strain),  # Use first strain as a proxy for all strains
            mouse_life_stage = first(mouse_life_stage),  # Use first life stage as a proxy for all stages
            parameter_name = first(parameter_name)) %>%  # Use first parameter name as a proxy for all phenotypes
  ungroup() %>%
  arrange(desc(log_fdr))  # Sort in descending order by log FDR
print(nrow(clustering_data))  # Check the number of unique gene_symbols (should be 123)

# Extract top 25 most significant gene_symbols based on FDR
top_25_genes <- clustering_data %>%
  top_n(25, log_fdr)

# Define user interface with number of clusters, gene selection & download heatmap button
ui <- fluidPage(
  
  titlePanel("Hierarchical Clustering of Mouse Gene_Symbols Based on Normalised Log FDR"),
  
  sidebarLayout(
    
    sidebarPanel(
      sliderInput("clusters", "Number of Clusters:", min = 2, max = 10, value = 3),
      
      selectizeInput("selected_genes", 
                     label = "Select Gene:", 
                     choices = clustering_data$gene_symbol, 
                     selected = top_25_genes$gene_symbol, 
                     multiple = TRUE, 
                     options = list(
                       create = TRUE, 
                       maxItems = 150, 
                       placeholder = 'Start typing a gene symbol')),
      
      textOutput("selected_genes_count"),
      actionButton("reset_button", "Reset to Top 25 gene_symbols"),
      actionButton("focus_button", "Focus on Specific gene_symbols (tuft1, sun5, mettl5)"),
      downloadButton("download_heatmap", "Download Heatmap"),
      downloadButton("download_data", "Download Cluster Data"),
      uiOutput("cluster_gene_symbols")
    ),
    
    mainPanel(
      plotOutput("heatmap", height = "700px", width = "750px")
    )
  )
)

# Server logic to handle reactive inputs and clustering
server <- function(input, output, session) {
  
  selected_genes_reactive <- reactiveVal(top_25_genes$gene_symbol)
  
  selected_genes_label <- reactive({
    num_selected <- length(input$selected_genes)
    paste("Select Gene (", num_selected, " selected) :", sep = "")
  })
  
  observe({
    updateSelectInput(session, "selected_genes", label = selected_genes_label())
  })
  
  observeEvent(input$reset_button, {
    selected_genes_reactive(top_25_genes$gene_symbol)
    updateSelectizeInput(session, "selected_genes", selected = top_25_genes$gene_symbol)
  })
  
  observeEvent(input$focus_button, {
    selected_genes_reactive(c("tuft1", "sun5", "mettl5"))
    updateSelectizeInput(session, "selected_genes", selected = c("tuft1", "sun5", "mettl5"))
  })
  
  reactive_heatmap <- reactive({
    selected_genes <- selected_genes_reactive()  # Get selected genes
    selected_genes_data <- clustering_data %>%
      filter(gene_symbol %in% selected_genes)
    
    clustering_matrix <- as.matrix(selected_genes_data[, c("log_fdr")])  # Use log-transformed FDR values for clustering
    rownames(clustering_matrix) <- selected_genes_data$gene_symbol
    
    clustering_matrix <- scale(clustering_matrix)  # Normalise log_fdr values
    
    colnames(clustering_matrix) <- "Normalised Log FDR"
    
    distance_matrix <- dist(clustering_matrix, method = "euclidean")
    hc <- hclust(distance_matrix, method = "complete")
    
    kmeans_result <- cutree(hc, k = input$clusters)
    
    annotation_row <- data.frame(Cluster = factor(kmeans_result),
                                 Mouse_Strain = selected_genes_data$mouse_strain,
                                 Mouse_Life_Stage = selected_genes_data$mouse_life_stage,
                                 Parameter_Name = selected_genes_data$parameter_name)
    rownames(annotation_row) <- rownames(clustering_matrix)
    
    cluster_levels <- unique(kmeans_result)
    cluster_colors <- brewer.pal(min(length(cluster_levels), 12), "Set3")
    names(cluster_colors) <- as.character(cluster_levels)
    annotation_colors <- list(Cluster = cluster_colors)
    
    cluster_genes <- lapply(cluster_levels, function(cluster_id) {
      genes_in_cluster <- rownames(clustering_matrix)[kmeans_result == cluster_id]
      paste(genes_in_cluster, collapse = ", ")
    })
    names(cluster_genes) <- paste("Cluster", cluster_levels)
    
    output$cluster_gene_symbols <- renderUI({
      tagList(
        lapply(names(cluster_genes), function(cluster_name) {
          div(
            strong(cluster_name),
            ": ",
            cluster_genes[[cluster_name]]
          )
        })
      )
    })
    
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
      color = colorRampPalette(c("cadetblue1", "plum1", "maroon2"))(50),
      fontsize_row = fontsize_row,
      fontsize_col = 10,
      display_numbers = TRUE,
      fontsize_number = fontsize_number,
      treeheight_row = 100,
      treeheight_col = 100,
      main = paste("Clustering of Genes by ", input$clusters, " Clusters", sep = ""),
      height = plot_height,
      width = plot_width,
      angle_col = 90,
      legend = TRUE,
      legend_title = "Significance Level",
      legend_position = "bottomright",
      legend_direction = "horizontal",
      annotation_legend = TRUE,
      annotation_legend_side = "right"
    )
  })
  
  output$heatmap <- renderPlot({
    reactive_heatmap()
  })
  
  observe({
    selected_genes_reactive(input$selected_genes)
  })
  
  output$selected_genes_count <- renderText({
    num_genes <- length(input$selected_genes)
    paste("Number of Selected Genes: ", num_genes)
  })
  
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
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("cluster_data.csv")
    },
    content = function(file) {
      selected_genes <- selected_genes_reactive()
      selected_genes_data <- clustering_data %>%
        filter(gene_symbol %in% selected_genes)
      
      kmeans_result <- cutree(hclust(dist(scale(as.matrix(selected_genes_data[, c("log_fdr")])))), k = input$clusters)
      
      selected_genes_data$Cluster <- kmeans_result
      
      write.csv(selected_genes_data[, c("gene_symbol", "mouse_strain", "mouse_life_stage", "parameter_name", "log_fdr", "Cluster")], file, row.names = FALSE)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)