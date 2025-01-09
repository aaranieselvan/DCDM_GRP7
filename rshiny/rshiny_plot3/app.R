#Heatmap3.3 

library(shiny)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Import cleaned data file
data <- read.csv("~/Desktop/DCDM_GRP7/outputs/clean_final_data.csv")
print(nrow(data))  # Check the total number of rows (should be 140931)

# Filter to retain statistically significant pvalues (p <= 0.05), by performing a series of transformation steps
filtered_p_value_data <- data %>%
  mutate(pvalue = as.numeric(pvalue)) %>% # Convert the 'pvalue' column to numeric 
  filter(pvalue > 0 & pvalue <= 0.05) %>% # Retain rows where the 'pvalue' is greater than 0 and less than or equal to 0.05
  mutate(log_p_value = -log10(pvalue)) %>% # Create a new column, 'log_p_value', where p values have been negatively log10 transformed
  filter(!is.na(log_p_value) & is.finite(log_p_value)) # Remove rows where log_p_value is NA or non-finite
print(nrow(filtered_p_value_data))  # Check the filtered number of rows (should be 910)

# Combine multiple identical gene_symbols (who have different p values) into a single value for clustering
clustering_data <- filtered_p_value_data %>%
  group_by(gene_symbol) %>% # Group data by unique 'gene_symbol'
  summarise(log_p_value = mean(log_p_value, na.rm = TRUE), # Calculate the mean log p value for each gene, ignoring missing values
            mouse_strain = first(mouse_strain),  # Use first strain as a proxy for all strains
            mouse_life_stage = first(mouse_life_stage),  # Use first life stage as a proxy for all stages
            parameter_name = first(parameter_name)) %>% # Use first parameter name as a proxy for all phenotypes
  ungroup() %>%
  arrange(desc(log_p_value))  # Sort in descending order by log_p_value
print(nrow(clustering_data))  # Check the number of unique gene_symbols (should be 123)
View(clustering_data) # Display resulting data frame 

# Extract top 25 most significant gene_symbols
top_25_genes <- clustering_data %>%
  top_n(25, log_p_value)

# Define user interface with number of clusters, gene selection & download heatmap button
ui <- fluidPage(
  
  # Application title displayed at the top 
  titlePanel("Hierarchical Clustering of Mouse Genes Based on Phenotypic Significance Associations"),
  
  # Define the layout with a sidebar and a main panel 
  sidebarLayout(
    
    # Sidebar panel for user input 
    sidebarPanel(
      
      # Slider for selecting the number of clusters to generate in the heatmap
      sliderInput("clusters", "Number of Clusters:", min = 2, max = 10, value = 3),
      
      # Dropdown menu for selecting genes, with multiple selection enabled
      selectizeInput("selected_genes", # Input data
                     label = "Select Gene:", # Label for the dropdown menu 
                     choices = clustering_data$gene_symbol, # Provide available gene symbols from the dataset
                     selected = top_25_genes$gene_symbol, # Pre-select the top 25 genes
                     multiple = TRUE, # Allow multiple genes to be selected 
                     options = list( # Additional options for the dropdown
                       create = TRUE, # Allow users to add custom gene symbols not in the list
                       maxItems = 150, # Limit the maximum number of selections to 150
                       placeholder = 'Start typing a gene symbol')), # Placeholder text for the dropdown
      
      # Display the count of selected genes dynamically 
      textOutput("selected_genes_count"),
      
      # Button to reset the gene selection to the top 25 gene symbols
      actionButton("reset_button", "Reset to Top 25 gene_symbols"),
      
      # Button to focus on specific genes of interest (tuft1, sun5, mettl5)
      actionButton("focus_button", "Focus on Specific Genes (tuft1, sun5, mettl5)"),
      
      # Button to download the heatmap as an image
      downloadButton("download_heatmap", "Download Heatmap"),
      
      # Button to download the clustering data as a file
      downloadButton("download_data", "Download Cluster Data"),
      
      # Dynamic UI output to display gene symbols in each cluster
      uiOutput("cluster_gene_symbols")
    ),
    
    # Main panel for displaying the heatmap 
    mainPanel(
      
      # Format the heatmap plot with specified dimensions
      plotOutput("heatmap", height = "700px", width = "750px")
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
    paste("Select Gene (", num_selected, " selected) :", sep = "")
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
  
  # Observe the focus button to pre-select tuft1, sun5, and mettl5
  observeEvent(input$focus_button, {
    selected_genes_reactive(c("tuft1", "sun5", "mettl5"))
    updateSelectizeInput(session, "selected_genes", selected = c("tuft1", "sun5", "mettl5"))
  })
  
  # Reactive to update heatmap based on selected genes
  reactive_heatmap <- reactive({
    selected_genes <- selected_genes_reactive()  # Get selected genes
    
    # Ensure that only genes from clustering_data are used
    selected_genes_data <- clustering_data %>%
      filter(gene_symbol %in% selected_genes)
    
    # Convert dataframe into a matrix for clustering, set row names to gene_symbols
    clustering_matrix <- as.matrix(selected_genes_data[, c("log_p_value")])  # We use log_p_value for clustering
    rownames(clustering_matrix) <- selected_genes_data$gene_symbol
    
    # Normalise log_p_values in matrix into normalised log p values
    clustering_matrix <- scale(clustering_matrix)
    
    # Rename the column to normalised log p-value
    colnames(clustering_matrix) <- "Normalised Log P-Value"
    
    # Calculate pairwise Euclidean distances between rows of matrix
    distance_matrix <- dist(clustering_matrix, method = "euclidean")
    hc <- hclust(distance_matrix, method = "complete")
    
    # Dendrogram is cut into clusters "k" determined by the input slider
    kmeans_result <- cutree(hc, k = input$clusters)
    
    # Dataframe contains clusters with their assigned gene_symbols
    annotation_row <- data.frame(Cluster = factor(kmeans_result),
                                 Mouse_Strain = selected_genes_data$mouse_strain,
                                 Mouse_Life_Stage = selected_genes_data$mouse_life_stage,
                                 Parameter_Name = selected_genes_data$parameter_name)
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
      tagList( # Create a list of UI elements to display 
        lapply(names(cluster_genes), function(cluster_name) { # Iterate over each cluster name in the 'cluster_genes' object   
          div( # Create a container for each cluster's display
            strong(cluster_name), # Display the cluster name in bold
            ": ", # Add a separator
            cluster_genes[[cluster_name]] # Display the genes associated with the cluster
          )
        })
      )
    })
    
    # Dynamically adjust font size and plot dimensions based on the number of selected genes
    num_genes <- length(selected_genes)
    
    # Adjust row font size: Larger font for fewer genes, smaller font for more genes
    fontsize_row <- ifelse(num_genes < 20, 15, 10)
    
    # Adjust font size for numbers displayed in heatmap cells : Similar concept as above
    fontsize_number <- ifelse(num_genes < 20, 15, 10)
    
    # Adjust heatmap plot height: Smaller plot for fewer genes, dynamically scaled for larger numbers
    plot_height <- ifelse(num_genes < 20, 750, 100 + num_genes * 20)
    
    # Adjust heatmap plot width : Smaller plot for fewer genes, dynamically scaled for larger numbers 
    plot_width <- ifelse(num_genes < 20, 750, 100 + num_genes * 20)
    
    # Inside the heatmap code
    pheatmap(
      clustering_matrix, # Input data 
      cluster_rows = TRUE, # Enable clustering across rows 
      cluster_cols = FALSE, # Disable clustering across columns
      annotation_row = annotation_row, # Provide row annotations
      annotation_colors = annotation_colors, # Define colours for the annotations
      color = colorRampPalette(c("cadetblue1", "plum1", "maroon2"))(50), # Set colour gradient for the heatmap with 50 shades 
      fontsize_row = fontsize_row, # Set the font size for row labels
      fontsize_col = 10, # Set the font size for column labels 
      display_numbers = TRUE,  # Display numerical values in the cells of the heatmap
      fontsize_number = fontsize_number, # Set the font size of the numbers displayed in the cells of the heatmap
      treeheight_row = 100, # Set height of the row dendogram
      treeheight_col = 100, # Set height of the column dendogram
      main = paste("Clustering of Genes by ", input$clusters, " Clusters", sep = ""), # Set the title of the heatmap, dynamically including the number of clusters
      height = plot_height, # Specify the overall height of the heatmap plot
      width = plot_width, # Specify the overall width of the heatmap plot
      angle_col = 90,  # Rotate columns for better readability 
      legend = TRUE,  # Display a figure legend 
      legend_title = "Significance Level", # Set the title of the heatmap legend 
      legend_position = "bottomright", # Position the legend at the bottom right 
      legend_direction = "horizontal", # Horizontal figure legend 
      annotation_legend = TRUE,  # Include a legend for the row/column annotations
      annotation_legend_side = "right" # Position the annotation legend on the right side of the plot
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
    paste("Number of Selected Genes: ", num_genes)
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
  
  # Download the cluster data
  # Inside the server function
  output$download_data <- downloadHandler(
    filename = function() {
      paste("cluster_data.csv")
    },
    content = function(file) {
      selected_genes <- selected_genes_reactive()  # Get selected genes
      
      # Ensure that only genes from clustering_data are used
      selected_genes_data <- clustering_data %>%
        filter(gene_symbol %in% selected_genes)
      
      # Get the clustering results
      kmeans_result <- cutree(hclust(dist(scale(as.matrix(selected_genes_data[, c("log_p_value")])))), k = input$clusters)
      
      # Add the clustering results to the dataframe
      selected_genes_data$Cluster <- kmeans_result
      
      # Write the data to a CSV file
      write.csv(selected_genes_data[, c("gene_symbol", "mouse_strain", "mouse_life_stage", "parameter_name", "log_p_value", "Cluster")], file, row.names = FALSE)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
