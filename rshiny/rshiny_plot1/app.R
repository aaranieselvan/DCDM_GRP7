# RShiny Plot 1

# Install and load libraries 
# The packages were installed in the console prior to loading
library(shiny)
library(ggplot2)
library(dplyr)

# Load the clean mouse data 
data_rshiny <- read.csv("~/Desktop/DCDM_GRP7/outputs/clean_final_data.csv")

# Apply the significance threshold, calculate FDR, and then transform to log scale
data_rshiny <- data_rshiny %>%
  mutate(
    Significant = pvalue <= 0.05,  # Apply the significance threshold (0.05)
    FDR = p.adjust(pvalue, method = "BH"),  # Apply Benjamini-Hochberg FDR adjustment
    log_p_value = -log10(FDR)  # Transform FDR-adjusted p-value to -log10 scale
  )

# Define the User Interface (UI) of the Shiny App
ui <- fluidPage(
  
  # Add a title panel to the app
  titlePanel("Statistically Significant P-Values by Phenotype"), 
  
  # Layout the app into a sidebar and a main content area
  sidebarLayout(  
    
    # User inputs will be placed here
    sidebarPanel( 
      
      # Create a dropdown for the user to choose a knockout gene
      selectInput("gene_symbol", # Input ID used to refer to this control in the server logic
                  "Select Knockout Gene:", # Label displayed above the dropdown
                  choices = unique(data_rshiny$gene_symbol),  # Dynamically populate choices
                  selected = unique(data_rshiny$gene_symbol)[1])  # Set the default selection to the first gene in list
    ),
    # Define the main panel, where the output plot will be displayed 
    mainPanel(
      
      # Display a plot as the main output of the app
      plotOutput("pvalplot", # Output ID used to reference this plot in the server logic
                 height = "700px")  # Set the height to 700 pixels for a better fit
    )
  )
)

# Define Server Logic
server <- function(input, output) {
  
  # Create a reactive plot output that can update dynamically 
  output$pvalplot <- renderPlot({
    
    # Filter the data for the selected gene symbol
    filtered_mouse_data <- data_rshiny %>%
      filter(Significant == TRUE & gene_symbol == input$gene_symbol)
    
    # Generate the plot
    plot_volcano(filtered_mouse_data)
  })
}

# Define a function to plot significant p-values
plot_volcano <- function(filtered_mouse_data) { 
  ggplot(filtered_mouse_data, aes(x = parameter_name, y = log_p_value)) +
    
    # Add vertical segments representing the -log10(FDR) from y=0 to the log_p_value
    geom_segment(aes(x = parameter_name, y = 0, yend = log_p_value), 
                 size = 1, color = "lightpink") +
    
    # Add points on the plot representing the -log10(FDR) for each phenotype
    geom_point(size = 3, color = "deeppink") +
    
    # Add text labels for the log_p_value on top of each point and adjust positioning 
    geom_text(aes(label = round(log_p_value, 2)), vjust = -1, hjust = 0.5, color = "black") +  
    
    # Apply a minimal theme
    theme_minimal() +
    
    # Add title and axis labels
    labs(
      title = "Significant Phenotypic Changes by Knockout Gene",
      x = "Phenotype",
      y = "-log10(FDR)"
    ) +
    
    # Adjust plot theme for better fit
    theme(
      plot.title = element_text(face = "bold", size = 16), # Bold the title and adjust font size 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Rotate x-axis labels for better fitting
      axis.title.x = element_text(size = 14),  # Adjust font size for x-axis label
      axis.title.y = element_text(size = 14),  # Adjust font size for y-axis label
      plot.margin = margin(10, 10, 40, 60)  # Increase the bottom and left margin for x-axis labels
    ) +
    
    # Expands limits to stop point being cut off
    coord_cartesian(clip = "off") 
  
}

# Run the Shiny App
shinyApp(ui = ui, server = server)