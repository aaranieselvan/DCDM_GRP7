# RShiny Plot 2

# Install and load libraries
# Packages were installed in the console prior to loading
library(shiny)
library(ggplot2)
library(dplyr)

# Load the clean mouse data
data <- read.csv("~/Desktop/DCDM_GRP7/outputs/clean_final_data.csv")


# Use the pipe function to transform the p-value to -10log scale, and then filter by significance
# We apply the -log10 scale for better visualisation of small p-values, emphasising statistically significant p-values
data <- data %>%
  mutate(
    log_p_value = -log10(pvalue),  # Transform p-value to -log10 scale
    Significant = pvalue <= 0.05,  # Add significance threshold
    FDR = p.adjust(pvalue, method = "BH"),  # Apply Benjamini-Hochberg FDR adjustment
    log_p_value = -log10(FDR)  # Transform FDR-adjusted p-value to -log10 scale
  )


# Define UI of the Shiny App
ui <- fluidPage(
  
  #Add a title panel to the app
  titlePanel("Phenotype by Gene Symbol"),
  
  # Layout the app into a sidebar and a main content area
  sidebarLayout(
    
    # User inputs will be placed here
    sidebarPanel(
      
      # Create a dropdown for the user to select a phenotype
      selectInput("parameter_name", "Select Phenotype:", # Input ID used to refer to this control in the server logic
                  choices = unique(data$parameter_name),  # Dynamically populate choices for parameter names
                  selected = unique(data$parameter_name)[1])  # Default to the first phenotype
    ),
    
    # Define the main panel, where the output plot will be displayed 
    mainPanel(
      
      # Display a plot as the main output of the app
      plotOutput("pvalplot", # Output ID used to reference this plot in the server logic
                 height = "700px") # Adjust height for better fit
    )
  )
)

# Define Server Logic
server <- function(input, output) {
  
  # Create a reactive plot output that can update dynamically 
  output$pvalplot <- renderPlot({
    
    # Filter the data based on the selected parameter name
    # This will keep only rows where, significance is true (i.e pval <0.05), and where the phenotype matches the user's selection 
    filtered_data <- data %>%
      filter(Significant == TRUE & parameter_name == input$parameter_name)
    
    # Generate the plot
    plot_lollipop(filtered_data)
  })
}


# Define a function to plot significant p-values
plot_lollipop <- function(filtered_data) {
  ggplot(filtered_data, aes(x = gene_symbol, y = pvalue)) +
    
    # Add vertical segments representing the -log10(p-value) from y=0 to the log_p_value
    geom_segment(aes(x = gene_symbol, y = 0, yend = pvalue), 
                 size = 1, color = "lightpink") +
    
    # Add points on the plot representing the -log10(FDR) for each phenotype
    geom_point(size = 4, color = "deeppink") +
    
    # Add text labels for the log_p_value on top of each point and adjust positioning 
    geom_text(aes(label = round(log_p_value, 2)), vjust = -1, hjust = 0.5, color = "black") +  
    
    # Apply a minimal theme
    theme_minimal() +
    
    # Add title and axis labels
    labs(
      title = "Statistical Scores of Knockout Mice for a Selected Phenotype",
      x = "Gene Symbol",
      y = "Log-Transformed Adjusted P-Values (FDR)"
    ) +
    
    # Adjust plot theme for better fit
    theme(
      plot.title = element_text(face = "bold", size = 16), # Bold the title and adjust font size 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Rotate x-axis labels for better fitting
      axis.title.x = element_text(size = 14),  # Adjust font size for x-axis label
      axis.title.y = element_text(size = 14),  # Adjust font size for y-axis label
      plot.margin = margin(10, 10, 30, 50)  # Increase the bottom and left margin for x-axis labels
      
    ) +
    
    # Expands limits to stop point being cut off
    coord_cartesian(clip = "off")
  
}

# Run the Shiny App
shinyApp(ui = ui, server=server)