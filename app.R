library(shiny)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(forestplot)
library(RColorBrewer)
library(grid)
library(ggiraph)
library(ggplotify)
library(here)


############## load in data and format for visualization ######################
here::i_am("app.R")
inflam_genes <- read_csv("proteomics_base_results.csv")
adjusted_p <- (0.05/(10*365))

inflam_genes2 <- inflam_genes %>%
  mutate(
    Significant = P < adjusted_p,
    statistic = BETA / SE,
    Phenotype = toupper(Phenotype)
  )

sig_proteins <- inflam_genes2 %>%
  filter(Significant) %>%
  pull(Phenotype)

sig_inflam_genes <- inflam_genes2 %>%
  filter(Phenotype %in% sig_proteins)

genes_2 <- sig_inflam_genes %>%
  select(Gene, Phenotype, statistic) %>%
  arrange(Phenotype) %>%
  pivot_wider(names_from = Phenotype, values_from = statistic) %>%
  column_to_rownames("Gene")

genes_matrix <- data.matrix(genes_2)
row_names <- rownames(genes_matrix)
col_names <- colnames(genes_matrix)

tooltip_df <- expand.grid(
  Gene = row_names,
  Protein = col_names
) %>%
  mutate(
    statistic = as.vector(genes_matrix),
    Significant = inflam_genes2$Significant[
      match(paste(Gene, Protein),
            paste(inflam_genes2$Gene, inflam_genes2$Phenotype))
    ],
    tooltip = paste0(
      "Gene: ", Gene, "\n",
      "Protein: ", Protein, "\n",
      "T-statistic: ", round(statistic, 2), "\n",
      "Significant: ", Significant
    ),
    data_id = paste(Gene, Protein, sep = "||")
  )


############# creating function for forest plot #######################

makeForestPlot <- function(selected_genes) {
  if (length(selected_genes) == 0) return(NULL)
  
  table_data <- inflam_genes2 %>%
    filter(Gene %in% selected_genes & Significant) %>%
    mutate(
      `Estimate (95% CI)` = paste0(
        round(BETA, 2), " (",
        round(BETA - 1.96 * SE, 2), "–",
        round(BETA + 1.96 * SE, 2), ")"
      ),
      mean = BETA,
      lower = BETA - 1.96 * SE,
      upper = BETA + 1.96 * SE
    )
  
  if (nrow(table_data) == 0) return(NULL)
  
  unique_colors <- RColorBrewer::brewer.pal(max(3, length(selected_genes)), "Set2")
  color_map <- setNames(unique_colors, selected_genes)
  table_data$color <- unname(color_map[table_data$Gene])
  
  fp_col <- fpShapesGp(
    lines = lapply(table_data$color, \(x) gpar(col = x)),
    box   = lapply(table_data$color, \(x) gpar(fill = x, col = x))
  )
  
  label_matrix <- as.matrix(table_data[, c("Gene", "Phenotype", "Estimate (95% CI)")])
  
  forestplot(
    labeltext = label_matrix,
    mean = table_data$mean,
    lower = table_data$lower,
    upper = table_data$upper,
    shapes_gp = fp_col,
    txt_gp = fpTxtGp(label = gpar(fontsize = 12)),
    title = "Effect of Selected CHIP Mutations on Proteins"
  )
}

################# UI code for Shiny #################################

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', sans-serif; }

      .chip-button {
  width: 100%;
  height: 130px;
  background-color: #C8A2C8;
  color: black;
  border: none;
  border-radius: 10px;
  font-size: 26px;
  font-weight: bold;
  transition: 0.2s;
}
.chip-button:hover {
  background-color: #A98DA9;
  cursor: pointer;
}"))
  ),
  
  ###Creating table####
  
  titlePanel("CHIP Mutations and the Inflammatory Proteome"),
  
  tabsetPanel(
    id = "mainTabs",
    
    ###### first page overview layout#########
    tabPanel(
      title = "Overview",
      
      br(),
      
      # CHIP button
      div(
        style = "width:100%;",
        actionButton(
          "chip_def_btn",
          label = "What is CHIP?",
          class = "chip-button",
          style = "width:100%;"
        )
      ),
      
      br(),
      
      # CHIP definition dropdown
      conditionalPanel(
        condition = "input.chip_def_btn % 2 == 1",
        wellPanel(
          h4("CHIP: Clonal Hematopoiesis of Indeterminate Potential"),
          p("CHIP mutations are somatic mutations in hematopoietic stem cells that lead to clonal expansion and are associated with increased cardiovascular disease risk")
        )
      ),
      
      br(),
      
      # PNG of study overview
      div(
        class = "overview-panel",
        h3("Study Overview"),
        img(
          src = "proteomics_goal.png",
          style = "width:100%; border-radius:10px; margin-bottom:15px;"
        )),
      
      br(),
      
      div(style = "display:flex;gap:20px;width:100%;",
        
        # Number of participants tile
        div(style = "
      flex:1;
      background-color:#BFDBFE;
      border-radius:20px;
      padding:20px;
      text-align:center;
      font-size:28px;
      font-weight:700;
      color:black;",
          "~ 48,000\nUKB Participants"),
        
        # Number of proteins tile
        div(
          style = "
      flex:1;
      background-color:#BFDBFE;
      border-radius:20px;
      padding:20px;
      text-align:center;
      font-size:28px;
      font-weight:700;
      color:black;
    ",
          "365\nInflammatory Proteins"
        )
      ),
      
      
      br()),
    
    ############ Results Tab #######################
    
    tabPanel(
      title = "Results",
      br(),
      fluidRow(
        column(
          width = 12,
          girafeOutput("heatmap", width = "100%", height = "500px"),
          tags$p(
            style = "text-align:center; font-size:16px; margin-top:10px;",
            "Figure 1. Heatmap showing gene–protein associations in genes with CHIP mutations."
          )
        )
      ),
      br(),
      fluidRow(
        column(
          width = 3,
          wellPanel(
            h4("Selected Genes"),
            selectInput(
              "selectedGenes", NULL,
              choices = sort(unique(inflam_genes2$Gene[inflam_genes2$Significant])),
              multiple = TRUE
            ),
            helpText("Click cells in the heatmap to populate this list.")
          )
        ),
        column(
          width = 9,
          plotOutput("forestPlot", height = "800px"),
          tags$p(
            style = "text-align:center; font-size:16px; margin-top:10px;",
            "Figure 2. Forest plot of estimated effect sizes for CHIP mutations."
          )
        )
      )
    )
  )
)

############## Server code for shiny ##################

server <- function(input, output, session) {
  
  observeEvent(input$go_to_results, {
    updateTabsetPanel(session, "mainTabs", selected = "CHIP Heatmap & Forest Plot")
  })
  
  output$heatmap <- renderGirafe({
    
    df <- tooltip_df %>%
      mutate(
        x = as.numeric(factor(Protein, levels = col_names)),
        y = as.numeric(factor(Gene, levels = rev(row_names)))
      )
    
    p <- ggplot(df, aes(
      x, y,
      fill = statistic,
      tooltip = tooltip,
      data_id = data_id
    )) +
      geom_tile_interactive() +
      scale_fill_gradient2(low = "#00008b", mid = "white", high = "darkgreen") +
      scale_y_continuous(breaks = 1:length(row_names), labels = rev(row_names)) +
      scale_x_continuous(breaks = 1:length(col_names), labels = col_names) +
      labs(
        x = "Protein",
        y = "Gene",
        title = "Strength of Associations\nBetween CHIP and Inflammatory Proteins",
        fill = "T-statistic"
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, size = 18, hjust = 1),
        axis.text.y = element_text(size = 18, face = "italic"),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title = element_text(size = 30, hjust = 0.5),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)
      )
    
    girafe(
      ggobj = p,
      height_svg = 10,
      width_svg = 20,
      options = list(
        opts_hover(css = "fill-opacity:0.7;cursor:pointer;"),
        opts_selection(type = "single", css = "stroke:black;stroke-width:2px;"),
        opts_toolbar(saveaspng = TRUE),
        opts_zoom(max = 5),
        opts_tooltip(css = "background-color:white;padding:6px;border-radius:4px;")
      )
    )
  })
  
  observeEvent(input$heatmap_selected, {
    req(input$heatmap_selected)
    
    ids <- input$heatmap_selected
    genes <- sapply(strsplit(ids, "\\|\\|"), `[`, 1)
    
    updateSelectInput(
      session,
      "selectedGenes",
      selected = unique(genes)
    )
  })
  
  output$forestPlot <- renderPlot({
    makeForestPlot(input$selectedGenes)
  })
}

shinyApp(ui, server)

