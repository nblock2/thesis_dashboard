##### load libraries #####

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

##### format data #####

here::i_am("app.R")

# Load base results
inflam_genes <- read_csv("proteomics_base_results.csv")

# Bonferroni-style correction
adjusted_p <- 0.05 / (11 * 365)

# Add significance + t-statistic + uppercase phenotype
inflam_genes2 <- inflam_genes %>%
  mutate(
    Significant = P < adjusted_p,
    statistic   = BETA / SE,
    Phenotype   = toupper(Phenotype)
  )

# getting sig only dataset
sig_only<-inflam_genes2 %>% filter(Significant==TRUE)

# Identify proteins with at least 1 significant association
sig_proteins <- inflam_genes2 %>%
  filter(Significant) %>%
  pull(Phenotype)

# Filter to significant protein–gene pairs
sig_inflam_genes <- inflam_genes2 %>%
  filter(Phenotype %in% sig_proteins)

# Wide matrix for heatmap
genes_2 <- sig_inflam_genes %>%
  select(Gene, Phenotype, statistic) %>%
  arrange(Phenotype) %>%
  pivot_wider(names_from = Phenotype, values_from = statistic) %>%
  column_to_rownames("Gene")

genes_matrix <- data.matrix(genes_2)
row_names <- rownames(genes_matrix)
col_names <- colnames(genes_matrix)

# Tooltip dataframe for interactive heatmap
tooltip_df <- expand.grid(Gene = row_names, Protein = col_names) %>%
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

##### creating function for forest plot #####

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
      mean  = BETA,
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
    mean      = table_data$mean,
    lower     = table_data$lower,
    upper     = table_data$upper,
    xlab = "Normalized Protein Expression (NPX)",
    shapes_gp = fp_col,
    txt_gp    = fpTxtGp(label = gpar(fontsize = 12)),
    title     = "Effect (NPX) of Selected CHIP Mutations\non Proteins"
  )
}

##### creating UI #####

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body { 
        font-family: 'Segoe UI', sans-serif; 
        background-color: #F5F6FA;
      }

      /* CHIP button */
      .chip-button {
        width: 100%;
        height: 130px;
        background-color: #C8A2C8;
        color: black;
        border: none;
        border-radius: 12px;
        font-size: 26px;
        font-weight: bold;
        transition: 0.2s;
        box-shadow: 0 3px 6px rgba(0,0,0,0.15);
      }
      .chip-button:hover {
        background-color: #A98DA9;
        cursor: pointer;
      }

      .dashboard-card {
        background: white;
        border-radius: 14px;
        padding: 25px;
        box-shadow: 0 3px 10px rgba(0,0,0,0.12);
        margin-bottom: 25px;
      }
      
      .info-tile {
        flex: 1;
        background-color: #BFDBFE;
        border-radius: 16px;
        padding: 25px;
        text-align: center;
        font-size: 24px;   /* Overview tiles */
        font-weight: 700;
        color: black;
        box-shadow: 0 3px 8px rgba(0,0,0,0.12);
        white-space: pre-line;
      }

      .results-tile {
        font-size: 20px !important;
        font-weight: 600;
        text-align: left !important;
        padding: 20px !important;
      }

      .insight-title {
        font-weight: 700;
        font-size: 22px;
        margin-bottom: 8px;
      }

      .insight-text {
        font-size: 18px;
        margin: 0;
      }

      .info-tile-container {
        display: flex;
        flex-direction: row;
        gap: 20px;
        width: 100%;
      }

      @media (max-width: 600px) {
        .info-tile-container {
          flex-direction: column;
          gap: 15px;
        }
      }

      .section-title {
        font-size: 28px;
        font-weight: 700;
        margin-bottom: 15px;
      }

      .footer-button {
        display: inline-block;
        padding: 12px 22px;
        background-color: #C8A2C8;
        color: black;
        font-size: 18px;
        font-weight: 600;
        border-radius: 10px;
        text-decoration: none;
        box-shadow: 0 3px 6px rgba(0,0,0,0.15);
        transition: 0.25s;
      }
      .footer-button:hover {
        background-color: #A98DA9;
        color: white !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0,0,0,0.18);
      }
    "))
  ),
  
  tabsetPanel(
    
    # Creating Overview tab
    tabPanel("Overview",
             div(
               class = "dashboard-card",
               h1("CHIP Mutations and the Inflammatory Proteome",
                  style = "text-align:center; font-weight:800; margin-bottom:10px;"),
               p("Interactive dashboard exploring gene–protein associations in CHIP carriers.",
                 style = "text-align:center; font-size:18px; margin-top:-10px;")
             ),
             
             div(
               class = "dashboard-card",
               actionButton("chip_def_btn", "What is CHIP?", class = "chip-button"),
               conditionalPanel(
                 condition = "input.chip_def_btn % 2 == 1",
                 div(style="margin-top:20px;",
                     h4("CHIP: Clonal Hematopoiesis of Indeterminate Potential"),
                     p("CHIP mutations are somatic mutations in hematopoietic stem cells that lead to clonal expansion. They are associated with increased cardiovascular disease and acute myeloid leukemia risk.")
                 )
               )
             ),
             
             div(
               class = "dashboard-card",
               div(class="section-title", "Study Overview"),
               img(src = "proteomics_goal.png",
                   style = "width:100%; border-radius:12px; margin-bottom:20px;"),
               div(
                 class = "info-tile-container",
                 div(class = "info-tile", "~ 46,000\nUKB Participants"),
                 div(class = "info-tile", "365\nInflammatory Proteins")
               )
             )
    ),
    
    # Creating results tab
    tabPanel("Results",
             div(
               class = "dashboard-card",
               div(class="section-title", "Results"),
               
               girafeOutput("heatmap", width = "100%"),
               tags$p(
                 style = "text-align:center; font-size:16px; margin-top:10px;",
                 "Figure 1. Heatmap showing gene–protein associations in genes with CHIP mutations."
               ),
               
               div(
                 class = "dashboard-card",
                 h4("Selected Genes"),
                 selectInput(
                   "selectedGenes", NULL,
                   choices = sort(unique(inflam_genes2$Gene[inflam_genes2$Significant])),
                   multiple = TRUE
                 ),
                 helpText("Click cells in the heatmap to populate this list or select genes from dropdown menu.")
               ),
               
               plotOutput("forestPlot", height = "900px"),
               tags$p(
                 style = "text-align:center; font-size:16px; margin-top:10px;",
                 "Figure 2. Forest plot of estimated effect sizes for CHIP mutations."
               )
             ),
             
             # adding main takeaways cards
             div(
               class = "dashboard-card",
               div(class="section-title", "Main Takeaways"),
               
               div(
                 class = "info-tile-container",
                 
                 # takeaway card 1
                 div(
                   class = "info-tile results-tile",
                   tags$div(class="insight-title",
                            "JAK2 drives the strongest inflammatory signature"),
                   tags$p(
                     class="insight-text",
                     "Out of CHIP mutations across all genes, those in JAK2 are estimated to upregulate the greatest number of circulating inflammation proteins."
                   )
                 ),
                 
                 # takeaway card 2
                 div(
                   class = "info-tile results-tile",
                   tags$div(class="insight-title",
                            "Inflammation may link CHIP to cardiovascular risk"),
                   tags$p(
                     class="insight-text",
                     "CHIP mutation–associated increases in circulating inflammation proteins could be contributing to elevated cardiovascular disease risk."
                   )
                 )
               )
             )
    ),
    
    # Dataset information tab
    tabPanel("Dataset Information",
             div(
               class = "dashboard-card",
               h2("Dataset Information"),
               
               h3("Dataset Description"),
               p("This dashboard analyzes associations between CHIP mutations and inflammatory protein levels in UK Biobank (UKB) participants. All UKB data is confidential, therefore all regression analyses were conducted directly on their platform. Additionally, our regression data is not publicly available as we are awaiting publication.
       The proteomics data used for this analysis are from Olink assays of blood plasma samples (Sun et al., 2023). The dataset used to generate this dashboard includes regression results for gene–protein pairs, significance indicators, and derived statistics used for visualization. Sources for the cohort and proteomics data are linked below."),
               
               br(),
               
               h3("Data Overview"),
               tableOutput("data_overview"),
               
               br(),
               
               h3("Significant Associations by Gene"),
               tableOutput("sig_by_gene"),
               
               br(),
               
               h3("Sources"),
               div(style = "display:flex; gap:20px; margin-top:10px;",
                   tags$a(
                     href = "https://www.nature.com/articles/s41586-018-0579-z",
                     target = "_blank",
                     class = "footer-button",
                     "UK Biobank Cohort Project"
                   ),
                   tags$a(
                     href = "https://www.nature.com/articles/s41586-023-06592-6",
                     target = "_blank",
                     class = "footer-button",
                     "UK Biobank Proteomics Project"
                   )
               )
             )
    ),
    
    div(
      style = "text-align:center; margin-top:40px; margin-bottom:30px;",
      tags$a(
        href = "https://github.com/nblock2/thesis_dashboard",
        class = "footer-button",
        target = "_blank",
        "View this Project on GitHub"
      )
    )
  ))


##### Creating server #####

server <- function(input, output, session) {
  
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
        legend.text = element_text(size = 14),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    girafe(
      ggobj = p,
      width_svg = 20,
      height_svg = 14,
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
  
  output$data_overview <- renderTable({
    tibble::tibble(
      Metric = c("Number of Participants", "Number of Significant Associations"),
      Value = c(
        as.integer(max(inflam_genes$N, na.rm = TRUE)),
        as.integer(nrow(sig_only))
      )
    )
  })
  
  output$sig_by_gene <- renderTable({
    sig_only %>%
      group_by(Gene) %>%
      summarise(
        Significant_Associations = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(Significant_Associations))
  })
}

##### run app #####

shinyApp(ui, server)
