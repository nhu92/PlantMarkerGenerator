library(shiny)
library(shinyFiles)

options(repos = BiocManager::repositories())
# Source original functions
source("chlgene_fingerprintable.R")

ui <- fluidPage(
  titlePanel("\U0001F9EC Gene Fingerprinting Marker Extractor"),
  sidebarLayout(
    sidebarPanel(
      # shared folder input
      textInput("dir", "Choose Working Directory:", value = "", placeholder = "Click 'Browse' to select a folder"),
      actionButton("browse_dir", "📁 Browse Folder"),
      br(),
      
      # function 1: homologous region search
      textInput("gname1", "Gene Name 1 (required):", value = "", placeholder = "e.g., matK"),
      textInput("gname2", "Gene Name 2 (optional):", value = "", placeholder = "e.g., psbA"),
      
      actionButton("run", "\U0001F680 Run Analysis", class = "btn btn-success"),
      br(), 
      sliderInput("t_offhit", "Minimum distance for the Off-target",
                  min = 0, max = 1, value = 0.6
      ),
      sliderInput("r_offhit", "Disagreement species ratio to consider as off-target",
                  min = 0, max = 1, value = 0.7
      ),
      sliderInput("t_indist", "Minimum distance for distinguishable ",
                  min = 0.0, max = 0.05, value = 0.01
      ),
      br(),
      actionButton("combine", "\U0001F4CA Combine & Evaluate Classifications", class = "btn btn-primary"),
      br(), 
      
      # function 2: full search
      selectInput("mode_common", "Mode for Common Region Processing:",
                  choices = c("distance", "length_indel"),
                  selected = "distance"),
      br(),
      actionButton("run_common_regions", "🔍 Process All Common Regions", class = "btn btn-info"),
      
      # function 3: Primer selection
      uiOutput("fasta_selector"),
      radioButtons("analysis_type", "Select Analysis Method:",
                   choices = c("Gappy" = "gappy", "Divergent" = "divergent")),
      br(),
      uiOutput("dynamic_params"),
      br(),
      actionButton("run_analysis", "🧬 Run Amplicon Analysis", class = "btn btn-warning")

    ),
    mainPanel(
      h4("\U0001F4AC Status Indicator"),
      verbatimTextOutput("status", placeholder = TRUE),
      br(),
      h4("\U0001F5C3 Final Classification Table"),
      tableOutput("final_table"),
      br(),
      h4("\U0001F4D1 Combined Summary Table"),
      tableOutput("combined_summary")
    )
  )
)

server <- function(input, output, session) {
  volumes <- c(Home = "~", getVolumes()())
  shinyDirChoose(input, "dir", roots = volumes)
  
  observeEvent(input$browse_dir, {
    dir_path <- rstudioapi::selectDirectory()
    if (!is.null(dir_path)) {
      updateTextInput(session, "dir", value = dir_path)
    }
  })
  
  observeEvent(input$run, {
    req(input$gname1)
    dir_path <- parseDirPath(volumes, input$dir)
    
    output$status <- renderText({ "Step 1: Extracting sequences..." })
    g1 <- input$gname1
    g2 <- input$gname2
    
    if (g2 == "") {
      out_fasta <- file.path(tempdir(), paste0(g1, ".fasta"))
      out_dist <- file.path(tempdir(), paste0(g1, "_distance.tsv"))
      out_class <- paste0(dir_path, "\\",  g1, "_classification.tsv")
      
      extract_shortest_interval_between_genes(g1, input_dir = dir_path, output_fasta = out_fasta)
      
      output$status <- renderText({ "Step 2: Aligning sequences and calculating distances..." })
      align_and_calc_distance(out_fasta, out_dist)
      
      output$status <- renderText({ "Step 3: Classifying species..." })
      result <- classify_species_by_distance(out_dist, input$t_offhit, input$r_offhit, input$t_indist)
      write.table(result, out_class, sep = "\t", row.names = FALSE, quote = FALSE)
      
      output$status <- renderText({ paste("✅ Analysis completed for", g1) })
    } else {
      combo_label <- paste0(g1, "_", g2)
      out_fasta <- file.path(dir_path, paste0(combo_label, ".fasta"))
      out_dist <- file.path(dir_path, paste0(combo_label, "_distance.tsv"))
      out_class <- paste0(dir_path, "\\",  combo_label, "_classification.tsv")
      
      extract_shortest_interval_between_genes(g1, g2, input_dir = dir_path, output_fasta = out_fasta)
      
      output$status <- renderText({ "Step 2: Aligning sequences and calculating distances..." })
      align_and_calc_distance(out_fasta, out_dist)
      
      output$status <- renderText({ "Step 3: Classifying species..." })
      result <- classify_species_by_distance(out_dist, input$t_offhit, input$r_offhit, input$t_indist)
      write.table(result, out_class, sep = "\t", row.names = FALSE, quote = FALSE)
      
      output$status <- renderText({ paste("✅ Analysis completed for", g1, "and", g2) })
    }
    
    output$final_table <- renderTable({
      if (input$gname2 == "") {
        read.delim(paste0(dir_path, "\\", input$gname1, "_classification.tsv"))
      } else {
        read.delim(paste0(dir_path, "\\", input$gname1, "_", input$gname2, "_classification.tsv"))
      }
    })
  })
  
  output$fasta_selector <- renderUI({
    req(input$dir)
    fasta_files <- list.files(input$dir, pattern = "\\.fasta$", full.names = TRUE)
    if (length(fasta_files) == 0) return(NULL)
    selectInput("selected_fasta", "Select a FASTA File:", choices = fasta_files)
  })
  
  output$dynamic_params <- renderUI({
    req(input$analysis_type)
    
    common_params <- tagList(
      numericInput("window_size", "Window Size:", value = 100),
      numericInput("window_step", "Step Size:", value = 20),
      numericInput("top_n_windows", "Top N Windows:", value = 3),
      numericInput("min_amplicon", "Minimum Amplicon Size:", value = 80),
      numericInput("max_amplicon", "Maximum Amplicon Size:", value = 600),
      numericInput("overlap_fraction", "Overlap Fraction:", value = 0.7),
      numericInput("max_per_amp_region", "Max per Amplicon Region:", value = 3),
      numericInput("output_top", "Max Amplicons in Output:", value = 20)
    )
    
    if (input$analysis_type == "divergent") {
      tagList(
        common_params,
        checkboxInput("save_distance_csv", "Save Distance Matrix (CSV)", value = FALSE),
        textInput("prefix", "Filename Prefix", value = "divergent")
      )
    } else {
      common_params
    }
  })
  
  observeEvent(input$run_analysis, {
    req(input$selected_fasta)
    
    if (input$analysis_type == "gappy") {
      result <- findAndDesignGappyAmplicons(
        fasta_file = input$selected_fasta,
        window_size = input$window_size,
        window_step = input$window_step,
        top_n_windows = input$top_n_windows,
        min_amplicon = input$min_amplicon,
        max_amplicon = input$max_amplicon,
        overlap_fraction = input$overlap_fraction,
        max_per_amp_region = input$max_per_amp_region,
        output_top = input$output_top
      )
    } else {
      result <- findAndDesignDivergentAmplicons(
        fasta_file = input$selected_fasta,
        window_size = input$window_size,
        window_step = input$window_step,
        top_n_windows = input$top_n_windows,
        min_amplicon = input$min_amplicon,
        max_amplicon = input$max_amplicon,
        overlap_fraction = input$overlap_fraction,
        max_per_amp_region = input$max_per_amp_region,
        output_top = input$output_top,
        save_distance_csv = input$save_distance_csv,
        prefix = input$prefix
      )
    }
    
    output$final_table <- renderTable(result)
    output$status <- renderText("✅ Amplicon analysis complete.")
  })
  
  
  observeEvent(input$combine, {
    dir_path <- parseDirPath(volumes, input$dir)
    output$status <- renderText({ "\U0001F50D Analyzing gene combinations..." })
    summary_result <- combine_and_evaluate_classification(dir_path)
    write.table(summary_result, file = "gene_combination_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    output$combined_summary <- renderTable(summary_result)
    output$status <- renderText({ "\U0001F389 Gene combination analysis complete!" })
  })
  
  observeEvent(input$run_common_regions, {
    req(input$dir)
    dir_path <- parseDirPath(volumes, input$dir)
    output_dir <- file.path(dir_path, "candidate")
    dir.create(output_dir, showWarnings = FALSE)
    
    output$status <- renderText({ "🔄 Processing common regions..." })
    
    log_file <- file.path(output_dir, "common_region_log.txt")
    con <- file(log_file, open = "wt")
    sink(con, type = "output")
    sink(con, type = "message")
    
    tryCatch({
      process_all_common_regions(
        input_dir = dir_path,
        output_dir = output_dir,
        t_offhit = input$t_offhit,
        r_offhit = input$r_offhit,
        t_indist = input$t_indist,
        mode = input$mode_common
      )
    }, error = function(e) {
      cat("❌ Error:", e$message, "\n")
    })
    
    sink(type = "message")
    sink()
    close(con)
    
    output$status <- renderText({
      paste("✅ Finished. Log saved to:", log_file)
    })
  })
}

shinyApp(ui, server)
