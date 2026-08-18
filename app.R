options(shiny.maxRequestSize = 250 * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(ggplot2)
  library(digest)
})

source("R/io.R", local = TRUE)
source("R/regions.R", local = TRUE)
source("R/analysis.R", local = TRUE)
source("R/primers.R", local = TRUE)

theme <- bs_theme(
  version = 5,
  bg = "#F7FAF8",
  fg = "#17352A",
  primary = "#167D50",
  secondary = "#5F7169"
)

status_badge <- function(text, type = "info") {
  colour <- switch(type, success = "#DDF5E7", warning = "#FFF2CC", danger = "#FFE2E2", "#E6F0FA")
  div(style = sprintf("background:%s;padding:.75rem 1rem;border-radius:.65rem;margin-bottom:1rem;", colour), text)
}

ui <- page_navbar(
  title = "ChloroMarker",
  theme = theme,
  header = tags$head(tags$style(HTML(
    "body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif}.navbar-brand{font-weight:750;letter-spacing:.02em}.card{border:0;box-shadow:0 8px 24px rgba(18,55,42,.07)}\
     .btn-primary,.btn-success{font-weight:650}.workflow-note{color:#53655D;font-size:.92rem}.dataTables_wrapper{font-size:.88rem}"
  ))),
  nav_panel(
    "1. Data",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        fileInput(
          "genomes", "Annotated chloroplast genomes",
          multiple = TRUE,
          accept = c(".gb", ".gbk", ".gbff", ".genbank", ".fsa", ".fasta", ".fna", ".tbl", ".csv")
        ),
        p(class = "workflow-note", "Upload one multi-record GenBank file, several GenBank files, or an NCBI .fsa + .tbl pair. An optional taxonomy-validation CSV can update names."),
        actionButton("parse", "Parse annotations", class = "btn-success w-100"),
        hr(),
        selectizeInput("taxa", "Taxa to compare", choices = NULL, multiple = TRUE),
        p(class = "workflow-note", "Marker discovery is most meaningful within a genus, family, or defined comparison set."),
        numericInput("min_coverage", "Minimum locus coverage", value = 0.8, min = 0.2, max = 1, step = 0.05),
        numericInput("max_candidates", "Candidates for precise alignment", value = 16, min = 3, max = 40, step = 1),
        numericInput("distance_threshold", "Minimum barcode gap", value = 0.01, min = 0, max = 0.1, step = 0.001),
        actionButton("analyze", "Discover markers", class = "btn-primary w-100")
      ),
      uiOutput("status"),
      card(
        card_header("Uploaded records"),
        DTOutput("records_table")
      )
    )
  ),
  nav_panel(
    "2. Markers",
    card(
      card_header("Ranked single markers"),
      p(class = "workflow-note", "Ranking prioritizes species distinguishability, coverage, alignment quality, and sequence divergence."),
      DTOutput("marker_table"),
      downloadButton("download_markers", "Download marker table")
    ),
    layout_columns(
      card(card_header("Two-marker combinations"), DTOutput("combination_table"), downloadButton("download_combinations", "Download combinations")),
      card(card_header("Selected-marker distance matrix"), selectInput("marker_key", "Marker", choices = NULL), plotOutput("distance_plot", height = 520), downloadButton("download_distance", "Download matrix")),
      col_widths = c(6, 6)
    )
  ),
  nav_panel(
    "3. Assay design",
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        selectInput("assay_marker", "Marker for assay design", choices = NULL),
        numericInput("min_amplicon", "Minimum amplicon (bp)", 120, min = 60, max = 1000),
        numericInput("max_amplicon", "Maximum amplicon (bp)", 800, min = 100, max = 3000),
        numericInput("max_degeneracy", "Maximum primer-pair degeneracy", 16, min = 1, max = 256),
        actionButton("design", "Design primers", class = "btn-primary w-100"),
        hr(),
        selectInput("assay_id", "Primer pair", choices = NULL),
        checkboxInput("use_digest", "Show PCR-RFLP gel", FALSE),
        selectInput("enzyme", "Restriction enzyme", choices = NULL)
      ),
      uiOutput("recommendation"),
      card(card_header("Candidate primer pairs"), DTOutput("primer_table"), downloadButton("download_primers", "Download primers")),
      layout_columns(
        card(card_header("Predicted gel"), plotOutput("gel_plot", height = 520)),
        card(card_header("Restriction-enzyme screen"), DTOutput("restriction_table"), downloadButton("download_restriction", "Download RFLP screen")),
        col_widths = c(7, 5)
      )
    )
  ),
  nav_panel(
    "Methods & interpretation",
    card(
      card_header("Decision framework"),
      h4("1. Fast screening, then precise analysis"),
      p("Annotated genes and adjacent intergenic spacers are extracted first. A 4-mer and length-based screen limits expensive multiple alignments to the most promising candidates."),
      h4("2. Marker performance"),
      p("Distances are uncorrected p-distances with pairwise deletion. A species is distinguishable when its minimum interspecific distance minus observed maximum intraspecific distance meets the selected barcode-gap threshold."),
      h4("3. Primer and gel prediction"),
      p("Primers are IUPAC consensus oligos screened for gap frequency, degeneracy, GC content, melting temperature, homopolymers, and compatible amplicon length. Gel calls use an approximate 2% size-resolution rule and are predictions, not laboratory validation."),
      h4("4. PCR-RFLP and escalation"),
      p("Common restriction enzymes are scanned in silico. If neither amplicon length nor restriction profiles fully resolve the taxa, the app recommends amplicon sequencing when sequence divergence is diagnostic. Plastid-unresolved groups are explicitly routed to nuclear markers or orthogonal evidence."),
      div(class = "alert alert-warning", "All primer, gel, and restriction outputs require wet-lab validation. Chloroplast inheritance, hybridization, introgression, and chloroplast capture can prevent species-level identification.")
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(
    records = NULL,
    result = NULL,
    combinations = NULL,
    primers = NULL,
    restrictions = NULL,
    message = "Upload annotated chloroplast genomes to begin.",
    message_type = "info"
  )
  cache <- new.env(parent = emptyenv())

  output$status <- renderUI(status_badge(values$message, values$message_type))

  observeEvent(input$parse, {
    req(input$genomes)
    values$message <- "Parsing sequences and annotations..."
    values$message_type <- "info"
    tryCatch({
      records <- read_chloromarker_uploads(input$genomes)
      values$records <- records
      values$result <- NULL
      taxa <- sort(unique(vapply(records, `[[`, character(1), "organism")))
      updateSelectizeInput(session, "taxa", choices = taxa, selected = taxa, server = TRUE)
      values$message <- sprintf("Parsed %d records and %d named gene annotations. Select the intended comparison set, then discover markers.",
        length(records), sum(vapply(records, function(x) nrow(representative_gene_features(x$features)), integer(1))))
      values$message_type <- "success"
    }, error = function(e) {
      values$message <- conditionMessage(e)
      values$message_type <- "danger"
    })
  })

  selected_records <- reactive({
    req(values$records)
    selected <- input$taxa
    if (is.null(selected) || !length(selected)) return(values$records)
    values$records[vapply(values$records, function(x) x$organism %in% selected, logical(1))]
  })

  observeEvent(input$analyze, {
    records <- selected_records()
    validate(need(length(records) >= 2L, "Select at least two records."))
    key <- digest(list(
      ids = names(records), coverage = input$min_coverage,
      candidates = input$max_candidates, threshold = input$distance_threshold
    ))
    values$message <- "Extracting candidate loci and aligning the pre-screened set..."
    values$message_type <- "info"
    tryCatch({
      if (exists(key, envir = cache, inherits = FALSE)) {
        bundle <- get(key, envir = cache)
      } else {
        bundle <- withProgress(message = "Marker discovery", value = 0, {
          catalog <- build_region_catalog(records, min_coverage = input$min_coverage)
          result <- analyze_catalog(
            catalog, length(records), max_candidates = input$max_candidates,
            threshold = input$distance_threshold, processors = 1L,
            progress = function(i, n, locus) incProgress(1 / n, detail = locus)
          )
          combinations <- evaluate_marker_combinations(result, max_k = 2L, threshold = input$distance_threshold)
          list(result = result, combinations = combinations)
        })
        assign(key, bundle, envir = cache)
      }
      values$result <- bundle$result
      values$combinations <- bundle$combinations
      values$primers <- NULL
      values$restrictions <- NULL
      choices <- stats::setNames(bundle$result$ranking$key, bundle$result$ranking$locus)
      updateSelectInput(session, "marker_key", choices = choices, selected = bundle$result$ranking$key[[1L]])
      updateSelectInput(session, "assay_marker", choices = choices, selected = bundle$result$ranking$key[[1L]])
      values$message <- sprintf("Finished: %d candidate regions were precisely aligned and ranked.", nrow(bundle$result$ranking))
      values$message_type <- "success"
    }, error = function(e) {
      values$message <- conditionMessage(e)
      values$message_type <- "danger"
    })
  })

  observeEvent(input$design, {
    req(values$result, input$assay_marker)
    analysis <- values$result$analyses[[input$assay_marker]]
    values$message <- sprintf("Designing shared primers for %s...", analysis$region$locus)
    values$message_type <- "info"
    tryCatch({
      primers <- design_primers(
        analysis,
        min_amplicon = input$min_amplicon,
        max_amplicon = input$max_amplicon,
        max_degeneracy = input$max_degeneracy,
        top_n = 12L
      )
      values$primers <- primers
      if (!nrow(primers$table)) stop("No primer pair met the current constraints. Increase degeneracy or widen the amplicon range.")
      updateSelectInput(session, "assay_id", choices = primers$table$assay_id, selected = primers$table$assay_id[[1L]])
      restrictions <- restriction_screen(primers$assays[[1L]], analysis$region$organisms)
      values$restrictions <- restrictions
      enzyme_choices <- if (nrow(restrictions$table)) restrictions$table$enzyme else character()
      updateSelectInput(session, "enzyme", choices = enzyme_choices, selected = if (length(enzyme_choices)) enzyme_choices[[1L]] else character())
      values$message <- sprintf("Designed %d primer pairs and screened %d restriction enzymes.", nrow(primers$table), nrow(restriction_enzymes))
      values$message_type <- "success"
    }, error = function(e) {
      values$primers <- NULL
      values$restrictions <- NULL
      values$message <- conditionMessage(e)
      values$message_type <- "warning"
    })
  })

  observeEvent(input$assay_id, {
    req(values$primers, values$result, input$assay_marker, input$assay_id)
    assay <- values$primers$assays[[input$assay_id]]
    analysis <- values$result$analyses[[input$assay_marker]]
    values$restrictions <- restriction_screen(assay, analysis$region$organisms)
    choices <- values$restrictions$table$enzyme
    updateSelectInput(session, "enzyme", choices = choices, selected = if (length(choices)) choices[[1L]] else character())
  }, ignoreInit = TRUE)

  output$records_table <- renderDT({
    req(values$records)
    datatable(dataset_summary(values$records), options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$marker_table <- renderDT({
    req(values$result)
    table <- values$result$ranking
    numeric <- vapply(table, is.numeric, logical(1))
    table[numeric] <- lapply(table[numeric], function(x) round(x, 4))
    datatable(table, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  output$combination_table <- renderDT({
    req(values$combinations)
    datatable(head(values$combinations, 100), options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$distance_plot <- renderPlot({
    req(values$result, input$marker_key)
    analysis <- values$result$analyses[[input$marker_key]]
    matrix <- analysis$distance
    labels <- make.unique(unname(analysis$region$organisms[rownames(matrix)]))
    rownames(matrix) <- colnames(matrix) <- labels
    long <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
    names(long) <- c("taxon_1", "taxon_2", "distance")
    ggplot(long, aes(taxon_1, taxon_2, fill = distance)) +
      geom_tile() +
      scale_fill_viridis_c(option = "C", na.value = "grey90") +
      coord_equal() +
      labs(x = NULL, y = NULL, fill = "p-distance", title = analysis$region$locus) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1), panel.grid = element_blank())
  })

  output$primer_table <- renderDT({
    req(values$primers)
    datatable(values$primers$table, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$restriction_table <- renderDT({
    req(values$restrictions)
    datatable(values$restrictions$table, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$recommendation <- renderUI({
    req(values$result, values$primers, input$assay_marker, input$assay_id)
    recommendation <- assay_recommendation(values$result$analyses[[input$assay_marker]], values$primers, input$assay_id)
    div(class = "alert alert-success", h4(recommendation$method), p(recommendation$detail))
  })

  output$gel_plot <- renderPlot({
    req(values$result, values$primers, input$assay_marker, input$assay_id)
    analysis <- values$result$analyses[[input$assay_marker]]
    assay <- values$primers$assays[[input$assay_id]]
    enzyme <- if (isTRUE(input$use_digest) && nzchar(input$enzyme %||% "")) input$enzyme else NULL
    plot_virtual_gel(gel_data(assay, analysis$region$organisms, restriction = enzyme))
  })

  csv_download <- function(filename, getter) {
    downloadHandler(filename = filename, content = function(file) utils::write.csv(getter(), file, row.names = FALSE))
  }
  output$download_markers <- csv_download("chloromarker_ranked_markers.csv", function() values$result$ranking)
  output$download_combinations <- csv_download("chloromarker_marker_combinations.csv", function() values$combinations)
  output$download_primers <- csv_download("chloromarker_primers.csv", function() values$primers$table)
  output$download_restriction <- csv_download("chloromarker_rflp_screen.csv", function() values$restrictions$table)
  output$download_distance <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_", input$marker_key), "_distance.csv"),
    content = function(file) utils::write.csv(values$result$analyses[[input$marker_key]]$distance, file)
  )
}

shinyApp(ui, server)
