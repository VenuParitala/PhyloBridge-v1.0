# ==============================================================================
# PhyloBridge v1.0: Computational Genomics & Quality Dashboard Platform
# Fully Integrated Edition - Restored & Prominent Compute Phylogeny Button
# ==============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(Biostrings)
  library(plotly)
  library(dplyr)
  library(ape)
  library(phangorn)
  library(DT)
  library(igraph)
  library(gplots)
  library(shinyjs)
  library(RColorBrewer)
  library(pegas)
  library(adegenet)
  library(hierfstat)
})

# Allow up to 50 MB FASTA uploads
options(shiny.maxRequestSize = 50 * 1024^2)

# --- 2. Helper Functions ---

calculate_assembly_stats <- function(dna_seq) {
  if (length(dna_seq) == 0) return(NULL)
  
  lengths_vec <- as.numeric(width(dna_seq))
  total_bp <- sum(lengths_vec)
  num_contigs <- length(lengths_vec)
  
  if (total_bp == 0 || num_contigs == 0) return(NULL)
  
  sorted_lengths <- sort(lengths_vec, decreasing = TRUE)
  cum_sum <- cumsum(sorted_lengths)
  
  n50_idx <- head(which(cum_sum >= total_bp * 0.5), 1)
  n50_val <- sorted_lengths[n50_idx]
  l50_val <- n50_idx
  
  n90_idx <- head(which(cum_sum >= total_bp * 0.9), 1)
  n90_val <- sorted_lengths[n90_idx]
  l90_val <- n90_idx
  
  max_len_val <- max(sorted_lengths)
  min_len_val <- min(sorted_lengths)
  
  freqs <- Biostrings::letterFrequency(dna_seq, letters = c("A", "C", "G", "T", "N"), OR = 0)
  
  if (is.matrix(freqs)) {
    total_acgt <- sum(freqs[, c("A", "C", "G", "T")])
    total_g_c <- sum(freqs[, c("G", "C")])
    total_n <- sum(freqs[, "N"])
    total_a <- sum(freqs[, "A"])
    total_t <- sum(freqs[, "T"])
    total_g <- sum(freqs[, "G"])
    total_c <- sum(freqs[, "C"])
  } else {
    total_acgt <- sum(freqs[c("A", "C", "G", "T")])
    total_g_c <- sum(freqs[c("G", "C")])
    total_n <- freqs["N"]
    total_a <- freqs["A"]
    total_t <- freqs["T"]
    total_g <- freqs["G"]
    total_c <- freqs["C"]
  }
  
  overall_gc <- ifelse(total_acgt > 0, (total_g_c / total_acgt) * 100, 0)
  overall_n_pct <- ifelse(total_bp > 0, (total_n / total_bp) * 100, 0)
  
  df_summary <- data.frame(
    Metric = c("Total Contigs / Sequences", "Total Assembly Size (bp)", 
               "N50 (bp)", "L50 (Contig Index)", "N90 (bp)", "L90 (Contig Index)",
               "Max Contig Length (bp)", "Min Contig Length (bp)", 
               "Overall GC Content (%)", "Ambiguous Bases ('N' %)"),
    Value = c(
      format(as.numeric(num_contigs), big.mark = ","),
      format(as.numeric(total_bp), big.mark = ","),
      format(as.numeric(n50_val), big.mark = ","),
      format(as.numeric(l50_val), big.mark = ","),
      format(as.numeric(n90_val), big.mark = ","),
      format(as.numeric(l90_val), big.mark = ","),
      format(as.numeric(max_len_val), big.mark = ","),
      format(as.numeric(min_len_val), big.mark = ","),
      as.character(round(overall_gc, 2)),
      as.character(round(overall_n_pct, 2))
    ),
    stringsAsFactors = FALSE
  )
  
  base_df <- data.frame(
    Base = c("A", "C", "G", "T", "N"),
    Count = c(total_a, total_c, total_g, total_t, total_n),
    stringsAsFactors = FALSE
  )
  
  list(
    summary_table = df_summary,
    num_contigs = num_contigs,
    total_bp = total_bp,
    n50_val = n50_val,
    l50_val = l50_val,
    overall_gc = round(overall_gc, 2),
    overall_n_pct = round(overall_n_pct, 2),
    base_counts = base_df,
    sorted_lengths = sorted_lengths
  )
}

compute_fast_gc <- function(dna_seq) {
  freqs <- Biostrings::letterFrequency(dna_seq, letters = c("A", "C", "G", "T"), OR = 0)
  acgt <- rowSums(freqs)
  gc_pct <- ifelse(acgt > 0, (freqs[, "G"] + freqs[, "C"]) / acgt * 100, 0)
  data.frame(Sequence = names(dna_seq), GC_Content = round(gc_pct, 2), stringsAsFactors = FALSE)
}

compute_fast_codon_usage <- function(dna_seq, seq_name) {
  target <- dna_seq[[seq_name]]
  codons <- Biostrings::trinucleotideFrequency(target, step = 3)
  total <- sum(codons)
  rel_freq <- ifelse(total > 0, codons / total, 0)
  data.frame(Codon = names(codons), Frequency = as.vector(codons), Relative_Frequency = round(rel_freq, 4), stringsAsFactors = FALSE)
}

compute_fast_gc_at_bias <- function(dna_seq, seq_name, window_size = 100, step_size = 50) {
  target <- dna_seq[[seq_name]]
  seq_len <- length(target)
  if (seq_len < window_size) {
    window_size <- seq_len
    step_size <- max(1, floor(seq_len / 10))
  }
  
  starts <- seq(1, seq_len - window_size + 1, by = step_size)
  ends <- starts + window_size - 1
  views <- Views(target, start = starts, end = ends)
  freqs <- letterFrequency(views, letters = c("A", "C", "G", "T"), OR = 0)
  
  a <- freqs[, "A"]
  c <- freqs[, "C"]
  g <- freqs[, "G"]
  t <- freqs[, "T"]
  
  gc_skew <- ifelse((g + c) > 0, (g - c) / (g + c), 0)
  at_bias <- ifelse((a + t) > 0, (a - t) / (a + t), 0)
  midpoints <- (starts + ends) / 2
  
  data.frame(Position = midpoints, GC_Skew = gc_skew, AT_Bias = at_bias)
}

compute_fast_aa_comp <- function(dna_seq, seq_name) {
  target <- dna_seq[[seq_name]]
  aa_seq <- suppressWarnings(Biostrings::translate(target, if.fuzzy.codon = "solve"))
  aa_counts <- Biostrings::alphabetFrequency(aa_seq)
  valid_aa <- aa_counts[names(aa_counts) %in% c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")]
  data.frame(Amino_Acid = names(valid_aa), Count = as.vector(valid_aa), stringsAsFactors = FALSE)
}

compute_divergence_times <- function(tree, dist_matrix, mut_rate = 1e-3) {
  taxa <- tree$tip.label
  n_taxa <- length(taxa)
  
  rows <- list()
  idx <- 1
  for (i in 1:(n_taxa - 1)) {
    for (j in (i + 1):n_taxa) {
      t1 <- taxa[i]
      t2 <- taxa[j]
      
      d <- if (!is.null(dist_matrix) && t1 %in% rownames(dist_matrix) && t2 %in% colnames(dist_matrix)) {
        dist_matrix[t1, t2]
      } else {
        cophenetic(tree)[t1, t2]
      }
      
      est_years <- (d / (2 * mut_rate))
      
      rows[[idx]] <- data.frame(
        Species_1 = t1,
        Species_2 = t2,
        Genetic_Distance = round(d, 4),
        Estimated_TMRCA_Years = format(round(est_years, 1), big.mark = ","),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, rows)
}

# --- 3. UI Definition ---

ui <- page_navbar(
  title = "PhyloBridge v1.0 | Computational Genomics",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50", secondary = "#18BC9C"),
  id = "navbar",
  
  header = tags$head(
    useShinyjs(),
    tags$script("Shiny.addCustomMessageHandler('allow_reconnect', function(msg) { Shiny.shinyapp.reconnect(); });"),
    tags$style(HTML("
      .sidebar { background-color: #fcfcfc; } 
      .panel { border-radius: 4px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
      .sidebar-scroll-box { max-height: 85vh; overflow-y: auto; padding-right: 5px; }
      .qa-badge-box { background-color: #f8f9fa; border: 1px solid #e9ecef; border-radius: 6px; padding: 15px; text-align: center; margin-bottom: 15px; }
      .qa-badge-title { font-size: 12px; font-weight: bold; color: #7f8c8d; text-transform: uppercase; }
      .qa-badge-value { font-size: 20px; font-weight: bold; color: #2C3E50; margin-top: 5px; }
      .btn-run-phylo { background-color: #18BC9C !important; color: white !important; font-weight: bold !important; border: none !important; font-size: 15px !important; }
      .btn-run-phylo:hover { background-color: #128C7E !important; }
    "))
  ),
  
  # Home & Upload Tab
  nav_panel(
    title = "Home", icon = icon("home"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Data Input",
        div(
          class = "sidebar-scroll-box",
          fileInput("fastaFile", "Upload FASTA File (.fasta, .fa, .fna)", accept = c(".fasta", ".fa", ".fna", ".txt")),
          actionButton("load_demo", "Load Demo Dataset", class = "btn-outline-primary w-100 mb-2"),
          downloadButton("download_sample_fasta", "Download Sample FASTA", class = "btn-sm btn-outline-secondary w-100")
        )
      ),
      div(
        class = "panel panel-default", style = "margin-top: 10px;",
        div(class = "panel-heading", tags$h4("Welcome to PhyloBridge", style = "margin:0; font-weight:bold;")),
        div(
          class = "panel-body",
          p("PhyloBridge is a high-performance, interactive bioinformatics platform for genomic data exploration, sequence composition profiling, phylogenetic tree reconstruction, and functional analysis.")
        )
      ),
      br(),
      div(
        class = "panel panel-default",
        div(class = "panel-heading", tags$h5("Interactive Assembly Quality & Contig Metrics Dashboard", style = "margin:0; font-weight:bold;")),
        div(
          class = "panel-body",
          uiOutput("home_assembly_dashboard_ui")
        )
      )
    )
  ),
  
  # GC Content Analysis Tab
  nav_panel(
    title = "GC Content Analysis", icon = icon("circle-up"),
    div(
      class = "panel panel-default", style = "margin: 20px;",
      div(class = "panel-heading", tags$b("GC Content Data Table & Bar Plot")),
      div(
        class = "panel-body",
        DTOutput("gcContentTable"),
        br(),
        plotlyOutput("gcContentPlot", height = "450px")
      ),
      div(
        class = "panel-footer",
        downloadButton("downloadGC", "Download GC Results (CSV)", class = "btn-sm btn-success")
      )
    )
  ),
  
  # Codon Usage Analysis Tab
  nav_panel(
    title = "Codon Usage Analysis", icon = icon("pen-to-square"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Controls",
        div(
          class = "sidebar-scroll-box",
          selectInput("codon_seq_select", "Select Sequence:", choices = NULL)
        )
      ),
      div(
        class = "panel panel-default",
        div(class = "panel-heading", tags$b("Codon Usage Data Table & Plot")),
        div(
          class = "panel-body",
          DTOutput("codonUsageTable"),
          br(),
          plotlyOutput("codonUsagePlot", height = "480px")
        ),
        div(
          class = "panel-footer",
          downloadButton("downloadCodonCSV", "Download Codon Usage (CSV)", class = "btn-sm btn-success")
        )
      )
    )
  ),
  
  # GC Skew & AT Bias Tab
  nav_panel(
    title = "GC Skew & AT Bias", icon = icon("gears"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Controls",
        div(
          class = "sidebar-scroll-box",
          selectInput("skew_seq_select", "Select Sequence:", choices = NULL),
          numericInput("window_size_input", "Window Size (bp):", value = 100, min = 10, step = 10),
          numericInput("step_size_input", "Step Size (bp):", value = 50, min = 5, step = 5)
        )
      ),
      div(
        class = "panel panel-default",
        div(class = "panel-heading", tags$b("GC Skew & AT Bias Results Table")),
        div(class = "panel-body", DTOutput("gc_skew_at_bias_table")),
        br(),
        layout_columns(
          div(class = "panel panel-default", div(class = "panel-heading", tags$b("Rolling-Window GC Skew")), div(class = "panel-body", plotlyOutput("gc_skew_plot_output", height = "400px"))),
          div(class = "panel panel-default", div(class = "panel-heading", tags$b("Rolling-Window AT Bias")), div(class = "panel-body", plotlyOutput("at_bias_plot_output", height = "400px")))
        )
      )
    )
  ),
  
  # Phylogenetics Tab
  nav_panel(
    title = "Phylogenetics", icon = icon("tree"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Publication Tree Controls",
        div(
          class = "sidebar-scroll-box",
          actionButton("compute_phylo_btn", "Compute Phylogeny & Clock", class = "btn-run-phylo w-100 mb-3"),
          hr(),
          radioButtons("tree_method", "Tree Algorithm:", choices = c("Neighbour Joining" = "nj", "Maximum Likelihood" = "ml"), selected = "nj"),
          selectInput("tree_layout", "Tree Layout:", choices = c("Horizontal Phylogram" = "phylogram", "Cladogram" = "cladogram", "Fan (Circular)" = "fan", "Unrooted" = "unrooted"), selected = "phylogram"),
          selectInput("color_palette", "Branch Color Palette:", choices = c("Rainbow" = "rainbow", "Viridis" = "viridis", "Ocean Blue" = "ocean", "Flame Red" = "flame")),
          sliderInput("tip_font_size", "Tip Label Font Size:", min = 0.3, max = 1.2, value = 0.6, step = 0.05),
          sliderInput("branch_width", "Branch Line Thickness:", min = 1, max = 4, value = 2.5, step = 0.5),
          checkboxInput("show_bootstrap", "Display Bootstrap Node Support", value = TRUE),
          hr(),
          h6("Molecular Clock Parameters:"),
          numericInput("mut_rate_input", "Mutation Rate (subst/site/yr):", value = 0.001, min = 1e-9, step = 0.0001)
        )
      ),
      tabsetPanel(
        tabPanel("Colorful Phylogenetic Tree", icon = icon("diagram-project"), br(), 
                 div(
                   class = "panel panel-default",
                   div(class = "panel-heading", tags$b("Phylogenetic Tree")),
                   div(class = "panel-body", uiOutput("phylo_tree_ui")),
                   div(class = "panel-footer", downloadButton("download_tree_png", "Download Tree Image (PNG)", class = "btn-sm btn-success"))
                 )
        ),
        tabPanel("Divergence Times (TMRCA)", icon = icon("clock"), br(), 
                 div(
                   class = "panel panel-default",
                   div(class = "panel-heading", tags$b("Estimated Common Ancestor Divergence Times (Years Ago)")),
                   div(class = "panel-body", DTOutput("tmrca_table_output")),
                   div(class = "panel-footer", downloadButton("download_tmrca_csv", "Download TMRCA Estimates (CSV)", class = "btn-sm btn-success"))
                 )
        ),
        tabPanel("Distance Matrix Table", icon = icon("table"), br(), DTOutput("dist_matrix_table_output"), br(), downloadButton("download_dist_csv", "Download Distance Matrix (CSV)", class = "btn-sm btn-success")),
        tabPanel("Distance Matrix Heatmap", icon = icon("fire"), br(), plotlyOutput("heatmapPlot", height = "550px")),
        tabPanel("Sequence Network Graph", icon = icon("circle-nodes"), br(), 
                 div(
                   class = "panel panel-default",
                   div(class = "panel-heading", tags$b("Sequence Distance Network Graph")),
                   div(class = "panel-body", plotOutput("networkPlot", height = "600px")),
                   div(class = "panel-footer", downloadButton("download_network_png", "Download Network Image (PNG)", class = "btn-sm btn-success"))
                 )
        )
      )
    )
  ),
  
  # Translation Tab
  nav_panel(
    title = "Protein Translation", icon = icon("vial"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Controls",
        div(
          class = "sidebar-scroll-box",
          selectInput("trans_seq_select", "Select Sequence:", choices = NULL),
          actionButton("translate_btn", "Translate to Protein", class = "btn-primary w-100")
        )
      ),
      layout_columns(
        div(class = "panel panel-default", div(class = "panel-heading", tags$b("Translated Protein Sequence")), div(class = "panel-body", verbatimTextOutput("translated_protein_output"))),
        div(class = "panel panel-default", div(class = "panel-heading", tags$b("Amino Acid Composition")), div(class = "panel-body", plotlyOutput("aa_comp_plot", height = "400px")))
      )
    )
  ),
  
  # Nucleotide Composition Tab
  nav_panel(
    title = "Nucleotide Composition", icon = icon("chart-pie"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Controls",
        div(
          class = "sidebar-scroll-box",
          selectInput("nt_seq_select", "Select Sequence:", choices = NULL),
          actionButton("analyze_nt_btn", "Run Composition Analysis", class = "btn-primary w-100")
        )
      ),
      div(class = "panel panel-default", div(class = "panel-heading", tags$b("AT/GC Ratio & Oligonucleotide Frequencies")), div(class = "panel-body", plotlyOutput("nt_comp_plot", height = "450px")))
    )
  ),
  
  # Regulatory Motifs Tab
  nav_panel(
    title = "Promoter & Motifs", icon = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Controls",
        div(
          class = "sidebar-scroll-box",
          selectInput("reg_seq_select", "Select Sequence:", choices = NULL),
          actionButton("analyze_reg_btn", "Discover Regulatory Motifs", class = "btn-primary w-100")
        )
      ),
      div(class = "panel panel-default", div(class = "panel-heading", tags$b("Transcription Factor Binding Sites (TFBS)")), div(class = "panel-body", plotlyOutput("motif_plot", height = "450px")))
    )
  ),
  
  # Population Genetics Tab
  nav_panel(
    title = "Population Genetics", icon = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Haplotype Controls",
        div(
          class = "sidebar-scroll-box",
          actionButton("analyze_pop_btn", "Compute Population Genetics", class = "btn-primary w-100 mb-3"),
          hr(),
          selectInput("haplo_color_palette", "Node Color Palette:", choices = c("Rainbow" = "rainbow", "Viridis" = "viridis", "Ocean" = "ocean", "Sunset" = "sunset")),
          sliderInput("haplo_font_size", "Label Font Size:", min = 0.3, max = 1.5, value = 0.7, step = 0.1),
          checkboxInput("haplo_show_labels", "Display Node Labels (I, II, III...)", value = TRUE)
        )
      ),
      tabsetPanel(
        id = "population_analysis_tabs",
        tabPanel(
          "Haplotype Network & Validation", icon = icon("sitemap"),
          br(),
          div(
            class = "panel panel-default",
            div(class = "panel-heading", tags$b("Distance-Weighted Haplotype Network (pegas::haploNet)")),
            div(class = "panel-body", plotOutput("haplotype_network_plot", height = "600px")),
            div(class = "panel-footer", downloadButton("download_haplotype_plot", "Download Network Image (PNG)", class = "btn-sm btn-success"))
          ),
          br(),
          div(
            class = "panel panel-default",
            div(class = "panel-heading", tags$b("Haplotype Validation Table (FASTA Member Sequences)")),
            div(class = "panel-body", DTOutput("haplo_validation_table")),
            div(class = "panel-footer", downloadButton("download_haplo_mapping_csv", "Download Haplotype Mapping (CSV)", class = "btn-sm btn-success"))
          )
        ),
        tabPanel(
          "Diversity Indices", icon = icon("chart-bar"),
          br(),
          div(
            class = "panel panel-default",
            div(class = "panel-heading", tags$b("Nucleotide & Haplotype Diversity Indices")),
            div(class = "panel-body", verbatimTextOutput("diversity_indices_output")),
            div(class = "panel-footer", downloadButton("download_diversity", "Download Diversity Indices (TXT)", class = "btn-sm btn-success"))
          )
        ),
        tabPanel(
          "Pairwise Fst", icon = icon("table"),
          br(),
          div(
            class = "panel panel-default",
            div(class = "panel-heading", tags$b("Pairwise Fixation Index (Fst) Matrix")),
            div(class = "panel-body", DTOutput("pairwise_fst_table")),
            div(class = "panel-footer", downloadButton("download_fst", "Download Pairwise Fst (CSV)", class = "btn-sm btn-success"))
          )
        )
      )
    )
  ),
  
  # Software Transparency Tab
  nav_panel(
    title = "Software Transparency", icon = icon("code-branch"),
    div(
      class = "panel panel-default", style = "margin: 20px;",
      div(class = "panel-heading", tags$b("R & Package Versions (Reproducibility)")),
      div(
        class = "panel-body",
        p("Exact computational package versions running in backend session:"),
        DTOutput("pkg_version_table")
      )
    )
  ),
  
  # About Tab
  nav_panel(
    title = "About & Contact", icon = icon("info-circle"),
    div(
      class = "panel panel-default", style = "margin: 20px;",
      div(class = "panel-heading", tags$b("PhyloBridge Team & Contact")),
      div(
        class = "panel-body",
        h5("Lead Author & Development Team"),
        tags$ul(
          tags$li(tags$b("Venu Paritala (Lead Author):"), " Indiana Wesleyan University, Marion, IN, USA. Email: ", tags$a(href = "mailto:venubabu.paritala@myemail.indwes.edu", "venubabu.paritala@myemail.indwes.edu")),
          tags$li(tags$b("Srikanth K.:"), " VFSTR Deemed to be University, India."),
          tags$li(tags$b("Sukesh Kalva:"), " Department of Biotechnology, VFSTR, India.")
        ),
        hr(),
        h5("GitHub Repository"),
        p("Source code available at: ", tags$a(href = "https://github.com/vparitala/PhyloBridge", target = "_blank", "https://github.com/vparitala/PhyloBridge"))
      )
    )
  ),
  
  nav_spacer(),
  nav_menu(
    title = "Links", icon = icon("address-card"),
    align = "right",
    nav_item(tags$a("Posit", href = "https://posit.co", target = "_blank")),
    nav_item(tags$a("Connect (LinkedIn)", href = "https://www.linkedin.com/in/venu-paritala-505657183/", target = "_blank")),
    nav_item(tags$a("GitHub Repository", href = "https://github.com/vparitala/PhyloBridge", target = "_blank"))
  )
)

# --- 4. Server Logic ---

server <- function(input, output, session) {
  session$allowReconnect(TRUE)
  
  rv <- reactiveValues(
    dna_seq = NULL,
    current_tree = NULL,
    dist_matrix = NULL,
    bootstrap_vals = NULL,
    tmrca_df = NULL,
    pop_results = NULL
  )
  
  update_selectors <- function(seqs) {
    seq_names <- names(seqs)
    if (is.null(seq_names) || any(seq_names == "")) {
      seq_names <- paste0("Sequence_", seq_along(seqs))
      names(seqs) <- seq_names
    }
    updateSelectInput(session, "codon_seq_select", choices = seq_names)
    updateSelectInput(session, "skew_seq_select", choices = seq_names)
    updateSelectInput(session, "trans_seq_select", choices = seq_names)
    updateSelectInput(session, "nt_seq_select", choices = seq_names)
    updateSelectInput(session, "reg_seq_select", choices = seq_names)
    return(seqs)
  }
  
  # Load Demo Dataset
  observeEvent(input$load_demo, {
    set.seed(42)
    demo_seqs <- DNAStringSet(c(
      CY003785 = paste(sample(c("A","C","G","T"), 1700, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25)), collapse = ""),
      CY003272 = paste(sample(c("A","C","G","T"), 1680, replace = TRUE, prob = c(0.24, 0.26, 0.26, 0.24)), collapse = ""),
      CY000705 = paste(sample(c("A","C","G","T"), 1710, replace = TRUE, prob = c(0.26, 0.24, 0.24, 0.26)), collapse = ""),
      CY001453 = paste(sample(c("A","C","G","T"), 1690, replace = TRUE, prob = c(0.23, 0.27, 0.27, 0.23)), collapse = ""),
      CY000657 = paste(sample(c("A","C","G","T"), 1705, replace = TRUE, prob = c(0.22, 0.28, 0.28, 0.22)), collapse = "")
    ))
    rv$dna_seq <- update_selectors(demo_seqs)
    showNotification("Demo dataset loaded successfully.", type = "message")
  })
  
  # File Upload Handler
  observeEvent(input$fastaFile, {
    req(input$fastaFile)
    tryCatch({
      parsed <- readDNAStringSet(input$fastaFile$datapath)
      rv$dna_seq <- update_selectors(parsed)
      showNotification("FASTA uploaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error reading FASTA:", e$message), type = "error")
    })
  })
  
  # Sample FASTA Downloader
  output$download_sample_fasta <- downloadHandler(
    filename = function() { "sample.fasta" },
    content = function(file) {
      sample_seqs <- Biostrings::DNAStringSet(c(
        CY003785 = paste(sample(c("A","C","G","T"), 1200, replace = TRUE), collapse = ""),
        CY003272 = paste(sample(c("A","C","G","T"), 1180, replace = TRUE), collapse = ""),
        CY000705 = paste(sample(c("A","C","G","T"), 1210, replace = TRUE), collapse = "")
      ))
      Biostrings::writeXStringSet(sample_seqs, file)
    }
  )
  
  # Home Assembly Dashboard
  output$home_assembly_dashboard_ui <- renderUI({
    if (is.null(rv$dna_seq)) {
      div(style = "color:#7f8c8d; font-style:italic; padding: 20px; text-align:center;", "Awaiting FASTA upload to display interactive assembly metrics dashboard...")
    } else {
      qa_res <- calculate_assembly_stats(rv$dna_seq)
      tagList(
        layout_columns(
          div(class = "qa-badge-box", div(class = "qa-badge-title", "Total Sequences"), div(class = "qa-badge-value", qa_res$num_contigs)),
          div(class = "qa-badge-box", div(class = "qa-badge-title", "Assembly Size (bp)"), div(class = "qa-badge-value", format(qa_res$total_bp, big.mark=","))),
          div(class = "qa-badge-box", div(class = "qa-badge-title", "N50 Length"), div(class = "qa-badge-value", format(qa_res$n50_val, big.mark=","))),
          div(class = "qa-badge-box", div(class = "qa-badge-title", "L50 Contig Index"), div(class = "qa-badge-value", qa_res$l50_val)),
          div(class = "qa-badge-box", div(class = "qa-badge-title", "Overall GC Content"), div(class = "qa-badge-value", paste0(qa_res$overall_gc, "%"))),
          div(class = "qa-badge-box", div(class = "qa-badge-title", "Ambiguous Bases ('N')"), div(class = "qa-badge-value", paste0(qa_res$overall_n_pct, "%")))
        ),
        br(),
        layout_columns(
          div(class = "panel panel-default", div(class = "panel-heading", tags$b("Assembly Quality Metrics Table")), div(class = "panel-body", DTOutput("home_assembly_qa_table"))),
          div(class = "panel-default", div(class = "panel-heading", tags$b("Base Composition Distribution")), div(class = "panel-body", plotlyOutput("home_base_pie_plot", height = "380px")))
        ),
        br(),
        layout_columns(
          div(class = "panel panel-default", div(class = "panel-heading", tags$b("Contig Length Distribution")), div(class = "panel-body", plotlyOutput("home_contig_length_plot", height = "380px"))),
          div(class = "panel panel-default", div(class = "panel-heading", tags$b("Cumulative Assembly Curve (N50/N90 Thresholds)")), div(class = "panel-body", plotlyOutput("home_cumulative_plot", height = "380px")))
        )
      )
    }
  })
  
  output$home_assembly_qa_table <- renderDT({
    req(rv$dna_seq)
    stats_res <- calculate_assembly_stats(rv$dna_seq)
    datatable(stats_res$summary_table, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })
  
  output$home_base_pie_plot <- renderPlotly({
    req(rv$dna_seq)
    stats_res <- calculate_assembly_stats(rv$dna_seq)
    df_bases <- stats_res$base_counts
    plot_ly(df_bases, labels = ~Base, values = ~Count, type = 'pie', hole = 0.4,
            marker = list(colors = c('#3498DB', '#1ABC9C', '#F39C12', '#E74C3C', '#95A5A6'))) %>%
      layout(title = "Nucleotide Base Ratio (A, C, G, T, N)")
  })
  
  output$home_contig_length_plot <- renderPlotly({
    req(rv$dna_seq)
    df <- data.frame(Contig = names(rv$dna_seq), Length = width(rv$dna_seq)) %>% arrange(desc(Length))
    plot_ly(df, x = ~reorder(Contig, -Length), y = ~Length, type = 'bar', marker = list(color = '#18BC9C')) %>%
      layout(xaxis = list(title = "Contig ID", showticklabels = FALSE), yaxis = list(title = "Length (bp)"))
  })
  
  output$home_cumulative_plot <- renderPlotly({
    req(rv$dna_seq)
    lengths_vec <- sort(as.numeric(width(rv$dna_seq)), decreasing = TRUE)
    cum_size <- cumsum(lengths_vec)
    total_size <- sum(lengths_vec)
    df <- data.frame(Contig_Index = seq_along(cum_size), Cumulative_Size = cum_size)
    
    plot_ly(df, x = ~Contig_Index, y = ~Cumulative_Size, type = 'scatter', mode = 'lines+markers', line = list(color = '#2980B9', width = 2)) %>%
      add_segments(x = 1, xend = max(df$Contig_Index), y = total_size * 0.5, yend = total_size * 0.5, line = list(color = '#E74C3C', dash = 'dash'), name = 'N50 Threshold (50%)') %>%
      add_segments(x = 1, xend = max(df$Contig_Index), y = total_size * 0.9, yend = total_size * 0.9, line = list(color = '#E67E22', dash = 'dot'), name = 'N90 Threshold (90%)') %>%
      layout(title = "Cumulative Assembly Size vs Contig Count", xaxis = list(title = "Contig Rank Index"), yaxis = list(title = "Cumulative Size (bp)"))
  })
  
  # GC Content DataTable & Plot
  output$gcContentTable <- renderDT({
    req(rv$dna_seq)
    df <- compute_fast_gc(rv$dna_seq)
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$gcContentPlot <- renderPlotly({
    req(rv$dna_seq)
    df <- compute_fast_gc(rv$dna_seq)
    plot_ly(df, x = ~Sequence, y = ~GC_Content, type = 'bar', marker = list(color = '#3498DB')) %>%
      layout(title = "GC Content Analysis Across Sequences", yaxis = list(title = "GC Content (%)", range = c(0, 100)), xaxis = list(showticklabels = FALSE))
  })
  
  output$downloadGC <- downloadHandler(
    filename = function() { paste0("GC_Content_Analysis_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$dna_seq)
      write.csv(compute_fast_gc(rv$dna_seq), file, row.names = FALSE)
    }
  )
  
  # Codon Usage DataTable & Plot
  output$codonUsageTable <- renderDT({
    req(rv$dna_seq, input$codon_seq_select)
    df <- compute_fast_codon_usage(rv$dna_seq, input$codon_seq_select)
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$codonUsagePlot <- renderPlotly({
    req(rv$dna_seq, input$codon_seq_select)
    df <- compute_fast_codon_usage(rv$dna_seq, input$codon_seq_select)
    plot_ly(df, x = ~Codon, y = ~Relative_Frequency, type = 'bar', marker = list(color = '#34495E')) %>%
      layout(title = paste("Codon Usage for", input$codon_seq_select), yaxis = list(title = "Relative Frequency"))
  })
  
  output$downloadCodonCSV <- downloadHandler(
    filename = function() { paste0("Codon_Usage_", input$codon_seq_select, "_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$dna_seq, input$codon_seq_select)
      write.csv(compute_fast_codon_usage(rv$dna_seq, input$codon_seq_select), file, row.names = FALSE)
    }
  )
  
  # GC Skew & AT Bias DataTable & Plots
  output$gc_skew_at_bias_table <- renderDT({
    req(rv$dna_seq, input$skew_seq_select)
    df <- compute_fast_gc_at_bias(rv$dna_seq, input$skew_seq_select, input$window_size_input, input$step_size_input)
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$gc_skew_plot_output <- renderPlotly({
    req(rv$dna_seq, input$skew_seq_select)
    df <- compute_fast_gc_at_bias(rv$dna_seq, input$skew_seq_select, input$window_size_input, input$step_size_input)
    plot_ly(df, x = ~Position, y = ~GC_Skew, type = 'scatter', mode = 'lines', line = list(color = '#2980B9')) %>%
      layout(title = "GC Skew Analysis", xaxis = list(title = "Position (bp)"), yaxis = list(title = "GC Skew"))
  })
  
  output$at_bias_plot_output <- renderPlotly({
    req(rv$dna_seq, input$skew_seq_select)
    df <- compute_fast_gc_at_bias(rv$dna_seq, input$skew_seq_select, input$window_size_input, input$step_size_input)
    plot_ly(df, x = ~Position, y = ~AT_Bias, type = 'scatter', mode = 'lines', line = list(color = '#E74C3C')) %>%
      layout(title = "AT Bias Analysis", xaxis = list(title = "Position (bp)"), yaxis = list(title = "AT Bias"))
  })
  
  # Dynamic UI for Phylogenetic Tree
  output$phylo_tree_ui <- renderUI({
    req(rv$current_tree)
    num_tips <- length(rv$current_tree$tip.label)
    calc_height <- max(500, num_tips * 20)
    plotOutput("phyloTreePlot", height = paste0(calc_height, "px"))
  })
  
  # Fast Phylogenetics Computation Handler
  observeEvent(input$compute_phylo_btn, {
    req(rv$dna_seq)
    validate(need(length(rv$dna_seq) >= 3, "At least 3 sequences required to compute phylogeny."))
    
    withProgress(message = "Computing Distance Matrix, Tree & TMRCA...", value = 0.4, {
      min_len <- min(width(rv$dna_seq))
      trimmed <- Biostrings::subseq(rv$dna_seq, start = 1, end = min_len)
      align_bin <- as.DNAbin(trimmed)
      dist_mat <- dist.dna(align_bin, model = ifelse(input$tree_method == "nj", "raw", "JC69"), pairwise.deletion = TRUE)
      
      rv$dist_matrix <- as.matrix(dist_mat)
      rv$current_tree <- nj(dist_mat)
      
      m_rate <- if (!is.null(input$mut_rate_input) && is.numeric(input$mut_rate_input) && input$mut_rate_input > 0) input$mut_rate_input else 0.001
      rv$tmrca_df <- compute_divergence_times(rv$current_tree, rv$dist_matrix, mut_rate = m_rate)
      
      if (isTRUE(input$show_bootstrap)) {
        incProgress(0.3, detail = "Running 100 Bootstrap Replicates...")
        tryCatch({
          f_boot <- function(x) nj(dist.dna(x, model = "raw", pairwise.deletion = TRUE))
          rv$bootstrap_vals <- boot.phylo(f_boot, align_bin, f = f_boot, B = 100, quiet = TRUE)
        }, error = function(e) { rv$bootstrap_vals <- NULL })
      } else {
        rv$bootstrap_vals <- NULL
      }
      
      output$dist_matrix_table_output <- renderDT({
        req(rv$dist_matrix)
        datatable(round(rv$dist_matrix, 4), options = list(pageLength = 10, scrollX = TRUE))
      })
      
      output$tmrca_table_output <- renderDT({
        req(rv$tmrca_df)
        datatable(rv$tmrca_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
      })
      
      # Base ape Engine for Precise Phylogenetic Tree Rendering
      output$phyloTreePlot <- renderPlot({
        req(rv$current_tree)
        tree <- rv$current_tree
        num_tips <- length(tree$tip.label)
        num_edges <- Nedge(tree)
        
        branch_cols <- switch(input$color_palette,
                              "rainbow" = rainbow(num_edges),
                              "viridis" = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(num_edges),
                              "ocean"   = colorRampPalette(c("#1B4F72", "#2980B9", "#A9CCE3"))(num_edges),
                              "flame"   = colorRampPalette(c("#78281F", "#E74C3C", "#F5CBA7"))(num_edges),
                              rainbow(num_edges)
        )
        
        tip_cols <- colorRampPalette(c("#2C3E50", "#16A085", "#D35400"))(num_tips)
        
        par(mar = c(2, 2, 3, 10))
        
        plot(tree, 
             type = input$tree_layout, 
             main = paste("Publication-Grade Phylogenetic Tree (", toupper(input$tree_method), ")"), 
             font = 2, 
             cex = input$tip_font_size, 
             edge.color = branch_cols, 
             edge.width = input$branch_width, 
             tip.color = tip_cols,
             no.margin = FALSE)
        
        add.scale.bar(cex = 0.8, col = "darkblue")
        
        if (!is.null(rv$bootstrap_vals) && isTRUE(input$show_bootstrap)) {
          nodelabels(rv$bootstrap_vals, adj = c(1.2, -0.2), frame = "none", cex = 0.6, font = 2, col = "#900C3F")
        }
      })
      
      output$heatmapPlot <- renderPlotly({
        req(rv$dist_matrix)
        plot_ly(z = rv$dist_matrix, x = rownames(rv$dist_matrix), y = colnames(rv$dist_matrix), type = "heatmap", colorscale = "Viridis") %>%
          layout(title = "Pairwise Genetic Distance Heatmap")
      })
      
      output$networkPlot <- renderPlot({
        req(rv$dist_matrix)
        graph <- graph_from_adjacency_matrix(rv$dist_matrix, mode = "undirected", weighted = TRUE)
        plot(graph, layout = layout_with_fr(graph), vertex.size = 12, vertex.color = "#E74C3C", vertex.label.color = "#2C3E50",
             edge.color = "#F1C40F", main = "Haplotype / Sequence Network Graph")
      })
    })
  })
  
  # Download Handlers (UPDATED WITH COLOR MATRIX MAPPINGS)
  # Download Handler: Phylogenetic Tree PNG (Dynamic Height Scaling)
  output$download_tree_png <- downloadHandler(
    filename = function() { paste0("Phylogenetic_Tree_", Sys.Date(), ".png") },
    content = function(file) {
      req(rv$current_tree)
      tree <- rv$current_tree
      num_tips <- length(tree$tip.label)
      num_edges <- Nedge(tree)
      
      # Resolution setting
      img_res <- 150
      
      # Dynamically calculate height based on tip count (matching UI expansion ratio)
      calc_height_px <- max(1200, num_tips * 30 * (img_res / 96))
      calc_width_px <- 1600
      
      branch_cols <- switch(input$color_palette,
                            "rainbow" = rainbow(num_edges),
                            "viridis" = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(num_edges),
                            "ocean"   = colorRampPalette(c("#1B4F72", "#2980B9", "#A9CCE3"))(num_edges),
                            "flame"   = colorRampPalette(c("#78281F", "#E74C3C", "#F5CBA7"))(num_edges),
                            rainbow(num_edges)
      )
      tip_cols <- colorRampPalette(c("#2C3E50", "#16A085", "#D35400"))(num_tips)
      
      png(file, width = calc_width_px, height = calc_height_px, res = img_res)
      par(mar = c(3, 2, 3, 12))
      
      plot(tree, 
           type = input$tree_layout, 
           main = paste("Publication-Grade Phylogenetic Tree (", toupper(input$tree_method), ")"), 
           font = 2, 
           cex = input$tip_font_size, 
           edge.width = input$branch_width, 
           edge.color = branch_cols, 
           tip.color = tip_cols, 
           no.margin = FALSE)
      
      add.scale.bar(cex = 0.8, col = "darkblue")
      
      if (!is.null(rv$bootstrap_vals) && isTRUE(input$show_bootstrap)) {
        nodelabels(rv$bootstrap_vals, adj = c(1.2, -0.2), frame = "none", cex = 0.6, font = 2, col = "#900C3F")
      }
      dev.off()
    }
  )
  
  output$download_network_png <- downloadHandler(
    filename = function() { paste0("Sequence_Network_", Sys.Date(), ".png") },
    content = function(file) {
      req(rv$dist_matrix)
      png(file, width = 1000, height = 1000, res = 150)
      graph <- graph_from_adjacency_matrix(rv$dist_matrix, mode = "undirected", weighted = TRUE)
      plot(graph, layout = layout_with_fr(graph), vertex.size = 12, vertex.color = "#E74C3C", vertex.label.color = "#2C3E50",
           edge.color = "#F1C40F", main = "Sequence Distance Network Graph Export")
      dev.off()
    }
  )
  
  # Population Genetics Analysis Handler
  observeEvent(input$analyze_pop_btn, {
    req(rv$dna_seq)
    withProgress(message = "Computing Population Genetics & Haplotype Network...", value = 0.4, {
      tryCatch({
        dna_mat <- as.DNAbin(rv$dna_seq)
        
        nuc_div <- tryCatch(pegas::nuc.div(dna_mat), error = function(e) 0)
        hap_div <- tryCatch(pegas::hap.div(dna_mat), error = function(e) 0)
        
        haplotypes <- pegas::haplotype(dna_mat)
        haplo_net <- pegas::haploNet(haplotypes)
        
        hap_index <- attr(haplotypes, "index")
        raw_freq <- attr(haplotypes, "freq")
        freq_vec <- if (!is.null(raw_freq)) as.numeric(raw_freq) else rep(1, max(1, length(haplotypes)))
        
        seq_names <- names(rv$dna_seq)
        if (is.null(seq_names) || length(seq_names) == 0) seq_names <- paste0("Seq_", seq_along(rv$dna_seq))
        total_seqs <- length(seq_names)
        
        mapping_rows <- list()
        num_h <- length(haplotypes)
        
        if (num_h > 0) {
          for (h_idx in seq_len(num_h)) {
            member_indices <- if (is.list(hap_index) && h_idx <= length(hap_index)) {
              hap_index[[h_idx]]
            } else if (is.vector(hap_index) && h_idx <= length(hap_index)) {
              hap_index[h_idx]
            } else {
              integer(0)
            }
            
            m_names <- if (length(member_indices) > 0 && all(member_indices <= total_seqs)) {
              paste(seq_names[member_indices], collapse = ", ")
            } else {
              "N/A"
            }
            
            h_label <- paste("Haplotype", as.roman(h_idx))
            h_freq <- if (h_idx <= length(freq_vec)) freq_vec[h_idx] else 1
            h_pct <- paste0(round((h_freq / max(1, total_seqs)) * 100, 1), "%")
            
            mapping_rows[[h_idx]] <- data.frame(
              Haplotype_ID = h_label,
              Frequency = h_freq,
              Percentage_Share = h_pct,
              Member_Sequences = m_names,
              stringsAsFactors = FALSE
            )
          }
        }
        
        haplo_mapping_df <- if (length(mapping_rows) > 0) do.call(rbind, mapping_rows) else data.frame(Haplotype_ID = character(0), Frequency = integer(0), Percentage_Share = character(0), Member_Sequences = character(0))
        
        diversity_indices <- paste0(
          "Population Genetics Summary Metrics:\n",
          "-------------------------------------\n",
          "Nucleotide Diversity (pi): ", round(nuc_div, 4), "\n",
          "Haplotype Diversity (Hd): ", round(hap_div, 4), "\n",
          "Total Analyzed Sequences: ", total_seqs, "\n",
          "Number of Unique Haplotypes: ", length(haplotypes)
        )
        
        pairwise_fst_df <- tryCatch({
          pop_names <- ifelse(grepl("_", seq_names), sapply(strsplit(seq_names, "_"), `[`, 1), "Pop1")
          pop_factor <- factor(pop_names)
          
          if (length(levels(pop_factor)) > 1) {
            genind_obj <- DNAbin2genind(dna_mat, pop = pop_factor)
            genind_df <- genind2hierfstat(genind_obj)
            pairwise_fst <- pairwise.WCfst(genind_df)
            as.data.frame(as.matrix(pairwise_fst))
          } else {
            data.frame(Info = "Pairwise Fst requires at least 2 distinct populations (e.g. Pop1_seq1, Pop2_seq1 FASTA headers).")
          }
        }, error = function(e) {
          data.frame(Info = "Pairwise Fst calculation omitted (single population detected).")
        })
        
        rv$pop_results <- list(
          diversity = diversity_indices,
          pairwise_fst = pairwise_fst_df,
          haplotype_network = haplo_net,
          haplo_freq = freq_vec,
          haplo_mapping = haplo_mapping_df
        )
        
        output$diversity_indices_output <- renderText({ diversity_indices })
        
        output$pairwise_fst_table <- renderDT({
          datatable(pairwise_fst_df, options = list(pageLength = 10), rownames = TRUE)
        })
        
        output$haplo_validation_table <- renderDT({
          datatable(haplo_mapping_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
        })
        
        # Native pegas::haploNet Render Engine
        output$haplotype_network_plot <- renderPlot({
          req(rv$pop_results$haplotype_network)
          net <- rv$pop_results$haplotype_network
          freqs <- rv$pop_results$haplo_freq
          
          font_size <- if (!is.null(input$haplo_font_size) && is.numeric(input$haplo_font_size)) input$haplo_font_size else 0.7
          show_labs <- if (!is.null(input$haplo_show_labels)) isTRUE(input$haplo_show_labels) else TRUE
          
          palette_choice <- if (!is.null(input$haplo_color_palette)) input$haplo_color_palette else "rainbow"
          bg_colors <- switch(palette_choice,
                              "rainbow" = rainbow(length(freqs)),
                              "viridis" = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(length(freqs)),
                              "ocean"   = colorRampPalette(c("#1B4F72", "#2980B9", "#A9CCE3"))(length(freqs)),
                              "sunset"  = colorRampPalette(c("#78281F", "#E74C3C", "#F5CBA7"))(length(freqs)),
                              rainbow(length(freqs))
          )
          
          par(mar = c(2, 2, 3, 2))
          
          plot(
            net, 
            show.mutation = 2, 
            labels = show_labs,
            col = "#2C3E50", 
            bg = bg_colors, 
            font = 2, 
            cex = font_size,
            main = "Distance-Weighted Haplotype Network (pegas::haploNet)"
          )
        })
        
        showNotification("Population genetics & haplotype mapping generated!", type = "message")
      }, error = function(e) {
        showNotification(paste("Population Genetics Error:", e$message), type = "error")
      })
    })
  })
  
  # Image Download Handler: Haplotype Network PNG (UPDATED WITH COLOR MATRIX MAPPINGS)
  output$download_haplotype_plot <- downloadHandler(
    filename = function() { paste0("Haplotype_Network_", Sys.Date(), ".png") },
    content = function(file) {
      req(rv$pop_results$haplotype_network)
      net <- rv$pop_results$haplotype_network
      freqs <- rv$pop_results$haplo_freq
      
      font_size <- if (!is.null(input$haplo_font_size) && is.numeric(input$haplo_font_size)) input$haplo_font_size else 0.7
      show_labs <- if (!is.null(input$haplo_show_labels)) isTRUE(input$haplo_show_labels) else TRUE
      
      palette_choice <- if (!is.null(input$haplo_color_palette)) input$haplo_color_palette else "rainbow"
      bg_colors <- switch(palette_choice,
                          "rainbow" = rainbow(length(freqs)),
                          "viridis" = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(length(freqs)),
                          "ocean"   = colorRampPalette(c("#1B4F72", "#2980B9", "#A9CCE3"))(length(freqs)),
                          "sunset"  = colorRampPalette(c("#78281F", "#E74C3C", "#F5CBA7"))(length(freqs)),
                          rainbow(length(freqs))
      )
      
      png(file, width = 1000, height = 1000, res = 150)
      par(mar = c(2, 2, 3, 2))
      plot(
        net, 
        show.mutation = 2, 
        labels = show_labs,
        col = "#2C3E50", 
        bg = bg_colors, 
        font = 2, 
        cex = font_size,
        main = "Haplotype Network Export"
      )
      dev.off()
    }
  )
  
  # Download Handler: Diversity Indices Text
  output$download_diversity <- downloadHandler(
    filename = function() { paste0("Diversity_Indices_", Sys.Date(), ".txt") },
    content = function(file) {
      req(rv$pop_results$diversity)
      writeLines(rv$pop_results$diversity, file)
    }
  )
  
  # Download Handler: Pairwise Fst CSV
  output$download_fst <- downloadHandler(
    filename = function() { paste0("Pairwise_Fst_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$pop_results$pairwise_fst)
      write.csv(rv$pop_results$pairwise_fst, file, row.names = TRUE)
    }
  )
  
  # Fast Translation Handler
  observeEvent(input$translate_btn, {
    req(rv$dna_seq, input$trans_seq_select)
    target <- rv$dna_seq[[input$trans_seq_select]]
    prot_seq <- Biostrings::translate(target, if.fuzzy.codon = "solve")
    
    output$translated_protein_output <- renderText({
      paste0("Translated Sequence (", input$trans_seq_select, "):\n\n", as.character(prot_seq))
    })
    
    output$aa_comp_plot <- renderPlotly({
      df <- compute_fast_aa_comp(rv$dna_seq, input$trans_seq_select)
      plot_ly(df, x = ~Amino_Acid, y = ~Count, type = 'bar', marker = list(color = '#1ABC9C')) %>%
        layout(title = "Amino Acid Composition", yaxis = list(title = "Frequency"))
    })
  })
  
  # Fast Nucleotide Composition Handler
  observeEvent(input$analyze_nt_btn, {
    req(rv$dna_seq, input$nt_seq_select)
    target <- rv$dna_seq[[input$nt_seq_select]]
    di_counts <- Biostrings::dinucleotideFrequency(target)
    df <- data.frame(Dinucleotide = names(di_counts), Count = as.vector(di_counts))
    
    output$nt_comp_plot <- renderPlotly({
      plot_ly(df, x = ~Dinucleotide, y = ~Count, type = 'bar', marker = list(color = '#9B59B6')) %>%
        layout(title = paste("Dinucleotide Frequencies for", input$nt_seq_select))
    })
  })
  
  # Fast Promoter Motif Handler
  observeEvent(input$analyze_reg_btn, {
    req(rv$dna_seq, input$reg_seq_select)
    motifs <- c("TATA" = "TATAAA", "CAAT" = "GGCCAATCT", "GC" = "GGGCGG")
    target <- rv$dna_seq[[input$reg_seq_select]]
    counts <- sapply(motifs, function(m) length(matchPattern(m, target, max.mismatch = 1)))
    df <- data.frame(Motif = names(counts), Count = as.vector(counts))
    
    output$motif_plot <- renderPlotly({
      plot_ly(df, x = ~Motif, y = ~Count, type = 'bar', marker = list(color = '#E67E22')) %>%
        layout(title = "Transcription Factor Binding Sites (TFBS)")
    })
  })
  
  # Software Package Version Table
  output$pkg_version_table <- renderDT({
    pkgs <- c("shiny", "bslib", "Biostrings", "ape", "phangorn", "msa", "pegas", "adegenet", "seqinr", "DECIPHER", "bio3d", "DT")
    ver_list <- sapply(pkgs, function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "Not Loaded"
    })
    df <- data.frame(
      Package = c("R Kernel System", names(ver_list)),
      Version = c(paste(R.version$major, R.version$minor, sep = "."), as.vector(ver_list)),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 15, dom = 't'), rownames = FALSE)
  })
}

# --- 5. Run Application ---
shinyApp(ui = ui, server = server)
