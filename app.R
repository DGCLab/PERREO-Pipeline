library(shiny)
library(tools)

# Allow uploading FASTQ files without size limits
options(shiny.maxRequestSize = Inf)

# Ensure a usable browser command is available when launching Shiny apps
set_browser_option <- function() {
  current_browser <- getOption("browser")
  if (is.null(current_browser) || identical(current_browser, "")) {
    fallback <- Sys.getenv("BROWSER")
    if (identical(fallback, "")) {
      fallback <- "xdg-open"
    }
    options(browser = fallback)
  }
}

set_browser_option()

# Always launch Shiny apps in the default web browser
options(shiny.launch.browser = TRUE)


app_dir <- local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  arg_match <- grep(file_flag, cmd_args)
  
  if (length(arg_match) > 0) {
    return(dirname(normalizePath(sub(file_flag, "", cmd_args[arg_match[1]]))))
  }
  
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  
  normalizePath(getwd())
})

addResourcePath("assets", file.path(app_dir, "www"))


theme <- bslib::bs_theme(
  version = 5,
  primary = "#BF5E49",
  secondary = "#54728C",
  bg = "#D9D9D9",
  fg = "#142240",
  base_font = "Inter, system-ui, -apple-system, Segoe UI, sans-serif"
)

tooltip_icon <- function(text) {
  tags$span(
    class = "info-icon",
    icon("circle-question"),
    tabindex = "0",
    `data-bs-toggle` = "tooltip",
    `data-bs-placement` = "right",
    title = text
  )
}

ui <- tagList(
  tags$head(
    tags$style(HTML(
      "
        .info-icon { color: #0b5fa5; font-size: 14px; margin-left: 6px; cursor: pointer; display: inline-flex; align-items: center; }
        .info-icon .fa { font-size: 14px; }
      "
    )),
    tags$script(src = "assets/tooltip.js")
    
  ),
  navbarPage(
    title = NULL,
    id = "main_tabs",
    inverse = FALSE,
    theme = theme,
    windowTitle = "Analysis Pipeline",
    tabPanel(
      title = "Home",
      value = "home",
      fluidPage(
        tags$head(
          tags$style(HTML(
            "
            body {
            background-color: #455F80;
            background-image:
              radial-gradient(circle at 20% 15%, rgba(0, 114, 178, 0.08), transparent 32%),
              radial-gradient(circle at 85% 10%, rgba(12, 132, 160, 0.07), transparent 30%),
              url('assets/background.svg');
            background-size: 420px 420px, 520px 520px, cover;
            background-repeat: no-repeat;
            background-attachment: fixed;
            background-position: 12% 12%, 88% 8%, center;
          }
          .hero {
            text-align: center;
            padding: 48px 24px;
            color: #fff;
            text-shadow: 0 2px 6px rgba(0,0,0,0.35);
            margin-bottom: 22px;
          }
          .hero h2 { font-weight: 800; letter-spacing: 0.2px; }
          .hero p, .hero .text-muted { font-weight: 600; color: #fff !important; }
          .hero-logo { max-height: 200px; margin-bottom: 15px; }
          .hero .badge { background: rgba(0,0,0,0.55); color: #fff; border: 1px solid rgba(255,255,255,0.3); }
          .hero hr { border-top: 1px solid rgba(255,255,255,0.4); }
          .section { background: white; border-radius: 10px; padding: 24px; margin-bottom: 24px; box-shadow: 0 4px 12px rgba(0,0,0,0.05); }
          .section h3 { margin-top: 0; }
          .contact a { color: #0072b2; font-weight: 600; }
          .contact-section { display: flex; gap: 18px; align-items: flex-start; justify-content: space-between; flex-wrap: wrap; }
          .contact-details { display: flex; flex-direction: column; gap: 12px; flex: 1; min-width: 220px; }
          .contact-list { display: flex; flex-direction: column; gap: 10px; margin-top: 4px; }
          .contact-row { display: flex; align-items: center; gap: 10px; }
          .contact-row .fa { color: #0b5fa5; font-size: 18px; width: 20px; text-align: center; }
          .contact-logo-large { max-width: 200px; width: 40%; height: auto; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.14); background: white; }
          .contact-logo-wrap { display: flex; justify-content: center; align-items: center; padding: 6px; }
          .analysis-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 18px; }
          .analysis-btn { width: 100%; text-align: left; padding: 16px; border: 2px solid #e3e3e3; border-radius: 10px; background: white; box-shadow: 0 2px 6px rgba(0,0,0,0.04); transition: all 0.2s ease; font-size: 16px; }
          .analysis-btn:hover { border-color: #0072b2; box-shadow: 0 6px 12px rgba(0,0,0,0.08); }
          .analysis-btn img { max-height: 60px; margin-right: 12px; }
          .analysis-title { font-weight: 700; font-size: 17px; margin-bottom: 4px; }
          .analysis-desc { color: #555; margin: 0; }
          .analysis-body { display: flex; align-items: center; gap: 12px; }
          .hero .badge { background: #e9f2fb; color: #0b5fa5; padding: 8px 12px; border-radius: 20px; display: inline-block; margin-bottom: 12px; font-weight: 600; }
          .info-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 14px; margin-top: 12px; }
          .info-card { border: 1px solid #e8e8e8; border-radius: 10px; padding: 14px; background: #fff; }
          .info-card strong { display: block; margin-bottom: 6px; }
          "
          ))
        ),
        div(
          class = "hero",
          img(src = "assets/pipeline_logo.svg", class = "hero-logo", alt = "Pipeline logo"),
          h2("Pipeline PERREO"),
          div(class = "badge", "SR-SE · SR-PE · LR"),
          p("One-shot pipeline designed for analyzing repeat RNAs expression in oncological contexts."),
          p(class = "text-muted", "PERREO provides three models depending on the sequencing technology used: SR-SE, SR-PE, and LR."),
          tags$hr(style = "max-width: 800px;")
        ),
        fluidRow(
          column(
            width = 6,
            div(
              class = "section",
              h3("Summary"),
              p("PERREO is applicable to data derived from different types of RNA-seq experiments and technologies, 
              as it can deal with single-end and paired-end data obtained with short-reads sequencing technologies,
              and also with long-reads sequencing data. In addition, the pipeline can analyze repeat RNAs in the contexts of 
              different organisms. The only requirement is that these species contain a reference genome, and genomic 
              and repeat elements annotations."),
              div(
                class = "info-grid",
                div(class = "info-card", strong("Input"), p("FASTQ data and library metadata (SR-SE, SR-PE, or LR).")),
                div(class = "info-card", strong("Output"), p("QC metrics, alignment, quantification, and downloadable reports.")),
                div(class = "info-card", strong("Audience"), p("Sequencing labs, bioinformatics teams, and clinical support."))
              )
            )
          ),
          column(
            width = 6,
            div(
              class = "section contact",
              h3("About Us"),
              div(
                class = "contact-section",
                div(
                  class = "contact-details",
                  div(
                    class = "contact-list",
                    div(
                      class = "contact-row",
                      icon("globe"),
                      tags$a(href = "https://grupo.us.es/dgclab/", target = "_blank", "https://grupo.us.es/dgclab/")
                    ),
                    div(
                      class = "contact-row",
                      icon("envelope"),
                      tags$a(href = "mailto:soporte@perreo.org", "soporte@perreo.org")
                    ),
                    div(
                      class = "contact-row",
                      icon("github"),
                      tags$a(href = "https://github.com/DGCLab/PERREO-Pipeline", target = "_blank", "https://github.com/DGCLab/PERREO-Pipeline")
                    )
                  )
                ),
                div(
                  class = "contact-logo-wrap",
                  img(src = "assets/dgc_logo.png", class = "contact-logo-large", alt = "Dinámica del Genoma en Cáncer Logo")
                )
              )
            )
          )
        ),
        div(
          class = "section",
          h3("Select a model"),
          p("PERREO offers three models that can be run based on the technology used."),
          div(
            class = "analysis-grid",
            actionButton(
              inputId = "go_analysis1",
              class = "analysis-btn",
              label = tagList(
                div(class = "analysis-body",
                    img(src = "assets/logo_srse.svg", alt = "SR-SE model"),
                    div(
                      div(class = "analysis-title", "PERREO SR-SE"),
                      p(class = "analysis-desc", "For RNA-seq with single-end short reads.")
                    )
                )
              )
            ),
            actionButton(
              inputId = "go_analysis2",
              class = "analysis-btn",
              label = tagList(
                div(class = "analysis-body",
                    img(src = "assets/logo_srpe.svg", alt = "SR-PE model"),
                    div(
                      div(class = "analysis-title", "PERREO SR-PE"),
                      p(class = "analysis-desc", "For RNA-seq with paired-end short reads.")
                    )
                )
              )
            ),
            actionButton(
              inputId = "go_analysis3",
              class = "analysis-btn",
              label = tagList(
                div(class = "analysis-body",
                    img(src = "assets/logo_lr.svg", alt = "LR model"),
                    div(
                      div(class = "analysis-title", "PERREO LR"),
                      p(class = "analysis-desc", "For direct RNA-seq with Nanopore long reads.")
                    )
                )
              )
            )
          )
        )
      )
    ),
    tabPanel(
      title = "Usage",
      value = "usage",
      fluidPage(
        tags$head(
          tags$style(HTML(
            "
            .usage-hero {
              background: linear-gradient(120deg, #e9f2fb 0%, #f8f0ec 100%);
              border-radius: 16px;
              padding: 28px;
              box-shadow: 0 10px 20px rgba(0,0,0,0.08);
              margin-bottom: 18px;
              color: #0f203b;
            }
            .usage-hero h2 { margin-top: 0; font-weight: 800; }
            .usage-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 14px; }
            .usage-card {
              background: #ffffff;
              border-radius: 14px;
              padding: 16px;
              border: 1px solid #e5e5e5;
              box-shadow: 0 8px 18px rgba(0,0,0,0.06);
            }
            .usage-card h4 { margin-top: 0; color: #0f203b; font-weight: 700; }
            .diagram {
              background: #ffffff;
              border-radius: 16px;
              padding: 18px;
              border: 1px solid #e5e5e5;
              box-shadow: 0 8px 18px rgba(0,0,0,0.07);
              margin-top: 16px;
            }
.usage-diagram-image {
              width: 100%;
              max-width: 980px;
              height: auto;
              display: block;
              margin: 14px auto 0;
            }
            .diagram-flow { display: flex; flex-wrap: wrap; align-items: center; gap: 12px; justify-content: center; }
            .diagram-node {
              background: #0b5fa5;
              color: #fff;
              padding: 10px 14px;
              border-radius: 12px;
              font-weight: 700;
              box-shadow: 0 6px 12px rgba(0,0,0,0.12);
              min-width: 140px;
              text-align: center;
            }
            .diagram-arrow {
              font-size: 22px;
              color: #0f203b;
            }
            .sample-pills { display: flex; flex-wrap: wrap; gap: 8px; justify-content: center; margin-top: 10px; }
            .sample-pill {
              background: #f0f4ff;
              color: #0f203b;
              padding: 8px 12px;
              border-radius: 12px;
              border: 1px solid #d7e2ff;
              font-weight: 600;
              box-shadow: 0 4px 10px rgba(0,0,0,0.04);
            }
          "
          ))
        ),
        div(
          class = "usage-hero",
          h2("Applications and utilities"),
          p("PERREO adapts to diverse projects: exploratory studies, biomarker validation, or clinical monitoring. 
          The pipeline integrates QC, trimming, alignment, repeat RNA quantification, coexpression analysis and 
          transcriptome assembly in a single flow."),
          p("It includes utilities that simplify traceability: reproducible execution, downloadable reports, 
          and compatibility with multiple sequencing technologies."),
          div(
            class = "usage-grid",
            div(
              class = "usage-card",
              h4("Applications"),
              tags$ul(
                tags$li("Repeat RNA profiling in cancer studies and therapy response."),
                tags$li("Comparison of samples from different tissues or biological origins."),
                tags$li("Novel repeat transcripts detection." ),
                tags$li("Centralized quality control with exportable results for reporting."),
                tags$li("Integration into bioinformatics pipelines or clinical environments." )
              )
            ),
            div(
              class = "usage-card",
              h4("Included utilities"),
              tags$ul(
                tags$li("Direct QC execution with FastQC/MultiQC."),
                tags$li("Three execution modes (SR-SE, SR-PE, LR) aligned with sequencing technology."),
                tags$li("Metadata and reference management to tailor analyses for each species."),
                tags$li("Downloadable results and logs for auditing." )
              )
            )
          )
        ),
        div(
          class = "diagram",
          h4("Usage diagram"),
          p("Samples can come from blood, plasma, tumors, cell lines, tissues, or other sources. All converge into PERREO modalities 
          to generate QC metrics, alignments, and quantifications."),
           img(
            src = "assets/variant_logo.png",
            class = "usage-diagram-image",
            alt = "PERREO pipeline usage diagram"
          ),
          div(
            class = "sample-pills",
            span(class = "sample-pill", "Blood"),
            span(class = "sample-pill", "Plasma"),
            span(class = "sample-pill", "Tumor"),
            span(class = "sample-pill", "Cell line"),
            span(class = "sample-pill", "Tissue"),
            span(class = "sample-pill", "Other sources")
          ),
          div(
            class = "diagram-flow",
            div(class = "diagram-node", "Samples"),
            span(class = "diagram-arrow", icon("arrow-right")),
            div(class = "diagram-node", "QC (FastQC / MultiQC)"),
            span(class = "diagram-arrow", icon("arrow-right")),
            div(class = "diagram-node", "PERREO execution"),
            span(class = "diagram-arrow", icon("arrow-right")),
            div(class = "diagram-node", "Downloadable results")
          )
        )
      )
    ),
    tabPanel(
      title = "QC - FastQC",
      value = "qc_fastqc",
      fluidPage(
        tags$head(
          tags$style(HTML(
            "
            .qc-wrapper {
              background: #ffffff;
              border-radius: 16px;
              padding: 18px;
              box-shadow: 0 8px 18px rgba(0, 0, 0, 0.1);
              border: 1px solid rgba(255, 255, 255, 0.65);
            }
            .qc-wrapper h3,
            .qc-wrapper h4 {
              color: #0f203b;
              font-weight: 700;
            }
            .qc-section + .qc-section {
              margin-top: 18px;
            }
          "
          ))
        ),
        div(
          class = "qc-wrapper",
          sidebarLayout(
            sidebarPanel(
              h3("FASTQ quality control"),
              fileInput(
                inputId = "fastq_files",
                label = tagList("Upload your FASTQ files", tooltip_icon("Select FASTQ files (.fastq/.fq/.gz) for quality control.")),
                multiple = TRUE,
                accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz")
              ),
              checkboxInput(
                inputId = "overwrite_qc",
                label = tagList("Overwrite previous results", tooltip_icon("Remove prior outputs before running FastQC again.")),
                value = FALSE
              ),
              actionButton(
                inputId = "run_fastqc",
                label = tagList("Run FastQC", tooltip_icon("Launch FastQC on the selected FASTQ files."))
              ),
              br(),
              br(),
              helpText(
                "Note: this tab calls the 'fastqc' program installed in the pipeline ",
                "environment (for example, the conda environment 'perreo')."
              )
            ),
            mainPanel(
              div(
                class = "qc-section",
                h3("FastQC results"),
                tableOutput("fastqc_table")
              ),
              div(
                class = "qc-section",
                h4("HTML reports"),
                uiOutput("fastqc_links")
              )
            )
          ),
          div(
            class = "qc-section",
            tags$hr(),
            h4("MultiQC report"),
            uiOutput("multiqc_report")
          )
        )
      )
    ),
    tabPanel(
      title = "PERREO SR-SE",
      value = "analysis1",
      fluidPage(
        tags$head(
          tags$style(HTML(
            "
            .analysis-card {
              background: #ffffff;
              border-radius: 16px;
              padding: 18px;
              box-shadow: 0 8px 18px rgba(0, 0, 0, 0.1);
              border: 1px solid rgba(255, 255, 255, 0.65);
              margin-bottom: 18px;
            }
.mode-hero {
              text-align: center;
              margin-bottom: 20px;
              color: #ffffff;
              text-shadow: 0 2px 6px rgba(0, 0, 0, 0.35);
            }
            .mode-hero h2 {
              color: #ffffff;
              font-weight: 800;
            }
            .mode-hero p {
              color: #eaf2ff;
              font-weight: 600;
            }
            .analysis-card h4 {
              margin-top: 0;
              color: #0f203b;
              font-weight: 700;
            }
            .arg-grid {
              display: grid;
              grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
              gap: 12px;
            }
          "
          ))
        ),
        div(
          class = "mode-hero",
          img(
            src = "assets/logo_srse.svg",
            alt = "PERREO SR-SE logo",
            style = "max-height: 120px; margin-bottom: 10px;"
          ),
          h2("PERREO SR-SE (single-end short reads)"),
          p("For RNA-seq generated with single-end short-read technology.")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Main files"),
          div(
            class = "arg-grid",
            fileInput(
              inputId = "srse_reference_genome",
              label = tagList("Reference genome", tooltip_icon("Upload the reference genome FASTA file (.fa/.fasta).")),
              multiple = FALSE,
              accept = c(".fa", ".fasta")
            ),
            fileInput(
              inputId = "srse_genome_annotations",
              label = tagList("Genome annotations", tooltip_icon("Provide genome annotation in GTF format.")),
              multiple = FALSE,
              accept = c(".gtf")
            ),
            fileInput(
              inputId = "srse_repeat_annotations",
              label = tagList("Repeat annotations", tooltip_icon("Upload repeat element annotations (GTF).")),
              multiple = FALSE,
              accept = c("gtf")
            ),
            fileInput(
              inputId = "srse_samplesheet",
              label = tagList("Samplesheet", tooltip_icon("Samples metadata file (.txt) describing your runs.")),
              multiple = FALSE,
              accept = c(".txt")
            )
          )
        ),
           div(
          class = "arg-section analysis-card",
          h4("Project settings"),
          div(
            class = "arg-grid",
            textInput(
              inputId = "srse_project_name",
              label = tagList("Project name", tooltip_icon("Name for the analysis project output.")),
              placeholder = "Example: PERREO_SRSE_Run1"
            )
          )
        ),
        
        div(
          class = "arg-section analysis-card",
          h4("Trimming"),
          tabsetPanel(
            tabPanel(
              title = "Adapters",
              textInput(
                inputId = "srse_adapter_sequence",
                label = tagList("Introduce your sequences adapters:", tooltip_icon("Adapter sequence to trim from single-end reads.")),
                placeholder = "Example: AGATCGGAAGAGCAC"
              )
            ),
            tabPanel(
              title = "Quality threshold trimming",
              sliderInput(
                inputId = "srse_trim_threshold",
                label = tagList("Select a quality threshold for trimming", tooltip_icon("Minimum base quality score to keep during trimming.")),
                min = 5,
                max = 40,
                value = 20,
                step = 1
              )
            ),
            tabPanel(
              title = "Trimming mode",
              radioButtons(
                inputId = "srse_trim_mode",
                label = tagList("Choose trimming mode", tooltip_icon("Select between standard or extra trimming settings.")),
                choices = c("simple", "extra"),
                selected = "simple"
              )
            ),
            tabPanel(
              title = "Minimum length",
              sliderInput(
                inputId = "srse_min_length",
                label = tagList("Minimum length after trimming", tooltip_icon("Discard reads shorter than this length after trimming.")),
                min = 10,
                max = 200,
                value = 16,
                step = 1
              )
            ),
            tabPanel(
              title = "Maximum length",
              sliderInput(
                inputId = "srse_max_length",
                label = tagList("Maximum length after trimming", tooltip_icon("Discard reads longer than this length after trimming.")),
                min = 10,
                max = 300,
                value = 200,
                step = 1
              )
            ),
            tabPanel(
              title = "Initial trim",
              numericInput(
                inputId = "srse_initial_trim",
                label = tagList("Initial trim (bases)", tooltip_icon("Trim a fixed number of bases from the start of each read.")),
                min = 0,
                max = 50,
                value = 0,
                step = 1
              )
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Alignment"),
          radioButtons(
            inputId = "srse_remove_duplicates",
            label = tagList("Remove duplicates", tooltip_icon("Choose whether to mark and remove duplicate reads.")),
            choices = c("TRUE", "FALSE"),
            selected = "FALSE",
            inline = TRUE
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("DEA"),
          div(
            class = "arg-grid",
            selectInput(
              inputId = "srse_dea_method",
              label = tagList("Method", tooltip_icon("Differential expression method to apply.")),
              choices = c("DESeq2", "edgeR"),
              selected = "DESeq2"
            ),
            numericInput(
              inputId = "srse_dea_log2fc",
              label = tagList("log2FC threshold", tooltip_icon("Minimum absolute log2 fold change to report.")),
              min = 0,
              max = 4,
              value = 1,
              step = 0.1
            ),
            numericInput(
              inputId = "srse_dea_fdr",
              label = tagList("FDR threshold", tooltip_icon("False discovery rate cutoff for significant results.")),
              min = 0,
              max = 0.2,
              value = 0.05,
              step = 0.01
            ),
            selectInput(
              inputId = "srse_dea_batch",
              label = tagList("Batch effect", tooltip_icon("Indicate whether to account for batch effects.")),
              choices = c("yes", "no"),
              selected = "no"
            )
          )
        ),
           div(
          class = "arg-section analysis-card",
          h4("Prediction model"),
          div(
            class = "arg-grid",
            selectInput(
              inputId = "srse_prediction_model",
              label = tagList("Prediction model", tooltip_icon("Enable or disable prediction model processing.")),
              choices = c("yes", "no"),
              selected = "no"
            ),
            textInput(
              inputId = "srse_positive_class",
              label = tagList("Positive class", tooltip_icon("Condition to use as the positive class.")),
              placeholder = "Example: Tumor"
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execution"),
          div(
            class = "arg-grid",
            numericInput(
              inputId = "srse_threads",
              label = tagList("Threads", tooltip_icon("Number of threads to use during execution.")),
              min = 1,
              max = 64,
              value = 8,
              step = 1
            ),
            numericInput(
              inputId = "srse_mismatch_align",
              label = tagList("Mismatch align", tooltip_icon("Maximum mismatches allowed during alignment.")),
              min = 0,
              max = 10,
              value = 2,
              step = 1
            ),
            
            checkboxInput(
              inputId = "srse_polya",
              label = tagList("PolyA", tooltip_icon("Enable polyA processing flag.")),
              value = FALSE
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Command preview"),
          verbatimTextOutput("srse_command")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execute pipeline"),
          p("Run the SR-SE pipeline with the parameters shown above."),
          actionButton("run_srse_pipeline", "Run SR-SE pipeline", class = "btn-primary")
        ),
                div(
          class = "arg-section analysis-card",
          h4("Execution output"),
          verbatimTextOutput("srse_run_output")
        ),
        actionButton("back_home1", tagList("Back to home", tooltip_icon("Return to the main landing page.")))
      )
    ),
              
    tabPanel(
      title = "PERREO SR-PE",
      value = "analysis2",
      fluidPage(
        tags$head(
          tags$style(HTML(
            "
            .analysis-card {
              background: #ffffff;
              border-radius: 16px;
              padding: 18px;
              box-shadow: 0 8px 18px rgba(0, 0, 0, 0.1);
              border: 1px solid rgba(255, 255, 255, 0.65);
              margin-bottom: 18px;
            }
 .mode-hero {
              text-align: center;
              margin-bottom: 20px;
              color: #ffffff;
              text-shadow: 0 2px 6px rgba(0, 0, 0, 0.35);
            }
            .mode-hero h2 {
              color: #ffffff;
              font-weight: 800;
            }
            .mode-hero p {
              color: #eaf2ff;
              font-weight: 600;
            }
            .analysis-card h4 {
              margin-top: 0;
              color: #0f203b;
              font-weight: 700;
            }
            .arg-grid {
              display: grid;
              grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
              gap: 12px;
            }
          "
          ))
        ),
        div(
          class = "mode-hero",
          img(
            src = "assets/logo_srpe.svg",
            alt = "PERREO SR-PE logo",
            style = "max-height: 120px; margin-bottom: 10px;"
          ),
          h2("PERREO SR-PE (paired-end short reads)"),
          p("For RNA-seq generated with paired-end short-read technology.")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Main files"),
          div(
            class = "arg-grid",
            fileInput(
              inputId = "srpe_reference_genome",
              label = tagList("Reference genome", tooltip_icon("Upload the reference genome FASTA file (.fa/.fasta).")),
              multiple = FALSE,
              accept = c(".fa", ".fasta")
            ),
            fileInput(
              inputId = "srpe_genome_annotations",
              label = tagList("Genome annotations", tooltip_icon("Provide genome annotation in GTF format.")),
              multiple = FALSE,
              accept = c(".gtf")
            ),
            fileInput(
              inputId = "srpe_repeat_annotations",
              label = tagList("Repeat annotations", tooltip_icon("Upload repeat element annotations (GTF).")),
              multiple = FALSE,
              accept = c("gtf")
            ),
            fileInput(
              inputId = "srpe_samplesheet",
              label = tagList("Samplesheet", tooltip_icon("Samples metadata file (.txt) describing your runs.")),
              multiple = FALSE,
              accept = c(".txt")
            )
          )
        ),
                 div(
          class = "arg-section analysis-card",
          h4("Project settings"),
          div(
            class = "arg-grid",
            textInput(
              inputId = "srpe_project_name",
              label = tagList("Project name", tooltip_icon("Name for the analysis project output.")),
              placeholder = "Example: PERREO_SRPE_Run1"
            )
          )
        ),
        
        div(
          class = "arg-section analysis-card",
          h4("Trimming"),
          tabsetPanel(
            tabPanel(
              title = "Adapters",
              textInput(
                inputId = "srpe_adapter_sequence_r1",
                label = tagList("Introduce the adapter sequence (read 1)", tooltip_icon("Adapter sequence to trim from read 1.")),
                placeholder = "Example: AGATCGGAAGAGCAC"
              ),
              textInput(
                inputId = "srpe_adapter_sequence_r2",
                label = tagList("Introduce the adapter sequence (read 2)", tooltip_icon("Adapter sequence to trim from read 2.")),
                placeholder = "Example: AGATCGGAAGAGCGT"
                 ),
              checkboxInput(
                inputId = "srpe_no_adapters",
label = tagList("No adapters (pass empty flags)", tooltip_icon("Send -adapt_r1= and -adapt_r2= to skip adapter trimming.")),
                value = FALSE
              )
            ),
            tabPanel(
              title = "Quality threshold trimming",
              sliderInput(
                inputId = "srpe_trim_threshold",
                label = tagList("Select a quality threshold for trimming", tooltip_icon("Minimum base quality score to keep during trimming.")),
                min = 5,
                max = 40,
                value = 20,
                step = 1
              )
            ),
            tabPanel(
              title = "Trimming mode",
              radioButtons(
                inputId = "srpe_trim_mode",
                label = tagList("Choose trimming mode", tooltip_icon("Select between standard or extra trimming settings.")),
                choices = c("simple", "extra"),
                selected = "simple"
              )
            ),
            tabPanel(
              title = "Minimum length",
              sliderInput(
                inputId = "srpe_min_length",
                label = tagList("Minimum length after trimming", tooltip_icon("Discard reads shorter than this length after trimming.")),
                min = 10,
                max = 200,
                value = 16,
                step = 1
              )
            ),
            tabPanel(
              title = "Maximum length",
              sliderInput(
                inputId = "srpe_max_length",
                label = tagList("Maximum length after trimming", tooltip_icon("Discard reads longer than this length after trimming.")),
                min = 10,
                max = 300,
                value = 200,
                step = 1
              )
            ),
            tabPanel(
              title = "Initial trim",
              numericInput(
                inputId = "srpe_initial_trim_r1",
                label = tagList("Initial trim read 1", tooltip_icon("Trim a fixed number of bases from the start of read 1.")),
                min = 0,
                max = 50,
                value = 0,
                step = 1
              ),
              numericInput(
                inputId = "srpe_initial_trim_r2",
                label = tagList("Initial trim read 2", tooltip_icon("Trim a fixed number of bases from the start of read 2.")),
                min = 0,
                max = 50,
                value = 0,
                step = 1
              )
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Alignment"),
          radioButtons(
            inputId = "srpe_remove_duplicates",
            label = tagList("Remove duplicates", tooltip_icon("Choose whether to mark and remove duplicate reads.")),
            choices = c("TRUE", "FALSE"),
            selected = "FALSE",
            inline = TRUE
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("DEA"),
          div(
            class = "arg-grid",
            selectInput(
              inputId = "srpe_dea_method",
              label = tagList("Method", tooltip_icon("Differential expression method to apply.")),
              choices = c("DESeq2", "edgeR"),
              selected = "DESeq2"
            ),
            numericInput(
              inputId = "srpe_dea_log2fc",
              label = tagList("log2FC threshold", tooltip_icon("Minimum absolute log2 fold change to report.")),
              min = 0,
              max = 4,
              value = 1,
              step = 0.1
            ),
            numericInput(
              inputId = "srpe_dea_fdr",
              label = tagList("FDR threshold", tooltip_icon("False discovery rate cutoff for significant results.")),
              min = 0,
              max = 0.2,
              value = 0.05,
              step = 0.01
            ),
            selectInput(
              inputId = "srpe_dea_batch",
              label = tagList("Batch effect", tooltip_icon("Indicate whether to account for batch effects.")),
              choices = c("yes", "no"),
              selected = "no"
            )
          )
        ),
          div(
          class = "arg-section analysis-card",
          h4("Prediction model"),
          div(
            class = "arg-grid",
            selectInput(
              inputId = "srpe_prediction_model",
              label = tagList("Prediction model", tooltip_icon("Enable or disable prediction model processing.")),
              choices = c("yes", "no"),
              selected = "no"
            ),
            textInput(
              inputId = "srpe_positive_class",
              label = tagList("Positive class", tooltip_icon("Condition to use as the positive class.")),
              placeholder = "Example: Tumor"
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execution"),
          div(
            class = "arg-grid",
            numericInput(
              inputId = "srpe_threads",
              label = tagList("Threads", tooltip_icon("Number of threads to use during execution.")),
              min = 1,
              max = 64,
              value = 8,
              step = 1
            ),
            numericInput(
              inputId = "srpe_mismatch_align",
              label = tagList("Mismatch align", tooltip_icon("Maximum mismatches allowed during alignment.")),
              min = 0,
              max = 10,
              value = 2,
              step = 1
            ),
          
            checkboxInput(
              inputId = "srpe_polya",
              label = tagList("PolyA", tooltip_icon("Enable polyA processing flag.")),
              value = FALSE
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Command preview"),
          verbatimTextOutput("srpe_command")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execute pipeline"),
          p("Run the SR-PE pipeline with the parameters shown above."),
          actionButton("run_srpe_pipeline", "Run SR-PE pipeline", class = "btn-primary")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execution output"),
          verbatimTextOutput("srpe_run_output")
        ),
        actionButton("back_home2", tagList("Back to home", tooltip_icon("Return to the main landing page.")))
      )
    ),
    tabPanel(
      title = "PERREO LR",
      value = "analysis3",
      fluidPage(
        tags$head(
          tags$style(HTML(
            "
            .analysis-card {
              background: #ffffff;
              border-radius: 16px;
              padding: 18px;
              box-shadow: 0 8px 18px rgba(0, 0, 0, 0.1);
              border: 1px solid rgba(255, 255, 255, 0.65);
              margin-bottom: 18px;
            }
 .mode-hero {
              text-align: center;
              margin-bottom: 20px;
              color: #ffffff;
              text-shadow: 0 2px 6px rgba(0, 0, 0, 0.35);
            }
            .mode-hero h2 {
              color: #ffffff;
              font-weight: 800;
            }
            .mode-hero p {
              color: #eaf2ff;
              font-weight: 600;
            }
            .analysis-card h4 {
              margin-top: 0;
              color: #0f203b;
              font-weight: 700;
            }
            .arg-grid {
              display: grid;
              grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
              gap: 12px;
            }
          "
          ))
        ),
        div(
  class = "mode-hero",
          img(
            src = "assets/logo_lr.svg",
            alt = "PERREO LR logo",
            style = "max-height: 120px; margin-bottom: 10px;"
          ),
          h2("PERREO LR (long-read Nanopore)"),
          p("For direct RNA-seq generated with Nanopore long-read technology.")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Main files"),
          div(
            class = "arg-grid",
            fileInput(
              inputId = "lr_reference_genome",
              label = tagList("Reference genome", tooltip_icon("Upload the reference genome FASTA file (.fa/.fasta).")),
              multiple = FALSE,
              accept = c(".fa", ".fasta")
            ),
            fileInput(
              inputId = "lr_genome_annotations",
              label = tagList("Genome annotations", tooltip_icon("Provide genome annotation in GTF format.")),
              multiple = FALSE,
              accept = c(".gtf")
            ),
            fileInput(
              inputId = "lr_repeat_annotations",
              label = tagList("Repeat annotations", tooltip_icon("Upload repeat element annotations (GTF).")),
              multiple = FALSE,
              accept = c("gtf")
            ),
            fileInput(
              inputId = "lr_samplesheet",
              label = tagList("Samplesheet", tooltip_icon("Samples metadata file (.txt) describing your runs.")),
              multiple = FALSE,
              accept = c(".txt")
            )
          )
        ),
         div(
          class = "arg-section analysis-card",
          h4("Project settings"),
          div(
            class = "arg-grid",
            textInput(
              inputId = "lr_project_name",
              label = tagList("Project name", tooltip_icon("Name for the analysis project output.")),
              placeholder = "Example: PERREO_LR_Run1"
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("DEA"),
          div(
            class = "arg-grid",
            selectInput(
              inputId = "lr_dea_method",
              label = tagList("Method", tooltip_icon("Differential expression method to apply.")),
              choices = c("DESeq2", "edgeR"),
              selected = "DESeq2"
            ),
            numericInput(
              inputId = "lr_dea_log2fc",
              label = tagList("log2FC threshold", tooltip_icon("Minimum absolute log2 fold change to report.")),
              min = 0,
              max = 4,
              value = 1,
              step = 0.1
            ),
            numericInput(
              inputId = "lr_dea_fdr",
              label = tagList("FDR threshold", tooltip_icon("False discovery rate cutoff for significant results.")),
              min = 0,
              max = 0.2,
              value = 0.05,
              step = 0.01
            ),
            selectInput(
              inputId = "lr_dea_batch",
              label = tagList("Batch effect", tooltip_icon("Indicate whether to account for batch effects.")),
              choices = c("yes", "no"),
              selected = "no"
            )
          )
        ),
         div(
          class = "arg-section analysis-card",
          h4("Prediction model"),
          div(
            class = "arg-grid",
            selectInput(
              inputId = "lr_prediction_model",
              label = tagList("Prediction model", tooltip_icon("Enable or disable prediction model processing.")),
              choices = c("yes", "no"),
              selected = "no"
            ),
            textInput(
              inputId = "lr_positive_class",
              label = tagList("Positive class", tooltip_icon("Condition to use as the positive class.")),
              placeholder = "Example: Tumor"
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execution"),
          div(
            class = "arg-grid",
            numericInput(
              inputId = "lr_threads",
              label = tagList("Threads", tooltip_icon("Number of threads to use during execution.")),
              min = 1,
              max = 64,
              value = 8,
              step = 1
            
            )
          )
        ),
        div(
          class = "arg-section analysis-card",
          h4("Command preview"),
          verbatimTextOutput("lr_command")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execute pipeline"),
          p("Run the LR pipeline with the parameters shown above."),
          actionButton("run_lr_pipeline", "Run LR pipeline", class = "btn-primary")
        ),
        div(
          class = "arg-section analysis-card",
          h4("Execution output"),
          verbatimTextOutput("lr_run_output")
        ),
        actionButton("back_home3", tagList("Back to home", tooltip_icon("Return to the main landing page.")))
      )
    )
  )
)

server <- function(input, output, session) {
  fastqc_outdir <- file.path(getwd(), "fastqc_results")
  dir.create(fastqc_outdir, showWarnings = FALSE, recursive = TRUE)
  
  addResourcePath("fastqc_reports", fastqc_outdir)
  
  fastqc_results <- reactiveVal(NULL)
  multiqc_report <- reactiveVal(NULL)
  srse_run_output <- reactiveVal("No pipeline execution yet.")
  srpe_run_output <- reactiveVal("No pipeline execution yet.")
  lr_run_output <- reactiveVal("No pipeline execution yet.")
  
  
  file_datapath <- function(file_input) {
    if (is.null(file_input) || is.null(file_input$datapath)) {
      return(NULL)
    }
    value <- file_input$datapath[1]
    if (is.null(value) || !nzchar(value)) {
      return(NULL)
    }
    value
  }
  
  file_displayname <- function(file_input) {
    if (is.null(file_input) || is.null(file_input$name)) {
      return(NULL)
    }
    value <- file_input$name[1]
    if (is.null(value) || !nzchar(value)) {
      return(NULL)
    }
    value
  }

   strip_leading_slashes <- function(path_value) {
    if (is.null(path_value)) {
      return(NULL)
    }
    sub("^/+", "", as.character(path_value))  }

  file_runtimepath <- function(file_input) {
    datapath <- file_datapath(file_input)
    display_name <- file_displayname(file_input)
    
    if (is.null(datapath) || is.null(display_name)) {
      return(NULL)
    }
    
    target_basename <- basename(display_name)
    target_path <- file.path(getwd(), target_basename)
    copied <- file.copy(datapath, target_path, overwrite = TRUE, copy.mode = TRUE)
    
    if (isTRUE(copied) || file.exists(target_path)) {
      return(target_basename)
      }
    strip_leading_slashes(normalizePath(datapath, mustWork = FALSE))
  }
  
  
  add_arg <- function(flag, value) {
    if (is.null(value) || length(value) == 0) {
      return(NULL)
    }
    if (is.character(value)) {
      value <- trimws(value)
      if (!nzchar(value)) {
        return(NULL)
      }
    }
    if (all(is.na(value))) {
       return(NULL)
    }
    c(flag, as.character(value))
  }
  
  add_empty_arg <- function(flag) {
c(flag, "")
  }  
  format_preview_token <- function(token) {
      if (identical(token, "")) {
      return("''")
    }
    if (grepl("\\s", token)) {
      return(shQuote(token))
    }
    token
  }
  
  build_command_preview <- function(script_name, args) {
    preview_args <- vapply(args, format_preview_token, character(1))
    paste(c("bash", format_preview_token(script_name), preview_args), collapse = " ")
  }
  
  build_run_command <- function(script_name, args) {
    list(
      script = file.path(app_dir, script_name),
      args = args
    )
  }
  
  srse_command <- reactive({
    args <- c(
       add_arg("-sample_list", file_runtimepath(input$srse_samplesheet)),
      add_arg("-reference_genome", file_runtimepath(input$srse_reference_genome)),
      add_arg("-genome_gtf", file_runtimepath(input$srse_genome_annotations)),
      add_arg("-repeat_gtf", file_runtimepath(input$srse_repeat_annotations)),
      add_arg("-threads", input$srse_threads),
      add_arg("-adapter", input$srse_adapter_sequence),
      add_arg("-trimming_quality_threshold", input$srse_trim_threshold),
      add_arg("-min_length_trim", input$srse_min_length),
      add_arg("-max_length_trim", input$srse_max_length),
      add_arg("-initial_trim_read", input$srse_initial_trim),
      add_arg("-mismatch_align", input$srse_mismatch_align),
      add_arg("-project_name", input$srse_project_name),
      add_arg("-remove_duplicates", input$srse_remove_duplicates),
      add_arg("-log2FC", input$srse_dea_log2fc),
      add_arg("-FDR", input$srse_dea_fdr),
      add_arg("-batch", input$srse_dea_batch),
      add_arg("-method", input$srse_dea_method),
      add_arg("-prediction_model", input$srse_prediction_model),
      add_arg("-positive_class", input$srse_positive_class)
    )
    if (isTRUE(input$srse_polya)) {
      args <- c(args, "-polya", "polya")
    }
    build_run_command("perreo_srse.sh", args)
  })
  
  srse_command_preview <- reactive({
    args <- c(
      add_arg("-sample_list", file_displayname(input$srse_samplesheet)),
      add_arg("-reference_genome", file_displayname(input$srse_reference_genome)),
      add_arg("-genome_gtf", file_displayname(input$srse_genome_annotations)),
      add_arg("-repeat_gtf", file_displayname(input$srse_repeat_annotations)),
      add_arg("-threads", input$srse_threads),
      add_arg("-adapter", input$srse_adapter_sequence),
      add_arg("-trimming_quality_threshold", input$srse_trim_threshold),
      add_arg("-min_length_trim", input$srse_min_length),
      add_arg("-max_length_trim", input$srse_max_length),
      add_arg("-initial_trim_read", input$srse_initial_trim),
      add_arg("-mismatch_align", input$srse_mismatch_align),
      add_arg("-project_name", input$srse_project_name),
      add_arg("-remove_duplicates", input$srse_remove_duplicates),
      add_arg("-log2FC", input$srse_dea_log2fc),
      add_arg("-FDR", input$srse_dea_fdr),
      add_arg("-batch", input$srse_dea_batch),
      add_arg("-method", input$srse_dea_method),
      add_arg("-prediction_model", input$srse_prediction_model),
      add_arg("-positive_class", input$srse_positive_class)
    )
    if (isTRUE(input$srse_polya)) {
      args <- c(args, "-polya", "polya")
    }
    build_command_preview("perreo_srse.sh", args)
  })
  
  srpe_command <- reactive({
       adapter_args <- if (isTRUE(input$srpe_no_adapters)) {
      c(add_empty_arg("-adapt_r1"), add_empty_arg("-adapt_r2"))
    } else {
      c(
        add_arg("-adapt_r1", input$srpe_adapter_sequence_r1),
        add_arg("-adapt_r2", input$srpe_adapter_sequence_r2)
      )
    }
    args <- c(
      add_arg("-sample_list", file_runtimepath(input$srpe_samplesheet)),
      add_arg("-reference_genome", file_runtimepath(input$srpe_reference_genome)),
      add_arg("-genome_gtf", file_runtimepath(input$srpe_genome_annotations)),
      add_arg("-repeat_gtf", file_runtimepath(input$srpe_repeat_annotations)),
      add_arg("-threads", input$srpe_threads),
      adapter_args,
      add_arg("-trimming_quality_threshold", input$srpe_trim_threshold),
      add_arg("-min_length_trim", input$srpe_min_length),
      add_arg("-max_length_trim", input$srpe_max_length),
      add_arg("-initial_trim_read1", input$srpe_initial_trim_r1),
      add_arg("-initial_trim_read2", input$srpe_initial_trim_r2),
      add_arg("-trimming", input$srpe_trim_mode),
      add_arg("-mismatch_align", input$srpe_mismatch_align),
      add_arg("-project_name", input$srpe_project_name),
      add_arg("-remove_duplicates", input$srpe_remove_duplicates),
      add_arg("-log2FC", input$srpe_dea_log2fc),
      add_arg("-FDR", input$srpe_dea_fdr),
      add_arg("-batch", input$srpe_dea_batch),
      add_arg("-method", input$srpe_dea_method),
      add_arg("-prediction_model", input$srpe_prediction_model),
      add_arg("-positive_class", input$srpe_positive_class)
    )
    if (isTRUE(input$srpe_polya)) {
      args <- c(args, "-polya", "polya")
    }
    build_run_command("perreo_srpe.sh", args)
  })
  srpe_command_preview <- reactive({
      adapter_args <- if (isTRUE(input$srpe_no_adapters)) {
      c(add_empty_arg("-adapt_r1"), add_empty_arg("-adapt_r2"))
    } else {
      c(
        add_arg("-adapt_r1", input$srpe_adapter_sequence_r1),
        add_arg("-adapt_r2", input$srpe_adapter_sequence_r2)
      )
    }
    args <- c(
      add_arg("-sample_list", file_displayname(input$srpe_samplesheet)),
      add_arg("-reference_genome", file_displayname(input$srpe_reference_genome)),
      add_arg("-genome_gtf", file_displayname(input$srpe_genome_annotations)),
      add_arg("-repeat_gtf", file_displayname(input$srpe_repeat_annotations)),
      add_arg("-threads", input$srpe_threads),
         adapter_args,
      add_arg("-trimming_quality_threshold", input$srpe_trim_threshold),
      add_arg("-min_length_trim", input$srpe_min_length),
      add_arg("-max_length_trim", input$srpe_max_length),
      add_arg("-initial_trim_read1", input$srpe_initial_trim_r1),
      add_arg("-initial_trim_read2", input$srpe_initial_trim_r2),
      add_arg("-trimming", input$srpe_trim_mode),
      add_arg("-mismatch_align", input$srpe_mismatch_align),
      add_arg("-project_name", input$srpe_project_name),
      add_arg("-remove_duplicates", input$srpe_remove_duplicates),
      add_arg("-log2FC", input$srpe_dea_log2fc),
      add_arg("-FDR", input$srpe_dea_fdr),
      add_arg("-batch", input$srpe_dea_batch),
      add_arg("-method", input$srpe_dea_method),
      add_arg("-prediction_model", input$srpe_prediction_model),
      add_arg("-positive_class", input$srpe_positive_class)
    )
    if (isTRUE(input$srpe_polya)) {
      args <- c(args, "-polya", "polya")
    }
    build_command_preview("perreo_srpe.sh", args)
  })
  lr_command <- reactive({
    args <- c(
      add_arg("-sample_list", file_runtimepath(input$lr_samplesheet)),
      add_arg("-reference_genome", file_runtimepath(input$lr_reference_genome)),
      add_arg("-repeat_gtf", file_runtimepath(input$lr_repeat_annotations)),
      add_arg("-threads", input$lr_threads),
      add_arg("-project_name", input$lr_project_name),
      add_arg("-log2FC", input$lr_dea_log2fc),
      add_arg("-FDR", input$lr_dea_fdr),
      add_arg("-batch_effect", input$lr_dea_batch),
      add_arg("-method", input$lr_dea_method),
      add_arg("-prediction_model", input$lr_prediction_model),
      add_arg("-positive_class", input$lr_positive_class)
    )
    build_run_command("perreo_lr.sh", args)
  })
  
  lr_command_preview <- reactive({
    args <- c(
      add_arg("-sample_list", file_displayname(input$lr_samplesheet)),
      add_arg("-reference_genome", file_displayname(input$lr_reference_genome)),
      add_arg("-repeat_gtf", file_displayname(input$lr_repeat_annotations)),
      add_arg("-threads", input$lr_threads),
      add_arg("-project_name", input$lr_project_name),
      add_arg("-log2FC", input$lr_dea_log2fc),
      add_arg("-FDR", input$lr_dea_fdr),
      add_arg("-batch_effect", input$lr_dea_batch),
      add_arg("-method", input$lr_dea_method),
      add_arg("-prediction_model", input$lr_prediction_model),
      add_arg("-positive_class", input$lr_positive_class)
    )
    build_command_preview("perreo_lr.sh", args)
  })
  run_pipeline <- function(command_spec, label, output_sink) {
    script_path <- command_spec$script
    if (!file.exists(script_path)) {
      showNotification(
        paste("Cannot find script:", script_path),
        type = "error",
        duration = NULL
      )
      return(invisible(NULL))
    }

    output_sink(paste0("[", Sys.time(), "] Starting ", label, " pipeline..."))
    showNotification(
      paste("Starting", label, "pipeline..."),
      type = "message"
    )
system_warning <- NULL
    command_output <- tryCatch(
      withCallingHandlers(
        system2(
          "bash",
          args = c(script_path, command_spec$args),
          wait = TRUE,
          stdout = TRUE,
          stderr = TRUE
        ),
        warning = function(w) {
          system_warning <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        paste("Execution error:", conditionMessage(e))
      }
    )

    exit_status <- attr(command_output, "status")
    
    if (!is.null(system_warning) && nzchar(system_warning)) {
      command_output <- c(command_output, paste("Warning:", system_warning))
    }
    if (!is.null(exit_status) && !identical(exit_status, 0L)) {
      command_output <- c(command_output, paste("Exit status:", exit_status))
    }
    if (length(command_output) == 0) {
      command_output <- "(No output was produced by the pipeline.)"
    }
    
    output_sink(paste(command_output, collapse = "\n"))
    showNotification(
      paste(label, "pipeline finished. Review the execution output panel."),
      type = "message"
    )
  }  
  observeEvent(input$run_srse_pipeline, {
    req(
      file_datapath(input$srse_samplesheet),
      file_datapath(input$srse_reference_genome),
      file_datapath(input$srse_genome_annotations),
      file_datapath(input$srse_repeat_annotations)
    )
    run_pipeline(srse_command(), "SR-SE", srse_run_output)  })
  
  observeEvent(input$run_srpe_pipeline, {
    req(
      file_datapath(input$srpe_samplesheet),
      file_datapath(input$srpe_reference_genome),
      file_datapath(input$srpe_genome_annotations),
      file_datapath(input$srpe_repeat_annotations)
    )
    run_pipeline(srpe_command(), "SR-PE", srpe_run_output)  })
  
  observeEvent(input$run_lr_pipeline, {
    req(
      file_datapath(input$lr_samplesheet),
      file_datapath(input$lr_reference_genome),
      file_datapath(input$lr_repeat_annotations)
    )
    run_pipeline(lr_command(), "LR", lr_run_output)  })
  observeEvent(input$run_fastqc, {
    req(input$fastq_files)
    
    if (isTRUE(input$overwrite_qc)) {
      unlink(list.files(fastqc_outdir, full.names = TRUE), recursive = TRUE)
    }
    
    fastq_df <- input$fastq_files
    
    withProgress(message = "Running FastQC...", value = 0, {
      n <- nrow(fastq_df)
      
      resultados <- lapply(seq_len(nrow(fastq_df)), function(i) {
        incProgress(1 / nrow(fastq_df), detail = fastq_df$name[i])
        
        f_in <- fastq_df$datapath[i]
        
        muestra_original <- fastq_df$name[i]
        
        base <- basename(f_in)
        
        base_noext <- base
        if (grepl("\\.gz$", base_noext)) {
          base_noext <- file_path_sans_ext(base_noext)
        }
        base_noext <- file_path_sans_ext(base_noext)
        
        cmd_args <- c(
          "-o", fastqc_outdir,
          "-f", "fastq",
          f_in
        )
        exit_code <- system2("fastqc", args = cmd_args)
        
        html_name <- paste0(base_noext, "_fastqc.html")
        html_path <- file.path(fastqc_outdir, html_name)
        html_exists <- file.exists(html_path)
        
        html_url <- if (html_exists) {
          paste0("/", file.path("fastqc_reports", html_name))
        } else {
          NA_character_
        }
        
        data.frame(
          sample    = muestra_original,
          html_file = html_name,
          html_path = html_path,
          html_url  = html_url,
          status    = ifelse(exit_code == 0 && html_exists,
                             "OK",
                             paste("ERROR", exit_code)),
          stringsAsFactors = FALSE
        )
      })
      
      
      resultados_df <- do.call(rbind, resultados)
      
      prev <- fastqc_results()
      if (!is.null(prev) && !isTRUE(input$overwrite_qc)) {
        resultados_df <- rbind(prev, resultados_df)
      }
      fastqc_results(resultados_df)
      
      multiqc_exit <- {
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)
        setwd(fastqc_outdir)
        system2("multiqc", args = c("--force", "."))
      }
      multiqc_html_name <- "multiqc_report.html"
      multiqc_html_path <- file.path(fastqc_outdir, multiqc_html_name)
      multiqc_html_url <- if (file.exists(multiqc_html_path)) {
        paste0("/", file.path("fastqc_reports", multiqc_html_name))
      } else {
        NA_character_
      }
      
      multiqc_report(list(
        status = ifelse(multiqc_exit == 0 && file.exists(multiqc_html_path), "OK", paste("ERROR", multiqc_exit)),
        html_path = multiqc_html_path,
        html_url = multiqc_html_url
      ))
    })
  })
  
  output$fastqc_table <- renderTable({
    res <- fastqc_results()
    if (is.null(res)) return(NULL)
    res[, c("sample", "status")]
  })
  
  output$fastqc_links <- renderUI({
    res <- fastqc_results()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    tagList(
      lapply(seq_len(nrow(res)), function(i) {
        tags$p(
          strong(res$sample[i]), "  -  ",
          if (res$status[i] == "OK" && !is.na(res$html_url[i])) {
            tags$a(
              "View FastQC report",
              href = res$html_url[i],
              target = "_blank"
            )
          } else {
            span("No report available (FastQC error)")
          }
        )
      })
    )
  })
  
  output$multiqc_report <- renderUI({
    res <- multiqc_report()
    if (is.null(res)) {
      return(tags$p("Run FastQC to generate the MultiQC report."))
    }
    
    if (!is.null(res$html_url) && !is.na(res$html_url) && res$status == "OK") {
      tags$iframe(
        src = res$html_url,
        width = "100%",
        height = "900px",
        style = "border: 1px solid #e0e0e0; border-radius: 8px;"
      )
    } else {
      tags$p(paste("Could not generate the MultiQC report (", res$status, ")."))
    }
  })
  
  observeEvent(input$go_analysis1, {
    updateTabsetPanel(session, "main_tabs", selected = "analysis1")
  })
  observeEvent(input$go_analysis2, {
    updateTabsetPanel(session, "main_tabs", selected = "analysis2")
  })
  observeEvent(input$go_analysis3, {
    updateTabsetPanel(session, "main_tabs", selected = "analysis3")
  })
  
  output$srse_command <- renderText({
    srse_command_preview()
  })
  
  output$srpe_command <- renderText({
    srpe_command_preview()
  })
  
  output$lr_command <- renderText({
    lr_command_preview()
  })

   output$srse_run_output <- renderText({
    srse_run_output()
  })
  
  output$srpe_run_output <- renderText({
    srpe_run_output()
  })
  
  output$lr_run_output <- renderText({
    lr_run_output()
  })
  observeEvent(input$back_home1, {
    updateTabsetPanel(session, "main_tabs", selected = "home")
  })
  observeEvent(input$back_home2, {
    updateTabsetPanel(session, "main_tabs", selected = "home")
  })
  observeEvent(input$back_home3, {
    updateTabsetPanel(session, "main_tabs", selected = "home")
  })
}

shinyApp(ui = ui, server = server)
