library(shiny)
ui <- navbarPage(title = "Protein Intensity Mapper",
                 tabPanel("Run",
                          flowLayout(
                            fileInput(inputId = "database_file", label = "Upload fasta database file",
                                      accept = c(".fasta", ".fa")),
                            fileInput(inputId = "peptide_file1", label = "Upload csv/tsv peptide output",
                                      accept = c(".csv", ".tsv", ".txt")),
                            textInput(inputId = "sample_name1", 
                                      label = "Input sample display name",
                                      value = "Sample 1"),
                            radioButtons(inputId = "combinedbool1",
                                         label = "Type of input file",
                                         choices = c("Individual Sample", "Combined"),
                                         selected = "Individual Sample"),
                            textInput(inputId = "sample_regex1",
                                      label = "For combined files, input sample name (RegEx)"),
                            selectInput(inputId = "file_type",
                                         label = "Search Software",
                                         choices = c("PEAKS","MSFragger", "MaxQuant", "Proteome Discover (Coming soon)"),
                                         selected = "PEAKS"),
                            textInput(inputId = "AccessionID", 
                                      label = "Input Protein Accession ID",
                                      value = ""),
                            radioButtons(inputId = "intensity_metric",
                                         label = "Intensity Metric",
                                         choices = c("PSM", "Area", "Intensity"),
                                         selected = "PSM")
                          ),
                          flowLayout(
                            textOutput("AccessionID"),
                            checkboxInput(inputId = "disp_origin",
                                          label = "Display Origin Peptides",
                                          value = F)
                          ),
                          
                          tabsetPanel(
                            tabPanel("Basic",
                                     withSpinner(
                                       plotlyOutput("plot_intensity1")
                                     )),
                            tabPanel("Two Samples",
                                     flowLayout(
                                       fileInput(inputId = "peptide_file2", label = "Upload 2nd .csv peptide output",
                                                 accept = c(".csv", ".tsv", ".txt")),
                                       radioButtons(inputId = "combinedbool2",
                                                    label = "Type of input file",
                                                    choices = c("Individual Sample", "Combined"),
                                                    selected = "Individual Sample"),
                                       textInput(inputId = "sample_regex2",
                                                 label = "For combined files, input sample name (RegEx)"),
                                       textInput(inputId = "sample_name2", 
                                                 label = "Input sample name",
                                                 value = "Sample 2"),
                                       radioButtons(inputId = "two_sample_comparison",
                                                    label = "Type of Comparison",
                                                    choices = c("Overlay", "Difference", "Fold Change"),
                                                    selected = "Overlay")),
                                     withSpinner(
                                       plotlyOutput("plot_intensity2")
                                     ),
                                    textOutput("sample1_label"),
                                    textOutput("sample2_label"), tableOutput("test_table")),
                            tabPanel("Annotation",
                                     flowLayout(
                                       selectInput(inputId = "annotation",
                                                   label = "Annotation",
                                                   choices = c("Trypsin", "LysC",
                                                               "N_glycosylation",
                                                               "CUSTOM"),
                                                   selected = "Trypsin"),
                                       textInput(inputId = "custom_annotation",
                                                 label = "Custom sequence annotation in RegEx"),
                                       textOutput("annotation_regex"),
                                       checkboxInput(inputId = "disp_overlay_annot",
                                                     label = "Display Both Samples",
                                                     value = F),
                                       selectInput(inputId = "annot_color",
                                                   label = "Color",
                                                   choices = c("red",
                                                               "yellow",
                                                               "green",
                                                               "pink",
                                                               "purple"),
                                                   selected = "green")
                                     ),
                                     withSpinner(
                                       plotlyOutput("annotation_plotly")
                                     ),
                                    ),
                            tabPanel("Unique Peptides", "Coming Soon"),
                            tabPanel("Stacked Peptide Plot",
                                      withSpinner(
                                        plotlyOutput("stacked_plotly") 
                                      )
                                     ),
                            tabPanel("Percent Coverage",
                                     textOutput("total_AA"),
                                     textOutput("covered_AA"),
                                     textOutput("percent_coverage"),
                                     textOutput("coverage_message")),
                            tabPanel("Export",
                                     flowLayout(
                                       
                                       downloadButton("download_intensity_df",
                                                      "Download Intensity CSV"),
                                       radioButtons(inputId = "download_sample_number",
                                                    label = "Download which samples?",
                                                    choices = c("One sample", "Both samples"),
                                                    selected = "One sample")
                                       )
                                      
                                     )
                          )
                        ),
                 tabPanel("Documentation", 
                          "By Simon D. Weaver",
                          tags$hr(),
                          "Champion Lab - University of Notre Dame")
  
)