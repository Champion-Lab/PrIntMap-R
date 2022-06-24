library(shiny)
ui <- navbarPage(title = "PrIntMap-R",
                 tabPanel("Documentation",
                          includeMarkdown("www/Documentation.md")),
                 tabPanel("Run",
                          flowLayout(
                            fileInput(inputId = "database_file", label = "Upload fasta database file",
                                      accept = c(".fasta", ".fa")),
                            fileInput(inputId = "peptide_file1", label = "Upload csv/tsv peptide output",
                                      accept = c(".csv", ".tsv", ".txt", ".psmtsv")),
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
                                         choices = c("PEAKS","MSFragger", "MaxQuant", "MetaMorpheus", "Proteome Discover"),
                                         selected = "PEAKS"),
                            textInput(inputId = "AccessionID", 
                                      label = "Input Protein Accession ID",
                                      value = ""),
                            uiOutput("intensity"),
                            radioButtons(inputId = "y_axis_scale",
                                         label = "y-Axis Scale",
                                         choices = c("Linear", "log"),
                                         selected = "Linear")
                          ),
                          flowLayout(
                            textOutput("AccessionID"),
                            textOutput("peptides1_sample_count"),
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
                                                 accept = c(".csv", ".tsv", ".txt", ".psmtsv")),
                                       checkboxInput( inputId = "duplicate_file2",
                                                                label = "Use same peptide output file?", 
                                                                value = F),
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
                                                    selected = "Overlay"),
                                       textOutput("peptides2_sample_count")),
                                     withSpinner(
                                       plotlyOutput("plot_intensity2")
                                     )),
                            tabPanel("Multiple Samples",
                                     fluidPage(
                                       fluidRow(column(3, numericInput(inputId = "number_sample", label = "Choose number of samples", 
                                                             value = 3,min =2, step = 1)),
                                                column(2, radioButtons(inputId = "mult_sample_comparison",
                                                             label = "Type of Comparison",
                                                             choices = c("Overlay", "Difference", "Fold Change"),
                                                             selected = "Overlay")),
                                                column(1, actionButton(inputId = "mult_go", label = "Go"))
                                                ),
                                       uiOutput("sample_numbers"),
                                     ),
                                     
                                     fluidRow(withSpinner(
                                       plotlyOutput("plot_intensity_mult")
                                     )),
                                     
                            ),
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
                                       radioButtons(inputId = "disp_overlay_annot",
                                                     label = "Type of Plot",
                                                     choices = c("One Sample", "Two Samples", 
                                                                 "Multiple Samples"),
                                                    selected = "One Sample"),
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
                            tabPanel("Unique Peptides",
                                     textOutput("unique_text"),
                                     actionButton(inputId = "run_unique",
                                                  label = "Run/Update"),
                                     withSpinner(
                                       plotlyOutput("unique_peps")
                                     )),
                            tabPanel("Stacked Peptide Plot",
                                     radioButtons(inputId = "stacked_peptides_yunits",
                                                  label = "Units for y-axis",
                                                  choices = c("AA Position", "Intensity"),
                                                  selected = "AA Position"),
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
                                                    choices = c("One sample", "Both samples", 
                                                                "Multiple Samples"),
                                                    selected = "One sample")
                                       )
                                      
                                     )
                          )
                        ),
                 tabPanel("Examples", 
                          "Coming Soon"                        
                         )
  
)