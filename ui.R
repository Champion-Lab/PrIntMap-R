library(shiny)
ui <- navbarPage(title = "PrIntMap-R",
                 tabPanel("Documentation",
                          includeMarkdown("www/Documentation.md"),
                          "Version 0.0.0"), #update version here for each push
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
                            uiOutput("combined_method_display"),
                            selectInput(inputId = "file_type",
                                         label = "Search Software",
                                         choices = c("PEAKS","MSFragger", "MaxQuant", "MetaMorpheus", "Proteome Discoverer", "Generic csv"),
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
                                     tabsetPanel(
                                       tabPanel("Traces",
                                                withSpinner(
                                                  plotlyOutput("plot_intensity2")
                                                )),
                                       tabPanel("Volcano Plot",
                                                flowLayout(
                                                  numericInput("p_cutoff", "p-value cutoff",
                                                               value = 0.05, min = 0, max = 1, step = 0.05),
                                                  numericInput("l2fc_cutoff", "log2 fold-change cutoff",
                                                               value = 1, min = 0, max = 100),
                                                  numericInput("min_valid_sample", "Miniumum number of observations per sample",
                                                               value = 2, min = 1),
                                                  checkboxInput("equal_var", "Equal Variance",
                                                                value = T),
                                                  checkboxInput("remove_na", "Remove NA values",
                                                                value = T),
                                                  numericInput("set_na_value", "Set NA values",
                                                               value = 0),
                                                  checkboxInput("BH_correction", "Apply Benjamini-Hochburg correction",
                                                                value = F),
                                                  checkboxInput("display_comp_vals", "Display compromised values",
                                                                value = T),
                                                  checkboxInput("display_infinite_vals", "Display Infinite Fold Changes",
                                                                value = T)
                                                ),
                                                withSpinner(
                                                  plotlyOutput("volcano", height = 600)
                                                ),
                                                textOutput("volcano_text"))),
                                     ),
                            tabPanel("Multiple Samples",
                                     fluidPage(
                                       fluidRow(column(3, numericInput(inputId = "number_sample", label = "Choose number of samples", 
                                                             value = 3,min =2, step = 1)),
                                                column(2, radioButtons(inputId = "mult_sample_comparison",
                                                             label = "Type of Comparison",
                                                             choices = c("Overlay", "Difference", "Fold Change"),
                                                             selected = "Overlay")),
                                                column(1, actionButton(inputId = "mult_go", label = "Run/Update"))
                                                ),
                                       uiOutput("sample_numbers"),
                                     ),
                                     
                                     fluidRow(withSpinner(
                                       plotlyOutput("plot_intensity_mult")
                                     ))
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
                                                                 "Multiple Samples", "Stacked Peptides"),
                                                    selected = "One Sample"),
                                       selectInput(inputId = "annot_color",
                                                   label = "Color",
                                                   choices = c("red",
                                                               "yellow",
                                                               "green",
                                                               "pink",
                                                               "purple"),
                                                   selected = "red")
                                     ),
                                     withSpinner(
                                       plotlyOutput("annotation_plotly")
                                     ),
                                    ),
                            tabPanel("PTMs",
                                     flowLayout(
                                       radioButtons(inputId = "PTM",
                                                          label = "PTMs",
                                                          choices = c("\\(\\+57.02\\)",
                                                                           "\\(\\+0.98\\)",
                                                                           "\\(\\+15.99\\)",
                                                                            "\\(\\-18.01\\)", 
                                                                           "\\(\\-17.03\\)",
                                                                          "\\(\\+21.98\\)")),
                                       textInput(inputId = "custom_PTM",
                                                 label = "Custom PTM annotation in RegEx"),
                                       radioButtons(inputId = "disp_overlay_PTM",
                                                    label = "Type of Plot",
                                                    choices = c("One Sample", "Two Samples", 
                                                                "Multiple Samples"),
                                                    selected = "One Sample")),
                                       withSpinner(
                                         plotlyOutput("PTM_plotly")
                                       
                                     )),
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