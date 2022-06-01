#load necessary functions
source("global.R")


#increase allowed upload size to 100 mb
options(shiny.maxRequestSize=100*1024^2)

server <- function(input, output) {
  
  database <- reactive({
    validate(
      need(!is.null(input$database_file), "No Database File Provided")
    )
    import_db(input$database_file$datapath)
  })
  
  
  protein_obj1 <- reactive({
    validate(
      need(input$AccessionID != "", "No Accession Number Provided")
    )
    
    protein_attributes(select_prot(db = database(),
                             Accession = toupper(input$AccessionID)))
  })
  
  display_origin_peps <- reactive(input$disp_origin)
  
  output$AccessionID <- renderText({
    paste0("Protein: ", protein_obj1()[2])
  })
  
  
  AA_df1 <- reactive({
    create_AA_df(protein_obj1()[1])
  })
  
  peptides1 <- reactive({
    validate(
      need(!is.null(input$peptide_file1), "no peptide file provided")
    )
    if (input$file_type == "PEAKS" && input$combinedbool1 == "Individual Sample") {
      read_peptide_csv_PEAKS_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "PEAKS" && input$combinedbool1 == "Combined"){
      read_peptide_csv_PEAKS_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1)
    } else if (input$file_type == "MSFragger" && input$combinedbool1 == "Individual Sample"){
      read_peptide_tsv_MSFragger_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "MSFragger" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_MSFragger_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1)
    } else if(input$file_type == "MaxQuant" && input$combinedbool1 == "Individual Sample") {
      read_peptide_tsv_MaxQuant_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1)
    } else if (input$file_type == "MaxQuant" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_MaxQuant_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1)
    } else if (input$file_type == "Proteome Discover" && input$combinedbool1 == "Individual Sample") {
      #coming soon
    } else if (input$file_type == "Proteome Discover" && input$combinedbool1 == "Combined"){
      #coming soon
    }
  })
  
  
  intensity_vec1 <- reactive({
    create_intensity_vec(peptides1(), protein_obj1()[1], intensity = input$intensity_metric)
  })
  
  AA_df_intensity1 <- reactive({
    combine_AA_intensity(AA_df1(), intensity_vec1())
  })
  

  origin_vec1 <- reactive({
    create_peptide_origin_vec(peptides1(), protein_obj1()[1], intensity = input$intensity_metric)
  })
  
  AA_df_origin1 <- reactive({
    add_origin_peptide_vector(AA_df_intensity1(), origin_vec1())
  })
  
  

  output$plot_intensity1 <- renderPlotly({
    validate(
      need(nrow(AA_df_origin1()) > 0, "Invalid peptide file. Check format.")
    )
    if (input$disp_origin) {
      create_plotly({
        plot_origin(AA_df_origin1(), protein_obj1()[2],
                    intensity_label = input$intensity_metric)
      })
    } else {
      create_plotly({
        plot_intensity(AA_df_intensity1(), protein_obj1()[2],
                         intensity_label = input$intensity_metric)
      })
    } 
  })

  AA_df2 <- reactive({
    create_AA_df(protein_obj1()[1])
  })

  peptides2 <- reactive({
    validate(
      need(!is.null(input$peptide_file2), "no peptide file provided")
    )
    if (input$file_type == "PEAKS" && input$combinedbool2 == "Individual Sample") {
      read_peptide_csv_PEAKS_bysamp(input$peptide_file2$datapath)
    } else if (input$file_type == "PEAKS" && input$combinedbool2 == "Combined"){
      read_peptide_csv_PEAKS_comb(input$peptide_file2$datapath, sample_pattern = input$sample_regex2)
    } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Individual Sample"){
      read_peptide_tsv_MSFragger_bysamp(input$peptide_file2$datapath)
    } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Combined") {
      read_peptide_tsv_MSFragger_comb(input$peptide_file2$datapath, sample_pattern = input$sample_regex2)
    } else if(input$file_type == "MaxQuant" && input$combinedbool2 == "Individual Sample") {
      read_peptide_tsv_MaxQuant_comb(input$peptide_file2$datapath, sample_pattern = input$sample_regex2)
    } else if (input$file_type == "MaxQuant" && input$combinedbool2 == "Combined") {
      read_peptide_tsv_MaxQuant_comb(input$peptide_file2$datapath, sample_pattern = input$sample_regex2)
    } else if (input$file_type == "Proteome Discover" && input$combinedbool2 == "Individual Sample") {
      #coming soon
    } else if (input$file_type == "Proteome Discover" && input$combinedbool2 == "Combined"){
      #coming soon
    }
  })

  intensity_vec2 <- reactive({
    create_intensity_vec(peptides2(), protein_obj1()[1], intensity = input$intensity_metric)
  })

  AA_df_intensity2 <- reactive({
    combine_AA_intensity(AA_df2(), intensity_vec2())
  })


  origin_vec2 <- reactive({
    create_peptide_origin_vec(peptides2(), protein_obj1()[1], intensity = input$intensity_metric)
  })

  AA_df_origin2 <- reactive({
    add_origin_peptide_vector(AA_df_intensity2(), origin_vec2())
  })


  sample1_df <- reactive(data.frame(sample = rep(input$sample_name1, nrow(AA_df_intensity1()))))
  sample2_df <- reactive(data.frame(sample = rep(input$sample_name2, nrow(AA_df_intensity1()))))

  AA_df_origin1b <- reactive(bind_cols(AA_df_origin1(), sample1_df()))
  AA_df_origin2b <- reactive(bind_cols(AA_df_origin2(), sample2_df()))
  AA_df_intensity1b <- reactive(bind_cols(AA_df_intensity1(), sample1_df()))
  AA_df_intensity2b <- reactive(bind_cols(AA_df_intensity2(), sample2_df()))

  AA_df_origin_comb <- reactive(bind_rows(AA_df_origin1b(), AA_df_origin2b()))
  AA_df_intensity_comb <- reactive(bind_rows(AA_df_origin1b(), AA_df_origin2b()))
  
  AA_df_origin_comb_wide <- reactive({
    cbind(AA_df_origin_comb()[AA_df_origin_comb()$sample == input$sample_name1,],
          AA_df_origin_comb()[AA_df_origin_comb()$sample == input$sample_name2,])
    })
  AA_df_intensity_comb_wide <- reactive({
    cbind(AA_df_intensity_comb()[AA_df_intensity_comb()$sample == input$sample_name1,],
          AA_df_intensity_comb()[AA_df_intensity_comb()$sample == input$sample_name2,])
  })
  
  wide_data_comb <- reactive({
    if (input$disp_origin) {
      data.frame(AA_df_origin_comb_wide())
    } else {
      data.frame(AA_df_intensity_comb_wide())
    }
  })
  
 difference_vec <- reactive({
   wide_data_comb()$intensity.1 - wide_data_comb()$intensity
 })
 
 fold_change_vec_nan <- reactive({
   wide_data_comb()$intensity.1 / wide_data_comb()$intensity
 })
 
 fold_change_label_vec <- reactive({
   return_vec <- rep("", length(fold_change_vec_nan()))
   for (i in 1:length(return_vec)) {
     if(is.infinite(wide_data_comb()$intensity.1[i] / wide_data_comb()$intensity[i])) {
       return_vec[i] <- "Infinite"
     } else if (is.na(wide_data_comb()$intensity.1[i] / wide_data_comb()$intensity[i])) {
       return_vec[i] <- "Zero"
     } else {
       return_vec[i] <- "Normal"
     }
   }
   return(return_vec)
 })
 
 fold_change_vec <- reactive({
   vec <- fold_change_vec_nan()
   vec[is.infinite(vec)] <- max(vec[is.finite(vec)], na.rm = T)*1.25
   vec[is.na(vec)] <- 0
   return(vec)
 })
 
 
 wide_data_comb_plot <- reactive(cbind(wide_data_comb(),
                                       data.frame(difference = difference_vec()),
                                       data.frame(fold_change = fold_change_vec()),
                                       data.frame(fold_change_label = fold_change_label_vec())))
 
  
  ggplot_intensity2 <- reactive({
    if (input$two_sample_comparison == "Difference") {
      plot_difference_comb(wide_data_comb_plot(), protein_obj1()[2])
    } else if (input$two_sample_comparison == "Fold Change"){
      plot_foldchange_comb(wide_data_comb_plot(), protein_obj1()[2])
    } else if (input$two_sample_comparison == "Overlay"){
      if (input$disp_origin) {
        plot_origin_comb(AA_df_origin_comb(), protein_obj1()[2],
                         intensity_label = input$intensity_metric)
      } else {
        plot_intensity_comb(AA_df_intensity_comb(), protein_obj1()[2],
                            intensity_label = input$intensity_metric)
      }
    }
  })
  
  output$plot_intensity2 <- renderPlotly(create_plotly(ggplot_intensity2()))
  
  
  output$sample1_label <- renderText(paste0("Sample 1 (orange): ", input$sample_name1))
  output$sample2_label <- renderText(paste0("Sample 2 (blue): ", input$sample_name2))
  
  
  annotation_regex <- reactive({
    if(input$annotation != "CUSTOM") {
      annotations[[input$annotation]]
    } else {
      input$custom_annotation
    }
  })
  
  annotation_name <- reactive({
    if(input$annotation != "CUSTOM") {
      input$annotation
    } else {
      input$custom_annotation
    }
  })
  
  output$annotation_regex <- renderText(annotation_regex())
  
  annotation_vec <- reactive({
    create_annotation_vec(pattern = annotation_regex(),
                          protein = protein_obj1()[1])
  })
  
  annotation_df <- reactive({
    create_annotation_df(AA_df = AA_df_origin1(),
                         annotation_vec = annotation_vec(),
                         annotation_name = annotation_name())
  })
  
  annotation_plot <- reactive({
    if (input$disp_origin) {
        plot_origin(AA_df_origin1(), protein_obj1()[2],
                    intensity_label = input$intensity_metric)
    } else {
        plot_intensity(AA_df_intensity1(), protein_obj1()[2],
                       intensity_label = input$intensity_metric)
    } 
  })
  
  annotation_plot2 <- reactive({
    if (input$disp_overlay_annot) {
      add_annotation_layer(plot = ggplot_intensity2(),
                           annotation_df = annotation_df(),
                           color = input$annot_color)
    } else {
      add_annotation_layer(plot = annotation_plot(),
                           annotation_df = annotation_df(),
                           color = input$annot_color)
    }
  })
  
  output$annotation_plotly <- renderPlotly(annotation_plot2())
  
  
  stacked_plot_dataframe <- reactive({
    create_stack_line_df(peptide_df = peptides1(),
                         protein = protein_obj1()[1],
                         intensity = input$intensity_metric)
  })
  
  stacked_plot <- reactive({
    create_stacked_line_plot(stacked_plot_dataframe(),
                             protein_name = protein_obj1()[2],
                             protein_seq = protein_obj1()[1])
  })
  
  output$stacked_plotly <- renderPlotly({
    create_stack_line_plotly(stacked_plot())})
  
  percent_coverage_list <- reactive({
    calculate_percent_cov(AA_df_intensity1())
  })
  
  output$total_AA <- renderText(paste0("Total Amino Acids: ",
                                       percent_coverage_list()[[1]]))
  output$covered_AA <- renderText(paste0("Amino Acids Observed: ",
                                         percent_coverage_list()[[2]]))
  output$percent_coverage <- renderText(paste0("Percent Coverage: ",
                                               percent_coverage_list()[[3]]))
  output$coverage_message <- renderText("Use PSM for traditional % coverage calculation. NA values not included for Area/Intensity")
  
  download_peptide_csv <- reactive({
    if (input$download_sample_number == "One sample") {
      covert_to_download_df(AA_df_origin1(), input$intensity_metric)
    } else {
      covert_to_download_df(AA_df_origin_comb(), input$intensity_metric)
    }
    
  })
  
  download_peptide_filename <- reactive({
    if (input$download_sample_number == "One sample") {
      paste0(input$sample_name1, "_", input$intensity_metric, "_", str_replace_all(protein_obj1()[2], "\\|", "_"), "_intensity_by_aminoacid.csv")
    } else {
      paste0(input$intensity_metric, "_", str_replace_all(protein_obj1()[2], "\\|", "_"), "_intensity_by_aminoacid_two_sample.csv")
    }
    
  })
  
  
  output$download_intensity_df <- downloadHandler(
    filename = function() {
      download_peptide_filename()
    },
    content = function(file) {
      write.csv(download_peptide_csv(), file, row.names = FALSE)
    }
  )
  
}
  
 