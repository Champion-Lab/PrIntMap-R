#load necessary functions
source("global.R")


#increase allowed upload size to 100 mb
options(shiny.maxRequestSize=100*1024^2)

server <- function(input, output, session) {
  
  database <- reactive({
    validate(
      need(!is.null(input$database_file), "No Database File Provided")
    )
    import_db(input$database_file$datapath)
  })
  
  observeEvent(input$file_type, {
    if(input$file_type == "MaxQuant"){
      updateRadioButtons(session,"combinedbool1",
                         selected = "Combined")
      updateRadioButtons(session,"combinedbool2",
                         selected = "Combined")
    }
  })
  
  protein_obj1 <- reactive({
    validate(
      need(input$AccessionID != "", "No Accession Number Provided")
    )
    
    protein_attributes(select_prot(db = database(),
                             Accession = input$AccessionID))
  })
  
  display_origin_peps <- reactive(input$disp_origin)
  
  output$AccessionID <- renderText({
    paste0("Protein: ", protein_obj1()[2])
  })
  
  
  AA_df1 <- reactive({
    create_AA_df(protein_obj1()[1])
  })
  
  peptides1_list <- reactive({
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
    } else if (input$file_type == "MaxQuant" ) {
      read_peptide_tsv_MaxQuant_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1)
    } else if (input$file_type == "Proteome Discover" && input$combinedbool1 == "Individual Sample") {
      #coming soon
    } else if (input$file_type == "Proteome Discover" && input$combinedbool1 == "Combined"){
      #coming soon
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool1 == "Individual Sample"){
      read_peptide_tsv_Metamorpheus_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_Metamorpheus_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1)
    } 
  })
  
  peptides1 <- reactive(peptides1_list()[[1]])
  
  output$peptides1_sample_count <- renderText({
    if (input$combinedbool1 == "Combined") {
      paste0("Samples/Replicates Combined: ", peptides1_list()[[2]])
    }
  })

  output$intensity <- renderUI({
    mychoicesint <- intensity_metric_choices(peptides1())
    radioButtons(inputId = "intensity_metric",
                 label = "Intensity Metric",
                 choices = mychoicesint)
  })

  intensity_vec1 <- reactive({
    create_intensity_vec(peptides1(), protein_obj1()[1], intensity = input$intensity_metric)
  })
  
  AA_df_intensity1 <- reactive({
    combine_AA_intensity(AA_df1(), intensity_vec1())
  })
  

  origin_vec1 <- reactive({
    validate(need(input$intensity_metric, message = FALSE))
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
        return_plot <- plot_origin(AA_df_origin1(), protein_obj1()[2],
                    intensity_label = input$intensity_metric)
        if(input$y_axis_scale == "log"){
          return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                          breaks = base_breaks())
        }
        return(return_plot)
      })
    } else {
      create_plotly({
        return_plot <- plot_intensity(AA_df_intensity1(), protein_obj1()[2],
                         intensity_label = input$intensity_metric)
        if(input$y_axis_scale == "log"){
          return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                          breaks = base_breaks())
        }
        return(return_plot)
      })
    } 
  })

  AA_df2 <- reactive({
    create_AA_df(protein_obj1()[1])
  })

  peptide_file2<- reactive({
    if(input$duplicate_file2){
      input$peptide_file1$datapath
    }else{
      validate(
        need(!is.null(input$peptide_file2), "no peptide file provided")
      )
      input$peptide_file2$datapath}
  })
  
  
  
  peptides2_list <- reactive({

      if (input$file_type == "PEAKS" && input$combinedbool2 == "Individual Sample") {
        read_peptide_csv_PEAKS_bysamp(peptide_file2())
      } else if (input$file_type == "PEAKS" && input$combinedbool2 == "Combined"){
        read_peptide_csv_PEAKS_comb(peptide_file2(), sample_pattern = input$sample_regex2)
      } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Individual Sample"){
        read_peptide_tsv_MSFragger_bysamp(peptide_file2())
      } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Combined") {
        read_peptide_tsv_MSFragger_comb(peptide_file2(), sample_pattern = input$sample_regex2)
      } else if (input$file_type == "MaxQuant" && input$combinedbool2 == "Combined") {
        read_peptide_tsv_MaxQuant_comb(peptide_file2(), sample_pattern = input$sample_regex2)
      } else if (input$file_type == "Proteome Discover" && input$combinedbool2 == "Individual Sample") {
        #coming soon
      } else if (input$file_type == "Proteome Discover" && input$combinedbool2 == "Combined"){
        #coming soon
      } else if (input$file_type == "MetaMorpheus" && input$combinedbool2 == "Individual Sample"){
        read_peptide_tsv_Metamorpheus_bysamp(peptide_file2())
      } else if (input$file_type == "MetaMorpheus" && input$combinedbool2 == "Combined") {
        read_peptide_tsv_Metamorpheus_comb(peptide_file2(), sample_pattern = input$sample_regex2)
      } 
   
  })
  
  peptides2 <- reactive(peptides2_list()[[1]])
  
  output$peptides2_sample_count <- renderText({
    if (input$combinedbool2 == "Combined") {
      paste0("Samples/Replicates Combined: ", peptides2_list()[[2]])
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
  AA_df_intensity_comb <- reactive(bind_rows(AA_df_intensity1b(), AA_df_intensity2b()))
  
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
   vec <- wide_data_comb()$intensity.1 / wide_data_comb()$intensity
   return(vec)
 })
 
 fold_change_label_vec <- reactive({
   return_vec <- rep("", length(fold_change_vec_nan()))
   for (i in 1:length(return_vec)) {
     if(is.infinite(fold_change_vec_nan()[i])) {
       return_vec[i] <- "Infinite"
     } else if (is.na(fold_change_vec_nan()[i])) {
       return_vec[i] <- "Normal"
     } else if (fold_change_vec_nan()[i] == 0){
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
      plot_difference_comb(wide_data_comb_plot(), protein_obj1()[2],
                           intensity_label = input$intensity_metric)
    } else if (input$two_sample_comparison == "Fold Change"){
      plot_foldchange_comb(wide_data_comb_plot(), protein_obj1()[2],
                           intensity_label = input$intensity_metric)
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
  
  output$plot_intensity2 <- renderPlotly({
    if(input$y_axis_scale == "log"){
      return_plot <- ggplot_intensity2() + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                      breaks = base_breaks())
      create_plotly(return_plot)
    } else {
      create_plotly(ggplot_intensity2())
    }
      
  })
  
  number <- reactive({
    if(is.integer(input$number_sample)){ 
      as.integer(input$number_sample)}
    else{
      stop("Please choose a number of samples.")
    }
    
  })  
  output$sample_numbers <- renderUI({
    if(number() >= 2){
      lapply(2:number(), function(i){
        fluidRow(list(
        column(3,fileInput(inputId = paste0("peptide_file_mult", i), 
                  label = paste0("Upload .csv peptide output number ",i),
                  accept = c(".csv", ".tsv", ".txt", ".psmtsv"))),
        column(2, checkboxInput( inputId = paste0("duplicate_file_mult", i),
                                 label = "Use same peptide output file?", 
                                 value = F)),
        if(input$file_type =="MaxQuant"){
          column(2,radioButtons(inputId = paste0("combinedbool_mult",i),
                                label = "Type of input file",
                                choices = c("Combined"),
                                selected = "Combined"))}
          else {
            column(2,radioButtons(inputId = paste0("combinedbool_mult",i),
                                      label = "Type of input file",
                                      choices = c("Individual Sample", "Combined"),
                                      selected = "Individual Sample"))},
        column(2, textInput(inputId = paste0("sample_regex_mult",i),
                  label = "For combined files, input sample name (RegEx)")),
        column(2, textInput(inputId = paste0("sample_name_mult",i),
                  label = "Input sample display name",
                  value = paste("Sample", i))),
        uiOutput(paste0("mult_sample_count", i))
        )
        )}
      )}
    else{
      stop("Please choose a valid number of samples.")
    }
  })
  
    
  AA_df_list <- reactive({
    lapply(2:number(), function(i){
      create_AA_df(protein_obj1()[1])
    })
  })
  
  peptide_file_mult<- reactive({
    lapply(2:number(), function(i){
      if(input[[paste0("duplicate_file_mult", i)]]){
        input$peptide_file1$datapath
      }else {
      validate(
        need(!is.null(input[[paste0("peptide_file_mult", i)]]), "no peptide file provided")
      )
        input[[paste0("peptide_file_mult",i)]][["datapath"]]
        }
    })
  })
  
  peptidesmult_list <- reactive({
    lapply(2:number(), function(i){
      
      if (input$file_type == "PEAKS" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample") {
        read_peptide_csv_PEAKS_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "PEAKS" && input[[paste0("combinedbool_mult",i)]] == "Combined"){
        read_peptide_csv_PEAKS_comb(peptide_file_mult()[[i-1]], 
                                    sample_pattern = input[[paste0("sample_regex_mult", i)]])
      } else if (input$file_type == "MSFragger" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample"){
        read_peptide_tsv_MSFragger_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "MSFragger" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_tsv_MSFragger_comb(peptide_file_mult()[[i-1]],
                                        sample_pattern = input[[paste0("sample_regex_mult", i)]])
      } else if (input$file_type == "MaxQuant" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_tsv_MaxQuant_comb(peptide_file_mult()[[i-1]], 
                                       sample_pattern = input[[paste0("sample_regex_mult", i)]])
      } else if (input$file_type == "Proteome Discover" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample") {
        #coming soon
      } else if (input$file_type == "Proteome Discover" && input[[paste0("combinedbool_mult",i)]] == "Combined"){
        #coming soon
      } else if (input$file_type == "MetaMorpheus" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample"){
        read_peptide_tsv_Metamorpheus_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "MetaMorpheus" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_tsv_Metamorpheus_comb(peptide_file_mult()[[i-1]], sample_pattern = input[[paste0("sample_regex_mult", i)]])
      } 
    }
    )}
  )
  
  peptidesmult <- reactive({ map(peptidesmult_list(), 1)
    })
  peptidesmult_count <- eventReactive(
    input$mult_go, { map(peptidesmult_list(), 2)
    })

  intensity_vec_list <- reactive({
    lapply(2:number(), function(i){
      create_intensity_vec(peptidesmult()[[i-1]], protein_obj1()[1], intensity = input$intensity_metric)
    })
  })
  
  AA_df_intensity_list <- reactive({
    lapply(2:number(), function(i){
      combine_AA_intensity(AA_df_list()[[i-1]], intensity_vec_list()[[i-1]])
    })
  })
  
  origin_vec_list <- reactive({
    lapply(2:number(), function(i){
      create_peptide_origin_vec(peptidesmult()[[i-1]], protein_obj1()[1], intensity = input$intensity_metric)
    })
  })
  
  AA_df_origin_list <- reactive({
    lapply(2:number(), function(i){
      add_origin_peptide_vector(AA_df_intensity_list()[[i-1]], origin_vec_list()[[i-1]])
    })
  })
  
  sample_mult_df <- reactive({
    lapply(2:number(), function(i){
      data.frame(sample = rep(input[[paste0("sample_name_mult",i)]], nrow(AA_df_intensity1())))
    })
  })
  
  AA_df_origin_multb <- reactive({
    lapply(2:number(), function(i){
      bind_cols(AA_df_origin_list()[[i-1]], sample_mult_df()[[i-1]])
    })
  })
  
 

  AA_df_intensity_multb <- reactive({
    lapply(2:number(), function(i){
      bind_cols(AA_df_intensity_list()[[i-1]], sample_mult_df()[[i-1]])
    })
  })
  
  
  AA_df_origin_comb_mult <- reactive(bind_rows(AA_df_origin1b(), AA_df_origin_multb()))
  AA_df_intensity_comb_mult <- reactive(bind_rows(AA_df_intensity1b(), AA_df_intensity_multb()))
  
  AA_df_origin_comb_mult_compare <- reactive(cbind(AA_df_origin_comb_mult(),
                                                   data.frame(intensity_original = rep(AA_df_origin1b()$intensity,
                                                                                       number()))))
  
  
  difference_vec_mult <- reactive({
    AA_df_origin_comb_mult_compare()$intensity - AA_df_origin_comb_mult_compare()$intensity_original
  })
  
  fold_change_vec_nan_mult <- reactive({
    vec <- AA_df_origin_comb_mult_compare()$intensity / AA_df_origin_comb_mult_compare()$intensity_original
    return(vec)
  })
  
  fold_change_label_vec_mult <- reactive({
    return_vec <- rep("", length(fold_change_vec_nan_mult()))
    for (i in 1:length(return_vec)) {
      if(is.infinite(fold_change_vec_nan_mult()[i])) {
        return_vec[i] <- "Infinite"
      } else if (is.na(fold_change_vec_nan_mult()[i])) {
        return_vec[i] <- "Normal"
      } else if (fold_change_vec_nan_mult()[i] == 0){
        return_vec[i] <- "Zero"
      } else {
        return_vec[i] <- "Normal"
      }
    }
    return(return_vec)
  })
  
  fold_change_vec_mult <- reactive({
    vec <- fold_change_vec_nan_mult()
    vec[is.infinite(vec)] <- max(vec[is.finite(vec)], na.rm = T)*1.25
    vec[is.na(vec)] <- 0
    return(vec)
  })
  
  AA_df_origin_comb_mult_compare_plot <- reactive(cbind(AA_df_origin_comb_mult_compare(),
                                        data.frame(difference = difference_vec_mult()),
                                        data.frame(fold_change = fold_change_vec_mult()),
                                        data.frame(fold_change_label = fold_change_label_vec_mult())))
  
  AA_df_origin_comb_mult_compare_plot_woSamp1 <- reactive({
    AA_df_origin_comb_mult_compare_plot()[(nrow(AA_df_origin1b())+1):nrow(AA_df_origin_comb_mult_compare_plot()),]
    })
  
  ggplot_intensity_mult <- eventReactive(input$mult_go, {
    if (input$mult_sample_comparison == "Difference") {
      plot_difference_comb_mult(AA_df_origin_comb_mult_compare_plot_woSamp1(), protein_obj1()[2],
                                intensity_label = input$intensity_metric)
    } else if (input$mult_sample_comparison == "Fold Change"){
      plot_foldchange_comb_mult(AA_df_origin_comb_mult_compare_plot_woSamp1(), protein_obj1()[2],
                                intensity_label = input$intensity_metric)
    } else if (input$mult_sample_comparison == "Overlay"){
      if (input$disp_origin) {
        plot_origin_comb(AA_df_origin_comb_mult(), protein_obj1()[2],
                         intensity_label = input$intensity_metric)
      } else {
        plot_intensity_comb(AA_df_intensity_comb_mult(), protein_obj1()[2],
                            intensity_label = input$intensity_metric)
      }
    }
    
  })
  
  observeEvent(input$mult_go, { 
    lapply(2:number(), function(i){ 
      output[[paste0("mult_sample_count",i)]] <- renderText({
        paste0("Samples/Replicates Combined: ", peptidesmult_count()[[i-1]])
      })
      })

  })
  
  output$plot_intensity_mult <- renderPlotly({
    if(input$y_axis_scale == "log"){
      return_plot <- ggplot_intensity_mult() + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                                  breaks = base_breaks())
      create_plotly(return_plot)
    } else {
      create_plotly(ggplot_intensity_mult())
    }
  }) 
 
  
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
        return_plot <- plot_origin(AA_df_origin1(), protein_obj1()[2],
                    intensity_label = input$intensity_metric)
    } else {
        return_plot <- plot_intensity(AA_df_intensity1(), protein_obj1()[2],
                       intensity_label = input$intensity_metric)
    } 
    if(input$y_axis_scale == "log"){
      return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                      breaks = base_breaks())
    }
    return(return_plot)
  })
  
  annotation_plot2 <- reactive({
    if (input$disp_overlay_annot == "Two Samples") {
      return_plot <- add_annotation_layer(plot = ggplot_intensity2(),
                           annotation_df = annotation_df(),
                           color = input$annot_color)
    }else if(input$disp_overlay_annot == "Multiple Samples"){
      return_plot <- add_annotation_layer(plot = ggplot_intensity_mult(),
                                          annotation_df = annotation_df(),
                                          color = input$annot_color)
    } else {
      return_plot <- add_annotation_layer(plot = annotation_plot(),
                           annotation_df = annotation_df(),
                           color = input$annot_color)
    }
    if(input$y_axis_scale == "log"){
      return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                      breaks = base_breaks())
    }
    return(return_plot)
  })
  
  output$annotation_plotly <- renderPlotly(annotation_plot2())
  
  
  stacked_plot_dataframe <- reactive({
    create_stack_line_df(peptide_df = peptides1(),
                         protein = protein_obj1()[1],
                         intensity = input$intensity_metric)
  })
  
  stacked_plot <- reactive({
    if (input$stacked_peptides_yunits == "AA Position") {
      return_plot <- create_stacked_line_plot_yval(stacked_plot_dataframe(),
                                    protein_name = protein_obj1()[2],
                                    protein_seq = protein_obj1()[1])
    } else {
      stacked_intensity_df <- stacked_plot_dataframe()
      stacked_intensity_df$intensity_value[is.na(stacked_intensity_df$intensity_value)] <- 0
      return_plot <- create_stacked_line_plot_intensity(stacked_intensity_df,
                                    protein_name = protein_obj1()[2],
                                    protein_seq = protein_obj1()[1],
                                    intensity_label = input$intensity_metric)
      if(input$y_axis_scale == "log"){
        return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                        breaks = base_breaks())
      }
    }
    return(return_plot)
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
    } else if(input$download_sample_number == "Multiple Samples"){
      covert_to_download_df(AA_df_origin_comb_mult(), input$intensity_metric)
    } else {
      covert_to_download_df(AA_df_origin_comb(), input$intensity_metric)
    }
    
  })
  
  download_peptide_filename <- reactive({
    if (input$download_sample_number == "One sample") {
      paste0(input$sample_name1, "_", input$intensity_metric, "_", str_replace_all(protein_obj1()[2], "\\|", "_"), "_intensity_by_aminoacid.csv")
    } else if(input$download_sample_number == "Multiple Samples"){
      paste0(input$intensity_metric, "_", str_replace_all(protein_obj1()[2], "\\|", "_"), "_intensity_by_aminoacid_multiple_samples.csv")
    }else {
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
  
  
  repeat_vector <- reactive({
    create_repeat_vec(peptides1(),
                      protein_obj1()[1],
                      protein_obj1()[2],
                      database(),
                      intensity = input$intensity_metric)
  })
  
  AA_df_unique_intensity <- eventReactive(input$run_unique, {
    return_df <- AA_df_intensity1()
    return_df$repeated <- repeat_vector()
    return(return_df)
  })
  
  AA_df_unique_origin <- eventReactive(input$run_unique, {
    return_df <- AA_df_origin1()
    return_df$repeated <- repeat_vector()
    return(return_df)
  })
  
  unique_plot <- eventReactive(input$run_unique, {
    if (input$disp_origin) {
      plot <- create_unique_plot_origin(AA_df_unique_origin(),
                                 protein_obj1()[2],
                                 intensity_label = input$intensity_metric)
    } else {
      plot <- create_unique_plot_intensity(AA_df_unique_intensity(),
                                 protein_obj1()[2],
                                 intensity_label = input$intensity_metric)
    }
    if(input$y_axis_scale == "log"){
      plot <- plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                      breaks = base_breaks())
    }
    return(plot)
  })
  
  output$unique_peps <- renderPlotly(create_plotly(unique_plot()))
  
  output$unique_text <- renderText("WARNING: This search can take up to 10 minutes depending on the length of protein")
  
}
  
 