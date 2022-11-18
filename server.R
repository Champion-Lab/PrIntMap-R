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
      read_peptide_csv_PEAKS_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, comb_method = input$combined_method)
    } else if (input$file_type == "MSFragger" && input$combinedbool1 == "Individual Sample"){
      read_peptide_tsv_MSFragger_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "MSFragger" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_MSFragger_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, comb_method = input$combined_method)
    } else if (input$file_type == "MaxQuant" ) {
      read_peptide_tsv_MaxQuant_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, comb_method = input$combined_method)
    } else if (input$file_type == "Proteome Discoverer" && input$combinedbool1 == "Individual Sample") {
      read_peptide_tsv_ProteomeDiscover_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "Proteome Discoverer" && input$combinedbool1 == "Combined"){
      read_peptide_tsv_ProteomeDiscover_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool1 == "Individual Sample"){
      read_peptide_tsv_Metamorpheus_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_Metamorpheus_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, comb_method = input$combined_method)
    } else if (input$file_type == "Generic csv" && input$combinedbool1 == "Individual Sample") {
      read_peptide_csv_generic_bysamp(input$peptide_file1$datapath)
    } else if (input$file_type == "Generic csv" && input$combinedbool1 == "Combined") {
      read_peptide_csv_generic_comb(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, comb_method = input$combined_method)
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
    tipify(radioButtons(inputId = "intensity_metric",
                 label = "Intensity Metric",
                 choices = mychoicesint), "Select which metric should be used for the y-axis. This applies to all plots. Different input data may have different otpions.")
  })
  
  output$combined_method_display <- renderUI({
    if (input$combinedbool1 == "Combined") {
      tipify(radioButtons(inputId = "combined_method",
                   label = "Combination Method",
                   choices = c("Sum", "Average")), "This method will be used to combine data for all included columns. Options to sum data or average data.")
    }
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
      need(nrow(peptides1()) > 0, "Invalid peptide file. Check format.")
    )
    if (input$disp_origin) {
        return_plot <- plot_origin(AA_df_origin1(), protein_obj1()[2],
                                   intensity_label = input$intensity_metric)
        if(input$displayAnnotations) {
          return_plot <- add_annotation_layer(plot = return_plot,
                                              annotation_df = annotation_df(),
                                              color = input$annot_color)
        }
        
        if(input$displayPTMs) {
          return_plot <- add_PTM_layer_origin(plot = return_plot,
                                              PTM_df = PTM_df_plot_bound())
        }

        if(input$y_axis_scale == "log"){
          return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                          breaks = base_breaks())
        }
        
        
    } else {
        return_plot <- plot_intensity(AA_df_intensity1(), protein_obj1()[2],
                         intensity_label = input$intensity_metric)
        
        if(input$displayAnnotations) {
          return_plot <- add_annotation_layer(plot = return_plot,
                                              annotation_df = annotation_df(),
                                              color = input$annot_color)
        }
        
        if(input$displayPTMs) {
          return_plot <- add_PTM_layer(plot = return_plot,
                                              PTM_df = PTM_df_plot_bound())
        }
        
        if(input$y_axis_scale == "log"){
          return_plot <- return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                          breaks = base_breaks())
        }
        
        
    } 
    create_plotly(return_plot)
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
        read_peptide_csv_PEAKS_comb(peptide_file2(), sample_pattern = input$sample_regex2, comb_method = input$combined_method)
      } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Individual Sample"){
        read_peptide_tsv_MSFragger_bysamp(peptide_file2())
      } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Combined") {
        read_peptide_tsv_MSFragger_comb(peptide_file2(), sample_pattern = input$sample_regex2, comb_method = input$combined_method)
      } else if (input$file_type == "MaxQuant" && input$combinedbool2 == "Combined") {
        read_peptide_tsv_MaxQuant_comb(peptide_file2(), sample_pattern = input$sample_regex2, comb_method = input$combined_method)
      } else if (input$file_type == "Proteome Discoverer" && input$combinedbool2 == "Individual Sample") {
        read_peptide_tsv_ProteomeDiscover_bysamp(peptide_file2())
      } else if (input$file_type == "Proteome Discoverer" && input$combinedbool2 == "Combined"){
        read_peptide_tsv_ProteomeDiscover_bysamp(peptide_file2())
      } else if (input$file_type == "MetaMorpheus" && input$combinedbool2 == "Individual Sample"){
        read_peptide_tsv_Metamorpheus_bysamp(peptide_file2())
      } else if (input$file_type == "MetaMorpheus" && input$combinedbool2 == "Combined") {
        read_peptide_tsv_Metamorpheus_comb(peptide_file2(), sample_pattern = input$sample_regex2, comb_method = input$combined_method)
      } else if (input$file_type == "Generic csv" && input$combinedbool2 == "Individual Sample"){
        read_peptide_csv_generic_bysamp(peptide_file2())
      } else if (input$file_type == "Generic csv" && input$combinedbool2 == "Combined") {
        read_peptide_csv_generic_comb(peptide_file2(), sample_pattern = input$sample_regex2, comb_method = input$combined_method)
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
      
    } else {
      return_plot <- ggplot_intensity2()
      
    }
    
    if(input$displayAnnotations) {
      return_plot <- add_annotation_layer(plot = return_plot,
                                          annotation_df = annotation_df(),
                                          color = input$annot_color)
    }
    
    if(input$displayPTMs) {
      if(input$two_sample_comparison == "Difference" | 
         input$two_sample_comparison == "Fold Change"){
        stop("PTM only plotted on Overlay Graphs")
      }
      
      if(input$disp_origin) {
        return_plot <- add_PTM_layer_origin(plot = return_plot, 
                                            PTM_df = PTM_df_plot_bound())
        return_plot <- add_PTM_layer_origin(plot = return_plot, PTM_df = PTM_df_plot_bound2())
      } else {
        return_plot <- add_PTM_layer(plot = return_plot, 
                                            PTM_df = PTM_df_plot_bound())
        return_plot <- add_PTM_layer(plot = return_plot, PTM_df = PTM_df_plot_bound2())
      }
     
    }
    
    create_plotly(return_plot) 
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
                                    sample_pattern = input[[paste0("sample_regex_mult", i)]], comb_method = input$combined_method)
      } else if (input$file_type == "MSFragger" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample"){
        read_peptide_tsv_MSFragger_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "MSFragger" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_tsv_MSFragger_comb(peptide_file_mult()[[i-1]],
                                        sample_pattern = input[[paste0("sample_regex_mult", i)]], comb_method = input$combined_method)
      } else if (input$file_type == "MaxQuant" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_tsv_MaxQuant_comb(peptide_file_mult()[[i-1]], 
                                       sample_pattern = input[[paste0("sample_regex_mult", i)]], comb_method = input$combined_method)
      } else if (input$file_type == "Proteome Discoverer" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample") {
        read_peptide_tsv_ProteomeDiscover_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "Proteome Discoverer" && input[[paste0("combinedbool_mult",i)]] == "Combined"){
        read_peptide_tsv_ProteomeDiscover_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "MetaMorpheus" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample"){
        read_peptide_tsv_Metamorpheus_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "MetaMorpheus" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_tsv_Metamorpheus_comb(peptide_file_mult()[[i-1]], sample_pattern = input[[paste0("sample_regex_mult", i)]], comb_method = input$combined_method)
      } else if (input$file_type == "Generic csv" && input[[paste0("combinedbool_mult",i)]] == "Individual Sample"){
        read_peptide_csv_generic_bysamp(peptide_file_mult()[[i-1]])
      } else if (input$file_type == "Generic csv" && input[[paste0("combinedbool_mult",i)]] == "Combined") {
        read_peptide_csv_generic_comb(peptide_file_mult()[[i-1]], sample_pattern = input[[paste0("sample_regex_mult", i)]], comb_method = input$combined_method)
      }}
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
  
  ggplot_intensity_mult2 <- eventReactive(input$mult_go, {
    if(input$y_axis_scale == "log"){
      return_plot <- ggplot_intensity_mult() + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                                                  breaks = base_breaks())
      
    } else {
      return_plot <- ggplot_intensity_mult()
    }
    
    if(input$displayAnnotations) {
      return_plot <- add_annotation_layer(plot = return_plot,
                                          annotation_df = annotation_df(),
                                          color = input$annot_color)
    }
    
    if(input$displayPTMs) {
      if(input$mult_sample_comparison == "Difference" | 
         input$mult_sample_comparison == "Fold Change"){
        stop("PTM only plotted on Overlay Graphs")
      }
      return_plot <- add_PTM_layer_origin(plot = return_plot,
                                          PTM_df = PTM_df_plot_bound())
      for(i in 2:number()){
        return_plot <- add_PTM_layer_origin(plot = return_plot,
                                            PTM_df = PTM_df_plot_bound_mult()[[i-1]])}
    }
    return(return_plot)
  })
  
  output$plot_intensity_mult <- renderPlotly({
    create_plotly(ggplot_intensity_mult2())
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
    } else if (input$disp_overlay_annot == "Stacked Peptides") {
      return_plot <- add_annotation_layer(plot = stacked_plot(),
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
  
PTM_regex <- reactive({
  if(input$custom_PTM_check){
    regex_list <- c(input$PTM, input$custom_PTM)
  }else{
    regex_list <- input$PTM
  }
  return(regex_list)
  })

PTM_regex_length <- reactive(length(PTM_regex()))
  
  PTM_df <- reactive({
    validate(
      need(!is.null(PTM_regex()), "no PTM selected")
    )
    if(input$file_type == "PEAKS" | input$file_type == "Generic csv"){
      lapply(1:PTM_regex_length(), function(i){
        create_PTM_df_PEAKS(peptide_df = peptides1(), 
                            protein = protein_obj1()[1], 
                            regex_pattern = PTM_regex()[[i]],
                            intensity = input$intensity_metric)
      })
    }else{
      stop("PTMs are not supported for this search software.")
    }
    
})
  PTM_df2 <- reactive({
    lapply(1:PTM_regex_length(), function(i){
      create_PTM_df_PEAKS(peptide_df = peptides2(), 
                          protein = protein_obj1()[1], 
                          regex_pattern = PTM_regex()[[i]],
                          intensity = input$intensity_metric)
    })
  })
  
  PTM_stacked <- reactive({
    lapply(1:PTM_regex_length(), function(i){
    create_PTM_df_stacked(peptide_df = PTM_df()[[i]], 
                          protein = protein_obj1()[1],
                          regex_pattern = PTM_regex()[[i]],
                          stacked_df = stacked_plot_dataframe()
                          )
      })
    })
  
  PTM_stacked_bound <- reactive(bind_rows(PTM_stacked()))
  
  
  PTM_df_mult <- reactive({
    lapply(2:number(), function(i){
      lapply(1:PTM_regex_length(), function(j){
        create_PTM_df_PEAKS(peptide_df = peptidesmult()[[i-1]], 
                            protein = protein_obj1()[1], 
                            regex_pattern = PTM_regex()[[j]],
                            intensity = input$intensity_metric)
      })
    })
  }) 
  
 
  PTM_df_plot <- reactive({
    lapply(1:PTM_regex_length(), function(i){
    create_PTM_plot_df_PEAKS(peptide_df = PTM_df()[[i]], protein = protein_obj1()[1],
                               regex_pattern = PTM_regex()[[i]], 
                             intensity_vector = intensity_vec1(),
                             origin_pep_vector = origin_vec1())
    })
  })
  
  PTM_df_plot_bound <- reactive({
        bind_rows(PTM_df_plot())
    })
  
  output$test_df <- renderTable(PTM_df_plot_bound())
  
  PTM_df_plot2 <- reactive({
    lapply(1:PTM_regex_length(), function(i){
      create_PTM_plot_df_PEAKS(peptide_df = PTM_df2()[[i]], protein = protein_obj1()[1],
                               regex_pattern = PTM_regex()[[i]], 
                               intensity_vector = intensity_vec2(),
                               origin_pep_vector = origin_vec2())
    })
  })
  
  PTM_df_plot_bound2 <- reactive(bind_rows(PTM_df_plot2()))
  
  PTM_df_plot_mult <- reactive({
    lapply(2:number(), function(i){
      lapply(1:PTM_regex_length(), function(j){
        create_PTM_plot_df_PEAKS(peptide_df = PTM_df_mult()[[i-1]][[j]], protein = protein_obj1()[1],
                                 regex_pattern = PTM_regex()[[j]], 
                                 intensity_vector = intensity_vec_list()[[i-1]],
                                 origin_pep_vector = origin_vec_list()[[i-1]])
      })
    })
  })
  
  PTM_df_plot_bound_mult <- reactive({
    return_list <- list()
    for (i in 2:number()) {
      return_list[[i-1]] <- bind_rows(PTM_df_plot_mult()[i-1][[1]])
    }
    print(return_list)
    return(return_list)
  })
  
  PTM_plot <-reactive ({
    
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

  PTM_plot2 <- reactive ({
   if(input$disp_overlay_PTM == "Stacked Peptides"){
      if(input$stacked_peptides_yunits == "AA Position"){
        return_plot <- add_PTM_layer_stacked(plot = stacked_plot(),
                                             PTM_df = PTM_stacked_bound())
      }else{ 
        return_plot <- add_PTM_layer_stacked_inten(plot = stacked_plot(),
                                                   PTM_df = PTM_stacked_bound())}
   }else if(input$disp_origin){ 
      if(input$disp_overlay_PTM == "Two Samples"){
        if(input$two_sample_comparison == "Difference" | 
           input$two_sample_comparison == "Fold Change"){
          stop("PTM only plotted on Overlay Graphs")
        }
        return_plot <- add_PTM_layer_origin(plot = ggplot_intensity2(), 
                                     PTM_df = PTM_df_plot_bound())
        return_plot <- add_PTM_layer_origin(plot = return_plot, PTM_df = PTM_df_plot_bound2())
      }else if(input$disp_overlay_PTM == "Multiple Samples"){
        if(input$mult_sample_comparison == "Difference" | 
           input$mult_sample_comparison == "Fold Change"){
          stop("PTM only plotted on Overlay Graphs")
        }
        return_plot <- add_PTM_layer_origin(plot = ggplot_intensity_mult(),
                                     PTM_df = PTM_df_plot_bound())
        for(i in 2:number()){
          return_plot <- add_PTM_layer_origin(plot = return_plot,
                                       PTM_df = PTM_df_plot_bound_mult()[[i-1]])}
      }else if(input$disp_overlay_PTM == "Annotation"){
        if(input$disp_overlay_annot == "One Sample"){
          return_plot  <- add_PTM_layer_origin(plot = annotation_plot2(),
                                                             PTM_df = PTM_df_plot_bound())
        }else if(input$disp_overlay_annot == "Two Samples"){
          return_plot <- add_PTM_layer_origin(plot = annotation_plot2(), 
                                              PTM_df = PTM_df_plot_bound())
          return_plot <- add_PTM_layer_origin(plot = return_plot, PTM_df = PTM_df_plot_bound2()) 
        }else if(input$disp_overlay_annot == "Multiple Samples"){
          return_plot <- add_PTM_layer_origin(plot = annotation_plot2(),
                                              PTM_df = PTM_df_plot_bound())
          for(i in 2:number()){
            return_plot <- add_PTM_layer_origin(plot = return_plot,
                                                PTM_df = PTM_df_plot_bound_mult()[[i-1]])}
        }else{
          if(input$stacked_peptides_yunits == "AA Position"){
            return_plot <- add_PTM_layer_stacked(plot = annotation_plot2(),
                                                 PTM_df = PTM_stacked_bound())
          }
        else{ 
          return_plot <- add_PTM_layer_stacked_inten(plot = annotation_plot2(),
                                                     PTM_df = PTM_stacked_bound())}
        }
      }else{ 
        return_plot <- add_PTM_layer_origin(plot = PTM_plot(),
                                                PTM_df = PTM_df_plot_bound())}
    }else{
      if(input$disp_overlay_PTM == "Two Samples"){
        if(input$two_sample_comparison == "Difference" | 
           input$two_sample_comparison == "Fold Change"){
          stop("PTM only plotted on Overlay Graphs")
        }
        return_plot <- add_PTM_layer(plot = ggplot_intensity2(), 
                                     PTM_df = PTM_df_plot_bound())
        return_plot <- add_PTM_layer(plot = return_plot, PTM_df = PTM_df_plot_bound2())
      }else if(input$disp_overlay_PTM == "Multiple Samples"){
        if(input$mult_sample_comparison == "Difference" | 
           input$mult_sample_comparison == "Fold Change"){
          stop("PTM only plotted on Overlay Graphs")
        }
        return_plot <- add_PTM_layer(plot = ggplot_intensity_mult(),
                                     PTM_df = PTM_df_plot_bound())
       for(i in 2:number()){
          return_plot <- add_PTM_layer(plot = return_plot,
                                       PTM_df = PTM_df_plot_bound_mult()[[i-1]])}
      }else if(input$disp_overlay_PTM == "Annotation"){
        if(input$disp_overlay_annot == "One Sample"){
          return_plot  <- add_PTM_layer(plot = annotation_plot2(),
                                               PTM_df = PTM_df_plot_bound())
        }else if(input$disp_overlay_annot == "Two Samples"){
          return_plot <- add_PTM_layer(plot = annotation_plot2(), 
                                              PTM_df = PTM_df_plot_bound())
          return_plot <- add_PTM_layer(plot = return_plot, PTM_df = PTM_df_plot_bound2()) 
        }else if(input$disp_overlay_annot == "Multiple Samples"){
          return_plot <- add_PTM_layer(plot = annotation_plot2(),
                                              PTM_df = PTM_df_plot_bound())
          for(i in 2:number()){
            return_plot <- add_PTM_layer(plot = return_plot,
                                                PTM_df = PTM_df_plot_bound_mult()[[i-1]])}
        }else{
          if(input$stacked_peptides_yunits == "AA Position"){
            return_plot <- add_PTM_layer_stacked(plot = annotation_plot2(),
                                                 PTM_df = PTM_stacked_bound())
          }
          else{ 
            return_plot <- add_PTM_layer_stacked_inten(plot = annotation_plot2(),
                                                       PTM_df = PTM_stacked_bound())}}
      }else{
          return_plot <- add_PTM_layer(plot = PTM_plot(),
                                   PTM_df = PTM_df_plot_bound())}
    }
    if(input$y_axis_scale == "log"){
      return_plot <-return_plot + scale_y_continuous(trans = pseudo_log_trans(base = 2),
                                          breaks = base_breaks())
    }
    return(return_plot)
    })



  output$PTM_plotly <- renderPlotly(PTM_plot2())
    
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
    
    if(input$displayAnnotations) {
      return_plot <- add_annotation_layer(plot = return_plot,
                                          annotation_df = annotation_df(),
                                          color = input$annot_color)
    }
    
    if(input$displayPTMs) {
      if(input$stacked_peptides_yunits == "AA Position"){
        return_plot <- add_PTM_layer_stacked(plot = return_plot,
                                             PTM_df = PTM_stacked_bound())
      }else{ 
        return_plot <- add_PTM_layer_stacked_inten(plot = return_plot,
                                                   PTM_df = PTM_stacked_bound())
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
  
  download_volcano_filename <- reactive({
    paste0("Volcano_plot_data_", input$intensity_metric, "_", select_prot_volcano(db = database(),
                                                                                  Accession = input$AccessionID), ".csv")
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

  peptides1_volcano <- reactive({
    validate(
      need(!is.null(input$peptide_file1), "no peptide file provided")
    )
    if (input$file_type == "PEAKS" && input$combinedbool1 == "Individual Sample") {
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "PEAKS" && input$combinedbool1 == "Combined"){
      read_peptide_csv_PEAKS_volcano(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, min_valid_sample = input$min_valid_sample,
                                     intensity_metric = input$intensity_metric)
    } else if (input$file_type == "MSFragger" && input$combinedbool1 == "Individual Sample"){
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "MSFragger" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_MSfragger_volcano(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, min_valid_sample = input$min_valid_sample,
                                         intensity_metric = input$intensity_metric)
    } else if (input$file_type == "MaxQuant" ) {
      read_peptide_tsv_MaxQuant_volcano(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, min_valid_sample = input$min_valid_sample,
                                        intensity_metric = input$intensity_metric)
    } else if (input$file_type == "Proteome Discoverer" && input$combinedbool1 == "Individual Sample") {
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "Proteome Discoverer" && input$combinedbool1 == "Combined"){
      stop("Proteome Discoverer data currently not supported for volcano plots")
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool1 == "Individual Sample"){
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool1 == "Combined") {
      read_peptide_tsv_Metamorpheus_volcano(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, min_valid_sample = input$min_valid_sample,
                                            intensity_metric = input$intensity_metric)
    } else if (input$file_type == "Generic csv" && input$combinedbool1 == "Individual Sample"){
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "Generic csv" && input$combinedbool1 == "Combined") {
      read_peptide_csv_generic_volcano(input$peptide_file1$datapath, sample_pattern = input$sample_regex1, min_valid_sample = input$min_valid_sample,
                                            intensity_metric = input$intensity_metric)
    }
  })

  
  
  peptides2_volcano <- reactive({
    
    if (input$file_type == "PEAKS" && input$combinedbool2 == "Individual Sample") {
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "PEAKS" && input$combinedbool2 == "Combined"){
      read_peptide_csv_PEAKS_volcano(peptide_file2(), sample_pattern = input$sample_regex2, min_valid_sample = input$min_valid_sample,
                                     intensity_metric = input$intensity_metric)
    } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Individual Sample"){
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "MSFragger" && input$combinedbool2 == "Combined") {
      read_peptide_tsv_MSfragger_volcano(peptide_file2(), sample_pattern = input$sample_regex2, min_valid_sample = input$min_valid_sample,
                                     intensity_metric = input$intensity_metric)
    } else if (input$file_type == "MaxQuant" ) {
      read_peptide_tsv_MaxQuant_volcano(peptide_file2(), sample_pattern = input$sample_regex2, min_valid_sample = input$min_valid_sample,
                                        intensity_metric = input$intensity_metric)
    } else if (input$file_type == "Proteome Discoverer" && input$combinedbool2 == "Individual Sample") {
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "Proteome Discoverer" && input$combinedbool2 == "Combined"){
      stop("Proteome Discoverer data currently not supported for volcano plots")
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool2 == "Individual Sample"){
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "MetaMorpheus" && input$combinedbool2 == "Combined") {
      read_peptide_tsv_Metamorpheus_volcano(peptide_file2(), sample_pattern = input$sample_regex2, min_valid_sample = input$min_valid_sample,
                                            intensity_metric = input$intensity_metric)
    } else if (input$file_type == "Generic csv" && input$combinedbool2 == "Individual Sample"){
      stop("Volcano plots can only be created with combined files for variance calculation")
    } else if (input$file_type == "Generic csv" && input$combinedbool2 == "Combined") {
      read_peptide_csv_generic_volcano(peptide_file2(), sample_pattern = input$sample_regex2, min_valid_sample = input$min_valid_sample,
                                            intensity_metric = input$intensity_metric)
    }
  })
  
  combined_volcano <- reactive({
    df <- combine_two_volcano_dfs(peptides1_volcano(), peptides2_volcano(),
                                  min_valid_sample = input$min_valid_sample,
                                  fdr = input$p_cutoff,
                                  fold_change_cutoff_plot = input$l2fc_cutoff,
                                  equal_variance_bool = input$equal_var,
                                  remove_na = input$remove_na,
                                  set_na = input$set_na_value)
  })
  
  volcano_plot_list <- reactive({
    create_volcano_plot(combined_volcano(),
                        fdr = input$p_cutoff,
                        fold_change_cutoff_plot = input$l2fc_cutoff,
                        equal_variance_bool = input$equal_var,
                        intensity_metric = input$intensity_metric,
                        sample1 = input$sample_name1,
                        sample2 = input$sample_name2,
                        BH_correction = input$BH_correction,
                        protein_of_interest = select_prot_volcano(db = database(),
                                                                  Accession = input$AccessionID),
                        display_comp_vals = input$display_comp_vals,
                        display_infinites = input$display_infinite_vals)
  })
  
  volcano_plot <- reactive(volcano_plot_list()[[1]])
  volcano_plot_download_df <- reactive(volcano_plot_list()[[2]])
  
  output$volcano <- renderPlotly({
    if (!is.null(volcano_plot())) {
      create_volcano_plotly(volcano_plot())
    } else {
      NULL
    }
  })
  
  
  output$download_volcano_df <- downloadHandler(
    filename = function() {
      download_volcano_filename()
    },
    content = function(file) {
      write.csv(volcano_plot_download_df(), file, row.names = FALSE)
    }
  )
  
  output$volcano_text <- renderText({
    if (is.null(volcano_plot())) {
      "Not enough information for volcano plot"
    } else {
      NULL
    }
  })
  
  
  output$downloadHSdatabase <- downloadHandler(
    filename <- function() {
      "Human_database_Uniprot_downloaded_20220823.fasta"
    },
    
    content <- function(file) {
      file.copy("example_data/Human_database_Uniprot_downloaded_20220823.fasta", file)
    }
  )
  
  output$downloadPhosphoNoPTM <- downloadHandler(
    filename <- function() {
      "MSFragger.PO4ase.NoPTMsearch.combined_modified_peptide.tsv"
    },
    
    content <- function(file) {
      file.copy("example_data/MSFragger.PO4ase.NoPTMsearch.combined_modified_peptide.tsv", file)
    }
  )
  
  output$downloadPhosphoWithPTM <- downloadHandler(
    filename <- function() {
      "MSFragger.PO4ase.PhosphoPTMsearch.combined_modified_peptide.tsv"
    },
    
    content <- function(file) {
      file.copy("example_data/MSFragger.PO4ase.PhosphoPTMsearch.combined_modified_peptide.tsv", file)
    }
  )
  
  output$downloadHXMS <- downloadHandler(
    filename <- function() {
      "PEAKS.db.peptides.H4.OptimizationHXMS.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/PEAKS.db.peptides.H4.OptimizationHXMS.csv", file)
    }
  )
  
  output$downloadDeglycoLFQ <- downloadHandler(
    filename <- function() {
      "PEAKS.lfq.peptides.HumanSerum.PNGaseF.Deglycosylation.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/PEAKS.lfq.peptides.HumanSerum.PNGaseF.Deglycosylation.csv", file)
    }
  )
  
  output$downloadDeglycoDB <- downloadHandler(
    filename <- function() {
      "PEAKS.db.peptides.HumanSerum.PNGaseF.Deglycosylation.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/PEAKS.db.peptides.HumanSerum.PNGaseF.Deglycosylation.csv", file)
    }
  )
  
  output$downloadDeglycoCTRL <- downloadHandler(
    filename <- function() {
      "PEAKS.db.HumanSerum.PNGaseF_CTRL_1.peptides.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/PEAKS.db.HumanSerum.PNGaseF_CTRL_1.peptides.csv", file)
    }
  )
  
  output$downloadDeglycoPNGaseF <- downloadHandler(
    filename <- function() {
      "PEAKS.db.HumanSerum.PNGaseF_PNG_1.peptides.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/PEAKS.db.HumanSerum.PNGaseF_PNG_1.peptides.csv", file)
    }
  )
  
  output$downloadBGalBSAfusion <- downloadHandler(
    filename <- function() {
      "PEAKS.lfq.peptides.BGalBSAFusion.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/PEAKS.lfq.peptides.BGalBSAFusion.csv", file)
    }
  )
  
  output$downloadHSdatabaseBgalBSAfusion <- downloadHandler(
    filename <- function() {
      "Human_database_Uniprot_BGAL-BSAfusion_downloaded_20220823.fasta"
    },
    
    content <- function(file) {
      file.copy("example_data/Human_database_Uniprot_BGAL-BSAfusion_downloaded_20220823.fasta", file)
    }
  )
  
  output$downloadIndividualGeneric <- downloadHandler(
    filename <- function() {
      "individualGenericTemplate.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/individual.csv", file)
    }
  )
  
  output$downloadCombinedGeneric <- downloadHandler(
    filename <- function() {
      "combinedGenericTemplate.csv"
    },
    
    content <- function(file) {
      file.copy("example_data/combined.csv", file)
    }
  )

}
  
 