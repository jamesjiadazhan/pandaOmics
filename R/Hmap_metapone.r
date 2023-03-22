#' Hmap_metapone
#' Perform a summary of significant features from metapone pathway tables and generate a heatmap for p-value less than 0.05
#' @param path_heatmap A string representing the path of the input files.
#' @param pathway_focus A string representing the pathway focus. For example, V1_VITD_PC. Your final output will be: pathway_focus, "_metapone_overlap", ".csv"
#' @return A data.frame with the metabolite analysis results.
#' @export

Hmap_metapone <- function(path_heatmap, pathway_focus) {
  
  # Create a list of filenames that match the specified pattern
  list.filenames <- list.files(path = path_heatmap, 
                               pattern = paste0(".csv"),
                               full.names = TRUE, 
                               recursive = FALSE) 
  
  # Extract the unique identifier from the filenames
  nchar_path = nchar(path_heatmap) + 2
  res <- substring(list.filenames, nchar_path, nchar(list.filenames) - nchar(".csv"))
  
  # Initialize data frames
  myData1 <- read_csv(list.filenames[1])
  Mum_Pathway <- myData1[order(myData1$pathway),][1]
  Mum_Pathway_V1 <- Mum_Pathway_V2 <- Mum_Pathway
  
  # Loop through the list of filenames to read in and store the data
  for (i in 1:length(list.filenames)) {
    Path_Tem <- read_csv(list.filenames[i]) %>%
      arrange(pathway)
    
    # Store the p-value, overlap size, and pathway size for the current file in the appropriate data frames
    Mum_Pathway[1 + i] <- Path_Tem$p_value
    Mum_Pathway_V1[1 + i] <- Path_Tem$n_mapped_metabolites
    Mum_Pathway_V2[1 + i] <- Path_Tem$n_metabolites
    
    # Assign the unique identifier as the column name for each data frame
    colnames(Mum_Pathway)[1 + i] <- res[i]
    colnames(Mum_Pathway_V1)[1 + i] <- res[i]
    colnames(Mum_Pathway_V2)[1 + i] <- res[i]
  }
  
  # Determine the number of columns in the first data frame
  coln <- dim(Mum_Pathway)[2]
  
  # Calculate the number of pathways with adjusted p-value less than 0.05
  Mum_Pathway <- Mum_Pathway %>%
    mutate(Total0.05 = rowSums(.[-1] <= 0.05))
  
  # Calculate the average overlap size and average pathway size
  Mum_Pathway_V1 <- Mum_Pathway_V1 %>%
    mutate(overlap_average = rowMeans(.[-1], na.rm = TRUE) %>% ceiling())
  Mum_Pathway_V2 <- Mum_Pathway_V2 %>%
    mutate(total_average = rowMeans(.[-1], na.rm = TRUE) %>% ceiling())
  
  # Combine the three data frames into one final data frame
  Mum_Pathway_final <- cbind(Mum_Pathway, Mum_Pathway_V1["overlap_average"], Mum_Pathway_V2["total_average"]) %>%
    mutate(overlap_percent = round(overlap_average / total_average, 3)) %>%
    arrange(desc(Total0.05))
  
  # Write the final data frame to a CSV file
  write_csv(Mum_Pathway_final, file.path(path_heatmap, paste0(pathway_focus, "_metapone_overlap.csv")))
  
}

