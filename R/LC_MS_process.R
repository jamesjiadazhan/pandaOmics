#' LC_MS_process
#'
#' Process the raw feature table data from LC-MS and make it ready for Metabolome-wide Association Study (MWAS) analysis
#' @import dplyr
#' @import readr
#' @param raw_data The raw feature table file with columns representing samples and rows representing features (metabolites)
#' @param sample_id_file The sample id file that links the sample names with the participant names
#' @param metabolite_start_column The number of the column where metabolites are starting. If metabolite data is starting from column 10, the metabolite_start_column=10
#' @param replicates The number of replicates you performed for each sample. The default is NULL, meaning you don't have technical replicates. If you have 3 replicates for each sample, it should be replicates=3.
#' @param transformation How you would like transform the data for normalization. "log2" and "log10" are available as options
#' @param imputation Whether the missing data due to transformation should be imputated using the QRILC method. Default=NULL. To use imputation, use imputation=TRUE.
#' @param output_name The name of output file using this function. For example, "Met_meantri_log2_HILIC22"
#' @return The complete and clean feature table that is ready for MWAS analysis
#' @examples
#' LC_MS_process(raw_data, sample_id_file, metabolite_start_column=10, replicates=NULL, transformation="log2", imputation=TRUE, output_name)
#' @export

LC_MS_process = function(raw_data, sample_id_file, metabolite_start_column, replicates=NULL, transformation, imputation=NULL, output_name){

  #######Step 1.1 Import and Clean the Metabolomic Feature Table###########
  ##Read in the raw feature table and sample id file
  if (is.character(raw_data) == TRUE){
    ft_data = data.table::fread(raw_data)
  } else {
    ft_data = raw_data
  }
  
  if (is.character(sample_id_file) == TRUE){
    infofile = data.table::fread(sample_id_file)
  } else {
    infofile = sample_id_file
  }

  ##Ensure data quality - filter out features appear less than 5% of all samples & median_CV > 30
  sample_appear = ceiling((dim(ft_data)[2]-10)*0.05)

  if (!("NumPres.All.Samples" %in% colnames(ft_data))){
    NumPres.All.Samples = rowSums(ft_data != 0) - (metabolite_start_column-1)
    ft_data$NumPres.All.Samples = NumPres.All.Samples
    relocate(ft_data, NumPres.All.Samples)
    metabolite_start_column = metabolite_start_column + 1
  }

  mdat = subset(ft_data, ft_data$median_CV<=30 & ft_data$NumPres.All.Samples>sample_appear)

  print("Step 1.1 Import and Clean the Metabolomic Feature Table finished")

  #######Step 1.2 Merge sample id files with the feature table

  ##Import the sample sequence ID file
  colnames(infofile)[1]='File_Name'
  colnames(infofile)[2]='Sample_ID'

  ##Check the correctness of the sample id file
  if (!is.na(table(duplicated(infofile$Sample_ID))['TRUE'])){
    stop("Your sample ID column has duplicates. Please clean the column first")
  }

  ###Merge Info file with feature table sample column names
  ## create a vector of all of the column names (except the first "metabolite_start_column-1" columns)
  colnames.mdat<-colnames(mdat)[metabolite_start_column:dim(mdat)[2]]

  File.n <- as.data.frame(colnames.mdat)

  #Ensure the number of character of the mdat column names matches the number of character of the File_Name in the sample id file
  if (min(nchar(colnames.mdat)) != min(nchar(infofile$File_Name))){
    stop(paste("Your feature table column names (sample id) cannot match the sample id in the sample id file. For example, the first of your sample id in the feature table is", colnames.mdat[1], "but the first of your sample id in the sample id file is", infofile$File_Name[1], "Please make sure they are the same in the first place."))
  }

  #Cleaning up the column title - want to use the same name as the file you are merging to
  colnames(File.n)<-'File_Name'

  ##Merge two datasets and see if the sample info matches with the column names in the feature table
  #need to do this to make sure the data is in the correct order
  infofile.mdat <- merge.data.frame(File.n,infofile,by ='File_Name',all.x = T,sort = F)

  if (sum(is.na(infofile.mdat[,2])) != 0){
    stop("There is a discrepancy between the feature table and sample id file that the merging process failed. Please check those files.")
  }

  print("Step 1.2 Merge sample id files with the feature table finished")

  #######Step 1.3 Replace missing values (0) with NA

  ##Create a matrix with no first "metabolite_start_column-1" columns (only include raw sample ID)
  intensity <- as.matrix(mdat[,metabolite_start_column:dim(mdat)[2]])

  ##Replace missing values (0) with NA
  intensity[intensity == 0] <- NA

  print("Step 1.3 Replace missing values (0) with NA finished")

  ##Replacing the colnames in the feature table *aka, sequence ID) with sample ID
  colnames(intensity) <- infofile.mdat$Sample_ID #assigning sample ID names to the column names in intensity

  ##Add the first "metabolite_start_column-1" columns back to form a complete feature table data with cleaned sample ID
  mdat_comp <- cbind(mdat[,1:(metabolite_start_column-1)],intensity)

  # Optional, if mean intensity has not been calculated
  #######Step 1.4 Calculate Mean Intensity for each metabolic feature across the number of replicates in Metabolomic Feature Table
  if (is.null(replicates) != TRUE){
    column = 1:ncol(intensity)
    ind = as.data.frame(matrix(column,byrow = F,nrow = replicates))
    means = as.data.frame(lapply(ind, function(i) rowMeans(intensity[,i],na.rm = T)))
    colnames(means) = colnames(intensity)[c(T,rep(F,2))]
    mdat_comp = cbind(mdat[,1:(metabolite_start_column-1)],means)
  }

  print("Step 1.4 Calculate Mean Intensity for each metabolic feature across the number of replicates finished")

  #######Step 1.5 log-transformed the data for approximating normal distribution and impute data as needed
  mdat_comp_log = if (transformation == "log2"){
    log2(mdat_comp[,metabolite_start_column:dim(mdat_comp)[2]])
  } else if (transformation == "log10"){
    log10(mdat_comp[,metabolite_start_column:dim(mdat_comp)[2]])
  } else if (transformation != "log10" | transformation != "log2"){
    stop("Data transformation is mandatory. Choose transformation='log10' or transformation='log2'")
  }

  mdat_comp_log2 = cbind(mdat_comp[,1:(metabolite_start_column-1)], mdat_comp_log)

  ###Convert all NaN values to NA
  is.nan.data.frame = function(x)
    do.call(cbind, lapply(x, is.nan))

  mdat_comp_log2[is.nan(mdat_comp_log2)] = NA

  if (imputation == TRUE){
    mdat_comp_log2 = imputeLCMD::impute.QRILC(mdat_comp_log2)
    mdat_comp_log2=as.data.frame(mdat_comp_log2[1])
    print("Imputation completed. The imputation method is: QRILC")
  }

  print("Step 1.5 log-transformed the data for approximating normal distribution finished")

  ##The feature table is complete, clean and now ready for MWAS analysis
  ##save the dataset
  wd = getwd()

  output_name_csv = paste(wd, "/", output_name, ".csv", sep="")

  ##Save MWAS Result
  write.csv(mdat_comp_log2, output_name_csv)

  print(paste("Completed! Output is saved in your current dictionary,", wd, ",using your input output_name,", output_name))

}
