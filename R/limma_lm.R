#' limma_lm
#'
#' Perform Metabolome-wide Association Study (MWAS) analysis using multiple linear regression with FDR correction with the limma algorithm
#' @import dplyr
#' @import readr
#' @import limma
#' @param DATA_comp The merged, complete data with demographic data and feature table data path or already-read data. if the data is not read, it has to be csv.
#' @param DATA_met The cleaned feature table data ready for MWAS analysis. It should be data path or already-read data. if the data is not read, it has to be csv.
#' @param output_name The name of output file using this function. For example, "Met_meantri_log2_HILIC22"
#' @param met_start_colnum The number of the column where metabolites are starting. If metabolite data is starting from column 10, the met_start_colnum=10
#' @param met_end_colnum The number of the column where metabolites are ending If metabolite data is ending at column 13414, the met_end_colnum=13414
#' @param outcome The outcome of your interest. For example, if your interested outcome is T2D (type 2 diabetes) and the column name is T2D, outcome="T2D"
#' @param confounders The confounders of your interest. For example, if your interested confounders are "Age_at_enrollment", "Gender", "BMI_on_enrollment", outcome=c("Age_at_enrollment", "Gender", "BMI_on_enrollment"). Note: no space should be in the column names, not even with 'Age at enrollment', which would cause errors.
#' @return The multiple regression results with FDR correction
#' @examples
#' limma_lm(DATA_comp = finaldt1_cleaned, DATA_met = metlog, output_name = "ILD_HILIC_limma_test", met_start_colnum = 15, met_end_colnum = 25318, outcome = "ILD1", confounders=c("Age_at_enrollment", "Gender", "BMI_on_enrollment"))
#' @export

limma_lm = function(DATA_comp, DATA_met, output_name, met_start_colnum, met_end_colnum, outcome, confounders){

  if (is.character(DATA_comp) == TRUE){
    DATA_comp = read_csv(DATA_comp)
  } else {
    DATA_comp = DATA_comp
  }

  if (is.character(DATA_met) == TRUE){
    DATA_met = read_csv(DATA_met)
  } else {
    DATA_met = DATA_met
  }

  # Filter any confounders with NA
  DATA_comp_cleaned = DATA_comp %>%
    tidyr::drop_na(any_of(c(outcome, confounders)))

  if (is.null(confounders)==TRUE){
    limma_formula = as.formula(paste("~", outcome))
  } else(
    limma_formula = as.formula(paste("~", outcome, "+", paste(confounders, collapse = "+")))
  )

  if (is.null(confounders)==TRUE){
    print(paste("The formula you selected is ", paste("~", outcome))
    )
  } else(
    print(paste("The formula you selected is ", paste("~", outcome, "+",
                                                      paste(confounders, collapse = "+")))
    )
  )

  design = model.matrix(limma_formula,
                        data=DATA_comp_cleaned)

  metabolite = DATA_comp_cleaned[, met_start_colnum:met_end_colnum]
  metabolite_t = t(metabolite)

  fit <- lmFit(metabolite_t, design)
  fit2 <- eBayes(fit)

  limma_outcome = cbind(fit2$coefficients[,2], fit2$stdev.unscaled[,2], fit2$t[,2], fit2$p.value[,2])
  colnames(limma_outcome)<-c("b","std","t","p")

  last_colummn=grep("Max.Intensity", colnames(DATA_met))
  output = cbind(DATA_met[,1:last_colummn], limma_outcome)
  output$pfdr = p.adjust(output$p, method="BH")
  output=arrange(output, p)

  wd = getwd()

  output_name_csv = paste(wd, "/", output_name, ".csv", sep="")

  ##Save MWAS Result
  readr::write_csv(output, output_name_csv)


  print(paste("Completed! Output is saved in your current dictionary,", wd, ",using your input output_name,", output_name))

}
