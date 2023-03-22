#' Multi-Regression Analysis for 2-Omics Data with FDR Correction
#'
#' @param data_comp Data frame containing the data for the analysis.
#' @param sig_feature_list List of significant features for analysis from the 1st omics results. Make sure the list is a vector
#' @param outcome_start_column Integer indicating the starting column number for outcome variables, which is the 2nd omics.
#' @param outcome_end_column Integer indicating the ending column number for outcome variables, which is the 2nd omics.
#' @param confounders Character vector of confounder variable names. For example, if your interested confounders are "Age_at_enrollment", "Gender", "BMI_on_enrollment", outcome=c("Age_at_enrollment", "Gender", "BMI_on_enrollment"). Note: no space should be in the column names, not even with 'Age at enrollment', which would cause errors.
#' @param suffix Character string to be added as a suffix to the output file name.
#'
#' @return A list of data frames containing the results of the multi-regression analysis for each significant feature, saved as CSV files in the working directory.
#' @export

multi_reg_2omics_fdr = function(data_comp, sig_feature_list, 
                                outcome_start_column, outcome_end_column,
                                confounders, suffix){
  
  outcome_column_num = outcome_end_column-outcome_start_column+1
  
  nam_list = list()
  
  for (j in 1:length(sig_feature_list)) 
  {
    df = data.frame(matrix(nrow = outcome_column_num, ncol = 5))
    colnames(df) = c("feature","b","std","t","p")
    nam_list[j] = list(df)
    names(nam_list)[j] = paste(sig_feature_list[j], "_", suffix, sep="")
    
    
    current_sig_feature_list = sig_feature_list[j]
    
    
    for (i in outcome_start_column:outcome_end_column) 
    {
      current_outcome = colnames(data_comp)[i]
      
      if (is.null(confounders)==TRUE){
        glm_formula = as.formula(paste0(current_sig_feature_list, "~", current_outcome))
      } else(
        glm_formula = as.formula(paste0(current_sig_feature_list, "~", current_outcome, "+", paste(confounders, collapse = "+")))
      )
      
      
      if (is.null(confounders)==TRUE & i == 1){
        print(paste("The first formula you selected is ", paste(current_sig_feature_list, "~", current_outcome))
        )
      } else if (is.null(confounders)==FALSE & i == 1){
        print(paste("The first formula you selected is ", paste(current_sig_feature_list, "~", current_outcome, "+", paste(confounders, collapse = "+")))
        )
      }
      
      tryCatch(
        {
          glmfit<-glm(glm_formula, data=data_comp, family=gaussian())
          nam_list[[j]][i-1,1]<-current_outcome
          nam_list[[j]][i-1,2]<-summary(glmfit)$coefficients[2,1]
          nam_list[[j]][i-1,3]<-summary(glmfit)$coefficients[2,2]
          nam_list[[j]][i-1,4]<-summary(glmfit)$coefficients[2,3]
          nam_list[[j]][i-1,5]<-summary(glmfit)$coefficients[2,4]
        }, 
        error=function(e){}
      )
      
      nam_list[[j]]$pfdr<-p.adjust(nam_list[[j]]$p,method="BH")
      
      nam_list[[j]]=arrange(nam_list[[j]], p)
      
    }
    
    wd = getwd()
    
    output_name_csv = paste(wd, "/", paste(current_sig_feature_list, "_", suffix, sep=""), ".csv", sep="")
    
    ##Save MWAS Result
    readr::write_csv(nam_list[[j]], output_name_csv)
    
    
    print(paste("Completed! Output is saved in your current dictionary,", wd, ",using the output_name,", paste(current_sig_feature_list, "_", suffix, sep="")))
    
  }
}