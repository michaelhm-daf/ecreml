
# Write a function to randomly assign each environment to one of the cross-validation groups
# Input E is a character vector of all unique environments
#' @title Defines environment groups for k-fold cross validation
#' @description
#' Generates the cross validation groups for each environment to go into. These groups are then used for the environmental covariate selection procedure.
#' 
#' @param E A character vector of environments included in the multi-environment trial data
#' @param folds The number of fold in the k-fold cross validation scheme
#' 
#' @return A data frame with 2 columns. The first column is the environment term. The 2nd column is the cross validation group that each environment hsa been randomly allocated to.
#' @examples
#' 
#' @export
cv_groups <- function(E, folds=6) {
  env_cv_df <- data.frame(E = rlang::expr(!!E))
  #env_cv_df[[E]] <- env_cv_df$E
  n <- nrow(env_cv_df)
  # Generate the initial clustering groups for each environment
  env_cv_df$init_group <- rep(1:folds, times=ceiling(n/folds))[1:n]
  # Now randomly allocate environments to folds/groups
  env_cv_df$sample <- sample(1:n, replace=F)
  env_cv_df$cv_group <- env_cv_df$init_group[env_cv_df$sample]
  # Convert E and cv_group from character vectors to factors
  env_cv_df$E <- factor(env_cv_df$E)
  env_cv_df$cv_group <- factor(env_cv_df$cv_group)
  return(env_cv_df)
}





#Define %dopar% and %do% locally ----
#`%dopar%` <- foreach::`%dopar%`
#`%do%` <- foreach::`%do%`

#' @title Include an EC into the model and obtain the RMSE
#' @description
#' This function includes an environment covariate (EC) into the model and then performs k-fold cross validation to obtain predictions 
#' from both the baseline model and the model with an environmental covariate included.
#' If management practice \code{.M} and the environment covariates \code{.ec} are factors, then the predictions will be for each unique combination of \code{GxExM}.
#' If \code{.M} is continuous, then predictions will be very for every unique \code{.M} value observed in the multi-environment trial data.
#' 
#' @param .fm The baseline \code{asreml} model object that environmental covariates will be added to. 
#' Note that the data frame used in the baseline model will be the dataframe used to identify each of the corresponding terms in the model.  
#' @param .ec An expression with the environmental covariate to be included in the model.
#' @param .G The genotype term in the model as an expression
#' @param .E The environment term in the model as an expression
#' @param .M The management practice term in the model as an expression
#' @param .trial The trial term  in the model as an expression. \code{.trial} is only required when each environment is not uniquely define by a trial. 
#' An example is when each TOS within a trial constitutes an environment for the purposes of identifying important environmental covariates
#' @param .env_cv_df A data frame identifying which cross-validation group each environment belongs to. 
#' If not provided, the environments will be randomly allocated to cross validation by calling \code{cv_groups} internally
#' @param .cores The number of computer cores used during the cross validation scheme. Note that this should not be more than the number of cross validation groups in the model. 
#' Also note that by setting this greater than 2 (i.e. the default) will require additional \code{ASReml-R} licenses.
#' @param .kn The number of knot points for the spline component of the model for the environment covariate being tested.
#' Note that this is ignored of the environmental covariate is a factor
#' @param .ecs_in_bline_model An optional list of quosures where each element of the list pertains to an environmental covariate that is in the initial baseline model.
#' 
#' @return A list with components:
#' \itemize{
#'  \item \code{Cor}: The Pearson correlation of the predicted values from the baseline model 
#'  compared with the predictions obtained with the environment covariate included in an untested environment. 
#'  The predictions are obtained internally using \code{asreml::predict.asreml}.  
#'  \item \code{Rmse}: the full log-likelihood for each model
#'  \item \code{Resids_cv}: The squared differences between the baseline model and the cross-validation predictions for each \code{GxExM} combination
#'  }
#' @examples
#' 
#' @export
ec_cv_full <- function(.fm, .ec, .G, .E, .M, .trial=NULL, .env_cv_df=NULL,
                       .cores=2, .kn=6, .ecs_in_bline_model=rlang::quos(NULL)){
  # Obtain the data frame from the baseline model
  .df  <- base::eval(.fm$call$data)
  # Now also round all continuous variables in df to 4 decimal places to avoid errors later on due to merging of data frames
  .df <- .df %>% purrr::modify_if(is.numeric, round, digits=4)
  # Remove cv_group from .df if it exists to stop bugs from happening
  if(length(which(colnames(.df)=="cv_group"))>0){
    .df <- .df %>%  dplyr::select(!"cv_group") %>% as.data.frame()
  }
  
  # ADD AN ERROR MESSAGE IF AT() HAS LEVELS NUMBERED INSTEAD OF STATED!!!!!!!!!!!!!!!!!!!!!!!!
  
  # Set each of the inputs as expressions so that the user does not have to make them as expressions prior to input
  # .ec <- enexpr(.ec)
  # .G <- enexpr(.G)
  # .E <- enexpr(.E)
  # .M <- enexpr(.M)
  #.trial <- enexpr(.trial)
  
  # Identify the EC variables that are present in the baseline model
  vars_ec_bl <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist to the ECs that are in the baseline model
  baseline_ec_cols <- c()
  if(!rlang::quo_is_null(.ecs_in_bline_model[[1]])==TRUE){
    baseline_ec_cols <- unlist(purrr::map(.ecs_in_bline_model, rlang::eval_tidy, vars_ec_bl))
  }
  
  # Obtain response variable from the data frame as an expression
  response_term <- attr(.fm$formulae$fixed, "variables")[[2]] %>%
    as.character() %>%
    rlang::parse_expr()
  
  # Merge the data frame with cross validation groupings for environments to the phenotype data
  # First if statement to determine if environment groupings have been provided as an input
  # If not provided, then generate them randomly
  if(is.null(.env_cv_df)==TRUE){
    # Generate ec_cv dataframe
    E_char <- unique(.df[[.E]]) %>% as.character()
    # Now run ec cross-validation function
    .env_cv_df <- cv_groups(E=E_char)
    # Define the column heading for E to be same as it is in the input into ec_cv_full
    .env_cv_df[[.E]] <- .env_cv_df$E
  }
  # If there is no column heading with the same name as E in the original data-frame, then create one
  # Show an error message describing so
  if(is.null(.env_cv_df[[.E]])==TRUE){
    if(is.null(.env_cv_df$.E)==TRUE){
      #stop("The environment term is missing from the provided environment grouping for cross validation")
    } else{
      #If 'E' is specified in the cross-validation dataframe, make the environment column the same as the 'E' column
      .env_cv_df[[.E]] <- .env_cv_df$.E
    }
  }
  
  #make cv_group a factor
  .env_cv_df$cv_group <- factor(.env_cv_df$cv_group)
  .df <- dplyr::left_join(.df, .env_cv_df[, c("cv_group",rlang::as_string(.E))], by=rlang::as_string(.E))
  
  # Make cv_group a factor (should already be a factor though!)
  .df$cv_group <- factor(.df$cv_group)
  
  # Define a table that will be used to generate model predictions later on
  # Identify the classify list and levels of each factor based on whether a trial term is included in the model
  if(rlang::quo_is_null(.ecs_in_bline_model[[1]])==TRUE){
    expr_list <- rlang::exprs(!!.E, !!.G, !!.M, !!.ec)
    # Create list of values we want to predict for
    aux_parallel <- unique(.df[,purrr::map_chr(expr_list, rlang::as_string)]) %>%
      purrr::modify_if(is.numeric, round, digits=4) # Round all continuous variables to 4 decimal places
    # Remove rows where M is NA (need to check that the code still works when M is a factor!!)
    aux_parallel <- aux_parallel[!is.na(aux_parallel[[.M]]),]
    
    # Define the terms corresponding to G, M and ec to appear in the classify statement
    classify_terms <- rlang::expr(!!.G:!!.M:!!.ec)
    # Use expr_text() to add quotations around the expression
    classify_terms <- rlang::expr_text(classify_terms)
    
    # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
    levels_list <- list(aux_parallel[[.G]],
                        aux_parallel[[.M]],
                        aux_parallel[[.ec]])
    
    # Now give the headings of the levels_list names
    names(levels_list) <- c(rlang::expr_text(.G), rlang::expr_text(.M), rlang::expr_text(.ec))
  } else {
    expr_list <- rlang::exprs(!!.E, !!.G, !!.M, !!.ec)
    
    # define a tibble version of .df so that the next line works correctly for 1 baseline EC
    df_tib <- tibble::as_tibble(.df)
    
    # Create list of values we want to predict for
    aux_parallel <- unique(df_tib[, c(purrr::map_chr(expr_list, rlang::as_string),colnames(df_tib[,baseline_ec_cols]))]) %>%
      purrr::modify_if(is.numeric, round, digits=4) %>% # Round all continuous variables to 4 decimal places
      as.data.frame()
    # Remove rows where M is NA (need to check that the code still works when M is a factor!!)
    aux_parallel <- aux_parallel[!is.na(aux_parallel[[.M]]),]
    
    
    bl_ecs_colon <- rlang::parse_expr(paste(colnames(.df[baseline_ec_cols]), collapse=":"))
    # Define the terms corresponding to G, M and ec to appear in the classify statement
    classify_terms <- rlang::expr(!!.G:!!.M:!!bl_ecs_colon:!!.ec)
    # Use expr_text() to add quotations around the expression
    classify_terms <-  gsub("\\(|\\)", "", rlang::expr_text(classify_terms)) # Remove all parentheses from the character string
    # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
    levels_list <- purrr::map(aux_parallel[,-1], as.vector) # Remove the .E column which should always be column 1
  }
  # Used the cross-validation data frame to determine the total number of folds
  v <- length(levels(.df$cv_group))
  
  # Define number of cores in parallel
  total_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(min(c(total_cores[1]-1, .cores)))
  doParallel::registerDoParallel(cl)
  
  # Use foreach to perform parallel programming when implementing the cross validation scheme
  cv_pred <- foreach::foreach(i=c(1:v), .combine=rbind,
                              .packages=c("asreml", "ASExtras4", "tidyverse","rlang")) %dopar% {
                                
                                # Create a subsetted version of the data frame
                                subset_df <-  .df[!.df$cv_group%in%levels(.df$cv_group)[i],]
                                
                                # Identify all the terms that are factors in the model
                                # This grep removes anything that is not a factor including:
                                # (i) anything with a colon
                                # (ii) anything with at()
                                # (iii) anything with mv()
                                # (iv) anything with the pipe operator |
                                # (v) anything which says 'Intercept'
                                # (vi) anything with 'spl' (i.e. any spline terms)
                                which_continuous <- grep(":|at\\(|mv|\\(Intercept\\)|spl", .fm$factor.names)
                                # Also remove the M term from list of factors if M is continuous
                                if(is.double(subset_df[[.M]])==TRUE){
                                  remove_cont_M <- grep(rlang::expr_text(.M), .fm$factor.names)
                                  which_continuous <- unique(c(which_continuous,remove_cont_M))
                                }
                                # Also remove ECs in the model that are continuous
                                # Note: May need to add for loop if multiple ECs are included in the model
                                if(!rlang::quo_is_null(.ecs_in_bline_model[[1]])==TRUE){
                                  for(i in baseline_ec_cols){
                                    bl_ec <- subset_df[,i]
                                    if(is.double(subset_df[[i]])==TRUE){
                                      remove_cont_ec <- grep(colnames(subset_df)[i] , .fm$factor.names)
                                      which_continuous <- unique(c(which_continuous, remove_cont_ec))
                                    }
                                  }
                                }
                                
                                
                                # Probably good to check that the output of factor_terms is as expected
                                factor_terms <- .fm$factor.names[-which_continuous]
                                
                                # Re-factorise all the terms that are factors in the model for the subsetted data-frame
                                #Using mutate() to convert specific columns to factors
                                # Note the different between 'factor' and 'as.factor'
                                # Namely, 'as.factor()' does not update the number levels of each term
                                subset_df <- subset_df %>% dplyr::mutate(dplyr::across(factor_terms, factor))
                                
                                
                                # Gather terms from baseline fixed model
                                # Note that the '+' does not appear if the term is of length 1
                                fixed_bl_terms <- attr(.fm$formulae$fixed, "term.labels") %>% paste(collapse="+") %>%
                                  rlang::parse_expr()
                                # Find out which terms have str in them
                                which_str_bl <- grep("str", attr(.fm$formulae$random, "term.labels"))
                                
                                # CONTINUE HERE 18/04/2024
                                
                                # Remove random terms if the trial does not appear in the cross validation list
                                # First, identify environments that do not appear in the subsetted data frame
                                missing_env <- setdiff(as.character(.env_cv_df[[.E]]), as.character(unique(subset_df[[.E]])))
                                missing_env <- paste(missing_env, collapse="|")
                                which_missing_env <- grep(missing_env, attr(.fm$formulae$random, "term.labels"))
                                # Define an empty character vector to store random terms that need to be modified
                                random_bl_terms_modified <- c()
                                # If trial term included, also identify trials that do not appear in the subsetted data frame
                                if(is.null(.trial)==FALSE){
                                  missing_trial <- setdiff(as.character(unique(.df[[.trial]])), as.character(unique(subset_df[[.trial]])))
                                  if(length(missing_trial) > 0){
                                    missing_trial <- paste(missing_trial, collapse="|")
                                    which_missing_trial <- grep(missing_trial, attr(.fm$formulae$random, "term.labels"))
                                    which_missing_random <- c(which_str_bl, which_missing_env, which_missing_trial)
                                  }
                                } else {
                                  # If there are no missing trials, only remove terms associated with missing environments
                                  which_missing_random <- c(which_str_bl, which_missing_env)
                                }
                                # Gather terms from baseline random model that will be included in the ec model
                                if(length(which_missing_random)>0){
                                  random_bl_terms_removed <- attr(.fm$formulae$random, "term.labels")[which_missing_random]
                                  for(j in 1:length(random_bl_terms_removed)){
                                    # If the removed term contains at(.trial) and multiple levels, only remove the element consisting of the missing trial
                                    if(!is.null(.trial) & grepl(paste0("at\\(", rlang::expr(!!.trial)),random_bl_terms_removed[j]) &
                                       grepl("c\\(", random_bl_terms_removed[j])){
                                      # Split the expression into individual elements
                                      #elements <- strsplit(random_bl_terms_removed[j], ",")[[1]]
                                      elements <- stringr::str_split_1(random_bl_terms_removed[j], ",")
                                      n_ele <- length(elements)
                                      # Remove the desired element
                                      elements_rem <- elements[!grepl(missing_trial, elements)]
                                      # Reconstruct the expression
                                      output_expression <- paste(elements_rem, collapse = ",")
                                      # Add c() out the front if it is the first element
                                      if(grepl(missing_trial, elements)[2]==TRUE){
                                        output_elements <- stringr::str_split_fixed(output_expression, ",", 2)
                                        random_bl_terms_removed[j] <- paste0(output_elements[,1], ", c(", output_elements[,2])
                                      } else if(grepl(missing_trial, elements)[n_ele]==TRUE){
                                        # Add parentheses at the end if it is the last element
                                        random_bl_terms_removed[j] <- paste0(output_expression, "))", ":", stringr::str_split_fixed(random_bl_terms_removed[j], ":", 2)[,2])
                                      } else {
                                        # If the element is not the first nor last, it should be fine to just print the output_expression as is
                                        random_bl_terms_removed[j] <- paste0(output_expression)
                                      }
                                      #Finally, add the term to the character vector of random terms that were modified
                                      random_bl_terms_modified <- c(random_bl_terms_modified, random_bl_terms_removed[j])
                                      
                                      # If the removed term contains at(.E) and multiple levels, only remove the element consisting of the missing environment
                                    } else if(grepl(paste0("at\\(", rlang::expr(!!.E)) , random_bl_terms_removed[j]) & grepl("c\\(", random_bl_terms_removed[j])){
                                      # Split the expression into individual elements
                                      elements <- stringr::str_split_1(random_bl_terms_removed[j], ",")
                                      n_ele <- length(elements)
                                      # Remove the desired element
                                      elements_rem <- elements[!grepl(missing_env, elements)]
                                      # Reconstruct the expression
                                      output_expression <- paste(elements_rem, collapse = ",")
                                      # Add c() out the front if it is the first element
                                      if(grepl(missing_env, elements)[2]==TRUE){
                                        output_elements <- stringr::str_split_fixed(output_expression, ",", 2)
                                        random_bl_terms_removed[j] <- paste0(output_elements[,1], ", c(", output_elements[,2])
                                      } else if(grepl(missing_env, elements)[n_ele]==TRUE){
                                        # Add parentheses at the end if it is the last element
                                        random_bl_terms_removed[j] <- paste0(output_expression, "))", ":",
                                                                             stringr::str_split_fixed(random_bl_terms_removed[j], ":", 2)[,2])
                                      } else {
                                        random_bl_terms_removed[j] <- paste0(output_expression, "))")
                                      }
                                      #Finally, add the term to the character vector of random terms that were modified
                                      random_bl_terms_modified <- c(random_bl_terms_modified, random_bl_terms_removed[j])
                                    }
                                  }
                                }
                                # If the removed term contains at() and multiple levels, only remove the element consisting of the missing trial or environment
                                if(length(which_missing_random)>0){
                                  random_bl_terms <- c(attr(.fm$formulae$random, "term.labels")[-which_missing_random], random_bl_terms_modified)  %>%
                                    paste(collapse = "+") %>%
                                    rlang::parse_expr()
                                } else {
                                  # If there are no random terms to be removed, then keep it the same as the random terms in the baseline model
                                  random_bl_terms <- attr(.fm$formulae$random, "term.labels") %>%
                                    paste(collapse = "+") %>%
                                    rlang::parse_expr()
                                }
                                
                                # Gather terms from baseline residual model
                                residual_bl_terms <- attr(.fm$formulae$residual, "term.labels") %>% paste(collapse="+") %>%
                                  rlang::parse_expr()
                                # Obtain the response variable using the call from the baseline model
                                #response_term <- attr(.fm$formulae$fixed, "variables")[[2]]
                                
                                # Define an expression for the fixed and random effect EC terms to be added in the updated model
                                fixed_ec_terms <- rlang::expr(!!.ec + !!.ec:!!.G + !!.ec:!!.M + !!.ec:!!.G:!!.M)
                                # Change the model depending on whether M is a factor or a variate (i.e. continuous)
                                if(is.double(.df[[.M]])==TRUE){
                                  # Determine the correct str model for the subsetted data
                                  eg <- length(levels(subset_df[[.G]]))*length(levels(subset_df[[.E]]))
                                  if(is.null(.trial)==FALSE){
                                    # Calculate total number of random regression terms for the trial by genotype term
                                    tg <- length(levels(subset_df[[.G]]))*length(levels(subset_df[[.trial]]))
                                    # Define the random regression terms that are embedded within str()
                                    random_str_terms <- rlang::expr(str( ~ !!.trial:!!.G + !!.trial:!!.G:!!.M, ~corh(2):id(!!tg) ) +
                                                                      str(~ !!.E:!!.G + !!.E:!!.G:!!.M, ~corh(2):id(!!eg) ))
                                  } else {
                                    random_str_terms <- rlang::expr(str(~ !!.E:!!.G + !!.E:!!.G:!!.M, ~corh(2):id(!!eg) ))
                                  }
                                  if(is.double(.df[[.ec]])==TRUE){
                                    random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G +
                                                                     spl(!!.ec, k=!!.kn):!!.M + spl(!!.ec, k=!!.kn):!!.G:!!.M +
                                                                     !!.ec:spl(!!.M, k=!!.kn) + !!.ec:!!.G:spl(!!.M, k=!!.kn) +
                                                                     spl(!!.ec, k=!!.kn):spl(!!.M, k=!!.kn) +
                                                                     spl(!!.ec, k=!!.kn):!!.G:spl(!!.M, k=!!.kn))  #Remember to change str() and spatial effects as required
                                    
                                    subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                              random =~ !!random_bl_terms + !!random_ec_terms + !!random_str_terms,
                                                                              residual=~ !!residual_bl_terms,
                                                                              data=subset_df,
                                                                              na.action=asreml::na.method(x='include'),
                                                                              aom=T, maxit=300))
                                  } else if(is.double(.df[[.ec]])==FALSE) {
                                    random_ec_terms <- rlang::expr(!!.ec:spl(!!.M, k=!!.kn) + !!.ec:!!.G:spl(!!.M, k=!!.kn))
                                    subset_call <- rlang::expr(asreml::update.asreml(.fm,
                                                                                     fixed= . ~ . + !!fixed_ec_terms,
                                                                                     random =~ . + !!random_ec_terms,
                                                                                     residual=~ .,
                                                                                     data=subset_df,
                                                                                     #G.param=list(), R.param=list(), # to ensure the model does not get stuck in a local maxima
                                                                                     aom=T, maxit=300))
                                  } else {
                                    stop("Need to specify whether the environmental covariate is a continuous or categorical variable")
                                  }
                                } else if(is.factor(.df[[.M]])==TRUE) {
                                  if(is.double(.df[[.ec]])==TRUE){
                                    random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G +
                                                                     spl(!!.ec, k=!!.kn):!!.M + spl(!!.ec, k=!!.kn):!!.G:!!.M)
                                    subset_call <- rlang::expr(asreml(fixed= !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                      random =~ !!random_bl_terms + !!random_ec_terms,
                                                                      residual=~ !!residual_bl_terms,
                                                                      data=subset_df,
                                                                      #G.param=list(), R.param=list(), # to ensure the model does not get stuck in a local maxima
                                                                      aom=T, maxit=300))
                                  } else if(is.factor(.df[[.M]])==TRUE){
                                    # Generate the call to the updated asreml model
                                    subset_call <- rlang::expr(asreml::update.asreml(.fm,
                                                                                     fixed= . ~ . + !!fixed_ec_terms,
                                                                                     random =~ .,
                                                                                     residual=~ .,
                                                                                     data=subset_df,
                                                                                     #G.param=list(), R.param=list(), # to ensure the model does not get stuck in a local maxima
                                                                                     aom=T, maxit=300))
                                    # evaluate the updated asreml model
                                  } else {
                                    stop("Need to specify whether the management practice is a categorical or numerical variable")
                                  }
                                }
                                # Run the asreml model
                                subset_fm <- eval(subset_call)
                                
                                #missing_env
                                # Rep and Mainplot are equivalent for Emerald since there is only 1 level of TOS
                                # temp_fm <- asreml(fixed = Yield..total...t.ha. ~ Hybrid + TargetPop + Hybrid:TargetPop +
                                #                     (PrePAW + PrePAW:Hybrid + PrePAW:TargetPop + PrePAW:Hybrid:TargetPop),
                                #                   random = ~at(Trial):Rep + at(Trial):MainPlot +
                                #                     at(Trial, "Breeza Dry"):Column + Trial + Env + Trial:Hybrid +
                                #                     Trial:Hybrid:TargetPop + Env:Hybrid + Env:Hybrid:TargetPop +
                                #                     #at(Trial, c("Breeza Dry", "Breeza Irr")) +
                                #                     (spl(PrePAW, k = 6) + spl(PrePAW, k = 6):Hybrid +
                                #                        spl(PrePAW, k = 6):TargetPop + spl(PrePAW, k = 6):Hybrid:TargetPop),
                                #                   residual = ~dsum(~units | ResidualTOS), data = subset_df,
                                #                   aom = T, maxit = 300)
                                
                                # Obtain predictions for all environments (tested & untested) for subsetted model
                                subset_pred <- asreml::predict.asreml(object=subset_fm,
                                                                      classify= classify_terms,
                                                                      levels=levels_list,
                                                                      parallel=T, maxit=1, pworkspace="1gb")
                                # Incorporate environment information into the table of predictions
                                temp_subset.pred <- dplyr::left_join(subset_pred$pvals, aux_parallel,
                                                                     by=colnames(aux_parallel)[-1]) # Assumes the 1st column is always .E
                                
                                # Subset so that the data frame only contains predictions from untested environments
                                temp_subset_pred <- temp_subset.pred[temp_subset.pred[[.E]]%in%setdiff(as.character(.env_cv_df[[.E]]),
                                                                                                       as.character(unique(subset_df[[.E]]))), ]
                                
                                temp_subset_pred # Need this at the end of parallel loop to tell the loop what to perform the rbind over
                              }
  
  # Return the excess asreml licenses back to the RLM server
  parallel::stopCluster(cl)
  
  
  # Need to include trial in predict data frame to remove error message!!!
  # Note that instead I have removed the trial information
  cv_df <- dplyr::select(.df, .G, .M, .E, .ec, response_term) %>%
    dplyr::left_join(cv_pred, . ,
                     by=c(rlang::expr_text(.G), rlang::expr_text(.M), rlang::expr_text(.E), rlang::expr_text(.ec) )) %>%
    tidyr::drop_na(.M)
  
  # Calculate the correlation between the response variable and the CV predictions
  Cor <- stats::cor(cv_df$predicted.value, cv_df[[response_term]], use="pairwise.complete.obs")
  #plot(cv_df$predicted.value, cv_df[[response_term]])
  
  # Calculate the squared differences between the response variable and the CV predictions
  Resids_cv <- (cv_df[[response_term]] - cv_df$predicted.value)^2
  # Finally calculate the RMSE between the response variable and the CV predictions
  Rmse <- sqrt(mean(Resids_cv, na.rm=T))
  
  return_list <- list(Cor=Cor, Rmse = Rmse, Resids_cv = Resids_cv)
  return(return_list)
}





#' @title Include multiple ECs into the model one at a time
#' @description
#' This function implements the forward selection procedure for each environmental covariates (EC).
#' Each EC is incorporate into the model individually and then k-fold cross validation is performed.
#' Predcitions are then obtained in an untested environment, which are then compared back to the baseline model to obtain an RMSE for each environment.
#' The EC which minimises the RMSE is selected from the forward selection procedure. 
#' If management practice \code{.M} and the environment covariates \code{.ec} are factors, then the predictions will be for each unique combination of \code{GxExM}.
#' If \code{.M} is continuous, then predictions will be very for every unique \code{.M} value observed in the multi-environment trial data.
#' 
#' @param fm The baseline \code{asreml} model object that environmental covariates will be added to. 
#' Note that the data frame used in the baseline model will be the dataframe used to identify each of the corresponding terms in the model.  
#' @param ECs An expression with the environmental covariate to be included in the model.
#' @param G The genotype term in the model as an expression
#' @param E The environment term in the model as an expression
#' @param M The management practice term in the model as an expression
#' @param trial The trial term  in the model as an expression. \code{.trial} is only required when each environment is not uniquely define by a trial. 
#' An example is when each TOS within a trial constitutes an environment for the purposes of identifying important environmental covariates
#' @param env_cv_df A data frame identifying which cross-validation group each environment belongs to. 
#' If not provided, the environments will be randomly allocated to cross validation by calling \code{cv_groups} internally
#' @param cores The number of computer cores used during the cross validation scheme. Note that this should not be more than the number of cross validation groups in the model. 
#' Also note that by setting this greater than 2 (i.e. the default) will require additional \code{ASReml-R} licenses.
#' @param kn The number of knot points for the spline component of the model for the environment covariate being tested.
#' Note that this is ignored of the environmental covariate is a factor
#' @param ecs_in_bline_model An optional list of quosures where each element of the list pertains to an environmental covariate that is in the initial baseline model.
#' 
#' @return A list with components:
#' \itemize{
#'  \item \code{ec_selected}: The environmental covariate seected as minimising the RMSE (a scalar of type character).
#'  compared with the predictions obtained with the environment covariate included in an untested environment. 
#'  The predictions are obtained internally using \code{asreml::predict.asreml}.  
#'  \item \code{summary_ecs}: A data frame with the RMSE and Pearson correlation for each environmental covariate considered during the forward selection procedure.
#'  }
#' @examples
#' 
#' @export
ec_finder <- function(fm, ECs, G, E, M, env_cv_df=NULL,
                      cores=2, kn=6, trial=NULL, ecs_in_bline_model=rlang::quos(NULL))
{
  df  <- base::eval(fm$call$data)
  
  # Turn all the relevant inputs into expressions (except for trial)
  # G <- rlang::enexpr(G)
  # E <- rlang::enexpr(E)
  # M <- rlang::enexpr(M)
  # Only make ecs in baseline model an expression if it is not empty
  # if(is.null(ecs_in_bline_model)==FALSE){
  #   ecs_in_bline_model <- enexpr(ecs_in_bline_model)
  # }
  
  # If trial is not a NULL, change it from type NULL to an expression
  # if(trial=="NULLExpression"){
  #   trial <- NULL
  # } else {
  #   trial <- enexpr(trial)
  # }
  
  # Put ECs in a quosure
  # quo_ecs <- rlang::enquos(ECs)
  # quo_ecs_bl <- rlang::enquos(ecs_in_bline_model)
  
  # Don't need to if ECs are already a quosure
  quo_ecs <- ECs
  quo_ecs_bl <- ecs_in_bline_model
  # Send an error message if an EC appears both in the list of ECs to be searched as well as the list of ECs in the baseline model
  
  # quo_ecs <-  quos(PrePAW:ISW, PreFlwT:SeedSet)
  # quo_ecs <- quos(c(PrePAW, PreCumRad))
  # quo_ecs <-  quos(ECs)
  
  # Identify the EC variables
  vars <- as.list(magrittr::set_names(seq_along(df), names(df)))
  # Identify which columns in the data frame consist to the ECs to be assessed
  cols <- unlist(purrr::map(quo_ecs, rlang::eval_tidy, vars))
  
  # Create a data frame to store the RMSE for each EC
  #ec_df <- data.frame(expr(!!ECs) = unique(df[[E]]))
  ec_df <- data.frame(EC = unique(colnames(df[cols])))
  ec_df$Cor <- NA
  ec_df$RMSE <- NA
  
  for(i in 1:dim(ec_df)[1]){
    ec_curr <- rlang::parse_expr(ec_df$EC[i])
    # ec_cv_temp <- ec_cv_full(fm=YieldInit.fm, ec=ec_curr, G= expr(Hybrid), E = expr(Env),
    #            M = expr(density), trial=expr(Trial), env_cv_df=sorghum_env_cv_df,
    #            cores = 6, kn=6)
    ec_cv_temp <- ec_cv_full(.fm=fm, .ec=ec_curr, .G= G, .E = E,
                             .M = M, .trial=trial, .env_cv_df = env_cv_df,
                             .cores=cores, .kn=kn,
                             .ecs_in_bline_model=quo_ecs_bl)
    
    # Store correlation and RMSE in EC data frame ec_df
    ec_df$Cor[i] <- ec_cv_temp$Cor
    ec_df$RMSE[i] <- ec_cv_temp$Rmse
  }
  # Identify the EC with the lowest RMSE
  # Need to check that 'which' combined with 'min' and 'na.rm=TRUE' does not cause a bug
  ec_selected <- ec_df$EC[which(ec_df$RMSE == min(ec_df$RMSE,na.rm=TRUE))]
  return_ecs_list <- list(ec_selected=ec_selected, summary_ecs = ec_df)
  return(return_ecs_list)
}


