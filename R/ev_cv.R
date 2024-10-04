#' @title Calculate the RMSE for a given model
#' @description
#' This function calculates the RMSE for the current model by implementing the k-fold cross validation and then calculating the
#'  predictions obtained in an untested environment.
#'
#' @param .fm The baseline \code{asreml} model object that environmental covariates will be added to.
#' Note that the data frame used in the baseline model will be the dataframe used to identify each of the corresponding terms in the model.
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
rmse_calc <- function(.fm, .G, .E, .M, .trial=NULL, .env_cv_df=NULL,
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


  #.ecs_in_bline_model <- enexpr(.ecs_in_bline_model)

  # Identify the EC variables that are present in the baseline model
  vars_ec_bl <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist of the ECs that are in the baseline model
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
    } else {
      #If 'E' is specified in the cross-validation dataframe, make the environment column the same as the 'E' column
      .env_cv_df[[.E]] <- .env_cv_df$.E
    }
  }

  # Make cv_group an integer for merging (in case it is already a factor!) (Note, cv_group is not supposed to exist yet!!)
  #.df$cv_group <- as.integer(.df$cv_group)

  .df <- dplyr::left_join(.df, .env_cv_df[, c("cv_group",rlang::as_string(.E))], by=rlang::as_string(.E))

  # Make cv_group a factor (should already be a factor though!)
  .df$cv_group <- factor(.df$cv_group)


  # Define a table that will be used to generate model predictions later on
  # Identify the classify list and levels of each factor based on whether a trial term is included in the model
  if(rlang::quo_is_null(.ecs_in_bline_model[[1]])==TRUE){
    # Change the classify terms if M is missing from the baseline model
    if(is.null(.M)==TRUE){
      expr_list <- rlang::exprs(!!.E, !!.G)
      # Create list of values we want to predict for
      aux_parallel <- unique(.df[,purrr::map_chr(expr_list, rlang::as_string)]) %>%
        purrr::modify_if(is.numeric, round, digits=4) # Round all continuous variables to 4 decimal places

      # Define the terms corresponding to G, M and ec to appear in the classify statement
      classify_terms <- rlang::expr(!!.G)
      # Use expr_text() to add quotations around the expression
      classify_terms <- rlang::expr_text(classify_terms)

      # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
      levels_list <- list(aux_parallel[[.G]])

      # Now give the headings of the levels_list names
      names(levels_list) <- c(rlang::expr_text(.G))
    } else {
      expr_list <- rlang::exprs(!!.E, !!.G, !!.M)
      # Create list of values we want to predict for
      aux_parallel <- unique(.df[,purrr::map_chr(expr_list, rlang::as_string)]) %>%
        purrr::modify_if(is.numeric, round, digits=4) # Round all continuous variables to 4 decimal places
      # Remove rows where M is NA (need to check that the code still works when M is a factor!!)
      aux_parallel <- aux_parallel[!is.na(aux_parallel[[.M]]),]

      # Define the terms corresponding to G, M and ec to appear in the classify statement
      classify_terms <- rlang::expr(!!.G:!!.M)
      # Use expr_text() to add quotations around the expression
      classify_terms <- rlang::expr_text(classify_terms)

      # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
      levels_list <- list(aux_parallel[[.G]],
                          aux_parallel[[.M]])

      # Now give the headings of the levels_list names
      names(levels_list) <- c(rlang::expr_text(.G), rlang::expr_text(.M))
    }
  } else {
    if(is.null(.M)==TRUE){
      expr_list <- rlang::exprs(!!.E, !!.G)

      # define a tibble version of .df so that the next line works correctly for 1 baseline EC
      df_tib <- tibble::as_tibble(.df)

      # Create list of values we want to predict for
      aux_parallel <- unique(df_tib[, c(purrr::map_chr(expr_list, rlang::as_string),colnames(df_tib[,baseline_ec_cols]))]) %>%
        purrr::modify_if(is.numeric, round, digits=4) %>% # Round all continuous variables to 4 decimal places
        as.data.frame()

      bl_ecs_colon <- rlang::parse_expr(paste(colnames(.df[baseline_ec_cols]), collapse=":"))
      # Define the terms corresponding to G, M and ec to appear in the classify statement
      classify_terms <- rlang::expr(!!.G:!!bl_ecs_colon)
      # Use expr_text() to add quotations around the expression
      classify_terms <-  gsub("\\(|\\)", "", rlang::expr_text(classify_terms)) # Remove all parentheses from the character string
      # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
      levels_list <- purrr::map(aux_parallel[,-1], as.vector) # Remove the .E column which should always be column 1
    } else {
      expr_list <- rlang::exprs(!!.E, !!.G, !!.M)

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
      classify_terms <- rlang::expr(!!.G:!!.M:!!bl_ecs_colon)
      # Use expr_text() to add quotations around the expression
      classify_terms <-  gsub("\\(|\\)", "", rlang::expr_text(classify_terms)) # Remove all parentheses from the character string
      # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
      levels_list <- purrr::map(aux_parallel[,-1], as.vector) # Remove the .E column which should always be column 1
    }
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
                                #str(SorghumMET_subset.df)
                                #str(SorghumMET.df)
                                #require(gdata)

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
                                if(is.null(.M)==FALSE){
                                  if(is.double(subset_df[[.M]])==TRUE){
                                    remove_cont_M <- grep(rlang::expr_text(.M), .fm$factor.names)
                                    which_continuous <- unique(c(which_continuous,remove_cont_M))
                                  }
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
                                #str(subset_df)


                                # Gather terms from baseline fixed model
                                # Note that the '+' does not appear if the term is of length 1
                                fixed_bl_terms <- attr(.fm$formulae$fixed, "term.labels") %>% paste(collapse="+") %>%
                                  rlang::parse_expr()
                                # Find out which terms have str in them
                                which_str_bl <- grep("str", attr(.fm$formulae$random, "term.labels"))


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
                                # Change the model depending on whether M is a factor or a variate (i.e. continuous) or missing
                                if(is.null(.M)==TRUE){
                                  if(is.null(.trial)==FALSE){
                                    # Define the random regression terms
                                    random_str_terms <- rlang::expr(!!.trial:!!.G + !!.E:!!.G)
                                  } else {
                                    random_str_terms <- rlang::expr(!!.E:!!.G)
                                  }
                                  subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms,
                                                                            random =~ !!random_bl_terms + !!random_str_terms,
                                                                            residual=~ !!residual_bl_terms,
                                                                            data=subset_df,
                                                                            na.action=asreml::na.method(x='include'),
                                                                            aom=T, maxit=300))
                                 } else if(is.double(.df[[.M]])==TRUE){
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
                                    subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms,
                                                                              random =~ !!random_bl_terms + !!random_str_terms,
                                                                              residual=~ !!residual_bl_terms,
                                                                              data=subset_df,
                                                                              na.action=asreml::na.method(x='include'),
                                                                              aom=T, maxit=300))

                                 } else if(is.factor(.df[[.M]])==TRUE) {
                                    if(is.null(.trial)==FALSE){
                                      # Define the random regression terms
                                      random_str_terms <- rlang::expr(!!.trial:!!.G + !!.trial:!!.G:!!.M + !!.E:!!.G + !!.E:!!.G:!!.M)
                                    } else {
                                      random_str_terms <- rlang::expr(!!.E:!!.G + !!.E:!!.G:!!.M)
                                    }
                                    subset_call <- rlang::expr(asreml::asreml(fixed= !!response_term ~ !!fixed_bl_terms,
                                                                              random =~ !!random_bl_terms + !!random_str_terms,
                                                                              residual=~ !!residual_bl_terms,
                                                                              data=subset_df, aom=T, maxit=300))

                                  } else {
                                    stop("M needs to be continuous of categorical")
                                  }
                          # Run the asreml model
                          subset_fm <- eval(subset_call)

                          # Obtain predictions for all environments (tested & untested) for subsetted model
                          if(length(levels_list)>1){
                            subset_pred <- asreml::predict.asreml(object=subset_fm,
                                                                  classify= classify_terms,
                                                                  levels=levels_list,
                                                                  parallel=T, maxit=1, pworkspace="1gb")
                          } else {
                            subset_pred <- asreml::predict.asreml(object=subset_fm,
                                                                  classify= classify_terms,
                                                                  levels=levels_list,
                                                                  parallel=F, maxit=1, pworkspace="1gb")
                          }

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
  if(is.null(.M)==TRUE){
    # Note that this causes many-to-many warnings
    cv_df <- dplyr::select(.df, .G, .E, response_term) %>%
      dplyr::left_join(cv_pred, . ,
                       by=c(rlang::expr_text(.G), rlang::expr_text(.E) )) %>%
      unique()
  } else {
    cv_df <- dplyr::select(.df, .G, .M, .E, response_term) %>%
      dplyr::left_join(cv_pred, . ,
                       by=c(rlang::expr_text(.G), rlang::expr_text(.M), rlang::expr_text(.E))) %>%
      tidyr::drop_na(.M) %>%
      unique() # Make sure this doesn't cause any bugs
  }
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

#' @title Select the best EC and the significant terms
#' @description
#' This function bring together the entire subset selection process to (i) identify the most important environmental covariate
#' via the forward selection procedure and then drop the non-significant fixed and random terms from the model to achieve a parsimonious model.
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
#' @param ncores The number of computer cores used during the cross validation scheme. Note that this should not be more than the number of cross validation groups in the model.
#' Also note that by setting this greater than 2 (i.e. the default) will require additional \code{ASReml-R} licenses.
#' @param kn The number of knot points for the spline component of the model for the environment covariate being tested.
#' Note that this is ignored of the environmental covariate is a factor
#' @param ecs_in_bline_model An optional list of quosures where each element of the list pertains to an environmental covariate that is in the initial baseline model.
#'
#' @return A list with 2 components:
#' \itemize{
#'  \item \code{fm}: An updated model with the most important environmental covariate included in the model along with only the significant
#'  fixed and random effects pertaining to the environmental covariate identified.
#'  \item \code{rmse}: The value of the RMSE for the updated model.
#'  }
#' @examples
#'
#' @export
ec_iteration <- function(fm, ECs, G, E, M, env_cv_df=NULL, ncores=2, kn=6, trial=NULL, ecs_in_bline_model=rlang::maybe_missing()){
  curr_fm <- fm
  #quo_ecs_bl <- rlang::enquos(ecs_in_bline_model)

  # Set each of the inputs as expressions so that the user does not have to make them as expressions prior to input
  #ECs <- enexpr(ECs)
  # G <- enexpr(G)
  # E <- enexpr(E)
  # M <- enexpr(M)
  #trial <- enexpr(trial)


  # Identify the data frame from the model
  .df <<- base::eval(fm$call$data)

  # Put ECs in a quosure
  #quo_ecs <- rlang::enquos(ECs)

  #quo_ecs <- rlang::quos(PrePAW:PostPAW, PreCumRad)

  # Identify the EC variables
  vars <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist to the ECs to be assessed
  cols_ecs <- unlist(purrr::map(ECs, rlang::eval_tidy, vars))
  # Use the columns to obtain a vector of EC terms
  ec_terms_char <- colnames(.df[cols_ecs])


  if(rlang::is_missing(ecs_in_bline_model)==TRUE){
    ecs_in_bline_model <- rlang::quos(NULL)
    ec_bl_terms_char <- c()
  } else {
    # # Identify which columns in the data frame consist of ECs in the baseline model
    # # Use the columns to obtain a vector of EC terms
    cols_bl_ecs <- unlist(purrr::map(ecs_in_bline_model, rlang::eval_tidy, vars))
    ec_bl_terms_char <- colnames(.df[cols_bl_ecs])
  }


  rmse_curr <- rmse_calc(.fm= curr_fm, .G=G, .E=E, .M=M, .trial=trial, .env_cv_df = env_cv_df,
                         .cores=ncores, .kn=kn, .ecs_in_bline_model=ecs_in_bline_model)$Rmse

  # Remove cv_group from .df if it exists to stop bugs from happening
  if(length(which(colnames(.df)=="cv_group"))>0){
    .df <- .df %>%  dplyr::select(!"cv_group") %>% as.data.frame()
  }

  continue <- TRUE
  while(continue==TRUE){
    ec_search <- ec_finder(fm=curr_fm, ECs= ECs,
                           G= G, E = E,
                           M= M, trial = trial,
                           env_cv_df=env_cv_df,
                           ecs_in_bline_model=ecs_in_bline_model,
                           cores=ncores, kn=kn)

    # Determine the candidate EC model
    ec_summary <- ec_search$summary_ecs
    which_selected <- which(ec_summary$RMSE==min(ec_summary$RMSE, na.rm=TRUE))
    ec_candid <- ec_summary$EC[which_selected]
    rmse_candid <- ec_summary$RMSE[which_selected]
    # If model has improved (i.e. RMSE is lower) than include best EC as a candidate model, then find a parsimonious version of the model
    if(rmse_candid < rmse_curr){

      ec_candidate <- rlang::parse_expr(ec_candid)
      # First define the model with every EC term for the candidate EC
      candid_fm <- ec_full_model_constructor(.fm=curr_fm, .ec=ec_candidate, .G=G, .M=M, .kn=6)

      ec_candidate <- rlang::parse_quos(ec_candid, rlang::global_env())
      # Remove non-significant fixed and random effects for the ECs in the candidate model
      candid_fm <- simplify_ec_model(.fm=candid_fm, .ecs_in_model = ec_candidate, .G=rlang::expr_text(G), .M=rlang::expr_text(M))

      ec_candidate_bl <- rlang::parse_quos(c(ec_bl_terms_char, ec_candid), rlang::global_env())
      # Calculate the RMSE of the simplified model
      rmse_candid <- rmse_calc(.fm= candid_fm, .G=G, .E=E, .M=M, .trial=trial, .env_cv_df = env_cv_df,
                               .cores=ncores, .kn=kn, .ecs_in_bline_model=ec_candidate_bl)$Rmse

      # If the RMSE for the candidate model is better than the current model, then update the current model to be the candidate model
      if( rmse_candid < rmse_curr ){
        curr_fm <- candid_fm
        rmse_curr <- rmse_candid
        continue <- FALSE
      } else { # If the candidate EC model is not an improvement, remove that particular EC from the data frame containing the list of candidate ECs
        ec_summary <- ec_summary[-which_selected, ]
        # If there are no more candidate EC models to choose from, finish the while loop and return the initial model as the best model
        if(dim(ec_summary)[1]<=0){
          continue <- FALSE
        }
      }
    } else { # If the best full EC candidate model is worse then the current model, then there are no ECs that can further improve the model
      continue <- FALSE
    }
  }
  final_fm <- curr_fm
  return_list <- list(fm=final_fm, selected_ec=ec_candid, rmse=rmse_curr)
  return(return_list)
}




#' @title Include all important ECs into the model
#' @description
#' This function combines all of the component functions required to complete the environment covariate selection algorithm for a given baseline model.
#' The final output is an updated model containing each of the important ECs identified by the algorithm, along with the corresponding
#' environmental covariate terms that are statistically significant.
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
#' @param ncores The number of computer cores used during the cross validation scheme. Note that this should not be more than the number of cross validation groups in the model.
#' Also note that by setting this greater than 2 (i.e. the default) will require additional \code{ASReml-R} licenses.
#' @param kn The number of knot points for the spline component of the model for the environment covariate being tested.
#' Note that this is ignored of the environmental covariate is a factor
#' @param ecs_in_bline_model An optional list of quosures where each element of the list pertains to an environmental covariate that is in the initial baseline model.
#'
#' @return A list with 2 components:
#' \itemize{
#'  \item \code{fm}: An updated model with the all of the important environmental covariate included in the model along with only the significant
#'  fixed and random effects pertaining to each environmental covariate identified.
#'  \item \code{rmse}: The value of the RMSE for the final environmental covariate model
#'  }
#' @examples
#'
#' @export
ec_all <- function(fm, ECs, G, E, M, env_cv_df=NULL, ncores=2, kn=6, trial=NULL, ecs_in_bline_model=rlang::maybe_missing()){

  # Identify the data frame from the model
  .df <<- base::eval(fm$call$data)

  # Define the current model to be the initial model
  curr_fm <- fm

  # Put ECs in a quosure
  #quo_ecs <- rlang::enquos(ECs)

  #quo_ecs <- rlang::quos(PrePAW:PostPAW, PreCumRad)

  # Define a term that keeps track of whether to keep searching for ECs
  continue_ec_search <- TRUE

  # Identify the EC variables
  vars <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist to the ECs to be assessed
  cols_ecs <- cols_ecs <- unlist(purrr::map(ECs, rlang::eval_tidy, vars))
  # Use the columns to obtain a vector of EC terms
  ec_terms_char <- colnames(.df[cols_ecs])
  # Define a quosure for the ecs to select from
  ecs_to_select_from <- ECs


  if(rlang::is_missing(ecs_in_bline_model)==TRUE){
    ecs_in_bline_model <- rlang::quos(NULL)
    cols_bl_ecs <- c()
    ec_bl_terms_char <- c()
  } else {
    # # Identify which columns in the data frame consist of ECs in the baseline model
    # cols_ecs_bl <- unlist(purrr::map(quo_ecs_bl, rlang::eval_tidy, vars))
    # # Use the columns to obtain a vector of EC terms
    # ec_bl_terms_char <- colnames(.df[cols_ecs_bl])
    #ecs_in_bline_model <- enexpr(ecs_in_bline_model)

    cols_bl_ecs <- unlist(purrr::map(ecs_in_bline_model, rlang::eval_tidy, vars))
    ec_bl_terms_char <- colnames(.df[cols_bl_ecs])
  }

  # Do the same for ecs in current baseline model
  ecs_in_curr_bl_model <- ecs_in_bline_model


  # Run the algorithm for a single iteration to calculate rmse
  curr_rmse <- rmse_calc(.fm= curr_fm, .G=G, .E=E, .M=M, .trial=trial, .env_cv_df = env_cv_df,
                         .cores=ncores, .kn=kn, .ecs_in_bline_model=ecs_in_curr_bl_model)$Rmse


  while(continue_ec_search==TRUE){
    ec_iter <- ec_iteration(fm = curr_fm, ECs = ecs_to_select_from, G=G, E=E, M=M, env_cv_df=env_cv_df, ncores=ncores, kn=kn, trial=trial, ecs_in_bline_model=ecs_in_curr_bl_model)
    ec_candid_fm <- ec_iter$fm
    candid_rmse <- ec_iter$rmse
    # If the model has been has changed, set the current model to be the candidate model, otherwise finish
    if(candid_rmse < curr_rmse){
      curr_fm <- ec_candid_fm
      curr_rmse <- candid_rmse
      # Update the ec terms in the current model #CONTINUE HERE!
      col_candid_ec <- unlist(vars[ec_iter$selected_ec])
      # Update ECs in current (baseline?!) model
      cols_bl_ecs <- c(cols_bl_ecs, col_candid_ec)
      ecs_in_curr_bl_model <- rlang::parse_quos(names(vars[cols_bl_ecs]), env=rlang::global_env())
      # Remove selected EC from list of candidate ECs to choose from in the next iteration
      cols_ecs <- cols_ecs[cols_ecs!=col_candid_ec]
      # If statement to finish the model if there are no more ECs to choose from
      if(length(cols_ecs)==0){
        break
      }
      # Update the quosure of ECs to select from for the next iteration
      ecs_to_select_from <- rlang::parse_quos(names(vars[cols_ecs]), env=rlang::global_env())
    } else {
      continue_ec_search <- FALSE
    }
  }
  # Make final model equal to current model
  final_fm <- curr_fm
  final_rmse <- curr_rmse

  final_list <- list(fm=final_fm, rmse=final_rmse)
  return(final_list)
}












