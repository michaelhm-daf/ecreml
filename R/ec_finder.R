
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
#' library(agridat)
#' data("ortiz.tomato.yield")
#' # Subset to obtain a character vector of unique environments
#' unique_envs <- as.character(levels(ortiz.tomato.yield$env))
#' # Partition the environments into 6 groups/folds for k-folds cross validation
#' cv_groups(unique_envs, folds=6)
#'
#' @export
cv_groups <- function(E, folds = 6) {
  # Check if E is provided and is not NULL
  if (missing(E) || is.null(E)) {
    stop("Error: The argument 'E' must be provided and cannot be NULL.")
  }

  # Check if E is a valid input (e.g., a vector or data frame column)
  if (!is.vector(E) && !is.factor(E)) {
    stop("Error: The argument 'E' must be a vector or a factor.")
  }

  # Check if folds is a positive integer
  if (!is.numeric(folds) || folds <= 0 || folds != as.integer(folds)) {
    stop("Error: The argument 'folds' must be a positive integer.")
  }

  # Check if the number of folds is not greater than the number of elements in E
  if (length(E) < folds) {
    stop("Error: The number of folds cannot exceed the number of elements in 'E'.")
  }

  # Create the data frame
  env_cv_df <- data.frame(E = rlang::expr(!!E))

  # Check if the data frame was created successfully
  if (nrow(env_cv_df) == 0) {
    stop("Error: The input 'E' resulted in an empty data frame. Please check your input.")
  }

  n <- nrow(env_cv_df)

  # Generate the initial clustering groups for each environment
  env_cv_df$init_group <- rep(1:folds, times = ceiling(n / folds))[1:n]

  # Now randomly allocate environments to folds/groups
  env_cv_df$sample <- sample(1:n, replace = FALSE)
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
#' @param aliased Type logical with default set to FALSE. An option that enables aliased predictions from \code{predict.asreml} to be outputted.
#' We do not recommend setting this option to TRUE. As its inclusion in only for debugging purposes.
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
#' library(asreml)
#' library(foreach)
#' data(SorghumYield)
#' data(SorghumCvGroup)
#' # Run baseline model
#' baseline_asr <- asreml::asreml( Yld ~ Genotype + density + Genotype:density,
#'   random =~ at(Trial):Rep  + at(Trial):MainPlot +
#'    at(Trial,c('Breeza 1', 'Breeza 2', 'Emerald', 'Moree')):SubPlot +
#'    at(Trial,'Breeza 1'):Column +
#'    Trial + Env +
#'    spl(density, k=6) + spl(density, k=6):Genotype +
#'    str(~Trial:Genotype + Trial:Genotype:density,
#'        ~corh(2):id(48)) +
#'    str(~Env:Genotype + Env:Genotype:density,
#'        ~corh(2):id(136)),
#'   residual=~ dsum(~units|ResidualTOS),
#'   data = SorghumYield,
#'   na.action=na.method(x='include'),
#'   maxit=30, workspace="1Gb")
#'  ec_rmse_summary <- ec_single(.fm=baseline_asr, .ec=rlang::expr(PrePAW), .G =rlang::expr(Genotype), .E = rlang::expr(Env),
#'                    .M = rlang::expr(density), .trial=rlang::expr(Trial), .env_cv_df=SorghumCvGroup)
#'  ec_rmse_summary$Rmse
#'
#' @export
ec_single <- function(.fm, .ec, .G, .E, .M, .trial=NULL, .env_cv_df=NULL,
                       .cores=2, .kn=6, .ecs_in_bline_model=rlang::quos(NULL), aliased=FALSE){

  # Error handling for input arguments

  # Check if .fm is provided and is a valid asreml model object
  if (missing(.fm) || is.null(.fm)) {
    stop("Error: The argument '.fm' (baseline model) must be provided and cannot be NULL.")
  }
  if (!inherits(.fm, "asreml")) {
    stop("Error: The argument '.fm' must be a valid asreml model object.")
  }

  # Check if .ec is provided and is an expression
  if (missing(.ec) || is.null(.ec)) {
    stop("Error: The argument '.ec' (environmental covariate) must be provided and cannot be NULL.")
  }
  if (!rlang::is_expression(.ec)) {
    stop("Error: The argument '.ec' must be a valid R expression.")
  }

  # Check if .G, .E, and .M are provided and are expressions
  if (missing(.G) || is.null(.G)) {
    stop("Error: The argument '.G' (genotype term) must be provided and cannot be NULL.")
  }
  if (!rlang::is_expression(.G)) {
    stop("Error: The argument '.G' must be a valid R expression.")
  }

  if (missing(.E) || is.null(.E)) {
    stop("Error: The argument '.E' (environment term) must be provided and cannot be NULL.")
  }
  if (!rlang::is_expression(.E)) {
    stop("Error: The argument '.E' must be a valid R expression.")
  }

  if (!is.null(.M) && !rlang::is_expression(.M)) {
    stop("Error: The argument '.M' (management practice term) must be a valid R expression if provided.")
  }

  # Check if .trial is an expression if provided
  if (!is.null(.trial) && !rlang::is_expression(.trial)) {
    stop("Error: The argument '.trial' must be a valid R expression if provided.")
  }

  # Check if .env_cv_df is a data frame if provided
  if (!is.null(.env_cv_df) && !is.data.frame(.env_cv_df)) {
    stop("Error: The argument '.env_cv_df' must be a data frame if provided.")
  }

  # Check if .cores is a positive integer
  if (!is.numeric(.cores) || .cores <= 0 || .cores != as.integer(.cores)) {
    stop("Error: The argument '.cores' must be a positive integer.")
  }

  # Check if .kn is a positive integer
  if (!is.numeric(.kn) || .kn <= 0 || .kn != as.integer(.kn)) {
    stop("Error: The argument '.kn' must be a positive integer.")
  }

  # Check if .ecs_in_bline_model is a list of quosures
  if (!is.null(.ecs_in_bline_model) && !rlang::is_quosures(.ecs_in_bline_model)) {
    stop("Error: The argument '.ecs_in_bline_model' must be a list of quosures if provided.")
  }

  # Check if aliased is a logical value
  if (!is.logical(aliased)) {
    stop("Error: The argument 'aliased' must be a logical value (TRUE or FALSE).")
  }

  # Obtain the data frame from the baseline model
  .df <- base::eval(.fm$call$data)

  # Check if the data frame is valid
  if (is.null(.df) || !is.data.frame(.df)) {
    stop("Error: The data frame used in the baseline model could not be retrieved or is not valid.")
  }

  # Round all continuous variables in the data frame to 4 decimal places
  .df <- .df %>% purrr::modify_if(is.numeric, round, digits = 4)

  # Remove cv_group from .df if it exists to prevent bugs
  if ("cv_group" %in% colnames(.df)) {
    .df <- .df %>% dplyr::select(!"cv_group") %>% as.data.frame()
  }

  # Check if .env_cv_df contains the environment term
  if (!is.null(.env_cv_df) && is.null(.env_cv_df[[rlang::as_string(.E)]])) {
    stop("Error: The environment term is missing from the provided environment grouping for cross-validation.")
  }

  # Ensure the number of cores does not exceed the number of cross-validation groups
  if (!is.null(.env_cv_df) && .cores > length(unique(.env_cv_df$cv_group))) {
    stop("Error: The number of cores cannot exceed the number of cross-validation groups.")
  }

  # Ensure the number of cores does not exceed the available system cores
  total_cores <- parallel::detectCores()
  if (.cores > total_cores) {
    stop(paste("Error: The number of cores specified (", .cores, ") exceeds the available system cores (", total_cores, ").", sep = ""))
  }


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
  } else {
    baseline_ec_cols <- c()
  }

  # Obtain response variable from the data frame as an expression
  response_term <- attr(.fm$formulae$fixed, "variables")[[2]] %>%
    as.character() %>%
    rlang::parse_expr()

  # Merge the data frame with cross validation groupings for environments to the phenotype data
  # First if statement to determine if environment groupings have been provided as an input
  # If not provided, then generate them randomly
  if(is.null(.env_cv_df)==TRUE){
    tryCatch({
      # Generate ec_cv dataframe
      E_char <- unique(.df[[.E]]) %>% as.character()
      # Now run ec cross-validation function
      .env_cv_df <- cv_groups(E = E_char)
      # Define the column heading for E to be same as it is in the input into ec_single
      .env_cv_df[[.E]] <- .env_cv_df$E
    }, error = function(e) {
      stop("Error during cross-validation group generation: ", e$message)
    })
  }

  # If there is no column heading with the same name as E in the original data-frame, then create one
  # Show an error message describing so
  if(is.null(.env_cv_df[[.E]])==TRUE){
    if(is.null(.env_cv_df$.E)==TRUE){
      stop("The environment term is missing from the provided environment grouping for cross validation")
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

  # Used the cross-validation data frame to determine the total number of folds
  v <- length(levels(.df$cv_group))

  # Define number of cores in parallel
  total_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(min(c(total_cores[1]-1, .cores)))
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, "parallel_predict_list", envir=) # Should only be required when testing, nope still needed for the package it seems

  # Use foreach to perform parallel programming when implementing the cross validation scheme
  cv_pred <- foreach::foreach(i=c(1:v), .combine=rbind,
                              .packages=c("asreml", "tidyverse","rlang")) %dopar% {

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
                                #factor_names <- .fm$factor.names
                                factor_names <- names(.fm$noeff) # Need this for asreml 4.2
                                # The bottom part needed for the latest version of asreml (UPDATE IN THE FUTURE)
                                if(is.null(factor_names)==TRUE){
                                  factor_names <- c(attr(.fm$formulae$fixed, "term.labels"),
                                                    attr(.fm$formulae$random, "term.labels"))#,
                                                   # attr(.fm$formulae$residual, "term.labels"))
                                }
                                which_continuous <- grep(":|at\\(|mv|\\(Intercept\\)|spl", factor_names)
                                # Also remove the M term from list of factors if M is continuous
                                if(is.null(.M)==FALSE){
                                  if(is.double(subset_df[[.M]])==TRUE){
                                    #remove_cont_M <- grep(rlang::expr_text(.M), .fm$factor.names) #factor.names is in older version of R
                                    remove_cont_M <- grep(rlang::expr_text(.M), factor_names)
                                    which_continuous <- unique(c(which_continuous,remove_cont_M))
                                  }
                                }
                                # Also remove ECs in the model that are continuous
                                # Note: May need to add for loop if multiple ECs are included in the model
                                if(!rlang::quo_is_null(.ecs_in_bline_model[[1]])==TRUE){
                                  for(i in baseline_ec_cols){
                                    bl_ec <- subset_df[,i]
                                    if(is.double(subset_df[[i]])==TRUE){
                                      remove_cont_ec <- grep(colnames(subset_df)[i] , factor_names)
                                      which_continuous <- unique(c(which_continuous, remove_cont_ec))
                                    }
                                  }
                                }


                                # Probably good to check that the output of factor_terms is as expected
                                factor_terms <- factor_names[-which_continuous]

                                # Re-factorise all the terms that are factors in the model for the subsetted data-frame
                                #Using mutate() to convert specific columns to factors
                                # Note the different between 'factor' and 'as.factor'
                                # Namely, 'as.factor()' does not update the number levels of each term
                                # DEPRECATED NEEDS TO BE UPDATED
                                subset_df <- subset_df %>% dplyr::mutate(dplyr::across(factor_terms, factor))


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
                                  # Define an expression for the fixed and random effect EC terms to be added in the updated model
                                  fixed_ec_terms <- rlang::expr(!!.ec + !!.ec:!!.G)
                                  # if(is.null(.trial)==FALSE){               #This should be redundant because the random_str_terms should already appear in the random_bl_terms
                                  #   # Define the random regression terms
                                  #   random_str_terms <- rlang::expr(!!.trial:!!.G + !!.E:!!.G)
                                  # } else {
                                  #   random_str_terms <- rlang::expr(!!.E:!!.G)
                                  #}
                                  if(is.double(.df[[.ec]])==TRUE){
                                    random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G)

                                    subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                              random =~ !!random_bl_terms + !!random_ec_terms,# + !!random_str_terms,
                                                                              residual=~ !!residual_bl_terms,
                                                                              data=subset_df,
                                                                              na.action=asreml::na.method(x='include'),
                                                                              aom=T, maxit=300))
                                  } else if(is.double(.df[[.ec]])==FALSE) {
                                    subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                              random =~ !!random_bl_terms,# + !!random_str_terms,
                                                                              residual=~ !!residual_bl_terms,
                                                                              data=subset_df,
                                                                              na.action=asreml::na.method(x='include'),
                                                                              aom=T, maxit=300))
                                  } else {
                                    stop("Need to specify whether the environmental covariate is a continuous or categorical variable")
                                  }
                                } else {
                                  # Change the model depending on whether M is a factor or a variate (i.e. continuous) or missing
                                  if(is.double(.df[[.M]])==TRUE){
                                    # Define an expression for the fixed and random effect EC terms to be added in the updated model
                                    fixed_ec_terms <- rlang::expr(!!.ec + !!.ec:!!.G + !!.ec:!!.M + !!.ec:!!.G:!!.M)
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
                                      subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                                random =~ !!random_bl_terms + !!random_ec_terms + !!random_str_terms,
                                                                                residual=~ !!residual_bl_terms,
                                                                                data=subset_df,
                                                                                na.action=asreml::na.method(x='include'),
                                                                                aom=T, maxit=300))
                                    } else {
                                      stop("Need to specify whether the environmental covariate is a continuous or categorical variable")
                                    }
                                  } else if(is.factor(.df[[.M]])==TRUE) {
                                    # Define an expression for the fixed and random effect EC terms to be added in the updated model
                                    fixed_ec_terms <- rlang::expr(!!.ec + !!.ec:!!.G + !!.ec:!!.M + !!.ec:!!.G:!!.M)
                                    if(is.null(.trial)==FALSE){
                                      # Define the random regression terms
                                      random_str_terms <- rlang::expr(!!.trial:!!.G + !!.trial:!!.G:!!.M + !!.E:!!.G + !!.E:!!.G:!!.M)
                                    } else {
                                      random_str_terms <- rlang::expr(!!.E:!!.G + !!.E:!!.G:!!.M)
                                    }
                                      if(is.double(.df[[.ec]])==TRUE){
                                        random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G +
                                                                       spl(!!.ec, k=!!.kn):!!.M + spl(!!.ec, k=!!.kn):!!.G:!!.M)
                                        subset_call <- rlang::expr(asreml::asreml(fixed= !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                        random =~ !!random_bl_terms + !!random_ec_terms + !!random_str_terms,
                                                                        residual=~ !!residual_bl_terms,
                                                                        data=subset_df,
                                                                        #G.param=list(), R.param=list(), # to ensure the model does not get stuck in a local maxima
                                                                        aom=T, maxit=300))
                                      } else if(is.double(.df[[.ec]])==FALSE){
                                        # Generate the call to the updated asreml model
                                        subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                                                  random =~ !!random_bl_terms  + !!random_str_terms,
                                                                                  residual=~ !!residual_bl_terms,
                                                                                  data=subset_df,
                                                                                  na.action=asreml::na.method(x='include'),
                                                                                  aom=T, maxit=300))
                                      # evaluate the updated asreml model
                                      } else {
                                      stop("Need to specify whether the management practice is a categorical or numerical variable")
                                      }
                                    }
                                  }
                                # Run the asreml model
                                tryCatch({
                                  subset_fm <- eval(subset_call)
                                }, error = function(e) {
                                  stop("Error during asreml model evaluation: ", e$message)
                                })

                                # Generate table of predictions
                                tryCatch({
                                  predict_info <- parallel_predict_list(.df = .df, subset_df = subset_df, .ec = .ec, .G = .G, .E = .E,
                                                                        .M = .M, baseline_ec_cols = baseline_ec_cols)
                                }, error = function(e) {
                                  stop("Error during prediction list generation: ", e$message)
                                })

                                # Obtain predictions for all environments (tested & untested) for subsetted model
                                tryCatch({
                                  if (length(predict_info$levels_list) > 1) {
                                    subset_pred <- asreml::predict.asreml(object = subset_fm,
                                                                          classify = predict_info$classify_terms,
                                                                          levels = predict_info$levels_list,
                                                                          parallel = TRUE, maxit = 1, pworkspace = "1gb")
                                  } else {
                                    subset_pred <- asreml::predict.asreml(object = subset_fm,
                                                                          classify = predict_info$classify_terms,
                                                                          levels = predict_info$levels_list,
                                                                          parallel = FALSE, maxit = 1, pworkspace = "1gb")
                                  }
                                }, error = function(e) {
                                  stop("Error obtaining cross-validation predictions using predict.asreml: ", e$message)
                                })
                                # Incorporate environment information into the table of predictions
                                temp_subset.pred <- dplyr::left_join(subset_pred$pvals, predict_info$aux_parallel,
                                                                     by=colnames(predict_info$aux_parallel)[-1]) # Assumes the 1st column is always .E

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
    cv_df <- dplyr::select(.df, .G, .E, .ec, response_term) %>%
      dplyr::left_join(cv_pred, . ,
                       by=c(rlang::expr_text(.G), rlang::expr_text(.E), rlang::expr_text(.ec) )) %>%
      unique()
  } else {
    cv_df <- dplyr::select(.df, .G, .M, .E, .ec, response_term) %>%
      dplyr::left_join(cv_pred, . ,
                       by=c(rlang::expr_text(.G), rlang::expr_text(.M), rlang::expr_text(.E), rlang::expr_text(.ec) )) %>%
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
#'  \item \code{ec_selected}: The environmental covariate selected as minimising the RMSE (a scalar of type character).
#'  compared with the predictions obtained with the environment covariate included in an untested environment.
#'  The predictions are obtained internally using \code{asreml::predict.asreml}.
#'  \item \code{summary_ecs}: A data frame with the RMSE and Pearson correlation for each environmental covariate considered during the forward selection procedure.
#'  }
#' @examples
#' library(asreml)
#' library(foreach)
#' data(SorghumYield)
#' data(SorghumCvGroup)
#' # Run baseline model
#' baseline_asr <- asreml::asreml( Yld ~ Genotype + density + Genotype:density,
#'   random =~ at(Trial):Rep  + at(Trial):MainPlot +
#'    at(Trial,c('Breeza 1', 'Breeza 2', 'Emerald', 'Moree')):SubPlot +
#'    at(Trial,'Breeza 1'):Column +
#'    Trial + Env +
#'    spl(density, k=6) + spl(density, k=6):Genotype +
#'    str(~Trial:Genotype + Trial:Genotype:density,
#'        ~corh(2):id(48)) +
#'    str(~Env:Genotype + Env:Genotype:density,
#'        ~corh(2):id(136)),
#'   residual=~ dsum(~units|ResidualTOS),
#'   data = SorghumYield,
#'   na.action=na.method(x='include'),
#'   maxit=30, workspace="1Gb")
#' ec_search <- ec_finder(fm=baseline_asr, ECs= rlang::quos(PrePAW, PostPAW),
#'     G= rlang::expr(Genotype), E = rlang::expr(Env),
#'     M= rlang::expr(density), trial = rlang::expr(Trial),
#'      env_cv_df=SorghumCvGroup)
#' ec_search$ec_selected
#' ec_search$summary_ecs
#' @export
ec_finder <- function(fm, ECs, G, E, M, env_cv_df=NULL,
                      cores=2, kn=6, trial=NULL, ecs_in_bline_model=rlang::quos(NULL))
{
  # Error handling for input arguments

  # Check if fm is provided and is a valid asreml model object
  if (missing(fm) || is.null(fm)) {
    stop("Error: The argument 'fm' (baseline model) must be provided and cannot be NULL.")
  }
  if (!inherits(fm, "asreml")) {
    stop("Error: The argument 'fm' must be a valid asreml model object.")
  }

  # Check if ECs is provided and is a list of quosures
  if (missing(ECs) || is.null(ECs)) {
    stop("Error: The argument 'ECs' (environmental covariates) must be provided and cannot be NULL.")
  }
  if (!rlang::is_quosures(ECs)) {
    stop("Error: The argument 'ECs' must be a list of quosures.")
  }

  # Check if G, E, and M are provided and are expressions
  if (missing(G) || is.null(G)) {
    stop("Error: The argument 'G' (genotype term) must be provided and cannot be NULL.")
  }
  if (!rlang::is_expression(G)) {
    stop("Error: The argument 'G' must be a valid R expression.")
  }

  if (missing(E) || is.null(E)) {
    stop("Error: The argument 'E' (environment term) must be provided and cannot be NULL.")
  }
  if (!rlang::is_expression(E)) {
    stop("Error: The argument 'E' must be a valid R expression.")
  }

  if (missing(M) || is.null(M)) {
    stop("Error: The argument 'M' (management practice term) must be provided and cannot be NULL.")
  }
  if (!rlang::is_expression(M)) {
    stop("Error: The argument 'M' must be a valid R expression.")
  }

  # Check if trial is an expression if provided
  if (!is.null(trial) && !rlang::is_expression(trial)) {
    stop("Error: The argument 'trial' must be a valid R expression if provided.")
  }

  # Check if env_cv_df is a data frame if provided
  if (!is.null(env_cv_df) && !is.data.frame(env_cv_df)) {
    stop("Error: The argument 'env_cv_df' must be a data frame if provided.")
  }

  # Check if cores is a positive integer
  if (!is.numeric(cores) || cores <= 0 || cores != as.integer(cores)) {
    stop("Error: The argument 'cores' must be a positive integer.")
  }

  # Check if kn is a positive integer
  if (!is.numeric(kn) || kn <= 0 || kn != as.integer(kn)) {
    stop("Error: The argument 'kn' must be a positive integer.")
  }

  # Check if ecs_in_bline_model is a list of quosures
  if (!is.null(ecs_in_bline_model) && !rlang::is_quosures(ecs_in_bline_model)) {
    stop("Error: The argument 'ecs_in_bline_model' must be a list of quosures if provided.")
  }

  # Check if any environmental covariates in ECs are also in ecs_in_bline_model
  if (!is.null(ecs_in_bline_model)) {
    ecs_in_bline_model_names <- purrr::map_chr(ecs_in_bline_model, rlang::as_label)
    ECs_names <- purrr::map_chr(ECs, rlang::as_label)
    common_ecs <- base::intersect(ECs_names, ecs_in_bline_model_names)
    if (length(common_ecs) > 0) {
      stop(paste("Error: The following environmental covariates appear in both 'ECs' and 'ecs_in_bline_model':",
                 paste(common_ecs, collapse = ", ")))
    }
  }

  df  <- base::eval(fm$call$data)

  # Turn all the relevant inputs into expressions (except for trial)
  # G <- rlang::enexpr(G)
  # E <- rlang::enexpr(E)
  # M <- rlang::enexpr(M)
  # Only make ecs in baseline model an expression if it is not empty
  # if(is.null(ecs_in_bline_model)==FALSE){
  #   ecs_in_bline_model <- enexpr(ecs_in_bline_model)
  # }

  # Check if the data frame is valid
  if (is.null(df) || !is.data.frame(df)) {
    stop("Error: The data frame used in the baseline model could not be retrieved or is not valid.")
  }

  # Check if the columns corresponding to ECs exist in the data frame
  vars <- as.list(magrittr::set_names(seq_along(df), names(df)))
  cols <- unlist(purrr::map(ECs, rlang::eval_tidy, vars))
  if (any(is.na(cols))) {
    stop("Error: One or more environmental covariates in 'ECs' do not exist in the data frame used in the baseline model.")
  }

  # Ensure the number of cores does not exceed the number of cross-validation groups
  if (!is.null(env_cv_df) && cores > length(unique(env_cv_df$cv_group))) {
    stop("Error: The number of cores cannot exceed the number of cross-validation groups.")
  }

  # Ensure the number of cores does not exceed the available system cores
  total_cores <- parallel::detectCores()
  if (cores > total_cores) {
    stop(paste("Error: The number of cores specified (", cores, ") exceeds the available system cores (", total_cores, ").", sep = ""))
  }



  # Don't need to if ECs are already a quosure
  quo_ecs <- ECs
  quo_ecs_bl <- ecs_in_bline_model
  # Send an error message if an EC appears both in the list of ECs to be searched as well as the list of ECs in the baseline model

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
    # ec_cv_temp <- ec_single(fm=YieldInit.fm, ec=ec_curr, G= expr(Hybrid), E = expr(Env),
    #            M = expr(density), trial=expr(Trial), env_cv_df=sorghum_env_cv_df,
    #            cores = 6, kn=6)
    tryCatch({
      ec_cv_temp <- ec_single(.fm = fm, .ec = ec_curr, .G = G, .E = E,
                              .M = M, .trial = trial, .env_cv_df = env_cv_df,
                              .cores = cores, .kn = kn, .ecs_in_bline_model = quo_ecs_bl)
    }, error = function(e) {
      stop(paste("Error during cross-validation for environmental covariate '", ec_curr, "': ", e$message, sep = ""))
    })

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











