# library(rlang)
# library(tidyverse)

# str(SorghumMET.df)
# # Create CV dataframe
# #df <- SorghumMET.df
# sorghum_env_cv_df <- data.frame(Env = unique(SorghumMET.df$Env))
# # Incorporate Trial information into the ec data frame
# sorghum_env_cv_df <- left_join(sorghum_env_cv_df, unique(SorghumMET.df[,c("Env","Trial")]), by="Env")
# # Now define cv_group
# sorghum_env_cv_df$cv_group <- as.integer(sorghum_env_cv_df$Trial)
#
#
# # .fm <- YieldInit.fm
# # ec <- expr(PrePAW)
# # G <- expr(Hybrid)
# # E <- expr(Env)
# # M <- expr(density)
# # trial <- expr(Trial)
# .env_cv_df <- sorghum_env_cv_df
#
# # Run with M as a factor
# .fm <- YieldMFactor.fm
# .M <- expr(TargetPop)
#
# # Run with M as continuous with 1 EC
# .fm <- YieldC1.fm
#
# # Define 2 ECs in the baseline model as a quosure
# .ecs_in_bline_model <- quos(PostPAW)
# #.fm <- YieldC1.fm
# .M <- expr(density)
# .ec <- expr(PrePAW)
# .G <- expr(Hybrid)
# .E <- expr(Env)
# .trial <- expr(Trial)
# .kn <- 6
# .cores <- 6
#
# # Test for NULL quosures
# x <- quos(NULL)
# y <- quo(PostPAW)
# quo_is_null(x[[1]])
# quo_is_missing(x[[1]])
# quo_is_null(y)
# quo_is_missing(y)

# Run the baseline model with established plant population fitted as a factor
# YieldMFactor.fm <- asreml( Yield..total...t.ha. ~ Hybrid + TargetPop + Hybrid:TargetPop,
#                         random =~ at(Trial):Rep  + at(Trial):MainPlot +
#                           at(Trial, c('Breeza Dry', 'Breeza Irr', 'Emerald')):SubPlot +
#                           at(Trial,'Breeza Dry'):Column +
#                           Trial + Env +
#                           Trial:Hybrid + Trial:Hybrid:TargetPop +
#                           Env:Hybrid + Env:Hybrid:TargetPop,
#                         residual=~ dsum(~units|ResidualTOS),
#                         data = SorghumMET.df,
#                         na.action=na.method(x='include'),
#                         maxit=30, workspace=1e8, extra=10)

#Define %dopar% and %do% locally
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`

ec_cv_full <- function(.fm, .ec, .G, .E, .M, .trial=NULL, .env_cv_df=NULL,
                       .cores=2, .kn=6, .ecs_in_bline_model=rlang::quos(NULL)){
  # Obtain the data frame from the baseline model
  .df  <- base::eval(.fm$call$data)
  # Now also round all continuous variables in df to 4 decimal places to avoid errors later on due to merging of data frames
  .df <- .df %>% purrr::modify_if(is.numeric, round, digits=4)


  # ADD AN ERROR MESSAGE IF AT() HAS LEVELS NUMBERED INSTEAD OF STATED!!!!!!!!!!!!!!!!!!!!!!!!

  # Set each of the inputs as expressions so that the user does not have to make them as expressions prior to input
  # .ec <- enexpr(.ec)
  # .G <- enexpr(.G)
  # .E <- enexpr(.E)
  # .M <- enexpr(.M)
  #.trial <- enexpr(.trial)

  # Identify the EC variables that are present in the baseline model
  vars_ec_bl <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist to the ECs that are in the basline model
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
  ## CONTINUE HERE WITH 0 BASELINE COVARIATES 08/05/2024 ----

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
            #str(subset_df)


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
                random_str_terms <- rlang::expr(asreml::str( ~ !!.trial:!!.G + !!.trial:!!.G:!!.M, ~corh(2):id(!!tg) ) +
                                           asreml::str(~ !!.E:!!.G + !!.E:!!.G:!!.M, ~corh(2):id(!!eg) ))
              } else {
                random_str_terms <- rlang::expr(asreml::str(~ !!.E:!!.G + !!.E:!!.G:!!.M, ~corh(2):id(!!eg) ))
              }
              if(is.double(.df[[.ec]])==TRUE){
                random_ec_terms <- rlang::expr(asreml::spl(!!.ec, k=!!.kn) + asreml::spl(!!.ec, k=!!.kn):!!.G +
                asreml::spl(!!.ec, k=!!.kn):!!.M + asreml::spl(!!.ec, k=!!.kn):!!.G:!!.M +
                !!.ec:asreml::spl(!!.M, k=!!.kn) + !!.ec:!!.G:asreml::spl(!!.M, k=!!.kn) +
                asreml::spl(!!.ec, k=!!.kn):asreml::spl(!!.M, k=!!.kn) +
                asreml::spl(!!.ec, k=!!.kn):!!.G:asreml::spl(!!.M, k=!!.kn))  #Remember to change str() and spatial effects as required

                subset_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_bl_terms + !!fixed_ec_terms,
                                                 random =~ !!random_bl_terms + !!random_ec_terms + !!random_str_terms,
                                                 residual=~ !!residual_bl_terms,
                                                 data=subset_df,
                                                 na.action=asreml::na.method(x='include'),
                                                 aom=T, maxit=300))
              } else if(is.double(.df[[.ec]])==FALSE) {
                  random_ec_terms <- rlang::expr(!!.ec:asreml::spl(!!.M, k=!!.kn) + !!.ec:!!.G:asreml::spl(!!.M, k=!!.kn))
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
                  random_ec_terms <- rlang::expr(asreml::spl(!!.ec, k=!!.kn) + asreml::spl(!!.ec, k=!!.kn):!!.G +
                                            asreml::spl(!!.ec, k=!!.kn):!!.M + asreml::spl(!!.ec, k=!!.kn):!!.G:!!.M)
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
          subset_pred <- asreml::predict(object=subset_fm,
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

# E <- unique(df$Env) %>% as.character()
# E <- expr(env)
# str(E)
# typeof(E)
#
#
#
# # Now test the function ----
# Rmse_df <- ec_cv_full(.fm=YieldInit.fm, .ec=expr(PrePAW), .G = expr(Hybrid), .E = expr(Env),
#                       .M = expr(density), .trial=quo(Trial), .env_cv_df=sorghum_env_cv_df,
#                       .cores=6, .kn=6)
#
# # Test with unspecified cv (come back to this later!!)
# Rmse_df <- ec_cv_full(.fm=YieldInit.fm, .ec=expr(PrePAW), .G = expr(Hybrid), .E = expr(Env),
#                       .M = expr(density), .trial=expr(Trial), .env_cv_df=sorghum_env_cv_df,
#                       .cores=6, .kn=6)
#
# # Test for when M is a factor
# Rmse_df <- ec_cv_full(.fm=YieldMFactor.fm, .ec=expr(PrePAW), .G = expr(Hybrid), .E = expr(Env),
#                       .M = expr(TargetPop), .trial=expr(Trial), .env_cv_df=sorghum_env_cv_df,
#                       .cores=6, .kn=6)
#
# # Update to include an EC
# Rmse_df <- ec_cv_full(.fm=YieldC1.fm, .ec=expr(PrePAW), .G = expr(Hybrid), .E = expr(Env),
#                       .M = expr(density), .trial=quo(Trial), .env_cv_df=sorghum_env_cv_df,
#                       .cores=6, .kn=6, .ecs_in_bline_model = quos(PostPAW))
#
# Rmse_df <- ec_cv_full(.fm=YieldC2.fm, .ec=expr(PrePAW), .G = expr(Hybrid), .E = expr(Env),
#                       .M = expr(density), .trial=expr(Trial), .env_cv_df=sorghum_env_cv_df,
#                       .cores=6, .kn=6, .ecs_in_bline_model = quos(c(PostPAW,ISW)))
#
# Rmse_df$Cor
# Rmse_df$Rmse

# Write a function to randomly assign each environment to one of the cross-validation groups
# Input E is a character vector of all unique environments
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

# envs <- unique(SorghumMET.df$Env) %>% as.character()
# cv_sorghum <- cv_groups(E=envs)
# cv_sorghum
#
#
# test <- function(G){
#   G <- expr(!!G)
#   return(G)
# }
# test(Genotype)


# Update the function to accept multiple ECs
ec_finder <- function(fm, ECs, G, E, M, env_cv_df=NULL,
                        cores=2, kn=6, trial=NULL, ecs_in_bline_model=NULL)
{
  df  <- base::eval(fm$call$data)

  # Turn all the relevant inputs into expressions (except for trial)
  G <- rlang::enexpr(G)
  E <- rlang::enexpr(E)
  M <- rlang::enexpr(M)
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
  quo_ecs <- rlang::enquos(ECs)
  quo_ecs_bl <- rlang::enquos(ecs_in_bline_model)

  # Send an error message if an EC appears both in the list of ECs to be searched as well as the list of ECs in the baseline model

  # quo_ecs <-  quos(PrePAW:ISW, PreFlwT:SeedSet)
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
  ec_selected <- ec_df$EC[which(ec_df$RMSE == min(ec_df$RMSE,na.rm=TRUE) )]
  return_ecs_list <- list(ec_selected=ec_selected, summary_ecs = ec_df)
  return(return_ecs_list)
}

# # Test the ec_finder function
# ec_search <- ec_finder(fm=YieldInit.fm, ECs= c(PrePAW, PreCumRad),
#                        G= Hybrid, E = Env,
#                        M= density, trial = expr(Trial),
#                        env_cv_df=sorghum_env_cv_df,
#                        cores=6, kn=6)
#
# ec_search$ec_selected
# ec_search$summary_ecs
#
# # Now test if you give it a different set of groups for cross-validation
# ec_search <- ec_finder(fm=YieldInit.fm, ECs= c(PrePAW, PreCumRad),
#                        G= Hybrid, E = Env,
#                        M= density, trial = expr(Trial),
#                        env_cv_df=sorghum_env_cv_df,
#                        cores=6, kn=6)
#
#
# # Now test when M is a factor
# ec_search <- ec_finder(fm=YieldMFactor.fm, ECs= c(PrePAW, ISW, PreCumRad),
#                        G= Hybrid, E = Env,
#                        M= TargetPop, trial = expr(Trial),
#                        env_cv_df=sorghum_env_cv_df,
#                        cores=6, kn=6)
#
# # Now test when there is an EC already in the model
# ec_search <- ec_finder(fm=YieldC1.fm, ECs= c(PrePAW, PreCumRad),
#                        G= Hybrid, E = Env,
#                        M= density, trial = expr(Trial),
#                        env_cv_df=sorghum_env_cv_df,
#                        cores=6, kn=6,
#                        ecs_in_bline_model=PostPAW)
#
# # Test when there are 2 ECs already in the model
# ec_search <- ec_finder(fm=YieldC2.fm, ECs= c(PrePAW, PreCumRad),
#                        G= Hybrid, E = Env,
#                        M= density, trial = expr(Trial),
#                        env_cv_df=sorghum_env_cv_df,
#                        cores=6, kn=6,
#                        ecs_in_bline_model=c(PostPAW,ISW))
#
# ec_search$ec_selected
# ec_search$summary_ecs
#
# # Test when there are 4 ECs already in the model
# # Now test when there is an EC already in the model
# ec_search <- ec_finder(fm=YieldC4.fm, ECs=c(PostFlwT:PreMaxT, PostMaxT),
#                        G= Hybrid, E = Env,
#                        M= density, trial = expr(Trial),
#                        env_cv_df=sorghum_env_cv_df,
#                        cores=6, kn=6,
#                        ecs_in_bline_model=c(PrePAW:ISW, PreCumRad))
#
#
#
#
# # An alternate version of 'select'
# select2 <- function(data, ...) {
#     dots <- enquos(...)
#
#     vars <- as.list(set_names(seq_along(data), names(data)))
#     cols <- unlist(map(dots, eval_tidy, vars))
#
#     data[, cols, drop = FALSE]
# }
#
# dfB <- data.frame(a = 1, b = 2, c = 3, d = 4, e = 5)
# select2(dfB, -d)
#
#
#
# subset2 <- function(data, rows) {
#   rows <- enquo(rows)
#   rows_val <- eval_tidy(rows, data)
#   stopifnot(is.logical(rows_val))
#
#   data[rows_val, , drop = FALSE]
# }
# sample_df <- data.frame(a = 1:5, b = 5:1, c = c(5, 3, 1, 4, 1))
# subset2(sample_df, b == c)
#
# # Now try with expressions
# subset3 <- function(data, rows) {
#   rows <- quo(rows)
#   rows_val <- eval_tidy(rows, data)
#   stopifnot(is.logical(rows_val))
#
#   data[rows_val, , drop = FALSE]
# }
# subset2(sample_df, b == c)


#
# .fm <- YieldC2_full.fm
# .ecs_in_model <- quos(PostPAW, ISW)
#
# # Create a function to identify the most parsimonious EC random effects model by removing non-significant EC terms
# ec_random_model <- function(.fm, .ecs_in_model){
#
#   # Obtain the data frame from the baseline model
#   .df  <- base::eval(.fm$call$data)
#
#   # Identify each of the EC terms and place into a character vector
#   ec_terms_char <-  unlist(map(.ecs_in_model, quo_text))
#   # Define the term that will go into grep
#   # Specifically we want to search for terms that have 'spl' AND an EC term in it via grep
#   ec_terms_for_grep <- paste0("spl.*(",
#                              paste(ec_terms_char, collapse="|"),
#                              ")")
#
#   # First identify the random terms corresponding to one of the ECs in the model
#   which_random_ec_terms <- grep(ec_terms_for_grep, .fm$factor.names)
#   random_ec_terms <-  .fm$factor.names[which_random_ec_terms]
#
#   # Place the EC terms into a data frame
#   ec_terms_df <- data.frame(Term=random_ec_terms)
#   #str(ec_terms_df)
#   # Determine the margin (i.e. order of hierachy) by how many times the strings ":" and "spl" appear in the term
#   ec_terms_df$Margin <-  str_count(ec_terms_df$Term, ":") + str_count(ec_terms_df$Term, "spl")
#
#   #ec_terms_df$Margin <- factor(ec_terms_df$Margin)
#
#   # Denote current model
#   varcomp_init_df <- summary(.fm)$varcomp
#   # Set current varcomps to be equal to initial model for now
#   varcomp_curr_df <-varcomp_init_df
#   # Set current model to be initial model for now
#   curr_fm <- .fm
#
#   #Calculate the current AIC value
#   AIC_curr <- icREML(list(curr_fm))$AIC
#
#
#   # Set removed terms to be empty by default
#   removed_terms <- rlang::maybe_missing()
#
#   # Do a for loop for each EC to drop terms from the model
#   for(i in 1:length(ec_terms_char)){
#     # Perform a grep for each EC to subset the terms
#     which_ith_ec_terms <- grep(ec_terms_char[i], ec_terms_df$Term)
#     # CONTINUE HERE!! 15/05/2024
#     temp_ec_terms_df <- ec_terms_df[ which_ith_ec_terms, ]
#     for(j in rev(c(1:max(temp_ec_terms_df$Margin))) ) {
#       subset_ec_terms_df <- temp_ec_terms_df[temp_ec_terms_df$Margin==j, ]
#       # Do we need to include k for each term within each margin??
#       for(k in 1:dim(subset_ec_terms_df)[1]){
#         term_to_test <- expr(!!subset_ec_terms_df$Term[k]) #%>% parse_expr()
#         #Subset current varcomp to current term
#
#         # If the current term being tested is boundary is the current model,
#         # Update curr_fm by removing the considered term
#         if(varcomp_curr_df$bound[rownames(varcomp_curr_df)==term_to_test] =="B"){
#           # Add currently test term to remove_terms
#           # If remove_terms is empty, then define remove_terms
#           if(is_missing(removed_terms){
#             removed_terms <- parse_expr(expr_text(term_to_test))
#           } else {
#             removed_terms <- parse_expr( paste(expr_text(term_to_test),
#                                                expr_text(removed_terms),
#                                                sep=" + ") )
#           }
#           #Now make term_tos_test an expression
#           term_to_test <- parse_expr(term_to_test)
#           # Update the current call to asreml
#           curr_call <- rlang::expr(asreml::update.asreml(curr_fm,
#                                                          fixed= . ~ . ,
#                                                          random =~ . - !!removed_terms,
#                                                          residual=~ .,
#                                                          data=.df,
#                                                          G.param=curr_fm$G.param, R.param=curr_fm$R.param, # to ensure the model does not get stuck in a local maxima
#                                                          aom=T, maxit=300))
#         } else {
#           # If remove_terms is NOT empty, then update curr_fm
#           if(is_missing(removed_terms)==FALSE{
#             curr_fm <-  eval(curr_call)
#             #Update the AIC value for the model (should be the same as before??)
#             AIC_curr <- icREML(list(curr_fm))$AIC
#           }
#           test_call <- rlang::expr(asreml::update.asreml(curr_fm,
#                                                          fixed= . ~ . ,
#                                                          random =~ . - !!term_to_test,
#                                                          residual=~ .,
#                                                          data=.df,
#                                                          G.param=curr_fm$G.param, R.param=curr_fm$R.param, # to ensure the model does not get stuck in a local maxima
#                                                          aom=T, maxit=300))
#
#           # Run the test model
#           test_fm <- eval(test_call)
#
#           #Calculate the AIC value for the model being considered
#           AIC_test <- icREML(list(test_fm))$AIC
#           # If the AIC of the test model is better, then set all the current model object to the test object
#           if(AIC_test <= AIC_curr){
#             curr_call <- test_call
#             curr_fm <- test_fm
#             AIC_curr <- AIC_test
#             varcomp_curr_df <- summary(curr_fm)$varcomp
#           } else {
#             # Add code here that ensures that relevant lower-order terms are no longer tested
#
#             # temp_ec_terms_df$Margin['relevant terms based on condition' &
#             #                         (temp_ec_terms_df$Margin < k)] <- NA
#             temp_ec_terms_df$Margin[(temp_ec_terms_df$Margin < j)] <- NA
#           }
#         }
#       }
#     }
#   }
# }






