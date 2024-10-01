
# Algorithm to obtain final random model
# Run most complex model
# Use while loop with a logical value that determine whether to continue or not with the algorithm
# If terms in highest margin value are boundary, then drop them from the model and re-run without them
# If the terms in the highest margin value are still non-boundary, set k = total non-boundary max margin terms
# For each k, perform AIC test, if non significant, drop the term and re-run for other k
# For final model, re-evaluate Margin, if other terms are now equal to the max margin value, then assess them
# If the max margin value cannot be reduced further, then stop and finish the model


#' @title Identify the most parsimonious EC random effects model
#' @description
#' Performs the backwards selection procedure for the random component of the model to obtain a parsimonious random effects model 
#' for the all environmental covariates currently in the model. 
#' Testing for the significance of random effect terms is performed using the AIC criteria at the REML parameter estimates (Verbyla 2019). 
#' 
#' @param .fm The current \code{asreml} model object that random environmental covariate terms will be tested from.  
#' @param .ecs_in_model A list of expressions such that each expression is an environmental covariate that is present in the current model.
#' @param .G The genotype term in the model as an expression
#' @param .M The management practice term in the model as an expression
#' 
#' @return An updated asreml model such that the non-significant environmental covariate spline terms have been removed from the model.
#' 
#' @examples
#' 
#' @export
ec_random_model <- function(.fm, .ecs_in_model, .G, .M){
  
  # Obtain the data frame from the baseline model
  .df  <<- base::eval(.fm$call$data)
  
  # Identify each of the EC terms and place into a character vector
  #quo_ecs_in_model <- enquos(.ecs_in_model)
  quo_ecs_in_model <- .ecs_in_model
  
  # Identify the variables in the data frame
  vars <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist to the ECs to be assessed
  cols_ecs <- unlist(purrr::map(quo_ecs_in_model, rlang::eval_tidy, vars))
  # Use the columns to obtain a vector of EC terms
  ec_terms_char <- colnames(.df[cols_ecs])
  
  
  
  #ec_terms_char <-  purrr::map_chr(.ecs_in_model, rlang::quo_text)
  
  
  # Identify which columns in the data frame consist to the ECs to be assessed
  #cols <- unlist(purrr::map(.ecs_in_model, rlang::eval_tidy, vars))
  
  
  # Define the term that will go into grep
  # Specifically we want to search for terms that have 'spl' AND an EC term in it via grep
  ec_terms_for_grep <- paste0("spl.*(",
                              paste(ec_terms_char, collapse="|"),
                              ")")
  
  # First identify the random terms corresponding to one of the ECs in the model
  which_random_ec_terms <- grep(ec_terms_for_grep, .fm$factor.names)
  random_ec_terms <-  .fm$factor.names[which_random_ec_terms]
  
  # Place the EC terms into a data frame
  ec_terms_df <- data.frame(Term=random_ec_terms)
  #str(ec_terms_df)
  
  # Determine the margin (i.e. order of hierachy) by how many times the strings ":" and "spl" appear in the term
  #ec_terms_df$Margin <-  str_count(ec_terms_df$Term, ":") + str_count(ec_terms_df$Term, "spl")
  
  #ec_terms_df$Margin <- factor(ec_terms_df$Margin)
  
  # Denote current model
  varcomp_init_df <- summary(.fm)$varcomp
  # Set current varcomps to be equal to initial model for now
  varcomp_curr_df <-varcomp_init_df
  # Set current model to be initial model
  curr_fm <- .fm
  
  #Calculate the current AIC value
  AIC_curr <- icREML(list(curr_fm))$AIC
  
  # Set removed terms to be empty by default
  removed_terms <- rlang::maybe_missing()
  
  # Define margin term in data frame and set to 0 by default
  ec_terms_df$Margin <- 0
  
  # Define an expression for the fixed and random effect EC terms to be added in the updated model
  fixed_terms_curr <- curr_fm$call$fixed
  residual_terms_curr <- curr_fm$call$residual
  random_terms_curr <- curr_fm$call$random
  
  curr_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
                                          random =  !!random_terms_curr,
                                          residual= !!residual_terms_curr,
                                          data= .df,
                                          na.action = na.method(x = "include"),
                                          aom=T, maxit=30))
  # set test_call equal to curr_call initially
  test_call <- curr_call
  
  # Do a for loop for each EC to drop terms from the model
  for(i in 1:length(ec_terms_char)){
    # Perform a grep for each EC to subset the terms
    which_ith_ec_terms <- grep(ec_terms_char[i], ec_terms_df$Term)
    temp_ec_terms_df <- (ec_terms_df[ which_ith_ec_terms, ])
    # Define the ith ec term separately because I am getting an error otherwise for some reason
    # Calculate the margin for each term
    # Need to define ec as a character to resolve error message for some reason
    temp_ec_terms_df$Margin <- margin(terms=temp_ec_terms_df$Term,
                                      ec=as.character(ec_terms_char[i]), G=.G, M=.M)
    
    # Set j equal to the current max value of margin
    j <- max(temp_ec_terms_df$Margin)
    #Create a 2nd version of temp_ec_terms
    temp_ec_terms_df2 <- temp_ec_terms_df
    continue_model_search <- TRUE
    # Now run models using a while-loop
    while(continue_model_search==TRUE){
      subset_ec_terms_df <- temp_ec_terms_df[temp_ec_terms_df$Margin==j, ]
      # Do we need to include k for each term within each margin??
      for(k in 1:dim(subset_ec_terms_df)[1]){
        term_to_test <- rlang::expr(!!subset_ec_terms_df$Term[k]) #%>% parse_expr()
        #Subset current varcomp to current term
        
        # If the current term being tested is boundary is the current model,
        # Update curr_fm by removing the considered term
        if(varcomp_curr_df$bound[rownames(varcomp_curr_df)==term_to_test] =="B"){
          # Add currently test term to remove_terms
          # If remove_terms is empty, then define remove_terms
          if(rlang::is_missing(removed_terms)){
            removed_terms <- as.character(term_to_test) #parse_expr(expr_text(term_to_test))
          } else {
            # removed_terms <- parse_expr( paste(expr(!!term_to_test),
            #                                  expr_text(removed_terms), # Note that removed_terms is already an expression, whilst term_to_test is a character
            #                                  sep=" + "))
            
            # Make the removed terms a character vector
            removed_terms <- c(removed_terms,   term_to_test )
          }
          #Now make term_to_test an expression
          #term_to_test <- parse_expr(term_to_test)
          # Update the current call to asreml, removing the boundary term
          # change asreml.options
          #asreml.options(extra=2, random.order="user")
          
          random_terms_curr <- curr_fm$call$random
          
          # Use self define function to remove random terms from the model for testing
          random_terms_curr <- subtract_terms(main_expr = random_terms_curr,
                                              removed_char_vec = removed_terms,
                                              response=FALSE)
          
          curr_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
                                                  random =  !!random_terms_curr,
                                                  residual= !!residual_terms_curr,
                                                  data= .df,
                                                  #keep.order = TRUE,
                                                  na.action = na.method(x = "include"),
                                                  #G.param=curr_fm$G.param,
                                                  #R.param=curr_fm$R.param, # to ensure the model does not get stuck in a local maxima
                                                  aom=T, maxit=30))
          #random.order="R"))
          
          # Remove the non-significant term from the terms data frame
          temp_ec_terms_df <- temp_ec_terms_df[!temp_ec_terms_df$Term%in%rlang::expr(!!term_to_test),]
          next # Start the next (kth) iteration of the for-loop
        } else { # Now for the alternate where the term being consider is not boundary
          # Only rerun asreml if the current model has changed since the last time asreml was called (i.e. boundary terms have since been removed from the current model)
          if(curr_call!=test_call){
            print('Dropping a RANDOM effects term from the model')
            curr_fm <-  eval(curr_call) # Note that the order of terms in the random terms is changed based on the order of the teams in the fixed effects??
            #Update the AIC value for the model (should be the same as before??)
            AIC_curr <- icREML(list(curr_fm))$AIC
            # Update varcomp table
            varcomp_curr_df <- summary(curr_fm)$varcomp
          }
          #Now make term_to_test an expression
          #term_to_test <- expr_text(parse_expr(term_to_test))
          # Use self define function to remove random terms from the model for testing
          random_terms_test <- subtract_terms(main_expr = random_terms_curr,
                                              removed_char_vec = term_to_test)
          
          test_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
                                                  random =  !!random_terms_test,
                                                  residual= !!residual_terms_curr,
                                                  data= .df,
                                                  na.action = na.method(x = "include"),
                                                  aom=T, maxit=30))
          
          # Run the test model
          test_fm <- eval(test_call)
          
          #Calculate the AIC value for the model being considered
          AIC_test <- icREML(list(test_fm))$AIC
          # If the AIC of the test model is better, then set all the current model object to the test object
          if(AIC_test <= AIC_curr){
            curr_call <- test_call
            curr_fm <- test_fm
            AIC_curr <- AIC_test
            varcomp_curr_df <- summary(curr_fm)$varcomp
            random_terms_curr <- random_terms_test
            # Remove the non-significant term from the terms data frame
            temp_ec_terms_df <- temp_ec_terms_df[!temp_ec_terms_df$Term%in%rlang::expr(!!term_to_test),]
          } else{
            # If there is no significant difference, set to test_call equal to current call to avoid unnecessary calls to asreml
            test_call <- curr_call
          }
        }
      }
      # If the number of rows in the 2nd version differs from the original, then re-run the loop
      if(length(temp_ec_terms_df2$Term)!=length(temp_ec_terms_df$Term)){
        
        # Update the values of Margin
        temp_ec_terms_df$Margin <- margin(terms= temp_ec_terms_df$Term,
                                          ec=as.character(ec_terms_char[i]), G=.G, M=.M)
        # Set j equal to new_j and continue testing for lower-order interaction terms
        j <- max(temp_ec_terms_df$Margin)
        # Set temp2 equal to temp 1
        temp_ec_terms_df2 <- temp_ec_terms_df
      } else {
        #Define new version of temp
        temp_ec_terms_new_df <- temp_ec_terms_df
        #Update the values of margin
        temp_ec_terms_new_df$Margin <- margin(terms= temp_ec_terms_new_df$Term,
                                              ec=as.character(ec_terms_char[i]), G=.G, M=.M)
        new_j <- max(temp_ec_terms_df$Margin)
        # If the value of Margin is the same for all model terms after updating, then the model cannot be further improved
        if(sum((temp_ec_terms_new_df$Margin==temp_ec_terms_df$Margin)==FALSE)==0){
          # define the final model
          #final_fm <- curr_fm
          # End the while loop
          continue_model_search <- FALSE
        } else {
          # Set temp_terms equal to new_temp_terms
          temp_ec_terms_df <- temp_ec_terms_new_df
          # Else set j equal to new_j and continue testing for lower-order interaction terms
          j <- max(temp_ec_terms_df$Margin)
          # Set temp2 equal to temp 1
          temp_ec_terms_df2 <- temp_ec_terms_df
        }
      }
    }
  }
  final_fm <- curr_fm
  return(final_fm)
}




#' @title Identify the most parsimonious EC fixed effects model
#' @description
#' Performs the backwards selection procedure for the fixed component of the model to obtain a parsimonious fixed effects model,
#' removing all non-significant environmental covariate terms. 
#' Testing for the significance of fixed effect terms is performed using Wald tests with an approximate F-statistic (Kenward & Roger 1997). 
#' 
#' @param .fm The current \code{asreml} model object that the environmental covariate terms will be tested from.  
#' @param .ecs_in_model A list of expressions such that each expression is an environmental covariate that is present in the current model.
#' @param .G The genotype term in the model as an expression
#' @param .M The management practice term in the model as an expression
#' 
#' @return An updated asreml model such that the non-significant environmental covariate fixed effect terms have been removed from the model.
#' 
#' @examples
#' 
#' @export
ec_fixed_model <- function(.fm, .ecs_in_model, .G, .M){
  
  # Obtain the data frame from the baseline model
  # Note, use super assignment to modify .df in the global environment
  .df  <<- base::eval(.fm$call$data)
  
  # Identify each of the EC terms and place into a character vector
  #quo_ecs_in_model <- enquos(.ecs_in_model)
  quo_ecs_in_model <- .ecs_in_model
  
  # Identify the variables in the data frame
  vars <- as.list(magrittr::set_names(seq_along(.df), names(.df)))
  # Identify which columns in the data frame consist to the ECs to be assessed
  cols_ecs <- unlist(purrr::map(quo_ecs_in_model, rlang::eval_tidy, vars))
  # Use the columns to obtain a vector of EC terms
  ec_terms_char <- colnames(.df[cols_ecs])
  
  # Identify each of the EC terms and place into a character vector
  #ec_terms_char <-  unlist(map(.ecs_in_model, rlang::quo_text))
  
  # Define the term that will go into grep
  # Specifically we want to search for terms that have 'spl' AND an EC term in it via grep
  ec_terms_for_grep <- paste(ec_terms_char, collapse="|")
  
  # Identify the fixed terms corresponding to one of the ECs in the model
  which_fixed_ec_terms <- grep(ec_terms_for_grep, attr(.fm$formulae$fixed, "term.labels") )
  fixed_ec_terms <-  attr(.fm$formulae$fixed, "term.labels")[which_fixed_ec_terms]
  
  # Place the EC terms into a data frame
  #ec_terms_df <- data.frame(Term=fixed_ec_terms)
  #str(ec_terms_df)
  
  # Determine the margin (i.e. order of hierachy) by how many times the strings ":" and "spl" appear in the term
  #ec_terms_df$Margin <-  str_count(ec_terms_df$Term, ":") + str_count(ec_terms_df$Term, "spl")
  
  #ec_terms_df$Margin <- factor(ec_terms_df$Margin)
  
  
  
  # Define an expression for the fixed and random effect EC terms to be added in the updated model
  # fixed_terms_curr <- .fm$call$fixed
  # random_terms_curr <- .fm$call$random
  # residual_terms_curr <- .fm$call$residual
  #
  # curr_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
  #                                         random =  !!random_terms_curr,
  #                                         residual= !!residual_terms_curr,
  #                                         data= .df,
  #                                         na.action = na.method(x = "include"),
  #                                         aom=T, maxit=30,
  #                                         wald=list(denDF="numeric",
  #                                                   ssType="conditional")))
  
  .fm <- update_fixed_asr(.fm=.fm)
  
  #.df  <<- base::eval(.fm$call$data)
  
  # Denote current model
  #wald_init_df <- asreml::wald.asreml(curr_fm, denDF='numeric', ssType='conditional')$Wald
  wald_init_df <- as.data.frame(.fm$aov)
  # If the last conditional F-value is NA, rerun and use incremental F-value instead
  h <- dim(wald_init_df)[1]
  if(is.na(wald_init_df$Fcon[h])){
    wald_init_df$M[h] <- max(wald_init_df$M, na.rm=TRUE) + 1
    marg <- wald_init_df$M
    .fm <- update_fixed_asr(.fm=.fm, ssType="incremental")
    wald_init_df <- as.data.frame(.fm$aov)
    wald_init_df$M <- marg
  }
  
  # Set current wald table to be equal to initial wald table
  # Set removed terms to be empty by default
  removed_terms <- rlang::maybe_missing()
  
  # Define margin term as an integer by first defining it as a factor
  #wald_curr_df$Margin <- factor(wald_curr_df$Margin) %>% as.integer()
  
  # Identify if the corresponding spline term is in the model for each EC
  wald_curr_df <- dropFixedTerm(.fm, wald_init_df, quo_ecs_in_model, .M, randomTerms=random_terms_curr)
  
  # Do a for loop for each EC to drop terms from the model
  # Set j equal to the current max value of margin
  j <- max(wald_curr_df$M)
  #Create a 2nd version of temp_ec_terms
  #temp_ec_terms_df2 <- temp_ec_terms_df
  continue_model_search <- TRUE
  # Now run models using a while-loop
  while(continue_model_search==TRUE){
    subset_ec_terms_df <- rownames(wald_curr_df[wald_curr_df$M==j, ])
    # Count how many terms have been removed during this iteration
    num_removed <- 0L
    # Set removed terms to be empty by default
    removed_terms <- rlang::maybe_missing()
    # Do we need to include k for each term within each margin??
    for(k in 1:length(subset_ec_terms_df)){
      term_to_test <- rlang::expr(!!subset_ec_terms_df[k]) #%>% parse_expr()
      # identify which row in the wald data frame is of relevance
      term_row <- which(rownames(wald_curr_df)==term_to_test)
      
      # Remove Ec term if P-value greater than 0.05 and term has no spline equivalent, then remove the term from the model
      if(wald_curr_df$Fprob[term_row] > 0.05 & wald_curr_df$canDrop[term_row]==TRUE){
        if(rlang::is_missing(removed_terms)){
          removed_terms <- as.character(term_to_test) #parse_expr(expr_text(term_to_test))
          num_removed <- num_removed + 1
        } else {
          removed_terms <- c(removed_terms,   term_to_test )
          num_removed <- num_removed + 1
        }
      }
    }
    if(num_removed>0){
      # Use self define function to remove random terms from the model for testing
      # fixed_terms_curr <- subtract_terms(main_expr = fixed_terms_curr,
      #                                     removed_char_vec = removed_terms,
      #                                     response=TRUE)
      #
      # curr_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
      #                                         random =  !!random_terms_curr,
      #                                         residual= !!residual_terms_curr,
      #                                         data= .df,
      #                                         na.action = na.method(x = "include"),
      #                                         aom=T, maxit=30,
      #                                         wald=list(denDF="numeric",
      #                                           ssType="conditional")))
      
      .fm <- update_fixed_asr(.fm=.fm, term=removed_terms)
      
      #.fm <- eval(curr_call)
      wald_curr_df <- as.data.frame(.fm$aov)
      #.df  <- base::eval(.fm$call$data)
      # Define margin term as an integer by first defining it as a factor
      #wald_curr_df$Margin <- factor(wald_curr_df$Margin) %>% as.integer()
      # Identify if the corresponding spline term is in the model for each EC
      wald_curr_df <- dropFixedTerm(.fm,wald_curr_df, quo_ecs_in_model, .M, randomTerms=random_terms_curr)
      j <- max(wald_curr_df$M)
    } else {
      continue_model_search <- FALSE
    }
  }
  final_fm <- .fm
  return(final_fm)
}



#' @title Identify the most parsimonious EC model
#' @description
#' Performs the backwards selection procedure for all environmental covariates in the curent model. 
#' This is achieve by cyclcing between dropping fixed and random environmental covariate terms until all non-significant environmental
#' covariate terms have been dropped from the current model.    
#' 
#' @param .fm The current \code{asreml} model object that the environmental covariate terms will be tested from.  
#' @param .ecs_in_model A list of expressions such that each expression is an environmental covariate that is present in the current model.
#' @param .G The genotype term in the model as an expression
#' @param .M The management practice term in the model as an expression
#' 
#' @return An updated asreml model such that all non-significant environmental covariate terms have been removed from the model.
#' 
#' @examples
#' 
#' @export
simplify_ec_model <- function(.fm, .ecs_in_model, .G, .M){
  
  # Identify the data frame from the model
  .df <<- base::eval(.fm$call$data)
  
  # Define the current model
  curr_fm <- .fm
  KeepSimplifying <- TRUE
  # Keep rotating between simplifying the fixed effects and the random effects intil the model can be no longer simplified
  while(KeepSimplifying==TRUE){
    new_random_fm <- ec_random_model(.fm=curr_fm, .ecs_in_model= .ecs_in_model, .G=.G, .M=.M)
    new_fm <- ec_fixed_model(.fm=new_random_fm, .ecs_in_model= .ecs_in_model, .G=.G, .M=.M)
    # If the fixed and random terms are exactly the same, then finish simplfying the model, otherwise, keep running to try and simplify further
    if( (new_fm$call$fixed==curr_fm$call$fixed) && (new_fm$call$random==curr_fm$call$random) ){
      KeepSimplifying <- FALSE
    } else {
      curr_fm <- new_fm
    }
  }
  return(new_fm)
}



