#' @title Convert an expression into an asreml call term
#' @description
#' A function to converts an expression into the appropriate term to go into an asreml call
#' @param expr An expression
#'
#' @return A character vector
#'

expr_to_terms <- function(expr) {
  expr_char <- rlang::expr_text(expr)

  if(grepl(" + ",  expr_char)) {
    # If the expression string contains "+" operators, split terms by the "+"
    strsplit( expr_char, " \\+ ", fixed = TRUE)[[1]]
  } else {
    # If there is no "+", return the whole string as a single term
    return(expr_char)
  }
}




# Write a function to subtract one expression from another larger expression
# subtract_expr <- function(main_expr, remove_char_vec) {
#   # Convert the expressions to character vectors
#
#
#   main_terms <- as.character(main_expr)[-1] # Remove first element (tilde)
#   remove_terms <- as.character(remove_expr)[-1]
#
#   # Remove backticks from the term names
#   # main_terms <- gsub("`", "", main_terms)
#   # remove_terms <- gsub("`", "", remove_terms)
#
#   # Find unique terms not in the removal expression
#   unique_terms <- dplyr::setdiff(main_terms, remove_terms)
#
#   # Recombine the unique terms into an expression without using sym()
#   unique_terms_str <- paste(unique_terms, collapse=" + ")
#   result_expr <- parse_expr(paste0("~", unique_terms_str))
#
#   return(result_expr)
# }

#' @title Subtract one expression from a character vector
#' @description
#' A function that enables the subtraction of one character from a character vector. Used to ensure that the correct asreml model is called.
#' @param main_expr An expression
#' @param removed_char_vec A character vector? (Need to double check)
#' @param response A logical determining whether the response variable is included in main_expr. Default is FALSE.
#'
#' @return An expression to be included in an asreml call.
#'

# Write a function to subtract one expression from another larger expression
subtract_terms <- function(main_expr, removed_char_vec, response=FALSE) {
  # Convert the expressions to character vectors

  # convert expression to a character vector of length 1
  main_terms <- rlang::expr_text(main_expr) %>%
    stringr::str_replace_all("\n", "")# %>%  #remove any new line symbols
  #str_replace_all(" ", "") %>% #remove all spaces to stop any unexpected bugs from occurring
  # str_replace("~", "")         # Replace only the first tilde in the model
  # If the response variable is included then do not remove the tile
  if(response==TRUE){
    main_terms <- stringr::str_split_1(main_terms, pattern= '\\+')
    main_terms <- c(stringr::str_split_1(main_terms[1], pattern="\\~"), main_terms[-1])
  } else{
    main_terms <- main_terms %>% stringr::str_replace("~","")
    main_terms <- stringr::str_split_1(main_terms, pattern= '\\+')
  }

  # Create temporary version of main_terms without spaces
  main_temp <- main_terms %>% stringr::str_replace_all(" ","")

  removed_terms <- removed_char_vec %>%
    stringr::str_replace_all(" ", "") #%>%   #remove all spaces to stop any unexpected bugs from occurring
  #paste0("+", .) # include plus sign in the character vectors


  # Remove backticks from the term names
  # main_terms <- gsub("`", "", main_terms)
  # remove_terms <- gsub("`", "", remove_terms)

  # Use which to preserve terms that require a space
  which_removed_terms <- which(main_temp%in%removed_terms)
  # set unique terms equal to main term
  unique_terms <-  main_terms[-which_removed_terms]


  # Recombine the unique terms into an expression

  if(response==TRUE){
    # Combine the first 2 terms with a tilde
    unique_terms <- c(paste0(unique_terms[1], "~", unique_terms[2]), unique_terms[c(-1, -2)] )
    unique_terms_str <- paste(unique_terms, collapse=" + ")
    result_expr <- rlang::parse_expr(unique_terms_str)
  } else {
    unique_terms_str <- paste(unique_terms, collapse=" + ")
    result_expr <- rlang::parse_expr(paste0("~", unique_terms_str))
  }

  return(result_expr)
}





# # Write a function to determine the margin (i.e. order of hierarchy)
# temp_ec_terms_df
#
# terms <- str_replace_all(temp_ec_terms_df$Term, "PostPAW","x")
# terms <- str_replace_all(terms, "density","M")
# terms <- str_replace_all(terms, "Hybrid","G")

# # This works
# test <- margin(terms = temp_ec_terms_df$Term, ec="PostPAW", G="Hybrid", M="density")
# test
#
# # This does not work for some reason
# margin(terms = temp_ec_terms_df$Term,
#        ec=ec_terms_char[i], G=.G, M=.M)
#
# ec_terms_char[i] <- "PostPAW"
# # This does not work for some reason???
# margin(terms = temp_ec_terms_df$Term,
#        ec=ec_terms_char[i], G=.G, M=.M)
#
#
# margin(terms = temp_ec_terms_df$Term,
#        ec=ec_terms_char[i], G="Hybrid", M="density")
#
#
# # Compare with known results
# temp_ec_terms_df$Margin
#
# i <- 1
#
# ec_terms_char[i] <- "PostPAW"
# terms <- c("spl(x, k = 6)","M:spl(x, k = 6)","spl(x, k = 6):G","M:spl(x, k = 6):G","spl(M, k = 6):x","spl(M, k = 6):x:G","spl(M, k = 6):spl(x, k = 6)","spl(M, k = 6):spl(x, k = 6):G")
# ec=ec_terms_char[i]
# terms <- str_replace_all(terms, ec,"x")


# algorithm
#Start with the most complex term (set to i)
# If current margin value is the same, Assess whether every preceding term is nested within it,
# If there is a term which is nested within im repeat for for all i in reverse until i=2
# At the end subtract min(margin_value) from all terms so that the starting value is 1

# A function to determine which terms are nested within other terms
# Note that this function assumes that the similar lower-order terms come first and vice-versa

#' @title Determine the margin (i.e. order of hierarchy)
#' @description
#' This function calculates the equivalent of the margin terms (used in the wald tests for fixed effects) for the random terms in the model,
#' enabling ec_random_model to correctly determine which terms to test and when to finishing testing for significant interaction terms.
#' @param terms A character vector of terms in the asreml model
#' @param ec A character indicating which term in the model is the environmental covariate being tested
#' @param G A character indicating which term in the model is the genotype component
#' @param M A character indicating which term in the model is the management practice component
#' @param df A data frame containing the data used in the multi-environment trial analysis
#'
#' @return A vector of type double indicating the hierarchical structure of the random terms in the model for a particular environmental covariate
#'

margin <- function(terms, ec, G, M, df){
  if (is.na(ec) || ec == "") {
    stop("The 'ec' parameter is NA or an empty character.")
  }
  #putting together index of letters and symbols
  letter_array <- c(toupper(letters))
  # Define marginality matrix based upon whether M is categorical or continuous
  if(is.null(M)==TRUE){
    if(is.numeric(df[[ec]])==TRUE){
      # # Create a matrix determining the marginality
      terms_possible <- c("spl(x, k = 6)", "spl(x, k = 6):G")

      terms_matrix <- matrix(FALSE, nrow=length(terms_possible), ncol=length(terms_possible))
      rownames(terms_matrix) <- terms_possible
      colnames(terms_matrix) <- terms_possible

      # The matrix of type logical determines whether the term in the ith row is nested within the term in the jth column
      terms_matrix[1,] <- c(TRUE, TRUE)
      terms_matrix[2,] <- c(FALSE, TRUE)
      # Output error message if EC is not continuous
    } else {
      stop("The environmental covariate must be numeric when management practice is missing")
    }

  } else if(is.factor(df[[M]])==TRUE ){
    if(is.numeric(df[[ec]])==TRUE){
      # # Create a matrix determining the marginality
      terms_possible <- c("spl(x, k = 6)","M:spl(x, k = 6)","spl(x, k = 6):G","M:spl(x, k = 6):G")

      terms_matrix <- matrix(FALSE, nrow=length(terms_possible), ncol=length(terms_possible))
      rownames(terms_matrix) <- terms_possible
      colnames(terms_matrix) <- terms_possible

      # The matrix of type logical determines whether the term in the ith row is nested within the term in the jth column
      # Note THAT THIS MATRIX IS NON-SYMMETRIC
      terms_matrix[1,] <- c(TRUE, TRUE, TRUE, TRUE)
      terms_matrix[2,] <- c(FALSE, TRUE, FALSE, TRUE)
      terms_matrix[3,] <- c(FALSE, FALSE, TRUE, TRUE)
      terms_matrix[4,] <- c(FALSE, FALSE, FALSE, TRUE)
      # Output error message if EC is not continuous
    } else {
      stop("The environmental covariate must be numeric when management practice categorical")
    }
  } else {
    if(is.numeric(df[[ec]])==TRUE){
      # # Create a matrix determining the marginality
      terms_possible <- c("spl(x, k = 6)","M:spl(x, k = 6)","spl(x, k = 6):G","M:spl(x, k = 6):G","spl(M, k = 6):x","spl(M, k = 6):x:G","spl(M, k = 6):spl(x, k = 6)","spl(M, k = 6):spl(x, k = 6):G")

      terms_matrix <- matrix(FALSE, nrow=length(terms_possible), ncol=length(terms_possible))
      rownames(terms_matrix) <- terms_possible
      colnames(terms_matrix) <- terms_possible

      # The matrix of type logical determines whether the term in the ith row is nested within the term in the jth column
      terms_matrix[1,] <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
      terms_matrix[2,] <- c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
      terms_matrix[3,] <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
      terms_matrix[4,] <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE)
      terms_matrix[5,] <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
      terms_matrix[6,] <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE)
      terms_matrix[7,] <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
      terms_matrix[8,] <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
    } else {
      # # Create a matrix determining the marginality
      terms_possible <- c("spl(M, k = 6):x","spl(M, k = 6):x:G")

      terms_matrix <- matrix(FALSE, nrow=length(terms_possible), ncol=length(terms_possible))
      rownames(terms_matrix) <- terms_possible
      colnames(terms_matrix) <- terms_possible

      # The matrix of type logical determines whether the term in the ith row is nested within the term in the jth column
      terms_matrix[1,] <- c(TRUE, TRUE)
      terms_matrix[2,] <- c(FALSE, TRUE)
    }
  }
  # Replace ec, G, and M with their corresponding names
  terms <- stringr::str_replace_all(terms, ec,"x")
  if(is.null(M)==FALSE){
    terms <- stringr::str_replace_all(terms, M,"M")
  }
  terms <- stringr::str_replace_all(terms, G,"G")
  # Set the default margin values to 4 (i.e. all terms are nested by default)
  margin_value <- rep(4, length(terms))
  for(i in (length(terms)):2 ){
    term_i <- terms[i]
    which_col_i <- which(term_i==colnames(terms_matrix))[1]
    # Define a logical vector indicating if j is nested in i for each iteration of j
    is_j_nested_in_i <- TRUE
    # Determine which values less than i have equal margin values
    which_terms_with_equal_margin <- which(margin_value[1:(i-1)]==margin_value[i])
    for(j in which_terms_with_equal_margin) {
      term_j <- terms[j]
      which_row_j <- which(term_j==rownames(terms_matrix))[1]
      is_j_nested_in_i <- terms_matrix[which_row_j, which_col_i]
      # If j is nested within i, subtract the margin value of j by 1
      if(is_j_nested_in_i==TRUE){
        margin_value[j] <- margin_value[j]-1
      }
    }
    # if( (sum(is_j_nested_in_i==FALSE)>0)==FALSE ){
    #   margin_value[i:length(terms)] <- margin_value[i:length(terms)] + 1
    # }
  }
  # At the end subtract min(margin_value) from all terms so that the starting value is 1
  margin_value <- margin_value - min(margin_value) + 1
  return(margin_value)
}





# A function to determine if the corresponding spline term for a fixed effects EC term is currently in the model
# Input a wald table, returns the wald table with a extra column in the data frame of type logical

#' @title Determine if corresponding spline term is present for each linear term
#' @description
#' A function to determine if the corresponding spline term for a fixed effects EC term is currently in the model.
#' If the spline term is still in the model, then the linear term (fit as a fixed effect) cannot be dropped from the model, and hence will not be tested in \code{ec_fixed_model}.
#' @param .fm An \code{asreml} model object
#' @param wald_df A data frame containing the information pertaining to the wald test (see \code{asreml::wald.asreml})
#' @param .ecs_in_model Need to check (quosure?)
#' @param .M A character indicating which term in the model is the management practice component
#' @param randomTerms Check if this is still needed in the model
#'
#' @return An updated Wald table with an additional column of type logical indicating whether each fixed effect term in the model can be dropped.
#'
dropFixedTerm <- function(.fm, wald_df, .ecs_in_model, .M){
  randomTerms <- attr(.fm$formulae$random, "term.labels") %>%
    stringr::str_replace_all(" ", "") %>%  #remove all spaces to stop any unexpected bugs from occurring
    stringr::str_replace_all(",k=\\d+","") # Remove the term after the comma within each spl() function to help with matching with the fixed term
  # Identify each of the EC terms and place into a character vector
  ec_terms_char <-  unlist(purrr::map(.ecs_in_model, rlang::quo_text))
  # Generate EC spline term replacements
  spl_ec_terms_char <- stats::setNames(paste0("spl(", ec_terms_char, ")"), ec_terms_char)
  # Add spl_equiv_ec_terms to wald_df
  new_wald_df <- wald_df
  new_wald_df$spl_equiv_ec <- stringr::str_replace_all(rownames(wald_df), spl_ec_terms_char)
  # For each EC, check if the equivalent spline term is in the current model w.r.t spl(EC) & spl(M)
  new_wald_df$pres_spl_ec <- new_wald_df$spl_equiv_ec%in%randomTerms

  if(is.factor(.M)==TRUE){
    # Subset vector to include only elements with replacements (i.e., those containing 'spl(')
    new_wald_df$spl_equiv_M <- stringr::str_replace_all(rownames(wald_df), .M, paste0("spl(", .M, ")"))
    new_wald_df$pres_M_ec <- new_wald_df$spl_equiv_M%in%randomTerms

    # Define a logical variable indicating if each term can be dropped
    new_wald_df$canDrop <- !(new_wald_df$pres_spl_ec | new_wald_df$pres_M_ec)
  } else {
    # Define a spline equivalent with a spl(density) term
    new_wald_df$spl_equiv_m <- stringr::str_replace_all(rownames(wald_df), .M,
                                  paste0('spl(',.M,')') )
    # For each EC, check if the equivalent spline term is in the current model w.r.t spl(EC) & spl(M)
    new_wald_df$pres_spl_m <- new_wald_df$spl_equiv_m%in%randomTerms

    # Define a logical variable indicating if each term can be dropped
    # Note that it can only be dropped if both spl(.M) and spl(ec) are NOT in the current model
    new_wald_df$canDrop <- !(new_wald_df$pres_spl_ec | new_wald_df$pres_spl_m)
  }
  return(new_wald_df)
}


#' @title Run an asreml model with updated fixed effects component
#' @description
#' This function is used to runs an updated version of asreml with a modified fixed effects model within \code{ec_fixed_model} with higher order fixed terms dropped in the model.
#' This is necessary for regression modelling to ensure that the Wald tests are correct for lower order terms in the model.
#' @param .fm An \code{asreml} model object
#' @param term A character scalar indicating which term in the model is currently being tested
#' @param ssType The type of sums of squares used for testing. Default is "\code{conditional}" which is analogous to Type 3 sums of squares.
#' Use "\code{incremental}" for Type 1 sums of squares. See  \code{asreml::wald.asreml} for more information.
#'
#' @return An asreml model object with a new fixed effets model.
#'
update_fixed_asr <- function(.fm, term=rlang::maybe_missing(), denDF="numeric", ssType="conditional"){

  .df <<- base::eval(.fm$call$data)
  # Define an expression for the fixed and random effect EC terms to be added in the updated model

  curr_fm <- .fm

  fixed_terms_curr <- curr_fm$call$fixed
  random_terms_curr <- curr_fm$call$random
  residual_terms_curr <- curr_fm$call$residual
  curr_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
                                          random =  !!random_terms_curr,
                                          residual = !!residual_terms_curr,
                                          data= .df,
                                          na.action = na.method(x = "include"),
                                          aom=T, maxit=30,
                                          wald=list(denDF= !!denDF,
                                                    ssType= !!ssType)))


  #curr_fm <- eval(curr_call)
  if(rlang::is_missing(term)==FALSE){
    print('Simplifying the FIXED effects model')
    fixed_terms_curr <- subtract_terms(main_expr = fixed_terms_curr,
                                       removed_char_vec = term,
                                       response=TRUE)

    curr_call <- rlang::expr(asreml::asreml(fixed= !!fixed_terms_curr  ,
                                            random =  !!random_terms_curr,
                                            residual = !!residual_terms_curr,
                                            data= .df,
                                            na.action = na.method(x = "include"),
                                            aom=T, maxit=30,
                                            wald=list(denDF=!!denDF,
                                                      ssType="conditional")))

  }
  curr_fm <- eval(curr_call)
  curr_fm$aov
  if(denDF=="none"){
    wald_df <- asreml::wald.asreml(curr_fm, denDF="none", ssType="conditional")$Wald
    wald_df
    curr_fm$aov[,colnames(curr_fm$aov)=="Finc"] <- wald_df[,colnames(wald_df)=="Wald(inc)"]
    curr_fm$aov[,colnames(curr_fm$aov)=="Fcon"] <- wald_df[,colnames(wald_df)=="Wald(con)"]
    curr_fm$aov[,colnames(curr_fm$aov)=="Fprob"] <- wald_df[,colnames(wald_df)=="Pr(chisq)"][,2]
  }
  # If the last term has no conditional F value, use the incremental value instead
  h <- dim(curr_fm$aov)[1]
  if(is.na(curr_fm$aov[,colnames(curr_fm$aov)=="Fprob"][h])){
    curr_fm$aov[,colnames(curr_fm$aov)=="Fcon"][h] <- curr_fm$aov[,colnames(curr_fm$aov)=="Finc"][h]
    curr_fm$aov[,colnames(curr_fm$aov)=="Fprob"][h] <- 1-pchisq(curr_fm$aov[,colnames(curr_fm$aov)=="Fcon"][h],
                                                                curr_fm$aov[,colnames(curr_fm$aov)=="df"][h])
    # Make margin term 1 more than the previous term (assums the previous term is a lower order interaction; MAY NEED TO FIX)
    curr_fm$aov[,colnames(curr_fm$aov)=="M"][h] <- curr_fm$aov[,colnames(curr_fm$aov)=="M"][h-1] + 1
  }
  return(curr_fm)
}





#' @title Creates an updated model which is the same as the current model but adds one new EC into the model
#' @description
#' Runs an updated asreml model with the additional fixed and random effects term for the environmental covariate added to the model.
#'
#' @param .fm An \code{asreml} model object without the environmental covariate in the model.
#' @param .ec An expression with the environmental covariate to be included in the model.
#' @param .G An expression indicating which term in the model is the genotype component.
#' @param .M An expression indicating which term in the model is the management practice component.
#' @param .kn  The number of knot points to be included for the spline terms related to the environmental covariate \code{.ec} being added to the model.
#'
#' @return An asreml model object with all of the fixed and random terms related to the environmental covariate \code{.ec} being added to the model.
#' @examples
#' library(asreml)
#' data(SorghumYield)
#' data(SorghumCvGroup)
#' baseline_asr <- asreml( Yld ~ Genotype + density + Genotype:density,
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
#' postPAW_full_asr <-  ec_full_model_constructor(.fm=baseline_asr, .ec=rlang::expr(PostPAW), .G=rlang::expr(Genotype), .M=rlang::expr(density))
#' postPAW_full_asr$call
#'
#' @export
ec_full_model_constructor <- function(.fm, .ec, .G, .M, .kn=6){

  # Identify the data frame from the model
  .df <<- base::eval(.fm$call$data)

  # Allow these terms to be included in the model unquoted
  # .ec <- enexpr(.ec)
  # .G <- enexpr(.G)
  # .M <- enexpr(.M)

  # Identify the fixed, random and resiudal terms currently in the model
  response_term <- attr(.fm$formulae$fixed, "variables")[[2]] %>%
    as.character() %>%
    rlang::parse_expr()
  fixed_bl_terms <- attr(.fm$formulae$fixed, "term.labels") %>% paste(collapse="+") %>%
    rlang::parse_expr()

  #fixed_terms <- .fm$call$fixed
  random_bl_terms <- .fm$call$random
  residual_terms <- .fm$call$residual

  # CONTINUE HERE (NEED TO ACCOUNT FOR CONTINUOUS ECs)
  if(is.null(.M)==TRUE){
    fixed_ec_terms <-  rlang::expr(!!.ec + !!.ec:!!.G)
    if(is.numeric(.df[[.ec]])==TRUE){
      random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G)
      random_terms <- rlang::parse_expr(paste0( rlang::expr_text(random_bl_terms), "+", rlang::expr_text( random_ec_terms)))

    } else {
      random_terms <- rlang::parse_expr(paste0( rlang::expr_text(random_bl_terms)))
    }
  } else if(is.factor(.df[[.M]])==TRUE){ # Change the random terms based on whether M is categorical or continuous
      fixed_ec_terms <-  rlang::expr(!!.ec + !!.M:!!.ec + !!.ec:!!.G + !!.M:!!.ec:!!.G)
      if(is.numeric(.df[[.ec]])==TRUE){
        random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + !!.M:spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G +
                                         !!.M:spl(!!.ec, k=!!.kn):!!.G)
        random_terms <- rlang::parse_expr(paste0( rlang::expr_text(random_bl_terms), "+", rlang::expr_text( random_ec_terms)))

      } else {
        random_terms <- rlang::parse_expr(paste0(rlang::expr_text(random_bl_terms)))
      }
  } else { #For continuous M
      # Change the random terms based on whether M is categorical or continuous
      fixed_ec_terms <-  rlang::expr(!!.ec + !!.M:!!.ec + !!.ec:!!.G + !!.M:!!.ec:!!.G)
      if(is.numeric(.df[[.ec]])==TRUE){
        random_ec_terms <- rlang::expr(spl(!!.ec, k=!!.kn) + !!.M:spl(!!.ec, k=!!.kn) + spl(!!.ec, k=!!.kn):!!.G +
                                         !!.M:spl(!!.ec, k=!!.kn):!!.G +
                                         spl(!!.M, k=!!.kn):!!.ec + spl(!!.M, k=!!.kn):!!.ec:!!.G +
                                         spl(!!.M, k=!!.kn):spl(!!.ec, k=!!.kn) +
                                         spl(!!.M, k=!!.kn):spl(!!.ec, k=!!.kn):!!.G)
        random_terms <- rlang::parse_expr(paste0( rlang::expr_text(random_bl_terms), "+", rlang::expr_text( random_ec_terms)))
      } else {
        random_ec_terms <- rlang::expr(spl(!!.M, k=!!.kn):!!.ec + spl(!!.M, k=!!.kn):!!.ec:!!.G)
        random_terms <- rlang::parse_expr(paste0( rlang::expr_text(random_bl_terms), "+", rlang::expr_text( random_ec_terms)))
      }
  }
  # Combine all fixed terms into a single expression respectively
  fixed_terms <- rlang::parse_expr(paste0( rlang::expr_text(fixed_bl_terms), "+", rlang::expr_text( fixed_ec_terms)))

  asr_call <- rlang::expr(asreml::asreml(fixed = !!response_term ~ !!fixed_terms,
                                         random = !!random_terms, # removed tilde as the tilde should already be in the call
                                         residual= !!residual_terms, # removed tilde as the tilde should already be in the call
                                         data=.df,
                                         na.action=asreml::na.method(x='include'),
                                         aom=T, maxit=300))

  ec_full_fm <- eval(asr_call)
  return(ec_full_fm)
}


#' @title Creates the list of information required for asreml.predict with continuous variables
#' @description
#' Creates the list of information required for asreml.predict with continuous variables
#'
#' @param .df The data frame used in the core asreml model
#' @param subset_df The data frame used in the subbsetted version during the cross validation procedure
#' @param .ec An expression with the environmental covariate to be included in the model.
#' @param .G An expression indicating which term in the model is the genotype component.
#' @param .M An expression indicating which term in the model is the management practice component.
#' @param baseline_ec_cols  A vector of integers indicating which columns in both \code{.df} and \code{subset_df} contain
#' environmental covariates that are in the current model
#'
#' @return A list of information required for asreml.predict with continuous variables (TBC!)
#' @export
parallel_predict_list <- function(.df, subset_df, .ec=NULL, .G, .E, .M, baseline_ec_cols=NULL ){
  # Define a table that will be used to generate model predictions later on
  # Identify the classify list and levels of each factor based on whether a trial term is included in the model
  if(is.null(.ec)==TRUE){
    if(is.null(baseline_ec_cols)==TRUE){
      # Change the classify terms if M is missing from the baseline model
      if(is.null(.M)==TRUE){
        expr_list <- rlang::exprs(!!.E, !!.G) # E is kept in to ensure that if the same combination appears multiple times
        #across environments, then that weighting is preserved in the cross validation scheme
        # Create list of values we want to predict for
        aux_parallel <- unique(.df[,purrr::map_chr(expr_list, rlang::as_string)]) %>%
          purrr::modify_if(is.numeric, round, digits=4) # Round all continuous variables to 4 decimal places


        # Identify which .G are not in the subsetted data frame at all
        missing_g <- setdiff(levels(.df[[.G]]), levels(subset_df[[.G]]))
        # Remove .G that are not in the subsetted data frame
        aux_parallel <- aux_parallel[!aux_parallel[[.G]]%in%missing_g, ]
        # Update factor levels of genotype term
        aux_parallel[[.G]] <- factor(aux_parallel[[.G]])

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

        # Identify which .G are not in the subsetted data frame at all
        missing_g <- setdiff(levels(.df[[.G]]), levels(subset_df[[.G]]))
        # Remove .G that are not in the subsetted data frame
        aux_parallel <- aux_parallel[!aux_parallel[[.G]]%in%missing_g, ]
        # Update factor levels of genotype term
        aux_parallel[[.G]] <- factor(aux_parallel[[.G]])

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
  } else {
    if(is.null(baseline_ec_cols)==TRUE){
      # Change the classify terms if M is missing from the baseline model
      if(is.null(.M)==TRUE){
        expr_list <- rlang::exprs(!!.E, !!.G, !!.ec) # E is kept in to ensure that if the same combination appears multiple times
        #across environments, then that weighting is preserved in the cross validation scheme
        # Create list of values we want to predict for
        aux_parallel <- unique(.df[,purrr::map_chr(expr_list, rlang::as_string)]) %>%
          purrr::modify_if(is.numeric, round, digits=4) # Round all continuous variables to 4 decimal places

        # Identify which .G are not in the subsetted data frame at all
        missing_g <- setdiff(levels(.df[[.G]]), levels(subset_df[[.G]]))
        # Remove .G that are not in the subsetted data frame
        aux_parallel <- aux_parallel[!aux_parallel[[.G]]%in%missing_g, ]
        # Update factor levels of genotype term
        aux_parallel[[.G]] <- factor(aux_parallel[[.G]])

        # Define the terms corresponding to G, M and ec to appear in the classify statement
        classify_terms <- rlang::expr(!!.G:!!.ec)
        # Use expr_text() to add quotations around the expression
        classify_terms <- rlang::expr_text(classify_terms)

        # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
        levels_list <- list(aux_parallel[[.G]],
                            aux_parallel[[.ec]])

        # Now give the headings of the levels_list names
        names(levels_list) <- c(rlang::expr_text(.G), rlang::expr_text(.ec))
      } else {
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
      }
    } else {
      if(is.null(.M)==TRUE){
        expr_list <- rlang::exprs(!!.E, !!.G, !!.ec)

        # define a tibble version of .df so that the next line works correctly for 1 baseline EC
        df_tib <- tibble::as_tibble(.df)

        # Create list of values we want to predict for
        aux_parallel <- unique(df_tib[, c(purrr::map_chr(expr_list, rlang::as_string),colnames(df_tib[,baseline_ec_cols]))]) %>%
          purrr::modify_if(is.numeric, round, digits=4) %>% # Round all continuous variables to 4 decimal places
          as.data.frame()

        # Identify which .G are not in the subsetted data frame at all
        missing_g <- setdiff(levels(.df[[.G]]), levels(subset_df[[.G]]))
        # Remove .G that are not in the subsetted data frame
        aux_parallel <- aux_parallel[!aux_parallel[[.G]]%in%missing_g, ]
        # Update factor levels of genotype term
        aux_parallel[[.G]] <- factor(aux_parallel[[.G]])

        bl_ecs_colon <- rlang::parse_expr(paste(colnames(.df[baseline_ec_cols]), collapse=":"))
        # Define the terms corresponding to G, M and ec to appear in the classify statement
        classify_terms <- rlang::expr(!!.G:!!bl_ecs_colon:!!.ec)
        # Use expr_text() to add quotations around the expression
        classify_terms <-  gsub("\\(|\\)", "", rlang::expr_text(classify_terms)) # Remove all parentheses from the character string
        # Generate the list of levels that will be used to calculate the predictions during the cross validation scheme
        levels_list <- purrr::map(aux_parallel[,-1], as.vector) # Remove the .E column which should always be column 1
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
    }
  }
  return(list(levels_list=levels_list, aux_parallel=aux_parallel, classify_terms=classify_terms))
}



# if denDF is too difficult to calculate, automatically make it 30
wald_approx_pvalue <- function(aov_df){
  # Set the denDF tp be 30 for all terms as a large sample approximation
  aov_df$denDF[is.na(aov_df$denDF)] <- 30
  # Calculate p-values based on large same denDF
  aov_df$Fprob <- 1-pf(q=aov_df$Fcon, df1=aov_df$df, df2=aov_df$denDF)
  return(aov_df)
}


