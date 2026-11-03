capture_expr <- function(x) {

  # Handle expr(Genotype)
  if (rlang::is_call(x, "expr")) {
    return(x[[2]])
  }

  # Handle quosures
  if (inherits(x, "quosure")) {
    return(rlang::get_expr(x))
  }

  x
}

capture_expr_user <- function(x) {

  expr <- rlang::enexpr(x)

  if (rlang::is_call(expr, "expr")) {
    expr <- expr[[2]]
  }

  expr
}

normalise_expr <- function(x) {

  if (inherits(x, "quosure")) {
    return(rlang::get_expr(x))
  }

  if (rlang::is_call(x, "expr")) {
    return(x[[2]])
  }

  x
}


capture_ecs <- function(ECs, data) {

  vars <- as.list(
    setNames(
      seq_along(data),
      names(data)
    )
  )

  # Case 1:
  # already a quosure list
  if (rlang::is_quosures(ECs)) {

    cols <- unlist(
      purrr::map(
        ECs,
        rlang::eval_tidy,
        vars
      )
    )

    return(
      rlang::parse_quos(
        names(data)[cols],
        env = rlang::current_env()
      )
    )
  }

  # Case 2:
  # expression such as quos(EC1:EC2)
  if (rlang::is_call(ECs, "quos")) {

    quo_obj <- eval(
      ECs,
      envir = parent.frame()
    )

    cols <- unlist(
      purrr::map(
        quo_obj,
        rlang::eval_tidy,
        vars
      )
    )

    return(
      rlang::parse_quos(
        names(data)[cols],
        env = rlang::current_env()
      )
    )
  }

  # Case 3:
  # expression such as EC1:EC2
  cols <- rlang::eval_tidy(
    ECs,
    vars
  )

  rlang::parse_quos(
    names(data)[cols],
    env = rlang::current_env()
  )
}

normalise_name <- function(x){

  if(is.character(x)){
    return(x)
  }

  x <- normalise_expr(x)

  rlang::as_name(x)
}


capture_quos <- function(x) {

  if (rlang::is_quosures(x)) {
    return(x)
  }

  rlang::enquos(x)
}
