
## Definition of a step for filtering by prevalence (ie filter by fraction of zeros)

fraction_zeros <- function(x) {
  mean(x == 0, na.rm = TRUE)
}

step_filter_zero_fraction <- function(recipe, ..., 
                                      role = "predictor", 
                                      trained = FALSE, 
                                      min_zero_fraction = 0, 
                                      max_zero_fraction = 1, 
                                      columns = NULL, 
                                      skip = FALSE, 
                                      id = recipes::rand_id("filter_zero_fraction")) {
  
  recipes::add_step(
    recipe,
    step_filter_zero_fraction_new(
      # terms = enquos(...),
      terms = recipes::ellipse_check(...),
      role = role,
      trained = trained,
      min_zero_fraction = min_zero_fraction,
      max_zero_fraction = max_zero_fraction,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

step_filter_zero_fraction_new <- function(terms, role, trained, 
                                          min_zero_fraction, max_zero_fraction, 
                                          columns, skip, id) {
  recipes::step(
    subclass = "filter_zero_fraction",
    terms = terms,
    role = role,
    trained = trained,
    min_zero_fraction = min_zero_fraction,
    max_zero_fraction = max_zero_fraction,
    columns = columns,
    skip = skip,
    id = id
  )
}

prep.step_filter_zero_fraction <- function(x, training, info = NULL, ...) {
  # Use recipes features tools to evaluate columns with specified role
  col_names <- recipes::recipes_eval_select(x$terms, data = training, info = info)
  
  # Filter to only numeric columns
  num_col_names <- 
    training %>%
    select(all_of(col_names)) %>%
    select(where(is.numeric)) %>%
    names()
  
  zero_fractions <- 
    training %>% 
    select(all_of(num_col_names)) %>% 
    summarise(across(everything(), fraction_zeros)) %>% 
    pivot_longer(everything(), names_to = "feature", values_to = "zeros")
  
  keep_cols <- 
    zero_fractions %>% 
    filter(zeros >= x$min_zero_fraction,
           zeros <= x$max_zero_fraction) %>% 
    pull(feature)
  
  remove_cols <- num_col_names[ !num_col_names %in% keep_cols]
  
  step_filter_zero_fraction_new(
    terms = x$terms, 
    role = x$role,
    trained = TRUE,
    min_zero_fraction = x$min_zero_fraction,
    max_zero_fraction = x$max_zero_fraction,
    columns = remove_cols,
    skip = x$skip,
    id = x$id
  )
}

bake.step_filter_zero_fraction <- function(object, new_data, ...) {
  if (length(object$columns) > 0) {
    new_data <- 
      new_data %>%
      select(-any_of(object$columns))
  }
  new_data
}

print.step_filter_zero_fraction <- function(x, width = max(20, options()$width - 30), ...) {
  cat("Feature selection by number of zeros")
  
  if(recipes::is_trained(x)) {
    n <- length(x$columns)
    cat(paste0(" (", n, " excluded)"))
  }
  cat("\n")
  
  invisible(x)
  # cat("Step: Filter based on zero fraction retaining columns [", toString(x$columns), "]\n")
  # invisible(x)
}

tidy.step_filter_zero_fraction <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble(terms = x$columns)
  } else {
    term_names <- recipes::sel2char(x$terms)
    res <- tibble(terms = rlang::na_chr)
  }
  res$id <- x$id
  res
}
