## Definition of a step for filtering by predictive ability together with some set of defined values

lin_comb <- function(ds, predef_predictor_cols, new_predictor_cols, outcome_col) {
  
  baseline_formula <- paste("outcome_variable", "~", paste(ifelse(nchar(predef_predictor_cols[1])>0, predef_predictor_cols, 1), collapse = "+"), collapse = "")
  baseline_model <- glm(formula = as.formula(baseline_formula), data = ds %>% rename(outcome_variable = all_of(outcome_col)), family = "binomial")
  
  model_by_feature <-
    ds %>% 
    select(id, any_of(unique(c(predef_predictor_cols,
                               new_predictor_cols,
                               outcome_col)))) %>% 
    pivot_longer(-c(id, any_of(c(outcome_col, predef_predictor_cols)))) %>% 
    group_by(name) %>% 
    nest() %>% 
    mutate(models = map(.x = data, .f = function(x = .x) {
      
      if (!all(x$value == 0)) {
        ## Impute zeros and transform abundance
        x <- x %>% mutate(value = case_when(value == 0 ~ min(value[ value > 0])/2,
                                            TRUE ~ value) %>%
                            log10())
        
        ## Dummy outcome variable
        x <- x %>% rename(outcome_variable = target)
        
        x <- x %>% mutate(outcome_variable = as.integer(as.factor(outcome_variable))-1)
        if (all(predef_predictor_cols != "")) {
          tmp_formula <- paste("outcome_variable", "~", paste(c(predef_predictor_cols, "value"), collapse = "+"), collapse = "")
        } else {
          tmp_formula <- "outcome_variable ~ value"
        }
        
        
        glm(formula = tmp_formula, data = x, family = "binomial")
      } else {
        NA
      }
    })) %>% 
    ungroup() %>% 
    mutate(aic = sapply(models, function(mod_output) ifelse(inherits(mod_output, "glm"), summary(mod_output)$aic, NA))) %>% 
    mutate(auc = map2_dbl(.x = data, .y = models, .f = function(dat = .x, mod = .y) {
      if (inherits(mod, "glm")) {
        preds <- predict(mod, type = "response")
        roc(response = dat %>% pull(all_of(outcome_col)), levels = c("negative", "positive"), predictor = preds, direction = "<")$auc[[1]]  
      } else {
        NA
      }
    })) %>%
    mutate(lrt_p = map_dbl(.x = models, .f = function(mod = .x) {
      anova(baseline_model, mod, test = "Chisq") %>% tidy() %>% filter(str_detect(term, "(\\+|\\~) value$")) %>% pull(p.value)
    })) %>% 
    select(-c(data, models))
  
  model_by_feature
  # 
  
}

step_select_lin_comb <- function(recipe,
                           ...,
                           outcome,
                           role = "predictor",
                           trained = FALSE,
                           # threshold = NA,
                           n_top_contr = NA,
                           adjustment_dataset = "",
                           selection_criterion = NA,
                           min_lrt_p = NA,
                           exclude = NULL,
                           skip = FALSE,
                           id = recipes::rand_id("select_lin_comb")) {
  recipes::add_step(
    recipe,
    step_select_lin_comb_new(
      terms = recipes::ellipse_check(...),
      outcome = outcome,
      role = role,
      trained = trained,
      # threshold = threshold,
      n_top_contr = n_top_contr,
      adjustment_dataset = adjustment_dataset,
      selection_criterion = selection_criterion,
      min_lrt_p = min_lrt_p,
      exclude = exclude,
      skip = skip,
      id = id
    )
  )
}

step_select_lin_comb_new <-
  function(terms, outcome, role, trained, 
           n_top_contr, adjustment_dataset, selection_criterion, min_lrt_p,
           exclude, skip, id) {
    recipes::step(
      subclass = "select_lin_comb",
      terms = terms,
      outcome = outcome,
      role = role,
      trained = trained,
      # threshold = threshold,
      n_top_contr = n_top_contr,
      adjustment_dataset = adjustment_dataset,
      selection_criterion = selection_criterion,
      min_lrt_p = min_lrt_p,
      exclude = exclude,
      skip = skip,
      id = id
    )
  }

prep.step_select_lin_comb <- function(x, training, info = NULL, ...) {
  y_name <- recipes::recipes_eval_select(x$outcome, data = training, info = info)
  y_name <- x$outcome[1]
  recipes::check_type(training[, y_name], quant = FALSE)
  x_names <- recipes::recipes_eval_select(x$terms, data = training, info = info)

  if(length(x_names) > 0) {

    # recipes::check_type(training[, x_names], )
    
    n_top_contr <- x$n_top_contr
    selection_criterion <- x$selection_criterion
    if (x$adjustment_dataset != "") {
      adjustment_dataset <- read_tsv(x$adjustment_dataset, col_types = cols())  
      
      ## Define adjustment variables - these will not be filtered out if they are also in training data
      adjustment_vars <- names(adjustment_dataset %>% select(-id))
      
      adjustment_dataset <- 
        adjustment_dataset %>% 
        select(-which(names(adjustment_dataset) %in% x_names))
    } else {
      adjustment_vars <- ""
    }
    
    
    
    feature_evaluations <- 
      training %>% 
      (function(tmp_dat) {
        if (x$adjustment_dataset != "") {
          tmp_dat %>% 
            left_join(adjustment_dataset, by = join_by(id))
        } else {
          tmp_dat
        }
      }) %>% 
      lin_comb(predef_predictor_cols = adjustment_vars, 
               new_predictor_cols = x_names,
               outcome_col = y_name)
    
    if (selection_criterion %in% "auc") {
      predictors_to_keep <-
        feature_evaluations %>% 
        slice_max(order_by = auc, n = n_top_contr) %>% 
        (function(tmp) if(!is.na(x$min_lrt_p)) {
          tmp %>% filter(lrt_p <= x$min_lrt_p)
          } else {
            tmp
          }) %>% 
        # slice_max(order_by = !!sym(selection_criterion), n = n_top_contr) %>% 
        pull(name)
    } else if (selection_criterion %in% "aic")  {
      predictors_to_keep <-
        feature_evaluations %>% 
        slice_min(order_by = aic, n = n_top_contr) %>% 
        (function(tmp) if(!is.na(x$min_lrt_p)) {
          tmp %>% filter(lrt_p <= x$min_lrt_p)
          } else {
            tmp
          }) %>% 
        # slice_min(order_by = !!sym(selection_criterion), n = n_top_contr) %>% 
        pull(name)
    }
    
    
    exclude_chr <- feature_evaluations$name[ !feature_evaluations$name %in% predictors_to_keep]
    
  } else {
    exclude_chr <- character()
  }

  step_select_lin_comb_new(
    terms = x$terms,
    outcome = x$outcome,
    role = x$role,
    trained = TRUE,
    # threshold = x$threshold,
    n_top_contr = x$n_top_contr,
    selection_criterion = x$selection_criterion,
    adjustment_dataset = x$adjustment_dataset,
    min_lrt_p = x$min_lrt_p,
    exclude = exclude_chr,
    skip = x$skip,
    id = x$id
  )
}

bake.step_select_lin_comb <- function(object, new_data, ...) {
  
  if (length(object$exclude) > 0) {
    new_data <- new_data %>% dplyr::select(-dplyr::one_of(object$exclude))
  }
  new_data
}

print.step_select_lin_comb <- function(x, width = max(20, options()$width - 30), ...) {
  cat("Linear combination feature selection")

  if(recipes::is_trained(x)) {
    n <- length(x$exclude)
    cat(paste0(" (", n, " excluded)"))
  }
  cat("\n")

  invisible(x)
}

tidy.step_select_lin_comb <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble(terms = x$exclude)
  } else {
    term_names <- recipes::sel2char(x$terms)
    res <- tibble(terms = rlang::na_chr)
  }
  res$id <- x$id
  res
}

# tunable.step_select_lin_comb <- function(x, ...) {
#   tibble::tibble(
#     name = c("top_p", "threshold"),
#     call_info = list(
#       list(pkg = "recipeselectors", fun = "top_p"),
#       list(pkg = "dials", fun = "threshold", range = c(0, 1))
#     ),
#     source = "recipe",
#     component = "step_select_roc",
#     component_id = x$id
#   )
# }
