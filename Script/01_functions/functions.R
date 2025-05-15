# Generate Edgelists

## Dependent


generate_edgelists_parallel <- function(effects,
                                        effect_sizes,
                                        covar,
                                        num_events,
                                        num_cores,
                                        start,
                                        end,
                                        folder) {
  future::plan(future::multisession, workers = num_cores)
  environment(effects) <- environment()
  future.apply::future_lapply(start:end, function(i) {
    set.seed(i)
    edgelist <- remulate::remulateTie(effects,
                                      actors = 1:50,
                                      events = num_events,
                                      time   = 10000)
    edgelist <- as.data.frame(lapply(edgelist, function(x) {
      attributes(x) <- NULL
      x
    }))
    save(edgelist, file = paste0(folder, "edgelist_", i, ".RData"))
    rm(edgelist)
    gc()
  }, future.seed = NULL, future.globals  = list(
    effects      = effects,
    effect_sizes = effect_sizes,
    covar        = covar,
    num_events   = num_events,
    folder       = folder
  ), future.packages = "remulate")
  
}



## Independent


generate_edgelists_independent_parallel <- function(effects,
                                                    effect_sizes,
                                                    covar,
                                                    num_events,
                                                    num_cores,
                                                    start,
                                                    end,
                                                    folder) {
  future::plan(future::multisession, workers = num_cores)
  environment(effects) <- environment()
  future.apply::future_lapply(start:end, function(i) {
    set.seed(num_events + i)
    edgelist <- remulate::remulateTie(effects,
                                      actors = 1:50,
                                      events = num_events,
                                      time   = 10000)
    edgelist <- as.data.frame(lapply(edgelist, function(x) {
      attributes(x) <- NULL
      x
    }))
    save(edgelist, file = paste0(folder, "edgelist_", i, ".RData"))
    rm(edgelist)
    gc()
  }, future.seed = NULL, future.globals  = list(
    effects      = effects,
    effect_sizes = effect_sizes,
    covar        = covar,
    num_events   = num_events,
    folder       = folder
  ), future.packages = "remulate")
  
}



# Generate Function from Parameter List


generate_formula <- function(parameters) {
  terms <- c("1") # Start with intercept
  
  for (name in names(parameters)) {
    # Detect interaction terms
    if (grepl(":", name)) {
      # Split interaction terms
      components <- strsplit(name, ":")[[1]]
      
      # Process each part of interaction
      processed_components <- sapply(components, function(component) {
        if (grepl("^difference_", component)) {
          var <- sub("^difference_", "", component)
          return(paste0("difference('", var, "', scaling = 'std')"))
        } else if (grepl("^same_", component)) {
          var <- sub("^same_", "", component)
          return(paste0("same('", var, "')"))
        } else if (grepl("^send_", component)) {
          var <- sub("^send_", "", component)
          return(paste0("send('", var, "', scaling = 'std')"))
        } else if (grepl("^receive_", component)) {
          var <- sub("^receive_", "", component)
          return(paste0("receive('", var, "', scaling = 'std')"))
        } else if (grepl("^tie_", component)) {
          var <- sub("^tie_", "", component)
          return(paste0("tie('", var, "', scaling = 'std')"))
        } else if (grepl("^average_", component)) {
          var <- sub("^average_", "", component)
          return(paste0("average('", var, "', scaling = 'std')"))
        } else if (grepl("^minimum_", component)) {
          var <- sub("^minimum_", "", component)
          return(paste0("minimum('", var, "', scaling = 'std')"))
        } else if (grepl("^maximum_", component)) {
          var <- sub("^maximum_", "", component)
          return(paste0("maximum('", var, "', scaling = 'std')"))
        } else if (grepl("^event_", component)) {
          var <- sub("^event_", "", component)
          return(paste0("event('", var, "')"))
        } else if (grepl("^userStat_", component)) {
          var <- sub("^userStat_", "", component)
          return(paste0("userStat('", var, "')"))
        } else if (component == "inertia") {
          return("inertia(scaling = 'std')")
        } else if (component == "indegreeSender") {
          return("indegreeSender(scaling = 'std')")
        } else if (component == "indegreeReceiver") {
          return("indegreeReceiver(scaling = 'std')")
        } else if (component == "outdegreeSender") {
          return("outdegreeSender(scaling = 'std')")
        } else if (component == "outdegreeReceiver") {
          return("outdegreeReceiver(scaling = 'std')")
        } else if (component == "totaldegreeDyad") {
          return("totaldegreeDyad(scaling = 'std')")
        } else if (component == "totaldegreeSender") {
          return("totaldegreeSender(scaling = 'std')")
        } else if (component == "totaldegreeReceiver") {
          return("totaldegreeReceiver(scaling = 'std')")
        } else if (component == "degreeMin") {
          return("degreeMin(scaling = 'std')")
        } else if (component == "degreeMax") {
          return("degreeMax(scaling = 'std')")
        } else if (component == "degreeDiff") {
          return("degreeDiff(scaling = 'std')")
        } else if (component == "sp") {
          return("sp(scaling = 'std')")
        } else if (component == "reciprocity") {
          return("reciprocity(scaling = 'std')")
        } else if (component == "otp") {
          return("otp(scaling = 'std')")
        } else if (component == "itp") {
          return("itp(scaling = 'std')")
        } else if (component == "osp") {
          return("osp(scaling = 'std')")
        } else if (component == "isp") {
          return("isp(scaling = 'std')")
        } else if (component == "psABBA") {
          return("psABBA()")
        } else if (component == "psABBY") {
          return("psABBY()")
        } else if (component == "psABXA") {
          return("psABXA()")
        } else if (component == "psABXB") {
          return("psABXB()")
        } else if (component == "psABXY") {
          return("psABXY()")
        } else if (component == "psABAY") {
          return("psABAY()")
        } else if (component == "psABAB") {
          return("psABAB()")
        } else if (component == "rrankSend") {
          return("rrankSend()")
        } else if (component == "rrankReceive") {
          return("rrankReceive()")
        } else if (component == "recencySendSender") {
          return("recencySendSender()")
        } else if (component == "recencySendReceiver") {
          return("recencySendReceiver()")
        } else if (component == "recencyReceiveSender") {
          return("recencyReceiveSender()")
        } else if (component == "recencyReceiveReceiver") {
          return("recencyReceiveReceiver()")
        } else if (component == "recencyContinue") {
          return("recencyContinue()")
        } else if (component == "FEtype") {
          return("FEtype()")
        } else {
          return("") # Handle or skip any unrecognized components
        }
      })
      
      # Add interaction terms by joining effects with `:`
      terms <- c(terms, paste(processed_components, collapse = " : "))
      
    } else if (name == "inertia") {
      terms <- c(terms, "inertia(scaling = 'std')")
    } else if (grepl("^difference_", name)) {
      var <- sub("^difference_", "", name)
      terms <- c(terms, paste0("difference('", var, "', scaling = 'std')"))
    } else if (grepl("^same_", name)) {
      var <- sub("^same_", "", name)
      terms <- c(terms, paste0("same('", var, "')"))
    } else if (grepl("^send_", name)) {
      var <- sub("^send_", "", name)
      terms <- c(terms, paste0("send('", var, "', scaling = 'std')"))
    } else if (grepl("^receive_", name)) {
      var <- sub("^receive_", "", name)
      terms <- c(terms, paste0("receive('", var, "', scaling = 'std')"))
    } else if (grepl("^tie_", name)) {
      var <- sub("^tie_", "", name)
      terms <- c(terms, paste0("tie('", var, "', scaling = 'std')"))
    } else if (grepl("^average_", name)) {
      var <- sub("^average_", "", name)
      terms <- c(terms, paste0("average('", var, "', scaling = 'std')"))
    } else if (grepl("^minimum_", name)) {
      var <- sub("^minimum_", "", name)
      terms <- c(terms, paste0("minimum('", var, "', scaling = 'std')"))
    } else if (grepl("^maximum_", name)) {
      var <- sub("^maximum_", "", name)
      terms <- c(terms, paste0("maximum('", var, "', scaling = 'std')"))
    } else if (grepl("^event_", name)) {
      var <- sub("^event_", "", name)
      terms <- c(terms, paste0("event('", var, "')"))
    } else if (grepl("^userStat_", name)) {
      var <- sub("^userStat_", "", name)
      terms <- c(terms, paste0("userStat('", var, "')"))
    } else if (name == "indegreeSender") {
      terms <- c(terms, "indegreeSender(scaling = 'std')")
    } else if (name == "indegreeReceiver") {
      terms <- c(terms, "indegreeReceiver(scaling = 'std')")
    } else if (name == "outdegreeSender") {
      terms <- c(terms, "outdegreeSender(scaling = 'std')")
    } else if (name == "outdegreeReceiver") {
      terms <- c(terms, "outdegreeReceiver(scaling = 'std')")
    } else if (name == "totaldegreeDyad") {
      terms <- c(terms, "totaldegreeDyad(scaling = 'std')")
    } else if (name == "totaldegreeSender") {
      terms <- c(terms, "totaldegreeSender(scaling = 'std')")
    } else if (name == "totaldegreeReceiver") {
      terms <- c(terms, "totaldegreeReceiver(scaling = 'std')")
    } else if (name == "degreeMin") {
      terms <- c(terms, "degreeMin(scaling = 'std')")
    } else if (name == "degreeMax") {
      terms <- c(terms, "degreeMax(scaling = 'std')")
    } else if (name == "degreeDiff") {
      terms <- c(terms, "degreeDiff(scaling = 'std')")
    } else if (name == "sp") {
      terms <- c(terms, "sp(scaling = 'std')")
    } else if (name == "reciprocity") {
      terms <- c(terms, "reciprocity(scaling = 'std')")
    } else if (name == "otp") {
      terms <- c(terms, "otp(scaling = 'std')")
    } else if (name == "itp") {
      terms <- c(terms, "itp(scaling = 'std')")
    } else if (name == "osp") {
      terms <- c(terms, "osp(scaling = 'std')")
    } else if (name == "isp") {
      terms <- c(terms, "isp(scaling = 'std')")
    } else if (name == "psABBA") {
      terms <- c(terms, "psABBA()")
    } else if (name == "psABBY") {
      terms <- c(terms, "psABBY()")
    } else if (name == "psABXA") {
      terms <- c(terms, "psABXA()")
    } else if (name == "psABXB") {
      terms <- c(terms, "psABXB()")
    } else if (name == "psABXY") {
      terms <- c(terms, "psABXY()")
    } else if (name == "psABAY") {
      terms <- c(terms, "psABAY()")
    } else if (name == "psABAB") {
      terms <- c(terms, "psABAB()")
    } else if (name == "rrankSend") {
      terms <- c(terms, "rrankSend()")
    } else if (name == "rrankReceive") {
      terms <- c(terms, "rrankReceive()")
    } else if (name == "recencySendSender") {
      terms <- c(terms, "recencySendSender()")
    } else if (name == "recencySendReceiver") {
      terms <- c(terms, "recencySendReceiver()")
    } else if (name == "recencyReceiveSender") {
      terms <- c(terms, "recencyReceiveSender()")
    } else if (name == "recencyReceiveReceiver") {
      terms <- c(terms, "recencyReceiveReceiver()")
    } else if (name == "recencyContinue") {
      terms <- c(terms, "recencyContinue()")
    } else if (name == "FEtype") {
      terms <- c(terms, "FEtype()")
    }
  }
  
  # Combine terms to formula
  formula_text <- paste("~", paste(terms, collapse = " + "))
  formula <- as.formula(formula_text)
  return(formula)
}


# Subsetting Edgelists


subset_edgelist_single <- function(edgelist, m) {
  return(edgelist[1:m, ])
}


# Estimating MLE and ABR


estimate_edgelists_parallel <- function(edgelists,
                                        parameters,
                                        covar,
                                        m,
                                        folder,
                                        num_cores,
                                        seed = 123) {
  subset_edgelist_single <- function(edgelist, m) {
    return(edgelist[1:m, ])
  }
  
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  set.seed(seed)
  
  # Precompute shared formula
  effects <- generate_formula(parameters)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(
        m_val = m_val,
        edgelist_index = edgelist_index,
        edgelist = subset_edgelist_single(edgelists[[edgelist_index]], m_val)
      )
    })
  })
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      # Remify the edgelist
      reh <- remify(task$edgelist, directed = TRUE, model = "tie")
      
      # Compute statistics
      statistics <- remstats(reh = reh,
                             tie_effects = effects,
                             attr_actors = covar)
      
      # Estimate the model using MLE
      fit <- remstimate::remstimate(
        reh = reh,
        stats = statistics,
        method = "MLE",
        timing = "interval"
      )
      
      # Remove statistics (for memory)
      rm(statistics)
      
      # Extract coefficients and set up for shrinkage
      coefs <- summary(fit)$coefsTab
      estimates <- coef(fit)
      cov <- fit$vcov
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      
      # Apply shrinkage methods
      output <- list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = coefs,
        shrink_hs = shrinkem(estimates, cov, type = "horseshoe")$estimates,
        shrink_ridge = shrinkem(estimates, cov, type = "ridge")$estimates
      )
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "01a_estimates_dependent/estimates_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(output, file = path)
      
    }, error = function(e) {
      dir.create("errors", showWarnings = FALSE)
      err_path <- paste0(
        folder,
        "06_errors/error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = NULL,
        shrink_hs = NULL,
        shrink_ridge = NULL
      )
    })
  }, future.seed = TRUE)
  
}


# Extra Function for Missing Models


estimate_edgelists_parallel_extra <- function(edgelists,
                                              parameters,
                                              covar,
                                              m,
                                              folder,
                                              num_cores,
                                              seed = 123,
                                              task_list) {
  subset_edgelist_single <- function(edgelist, m) {
    return(edgelist[1:m, ])
  }
  
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  set.seed(seed)
  
  # Precompute shared formula
  effects <- generate_formula(parameters)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      # Remify the edgelist
      reh <- remify(task$edgelist, directed = TRUE, model = "tie")
      
      # Compute statistics
      statistics <- remstats(reh = reh,
                             tie_effects = effects,
                             attr_actors = covar)
      
      # Estimate the model using MLE
      fit <- remstimate::remstimate(
        reh = reh,
        stats = statistics,
        method = "MLE",
        timing = "interval"
      )
      
      # Remove statistics (for memory)
      rm(list = c("statistics", "reh"))
      
      # Extract coefficients and set up for shrinkage
      coefs <- summary(fit)$coefsTab
      estimates <- coef(fit)
      cov <- fit$vcov
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      
      # Apply shrinkage methods
      output <- list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = coefs,
        shrink_hs = shrinkem(estimates, cov, type = "horseshoe"),
        shrink_ridge = shrinkem(estimates, cov, type = "ridge")
      )
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "01a_estimates_dependent/estimates_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(output, file = path)
      
    }, error = function(e) {
      dir.create("errors", showWarnings = FALSE)
      err_path <- paste0(
        folder,
        "06_errors/error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = NULL,
        shrink_hs = NULL,
        shrink_ridge = NULL
      )
    })
  }, future.seed = TRUE)
  
}


# Function for Independent edgelists


estimate_edgelists_independent <- function(parameters,
                                           covar,
                                           folder,
                                           num_cores,
                                           seed = 123,
                                           task_list) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  set.seed(seed)
  
  # Precompute shared formula
  effects <- generate_formula(parameters)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      # Remify the edgelist
      reh <- remify(task$edgelist, directed = TRUE, model = "tie")
      
      # Compute statistics
      statistics <- remstats(reh = reh,
                             tie_effects = effects,
                             attr_actors = covar)
      
      # Set seed
      set.seed(seed)
      
      # Estimate the model using MLE
      fit <- remstimate::remstimate(
        reh = reh,
        stats = statistics,
        method = "MLE",
        timing = "interval"
      )
      
      # Remove statistics (for memory)
      rm(list = c("statistics", "reh"))
      
      # Extract coefficients and set up for shrinkage
      coefs <- summary(fit)$coefsTab
      estimates <- coef(fit)
      cov <- fit$vcov
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      
      # Apply shrinkage methods
      output <- list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = coefs,
        shrink_hs = shrinkem(estimates, cov, type = "horseshoe"),
        shrink_ridge = shrinkem(estimates, cov, type = "ridge")
      )
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "01c_estimates_independent/estimates_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(output, file = path)
      
    }, error = function(e) {
      dir.create("errors", showWarnings = FALSE)
      err_path <- paste0(
        folder,
        "06_errors/error_independent_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = NULL,
        shrink_hs = NULL,
        shrink_ridge = NULL
      )
    })
  }, future.seed = TRUE)
  
}


# Rewrite Data to Poisson Format


poisson_df <- function(events, tie_stats, tie_reh, t0 = 0) {
  # Creating risk set
  risk_set <- vector("list", dim(tie_stats)[2])
  for (i in 1:dim(tie_stats)[2]) {
    risk_set[[i]] <- getDyad(tie_reh, i)
  }
  
  risk_set <- do.call(rbind, risk_set)[, 2:3]
  
  M <- nrow(events) # number of events
  poisson_list <- vector("list", M) # initialize list to store data frames
  
  for (m in 1:M) {
    # Calculate time difference between current and previous event
    t_curr <- events$time[m]
    t_prev <- if (m == 1)
      t0
    else
      events$time[m - 1]
    delta_m <- t_curr - t_prev
    
    # Get statistics for current event m
    stats_m <- as.data.frame(tie_stats[m, , ])
    
    # Combine risk set with covariates
    df_m <- cbind(risk_set, stats_m)
    
    # Create binary outcome y
    df_m$y <- ifelse(df_m$actor1 == events$sender[m] &
                       df_m$actor2 == events$receiver[m],
                     1,
                     0)
    
    # Add offset
    df_m$logDelta <- log(delta_m)
    
    # Store data frame for event m
    poisson_list[[m]] <- df_m
  }
  
  # Combine all event data frames in list into one data frame
  df_poisson <- do.call(rbind, poisson_list)
  df_poisson <- df_poisson %>%
    dplyr::select(-actor1, -actor2, -baseline) %>%
    # turn all variables that include "same_" to factors
    mutate(across(contains("same_"), as.factor)) %>%
    # turn y to integer
    mutate(y = as.integer(y))
  
  return(df_poisson)
}


# Estimating brms Models


## Horseshoe


estimate_brms_hs_parallel <- function(edgelists,
                                      parameters,
                                      covar,
                                      folder,
                                      num_cores,
                                      task_list,
                                      seed = 123) {
  # Set up parallel backend
  future::plan(future::multisession, workers = floor(num_cores / 4))
  set.seed(seed)
  
  # Precompute shared formula
  effects <- generate_formula(parameters)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      # Remify the edgelist
      reh <- remify(
        edgelist = task$edgelist,
        directed = TRUE,
        model = "tie"
      )
      
      # Compute statistics
      statistics <- remstats(reh = reh,
                             tie_effects = effects,
                             attr_actors = covar)
      
      # Transform data frame
      df_poisson <- poisson_df(task$edgelist,
                               tie_stats = statistics,
                               tie_reh = reh)
      
      
      # Create Formula
      predictors <- names(df_poisson)[1:(which(colnames(df_poisson) == "y") - 1)]
      glm_formula <- as.formula(paste("y ~", paste(
        c(predictors, "offset(logDelta)"), collapse = " + "
      )))
      
      # Remove statistics (for memory)
      rm(statistics)
      
      # Estimate the model using brm
      model <- brm(
        formula = glm_formula,
        data = df_poisson,
        family = poisson(link = "log"),
        prior = set_prior(
          horseshoe(
            df = 3,
            scale_global = 1,
            df_global = 3,
            scale_slab = 2,
            df_slab = 4,
            par_ratio = NULL,
            autoscale = TRUE
          ),
          class = "b"
        ),
        backend = "cmdstanr",
        cores = 4
      )
      
      coefs_ebr_hs <- bayestestR::describe_posterior(model, centrality = "all")
      
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "01b_estimates_ebr/estimates_brms_hs_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(coefs_ebr_hs, file = path)
      
      fls <- list.files(tempdir(), full.names = T)
      file.remove(fls)
      
    }, error = function(e) {
      dir.create("errors", showWarnings = FALSE)
      err_path <- paste0(
        folder,
        "02b_errors_brms/error_brms_hs_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(edgelist_index = task$edgelist_index,
           m = task$m_val)
    })
  }, future.seed = TRUE)
  
}


# Extra Function for Missing Models


estimate_brms_parallel_extra <- function(edgelists,
                                         parameters,
                                         covar,
                                         m,
                                         folder,
                                         num_cores,
                                         seed = 123,
                                         task_list) {
  # Set up parallel backend with appropriate number of workers
  future::plan(future::multisession, workers = floor(num_cores / 4))
  set.seed(seed)
  options(cmdstanr_output_dir = file.path(getwd(), "tmp"))
  
  # Precompute shared formula effects
  effects <- generate_formula(parameters)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      # Remify the edgelist (creates a remify object)
      reh <- remify(
        edgelist = task$edgelist,
        directed = TRUE,
        model = "tie"
      )
      
      # Compute network statistics based on the effects and covariates
      statistics <- remstats(reh = reh,
                             tie_effects = effects,
                             attr_actors = covar)
      
      # Transform data frame for Poisson regression
      df_poisson <- poisson_df(task$edgelist,
                               tie_stats = statistics,
                               tie_reh = reh)
      
      # Create the regression formula based on the predictors
      predictors <- names(df_poisson)[1:(which(colnames(df_poisson) == "y") - 1)]
      glm_formula <- as.formula(paste("y ~", paste(
        c(predictors, "offset(logDelta)"), collapse = " + "
      )))
      
      # Remove statistics (for memory)
      rm(statistics)
      
      # Estimate the model using brms with cmdstanr backend
      model <- brm(
        formula = glm_formula,
        data = df_poisson,
        family = poisson(link = "log"),
        prior = set_prior(
          horseshoe(
            df = 3,
            scale_global = 1,
            df_global = 3,
            scale_slab = 2,
            df_slab = 4,
            par_ratio = NULL,
            autoscale = TRUE
          ),
          class = "b"
        ),
        backend = "cmdstanr",
        cores = 4
      )
      
      coefs_ebr_hs <- bayestestR::describe_posterior(model, centrality = "all")
      
      # Save the model estimates into a designated file
      estimates_path <- paste0(
        folder,
        "01b_estimates_ebr/estimates_brms_hs_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(coefs_ebr_hs, file = estimates_path)
      
      # Clean up temporary files
      fls <- list.files(tempdir(), full.names = TRUE)
      file.remove(fls)
      
    }, error = function(e) {
      # Ensure error folder exists (with recursion in case nested folders are needed)
      dir.create(
        paste0(folder, "02b_errors_brms"),
        showWarnings = FALSE,
        recursive = TRUE
      )
      err_path <- paste0(
        folder,
        "02b_errors_brms/error_brms_hs_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(edgelist_index = task$edgelist_index,
           m = task$m_val)
    })
  }, future.seed = TRUE)
}


# Selecting Variables


select_variables_parallel <- function(m, n_edgelists, folder, num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:n_edgelists, function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          folder,
          "01a_estimates_dependent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      # Selected variables using p<0.05 in MLE
      mle_05 <- rownames(output[["mle_coefs"]][output[["mle_coefs"]][, 4] < 0.05, ])
      
      # Selected variables from ABR with Horseshoe
      abr_hs_95 <- rownames(output[["shrink_hs"]][output[["shrink_hs"]]$nonzero == 1, ])
      abr_hs_mode_01 <- rownames(output[["shrink_hs"]][abs(output[["shrink_hs"]]$shrunk.mode) >= 0.1, ])
      abr_hs_median_01 <- rownames(output[["shrink_hs"]][abs(output[["shrink_hs"]]$shrunk.median) >= 0.1, ])
      abr_hs_mean_01 <- rownames(output[["shrink_hs"]][abs(output[["shrink_hs"]]$shrunk.mean) >= 0.1, ])
      
      # Selected variables from ABR with Ridge
      abr_ridge_95 <- rownames(output[["shrink_ridge"]][output[["shrink_ridge"]]$nonzero == 1, ])
      abr_ridge_mode_01 <- rownames(output[["shrink_ridge"]][abs(output[["shrink_ridge"]]$shrunk.mode) >= 0.1, ])
      abr_ridge_median_01 <- rownames(output[["shrink_ridge"]][abs(output[["shrink_ridge"]]$shrunk.median) >= 0.1, ])
      abr_ridge_mean_01 <- rownames(output[["shrink_ridge"]][abs(output[["shrink_ridge"]]$shrunk.mean) >= 0.1, ])
      
      # Selected variables from EBR with Horseshoe only for m = 100 and 200
      if (task$m_val %in% c(100, 200)) {
        load(
          paste0(
            folder,
            "01b_estimates_ebr/estimates_brms_hs_m_",
            task$m_val,
            "_edgelist_",
            task$edgelist_index,
            ".RData"
          )
        )

        coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter))
        ebr_hs_95 <- coefs_ebr_hs[!(coefs_ebr_hs$CI_low < 0 &
                                  coefs_ebr_hs$CI_high > 0), ]$Parameter
        ebr_hs_mode_01 <- coefs_ebr_hs[abs(coefs_ebr_hs$MAP) >= 0.1, ]$Parameter
        ebr_hs_median_01 <- coefs_ebr_hs[abs(coefs_ebr_hs$Median) >= 0.1, ]$Parameter
        ebr_hs_mean_01 <- coefs_ebr_hs[abs(coefs_ebr_hs$Mean) >= 0.1, ]$Parameter
      } else {
        ebr_hs_95 <- character(0)
        ebr_hs_mode_01 <- character(0)
        ebr_hs_median_01 <- character(0)
        ebr_hs_mean_01 <- character(0)
      }
      
      selected_vars <- list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_05 = mle_05,
        
        abr_hs_95 = abr_hs_95,
        abr_hs_mode_01 = abr_hs_mode_01,
        abr_hs_median_01 = abr_hs_median_01,
        abr_hs_mean_01 = abr_hs_mean_01,
        
        abr_ridge_95 = abr_ridge_95,
        abr_ridge_mode_01 = abr_ridge_mode_01,
        abr_ridge_median_01 = abr_ridge_median_01,
        abr_ridge_mean_01 = abr_ridge_mean_01,
        
        ebr_hs_95 = ebr_hs_95,
        ebr_hs_mode_01 = ebr_hs_mode_01,
        ebr_hs_median_01 = ebr_hs_median_01,
        ebr_hs_mean_01 = ebr_hs_mean_01
      )
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "03_selected_variables/selection_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(selected_vars, file = path)
      
    }, error = function(e) {
      dir.create("errors", showWarnings = FALSE)
      err_path <- paste0(
        folder,
        "06_errors/selection_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = NULL,
        shrink_hs = NULL,
        shrink_ridge = NULL
      )
    })
  }, future.seed = TRUE)
  
}


# Selecting Variables from brms


select_variables_brms_parallel <- function(m, n_edgelists, folder, num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:n_edgelists, function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          folder,
          "01b_estimates_ebr/estimates_brms_hs_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )

      coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter))
      
      
      # Selected variables
      hs_95 <- coefs_ebr_hs[!(coefs_ebr_hs$CI_low < 0 &
                         coefs_ebr_hs$CI_high > 0), ]$Parameter
      hs_map_01 <- coefs_ebr_hs[abs(coefs_ebr_hs$MAP) >= 0.1, ]$Parameter
      hs_median_01 <- coefs_ebr_hs[abs(coefs_ebr_hs$Median) >= 0.1, ]$Parameter
      hs_mean_01 <- coefs_ebr_hs[abs(coefs_ebr_hs$Mean) >= 0.1, ]$Parameter
      
      selected_vars <- list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        
        hs_95 = hs_95,
        hs_map_01 = hs_map_01,
        hs_median_01 = hs_median_01,
        hs_mean_01 = hs_mean_01
      )
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "03b_selected_variables_brms/selection_brms_hs_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(selected_vars, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        folder,
        "02b_errors_brms/selection_brms_hs_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = NULL,
        shrink_hs = NULL
      )
    })
  }, future.seed = TRUE)
  
}



select_variables_independent_parallel <- function(m, n_edgelists, folder, num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:n_edgelists, function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          folder,
          "01c_estimates_independent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      # Selected variables using p<0.05 in MLE
      mle_05 <- rownames(output[["mle_coefs"]][output[["mle_coefs"]][, 4] < 0.05, ])
      
      # Selected variables from ABR with Horseshoe
      abr_hs_95 <- rownames(output[["shrink_hs"]][output[["shrink_hs"]]$nonzero == 1, ])
      abr_hs_mode_01 <- rownames(output[["shrink_hs"]][abs(output[["shrink_hs"]]$shrunk.mode) >= 0.1, ])
      abr_hs_median_01 <- rownames(output[["shrink_hs"]][abs(output[["shrink_hs"]]$shrunk.median) >= 0.1, ])
      abr_hs_mean_01 <- rownames(output[["shrink_hs"]][abs(output[["shrink_hs"]]$shrunk.mean) >= 0.1, ])
      
      # Selected variables from ABR with Ridge
      abr_ridge_95 <- rownames(output[["shrink_ridge"]][output[["shrink_ridge"]]$nonzero == 1, ])
      abr_ridge_mode_01 <- rownames(output[["shrink_ridge"]][abs(output[["shrink_ridge"]]$shrunk.mode) >= 0.1, ])
      abr_ridge_median_01 <- rownames(output[["shrink_ridge"]][abs(output[["shrink_ridge"]]$shrunk.median) >= 0.1, ])
      abr_ridge_mean_01 <- rownames(output[["shrink_ridge"]][abs(output[["shrink_ridge"]]$shrunk.mean) >= 0.1, ])
      
      selected_vars <- list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_05 = mle_05,
        
        abr_hs_95 = abr_hs_95,
        abr_hs_mode_01 = abr_hs_mode_01,
        abr_hs_median_01 = abr_hs_median_01,
        abr_hs_mean_01 = abr_hs_mean_01,
        
        abr_ridge_95 = abr_ridge_95,
        abr_ridge_mode_01 = abr_ridge_mode_01,
        abr_ridge_median_01 = abr_ridge_median_01,
        abr_ridge_mean_01 = abr_ridge_mean_01
      )
      
      # Save estimates in separate files
      path <- paste0(
        folder,
        "03b_selected_variables_independent/selection_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(selected_vars, file = path)
      
    }, error = function(e) {
      dir.create("errors", showWarnings = FALSE)
      err_path <- paste0(
        folder,
        "06_errors/selection_error_independent_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        edgelist_index = task$edgelist_index,
        m = task$m_val,
        mle_coefs = NULL,
        shrink_hs = NULL,
        shrink_ridge = NULL
      )
    })
  }, future.seed = TRUE)
  
}




extract_estimates <- function(m, n_edgelists, results_folder, num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:n_edgelists, function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  estimates <- future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          results_folder,
          "01a_estimates_dependent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      ncoefs <- length(output[["mle_coefs"]][, 1])
      df <- data.frame(matrix(nrow = ncoefs, ncol = 7))
      colnames(df) <- c("m",
                        "iteration",
                        "coef",
                        "mle",
                        "abr_hs",
                        "abr_ridge",
                        "ebr_hs")
      
      df$m <- task$m_val
      df$iteration <- task$edgelist_index
      
      df$coef <- names(output[["mle_coefs"]][, 1])
      df$mle <- output[["mle_coefs"]][, 1]
      df$abr_hs <- output[["shrink_hs"]]$shrunk.mode
      df$abr_ridge <- output[["shrink_ridge"]]$shrunk.mode
      
      if (task$m_val %in% c(100, 200)) {
        load(
          paste0(
            results_folder,
            "01b_estimates_ebr/estimates_brms_hs_m_",
            task$m_val,
            "_edgelist_",
            task$edgelist_index,
            ".RData"
          )
        )
        
        df$ebr_hs <- coefs_ebr_hs$MAP
      } else {
        df$ebr_hs <- NA
      }
      
      
      return(df)
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  df_estimates <- do.call(rbind, estimates)
  
  return(df_estimates)
}




extract_estimates_independent <- function(m, n_edgelists, results_folder, num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:n_edgelists, function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  estimates <- future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          results_folder,
          "01c_estimates_independent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      ncoefs <- length(output[["mle_coefs"]][, 1])
      df <- data.frame(matrix(nrow = ncoefs, ncol = 6))
      colnames(df) <- c("m", "iteration", "coef", "mle", "abr_hs", "abr_ridge")
      
      df$m <- task$m_val
      df$iteration <- task$edgelist_index
      
      df$coef <- names(output[["mle_coefs"]][, 1])
      df$mle <- output[["mle_coefs"]][, 1]
      df$abr_hs <- output[["shrink_hs"]]$shrunk.mode
      df$abr_ridge <- output[["shrink_ridge"]]$shrunk.mode
      
      
      return(df)
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  df_estimates <- do.call(rbind, estimates)
  
  return(df_estimates)
}




# Plot Estimates


plot_estimates <- function(df_estimates,
                           parameters,
                           m_val,
                           estimate = c("mle", "abr_hs", "abr_ridge"),
                           uncertainty = c("ci", "iqr")) {
  uncertainty <- match.arg(uncertainty)
  param_vec <- unlist(parameters)
  
  df_true <- data.frame(coef = names(param_vec), true_param = param_vec) %>%
    filter(coef != "baseline")
  
  if (uncertainty == "ci") {
    df_summarized <- df_estimates %>%
      filter(m == m_val, coef != "baseline") %>%
      group_by(coef = factor(coef, levels = names(param_vec))) %>%
      summarise(
        coef_mean = mean(get(estimate)),
        coef_sd = sd(get(estimate)),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        se = coef_sd / sqrt(n),
        lower_ci = coef_mean - qnorm(0.975) * se,
        upper_ci = coef_mean + qnorm(0.975) * se
      )
  } else if (uncertainty == "iqr") {
    df_summarized <- df_estimates %>%
      filter(m == m_val, coef != "baseline") %>%
      group_by(coef = factor(coef, levels = names(param_vec))) %>%
      summarise(
        coef_mean = mean(get(estimate)),
        lower_ci = quantile(get(estimate), 0.25),
        upper_ci = quantile(get(estimate), 0.75),
        .groups = "drop"
      )
  }
  
  df_plot <- left_join(df_summarized, df_true, by = "coef")
  df_plot$coef <- factor(df_plot$coef, levels = rev(names(param_vec)[names(param_vec) != "baseline"]))
  
  ggplot(df_plot, aes(y = coef)) +
    geom_point(aes(x = true_param, color = "True"),
               alpha = 0.5,
               size = 1.5) +
    geom_point(aes(x = coef_mean, color = "Estimate"),
               alpha = 0.5,
               size = 1.5) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci),
                   height = 0.2,
                   color = "black") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
    xlim(c(-6, 6)) +
    labs(title = paste0("M=", m_val),
         x = "Value",
         y = "") +
    theme_linedraw() +
    scale_color_manual(values = c("Estimate" = "black", "True" = "red")) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.title = element_blank(),
      title = element_text(size = 10)
    )
}





## Assess Bias


bias_estimates <- function(parameters,
                           df_estimates,
                           method = c("mle", "abr_hs", "abr_ridge", "ebr_hs")) {
  # Extract coefficients and settings
  all_coefficients <- names(parameters)
  num_iterations <- max(df_estimates$iteration)
  num_coefficients <- length(all_coefficients)
  
  # Extract true coefficients
  true <- unlist(parameters)
  names(true) <- names(parameters)
  
  # Add true parameter to df_estimates and calculate bias
  data <- df_estimates %>%
    # add true parameter with name(true) == df_estimates$coef
    mutate(true = true[match(df_estimates$coef, names(true))]) %>%
    # calculate bias for all methods in method (each method has its own column)
    mutate(across(all_of(method), ~ .x - true)) %>%
    # add "bias_" to the method names
    rename_with(~ paste0("bias_", .x), all_of(method))
  
  summary <- data %>%
    group_by(m) %>%
    # mean and var of bias for each method
    summarise(across(starts_with("bias_"), list(mean = mean, var = var), na.rm = TRUE))
  return(list(summary = summary, data = data))
}




# Discovery Rates

## Compute Discovery Rates


discovery_rate <- function(folder,
                           parameters,
                           m,
                           num_edgelists,
                           endo = TRUE,
                           tdr = TRUE,
                           effect_size = NULL,
                           num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  if (tdr) {
    # Get non-zero coefficients
    # If effect_size is provided, filter the parameters
    if (!is.null(effect_size)) {
      param_names <- names(parameters[parameters != 0 &
                                        parameters == effect_size])
    } else {
      param_names <- names(parameters[parameters != 0])
    }
  } else {
    # Get zero coefficients if tdr == FALSE
    param_names <- names(parameters[parameters == 0])
  }
  
  if (endo) {
    # Get endogenous predictors if endo == TRUE
    coefficients <- param_names[!grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  } else {
    # Get exogenous predictors if endo == FALSE
    coefficients <- param_names[grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  }
  
  
  # All coefficients (both nonzero and zero)
  all_coefficients <- names(parameters)
  num_coefficients <- length(all_coefficients)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:num_edgelists, function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  selection <- future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          folder,
          "selection_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      
      
      # Create the data frame
      df <- data.frame(
        m = rep(task$m_val, num_coefficients),
        iteration = rep(task$edgelist_index, num_coefficients),
        coef = all_coefficients
      )
      
      # Define methods
      methods <- setdiff(names(selected_vars), c("edgelist_index", "m"))
      
      # Add columns to data for each method
      for (i in methods) {
        df[[paste0("selected_", i)]] <- "not selected"
        df[[paste0("selected_", i)]][df$coef %in% selected_vars[[i]]] <- "selected"
        df[[paste0("selected_", i)]] <- factor(df[[paste0("selected_", i)]], levels = c("not selected", "selected"))
      }
      
      return(df)
      
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  # Combine list to data frame
  df_selection <- do.call(rbind, selection)
  
  # remove list for memory
  rm(selection)
  
  
  # Compute discovery rate
  df_discovery <- df_selection %>%
    filter(coef %in% coefficients) %>%
    pivot_longer(
      cols = starts_with("selected_"),
      names_to = "method",
      values_to = "selected"
    ) %>%
    mutate(method = sub("^selected_", "", method)) %>%
    group_by(m, iteration, method) %>%
    summarise(
      selected = sum(selected == "selected"),
      total = n(),
      discovery_rate = selected / total,
      .groups = 'drop'
    ) %>%
    dplyr::select(-selected, -total)
  
  discovery_rate <- df_discovery %>%
    # calculate variance of discovery rate by m and method
    group_by(m, method) %>% summarise(
      mean_discovery_rate = mean(discovery_rate),
      se_discovery_rate = sd(discovery_rate) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(
      lower = mean_discovery_rate - qnorm(0.975) * se_discovery_rate,
      upper = mean_discovery_rate + qnorm(0.975) * se_discovery_rate
    )
  
  return(
    list(
      discovery_rate = discovery_rate,
      data = df_discovery,
      selection_data = df_selection
    )
  )
}


### Independent


discovery_rate_independent <- function(folder,
                                       parameters,
                                       m,
                                       num_edgelists,
                                       endo = TRUE,
                                       tdr = TRUE,
                                       num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  if (tdr) {
    # Get non-zero coefficients if tdr == TRUE
    param_names <- names(parameters[parameters != 0])
  } else {
    # Get zero coefficients if tdr == FALSE
    param_names <- names(parameters[parameters == 0])
  }
  
  if (endo) {
    # Get endogenous predictors if endo == TRUE
    coefficients <- param_names[!grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  } else {
    # Get exogenous predictors if endo == FALSE
    coefficients <- param_names[grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  }
  
  
  # All coefficients (both nonzero and zero)
  all_coefficients <- names(parameters)
  num_coefficients <- length(all_coefficients)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(1:num_edgelists, function(edgelist_index) {
      list(
        m_val = m_val,
        edgelist_index = edgelist_index,
        edgelist = get(load(
          paste0(
            "../GeneratedData/edgelists/independent/",
            m_val,
            "/edgelist_",
            edgelist_index,
            ".RData"
          )
        )) %>%
          head(m_val)
      )
    })
  })
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  selection <- future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          folder,
          "03_selected_variables/selection_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      
      
      # Create the data frame
      df <- data.frame(
        m = rep(task$m_val, num_coefficients),
        iteration = rep(task$edgelist_index, num_coefficients),
        coef = all_coefficients
      )
      
      # Define methods
      methods <- setdiff(names(selected_vars), c("edgelist_index", "m"))
      
      # Add columns to data for each method
      for (i in methods) {
        df[[paste0("selected_", i)]] <- "not selected"
        df[[paste0("selected_", i)]][df$coef %in% selected_vars[[i]]] <- "selected"
        df[[paste0("selected_", i)]] <- factor(df[[paste0("selected_", i)]], levels = c("not selected", "selected"))
      }
      
      return(df)
      
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  # Combine list to data frame
  df_selection <- do.call(rbind, selection)
  
  # remove list for memory
  rm(selection)
  
  
  # Compute discovery rate
  discovery_rate <- df_selection %>%
    filter(coef %in% coefficients) %>%
    pivot_longer(
      cols = starts_with("selected_"),
      names_to = "method",
      values_to = "selected"
    ) %>%
    mutate(method = sub("^selected_", "", method)) %>%
    group_by(m, iteration, method) %>% summarise(
      selected = sum(selected == "selected"),
      total = n(),
      discovery_rate = selected / total,
      .groups = 'drop'
    ) %>%
    # calculate variance of discovery rate by m and method
    group_by(m, method) %>% summarise(
      mean_discovery_rate = mean(discovery_rate),
      se_discovery_rate = sd(discovery_rate) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(
      lower = mean_discovery_rate - qnorm(0.975) * se_discovery_rate,
      upper = mean_discovery_rate + qnorm(0.975) * se_discovery_rate
    )
  
  return(list(discovery_rate = discovery_rate, data = df_selection))
}


### brms


discovery_rate_brms <- function(folder,
                                parameters,
                                m,
                                num_edgelists,
                                endo = TRUE,
                                tdr = TRUE,
                                effect_size = NULL,
                                num_cores) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  if (tdr) {
    # Get non-zero coefficients
    # If effect_size is provided, filter the parameters
    if (!is.null(effect_size)) {
      param_names <- names(parameters[parameters != 0 &
                                        parameters == effect_size])
    } else {
      param_names <- names(parameters[parameters != 0])
    }
  } else {
    # Get zero coefficients if tdr == FALSE
    param_names <- names(parameters[parameters == 0])
  }
  
  if (endo) {
    # Get endogenous predictors if endo == TRUE
    coefficients <- param_names[!grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  } else {
    # Get exogenous predictors if endo == FALSE
    coefficients <- param_names[grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  }
  
  
  # All coefficients (both nonzero and zero)
  all_coefficients <- names(parameters)
  num_coefficients <- length(all_coefficients)
  
  # Flatten list of tasks
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(
        m_val = m_val,
        edgelist_index = edgelist_index,
        edgelist = subset_edgelist_single(edgelists[[edgelist_index]], m_val)
      )
    })
  })
  task_list <- unlist(task_list, recursive = FALSE)
  
  # Run tasks in parallel
  selection <- future_lapply(task_list, function(task) {
    tryCatch({
      load(
        paste0(
          folder,
          "03b_selected_variables_brms/selection_brms_hs_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      
      
      # Create the data frame
      df <- data.frame(
        m = rep(task$m_val, num_coefficients),
        iteration = rep(task$edgelist_index, num_coefficients),
        coef = all_coefficients
      )
      
      # Define methods
      methods <- setdiff(names(selected_vars), c("edgelist_index", "m"))
      
      # Add columns to data for each method
      for (i in methods) {
        df[[paste0("selected_", i)]] <- "not selected"
        df[[paste0("selected_", i)]][df$coef %in% selected_vars[[i]]] <- "selected"
        df[[paste0("selected_", i)]] <- factor(df[[paste0("selected_", i)]], levels = c("not selected", "selected"))
      }
      
      return(df)
      
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  # Combine list to data frame
  df_selection <- do.call(rbind, selection)
  
  # remove list for memory
  rm(selection)
  
  
  # Compute discovery rate
  discovery_rate <- df_selection %>%
    filter(coef %in% coefficients) %>%
    pivot_longer(
      cols = starts_with("selected_"),
      names_to = "method",
      values_to = "selected"
    ) %>%
    mutate(method = sub("^selected_", "", method)) %>%
    group_by(m, iteration, method) %>% summarise(
      selected = sum(selected == "selected"),
      total = n(),
      discovery_rate = selected / total,
      .groups = 'drop'
    ) %>%
    # calculate variance of discovery rate by m and method
    group_by(m, method) %>% summarise(
      mean_discovery_rate = mean(discovery_rate),
      se_discovery_rate = sd(discovery_rate) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(
      lower = mean_discovery_rate - qnorm(0.975) * se_discovery_rate,
      upper = mean_discovery_rate + qnorm(0.975) * se_discovery_rate
    )
  
  return(list(discovery_rate = discovery_rate, data = df_selection))
}


## Plot Discovery Rates



plot_discovery_rate <- function(data, criteria, title = "") {
  legend_labels <- c(
    "mle_05"              = "\u03B1 = 0.05 (MLE)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95%-CI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)"
  )
  
  method_colors <- c(
    "mle_05"             = "#0072B2",
    
    # ABR HS methods
    "abr_hs_95"          = "#FFA500",
    "abr_hs_mode_01"     = "#FF4500",
    "abr_hs_median_01"   = "#FFA07A",
    "abr_hs_mean_01"     = "#FFD700",
    
    # ABR Ridge methods
    "abr_ridge_95"       = "#3CB371",
    "abr_ridge_mode_01"  = "#AFEEEE",
    "abr_ridge_median_01" = "#20B2AA",
    "abr_ridge_mean_01"  = "#40E0D0",
    
    # EBR HS
    "ebr_hs_95" = "#DF00FF",
    "ebr_hs_mode_01" = "#D88FD8",
    "ebr_hs_mean_01" = "#C54B8C",
    "ebr_hs_median_01" = "#5D3F6A"
  )
  
  
  
  # different line types for each method
  method_linetypes <- c(
    "mle_05"             = "solid",
    "abr_hs_95"          = "solid",
    "abr_hs_mode_01"     = "dotted",
    "abr_hs_median_01"   = "dashed",
    "abr_hs_mean_01"     = "dotdash",
    "abr_ridge_95"       = "solid",
    "abr_ridge_mode_01"  = "dotted",
    "abr_ridge_median_01" = "dashed",
    "abr_ridge_mean_01"  = "dotdash",
    "ebr_hs_95" = "solid",
    "ebr_hs_mode_01" = "dotted",
    "ebr_hs_mean_01" = "dashed",
    "ebr_hs_median_01" = "dotdash"
  )
  
  
  data %>%
    dplyr::filter(method %in% criteria) %>%
    ggplot(aes(x = m, y = mean_discovery_rate, group = method)) +
    
    # Confidence interval ribbon
    geom_ribbon(aes(
      ymin = lower,
      ymax = upper,
      fill = method
    ),
    alpha = 0.1,
    color = NA) +
    
    # Mean discovery rate as line
    geom_line(aes(color = method, linetype = method), size = 1) +
    geom_point(aes(color = method)) +
    
    # Logarithmic scaling
    scale_x_log10(breaks = c(100, 200, 400, 800, 1600, 3200, 6400, 12800)) +
    
    # Distinct color, fill, and linetype for each method
    scale_color_manual(values = method_colors, labels = legend_labels) +
    scale_fill_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    
    # Adjust the legend position and title
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    
    # Axis and title labels
    labs(
      title = title,
      x = "Events (M)",
      y = "Discovery Rate",
      color = "Selection Criterion",
      fill = "Selection Criterion",
      linetype = "Selection Criterion"
    )
}




## Plot Distance Metric


plot_distance_metric <- function(data, criteria, title = "") {
  legend_labels <- c(
    "mle_05"              = "\u03B1 = 0.05 (MLE)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95%-CI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)"
  )
  
  method_colors <- c(
    "mle_05"             = "#0072B2",
    
    # ABR HS methods
    "abr_hs_95"          = "#FFA500",
    "abr_hs_mode_01"     = "#FF4500",
    "abr_hs_median_01"   = "#FFA07A",
    "abr_hs_mean_01"     = "#FFD700",
    
    # ABR Ridge methods
    "abr_ridge_95"       = "#3CB371",
    "abr_ridge_mode_01"  = "#AFEEEE",
    "abr_ridge_median_01" = "#20B2AA",
    "abr_ridge_mean_01"  = "#40E0D0",
    
    # EBR HS
    "ebr_hs_95" = "#DF00FF",
    "ebr_hs_mode_01" = "#D88FD8",
    "ebr_hs_mean_01" = "#C54B8C",
    "ebr_hs_median_01" = "#5D3F6A"
  )
  
  
  
  # different line types for each method
  method_linetypes <- c(
    "mle_05"             = "solid",
    "abr_hs_95"          = "solid",
    "abr_hs_mode_01"     = "dotted",
    "abr_hs_median_01"   = "dashed",
    "abr_hs_mean_01"     = "dotdash",
    "abr_ridge_95"       = "solid",
    "abr_ridge_mode_01"  = "dotted",
    "abr_ridge_median_01" = "dashed",
    "abr_ridge_mean_01"  = "dotdash",
    "ebr_hs_95" = "solid",
    "ebr_hs_mode_01" = "dotted",
    "ebr_hs_mean_01" = "dashed",
    "ebr_hs_median_01" = "dotdash"
  )
  
  # Filter data and ensure the factor levels of 'method' match our legend labels.
  data %>%
    dplyr::filter(method %in% criteria) %>%
    dplyr::mutate(method = factor(method, levels = names(legend_labels))) %>%
    ggplot(aes(x = m, y = mean_distance, group = method)) +
    
    # Plot the distance metric as lines and add points for clarity.
    geom_line(aes(color = method, linetype = method), size = 1) +
    geom_point(aes(color = method), size = 1.5) +
    
    # Use a log10 scale for x-axis to minimize overlap in the event counts.
    scale_x_log10(breaks = c(100, 200, 400, 800, 1600, 3200, 6400, 12800)) +
    
    # Apply custom manual scales for color and linetype with our labels.
    scale_color_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    
    # Use minimal theme adjustments and position legend at the bottom.
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    
    # Label the axes and title.
    labs(
      title = title,
      x = "Events (M)",
      y = "Distance Metric",
      color = "Selection Criterion",
      linetype = "Selection Criterion"
    )
}



# Matthews' Correlation Coefficient


mcc <- function(data,
                parameters,
                endo = TRUE,
                selection_label = "selected") {
  param_names <- names(parameters)
  if (endo) {
    # Get endogenous predictors if endo == TRUE
    coefficients <- param_names[!grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  } else {
    # Get exogenous predictors if endo == FALSE
    coefficients <- param_names[grepl("_", param_names) &
                                  !(param_names %in% c("baseline", "Intercept"))]
  }
  
  selection_cols <- grep("^selected", names(data), value = TRUE)
  results <- data %>%
    dplyr::filter(coef %in% coefficients) %>%
    group_by(iteration, m) %>%
    summarize(across(all_of(selection_cols), ~ {
      pred_class <- ifelse(. == selection_label, 1, 0)
      true_class <- ifelse(parameters[as.character(coef)] != 0, 1, 0)
      
      # Calculate the confusion matrix counts.
      tp <- sum(pred_class == 1 &
                  true_class == 1)  # true positives
      tn <- sum(pred_class == 0 &
                  true_class == 0)  # true negatives
      fp <- sum(pred_class == 1 &
                  true_class == 0)  # false positives
      fn <- sum(pred_class == 0 &
                  true_class == 1)  # false negatives
      
      # Compute the numerator and denominator for MCC.
      numerator <- tp * tn - fp * fn
      denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      
      # Check if denominator is NA or zero, and return NA if so.
      if (is.na(denominator) || denominator == 0) {
        return(NA)
      } else {
        return(numerator / denominator)
      }
    }, .names = "mcc_{.col}")) %>%
    ungroup() %>%
    group_by(m) %>%
    # calculate mean and sd of mcc with na.rm = TRUE
    summarise(across(starts_with("mcc_"), list(
      mean = ~ mean(.x, na.rm = TRUE),
      se = ~ sd(.x, na.rm = TRUE) / sqrt(n())
    ))) %>%
    pivot_longer(
      cols = starts_with("mcc_selected_"),
      names_to = c("method", "stat"),
      names_pattern = "mcc_selected_(.*)_(mean|se)",
      values_to = "value"
    ) %>%
    # Reformat the data so that mean and se are in separate columns.
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(lower = mean - qnorm(0.975) * se,
           upper = mean + qnorm(0.975) * se) %>%
    # Select and order the desired columns: m, method, mean, sd, lower, upper.
    dplyr::select(m, method, mean, se, lower, upper)
  
  return(results)
}


## Plot MCC


plot_mcc <- function(data, criteria, title = "") {
  legend_labels <- c(
    "mle_05"              = "\u03B1 = 0.05 (MLE)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95%-CI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)"
  )
  
  method_colors <- c(
    "mle_05"             = "#0072B2",
    
    # ABR HS methods
    "abr_hs_95"          = "#FFA500",
    "abr_hs_mode_01"     = "#FF4500",
    "abr_hs_median_01"   = "#FFA07A",
    "abr_hs_mean_01"     = "#FFD700",
    
    # ABR Ridge methods
    "abr_ridge_95"       = "#3CB371",
    "abr_ridge_mode_01"  = "#AFEEEE",
    "abr_ridge_median_01" = "#20B2AA",
    "abr_ridge_mean_01"  = "#40E0D0",
    
    # EBR HS
    "ebr_hs_95" = "#DF00FF",
    "ebr_hs_mode_01" = "#D88FD8",
    "ebr_hs_mean_01" = "#C54B8C",
    "ebr_hs_median_01" = "#5D3F6A"
  )
  
  
  
  # different line types for each method
  method_linetypes <- c(
    "mle_05"             = "solid",
    "abr_hs_95"          = "solid",
    "abr_hs_mode_01"     = "dotted",
    "abr_hs_median_01"   = "dashed",
    "abr_hs_mean_01"     = "dotdash",
    "abr_ridge_95"       = "solid",
    "abr_ridge_mode_01"  = "dotted",
    "abr_ridge_median_01" = "dashed",
    "abr_ridge_mean_01"  = "dotdash",
    "ebr_hs_95" = "solid",
    "ebr_hs_mode_01" = "dotted",
    "ebr_hs_mean_01" = "dashed",
    "ebr_hs_median_01" = "dotdash"
  )
  
  data %>%
    dplyr::filter(method %in% criteria) %>%
    ggplot(aes(x = m, y = mean, group = method)) +
    
    # Mean discovery rate as line
    geom_line(aes(color = method, linetype = method), size = 1) +
    
    geom_point(aes(color = method), size = 1.5) +
    
    # Logarithmic scaling
    scale_x_log10(breaks = c(100, 200, 400, 800, 1600, 3200, 6400, 12800)) +
    
    # Distinct color, fill, and linetype for each method
    scale_color_manual(values = method_colors, labels = legend_labels) +
    scale_fill_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    
    # Adjust the legend position and title
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    
    # Axis and title labels
    labs(
      title = title,
      x = "Events (M)",
      y = "MCC",
      color = "Selection Criterion",
      fill = "Selection Criterion",
      linetype = "Selection Criterion"
    )
}


# Predictive Performance


predictive_performance_all <- function(edgelist,
                                       covar,
                                       coefficients,
                                       reh,
                                       statistics,
                                       quantile = 0.95,
                                       warnings = TRUE) {
  top5 <- rep(FALSE, nrow(edgelist))
  
  if (!warnings) {
    suppressWarnings({
      for (i in 1:nrow(edgelist)) {
        # Calculate event rates for every actor in R
        lambda <- exp(as.numeric(statistics[i, , ] %*% coefficients))
        
        # Extract top 5% of lambdas row indices
        top_5_lambda <- which(lambda > quantile(lambda, quantile))
        # extract specific dyad that belongs to this row
        top_5_dyads <- getDyad(x = reh, dyadID = top_5_lambda)[, 2:3]
        
        actual <- edgelist[i, c(2, 3)]
        
        # Check if actual is part of top_5_dyads
        top5[i] <- any(top_5_dyads$actor1 == actual$sender &
                         top_5_dyads$actor2 == actual$receiver)
        
      }
    })
  } else {
    for (i in 1:nrow(edgelist)) {
      lambda <- exp(as.numeric(statistics[i, , ] %*% coefficients))
      
      # Extract top 5% of lambdas row indices
      top_5_lambda <- which(lambda > quantile(lambda, quantile))
      # extract specific dyad that belongs to this row
      top_5_dyads <- getDyad(x = reh, dyadID = top_5_lambda)[, 2:3]
      
      actual <- edgelist[i, c(2, 3)]
      
      # Check if actual is part of top_5_dyads
      top5[i] <- any(top_5_dyads$actor1 == actual$sender &
                       top_5_dyads$actor2 == actual$receiver)
    }
  }
  
  return(top5 = sum(top5) / nrow(edgelist))
}


### 3.3.2 In-Sample Predictive Performance for All Edgelists

### All Coefficients


pp_is_parallel <- function(edgelists,
                           m,
                           parameters,
                           covar,
                           results_folder,
                           statistics_folder,
                           output_folder,
                           num_cores,
                           quantile = 0.95) {
  future::plan(future::multisession, workers = num_cores)
  
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        edgelists[[task$edgelist_index]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01a_estimates_dependent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          statistics_folder,
          "statistics_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      mle_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile
      )
      
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      abr_hs_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_hs$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      if (task$m_val %in% c(100, 200)) {
        load(
          paste0(
            results_folder,
            "01b_estimates_ebr/estimates_brms_hs_m_",
            task$m_val,
            "_edgelist_",
            task$edgelist_index,
            ".RData"
          )
        )
        
        ebr_hs_pp <- predictive_performance_all(
          edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
          covar = covar,
          coefficients = coefs_ebr_hs$MAP,
          reh = reh,
          statistics = statistics[1:task$m_val, , ],
          quantile = quantile,
          warnings = FALSE
        )
      } else {
        ebr_hs_pp <- NA
      }
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp,
        ebr_hs = ebr_hs_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_ridge = NA,
        abr_hs = NA,
        ebr_hs = NA
      )
    })
  }, future.seed = TRUE)
  
}



### Sparse Models


pp_is_sparse_parallel <- function(edgelists,
                                  m,
                                  parameters,
                                  covar,
                                  results_folder,
                                  statistics_folder,
                                  selected_folder,
                                  output_folder,
                                  num_cores,
                                  quantile = 0.95) {
  future::plan(future::multisession, workers = num_cores)
  
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        edgelists[[task$edgelist_index]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01a_estimates_dependent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          statistics_folder,
          "statistics_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          selected_folder,
          "selection_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      statistics_is <- statistics[1:task$m_val, , ]
      
      rm(statistics)
      
      mle_05_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output[["mle_coefs"]][, 1],!rownames(output[["mle_coefs"]]) %in% c("baseline", selected_vars$mle_05),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_95_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_95),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_mode_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_median_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_mean_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      if (task$m_val %in% c(100, 200)) {
        load(
          paste0(
            results_folder,
            "01b_estimates_ebr/estimates_brms_hs_m_",
            task$m_val,
            "_edgelist_",
            task$edgelist_index,
            ".RData"
          )
        )

        coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter))
        
        ebr_hs_95_pp <- predictive_performance_all(
          edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_95),
            0
          ),
          reh = reh,
          statistics = statistics_is,
          quantile = quantile
        )
        
        ebr_hs_mode_01_pp <- predictive_performance_all(
          edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_mode_01),
            0
          ),
          reh = reh,
          statistics = statistics_is,
          quantile = quantile
        )
        
        ebr_hs_median_01_pp <- predictive_performance_all(
          edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_median_01),
            0
          ),
          reh = reh,
          statistics = statistics_is,
          quantile = quantile
        )
        
        ebr_hs_mean_01_pp <- predictive_performance_all(
          edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_mean_01),
            0
          ),
          reh = reh,
          statistics = statistics_is,
          quantile = quantile
        )
      } else {
        ebr_hs_95_pp <- NA
        ebr_hs_mode_01_pp <- NA
        ebr_hs_median_01_pp <- NA
        ebr_hs_mean_01_pp <- NA
      }
      
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_is,
        quantile = quantile,
        warnings = FALSE
      )
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle_05 = mle_05_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs_95 = abr_hs_95_pp,
        abr_hs_mode_01 = abr_hs_mode_01_pp,
        abr_hs_median_01 = abr_hs_median_01_pp,
        abr_hs_mean_01 = abr_hs_mean_01_pp,
        abr_ridge_95 = abr_ridge_95_pp,
        abr_ridge_mode_01 = abr_ridge_mode_01_pp,
        abr_ridge_median_01 = abr_ridge_median_01_pp,
        abr_ridge_mean_01 = abr_ridge_mean_01_pp,
        ebr_hs_95 = ebr_hs_95_pp,
        ebr_hs_mode_01 = ebr_hs_mode_01_pp,
        ebr_hs_median_01 = ebr_hs_median_01_pp,
        ebr_hs_mean_01 = ebr_hs_mean_01_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_sparse_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_ridge = NA,
        abr_hs_95 = NA,
        abr_hs_mode_01 = NA,
        abr_hs_median_01 = NA,
        abr_hs_mean_01 = NA,
        abr_ridge_95 = NA,
        abr_ridge_mode_01 = NA,
        abr_ridge_median_01 = NA,
        abr_ridge_mean_01 = NA,
        ebr_hs_95 = NA,
        ebr_hs_mode_01 = NA,
        ebr_hs_median_01 = NA,
        ebr_hs_mean_01 = NA
      )
    })
  }, future.seed = TRUE)
  
}



#### Independent


pp_is_independent_parallel <- function(parameters,
                                       covar,
                                       results_folder,
                                       output_folder,
                                       num_cores,
                                       quantile = 0.95,
                                       task_list) {
  future::plan(future::multisession, workers = num_cores)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01c_estimates_independent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      
      mle_pp <- predictive_performance_all(
        edgelist = task$edgelist[1:task$m_val, ],
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile
      )
      
      abr_ridge_pp <- predictive_performance_all(
        edgelist = task$edgelist[1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      abr_hs_pp <- predictive_performance_all(
        edgelist = task$edgelist[1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_hs$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_independent_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_ridge = NA,
        abr_hs = NA
      )
    })
  }, future.seed = TRUE)
  
}


### Sparse Models


pp_is_independent_sparse_parallel <- function(parameters,
                                              covar,
                                              results_folder,
                                              selected_folder,
                                              output_folder,
                                              num_cores,
                                              quantile = 0.95,
                                              task_list) {
  future::plan(future::multisession, workers = num_cores)
  
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      statistics_is <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      
      load(
        paste0(
          results_folder,
          "01c_estimates_independent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          selected_folder,
          "selection_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      mle_05_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output[["mle_coefs"]][, 1],!rownames(output[["mle_coefs"]]) %in% c("baseline", selected_vars$mle_05),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_95_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_95),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_mode_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_median_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      abr_ridge_mean_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle_05 = mle_05_pp,
        abr_hs_95 = abr_hs_95_pp,
        abr_hs_mode_01 = abr_hs_mode_01_pp,
        abr_hs_median_01 = abr_hs_median_01_pp,
        abr_hs_mean_01 = abr_hs_mean_01_pp,
        abr_ridge_95 = abr_ridge_95_pp,
        abr_ridge_mode_01 = abr_ridge_mode_01_pp,
        abr_ridge_median_01 = abr_ridge_median_01_pp,
        abr_ridge_mean_01 = abr_ridge_mean_01_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_independent_sparse_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_hs_95 = NA,
        abr_hs_mode_01 = NA,
        abr_hs_median_01 = NA,
        abr_hs_mean_01 = NA,
        abr_ridge_95 = NA,
        abr_ridge_mode_01 = NA,
        abr_ridge_median_01 = NA,
        abr_ridge_mean_01 = NA
      )
    })
  }, future.seed = TRUE)
  
}



### 3.3.3 Out-of-Sample Predictive Performance for All Edgelists


pp_oos_parallel <- function(edgelists,
                            m,
                            parameters,
                            covar,
                            results_folder,
                            statistics_folder,
                            output_folder,
                            num_cores,
                            quantile = 0.95,
                            new = 1000) {
  future::plan(future::multisession, workers = num_cores)
  
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        edgelists[[task$edgelist_index]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01a_estimates_dependent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          statistics_folder,
          "statistics_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      edgelist_oos <- edgelists[[task$edgelist_index]][(task$m_val + 1):(task$m_val + new), ]
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ]
      
      mle_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      abr_hs_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_hs$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      if (task$m_val %in% c(100, 200)) {
        load(
          paste0(
            results_folder,
            "01b_estimates_ebr/estimates_brms_hs_m_",
            task$m_val,
            "_edgelist_",
            task$edgelist_index,
            ".RData"
          )
        )

        
        ebr_hs_pp <- predictive_performance_all(
          edgelist = edgelist_oos,
          covar = covar,
          coefficients = coefs_ebr_hs$MAP,
          reh = reh,
          statistics = statistics_oos,
          quantile = quantile,
          warnings = FALSE
        )
      } else {
        ebr_hs_pp <- NA
      }
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp,
        ebr_hs = ebr_hs_pp
      )
      
      path <- paste0(
        output_folder,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_oos_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_ridge = NA,
        abr_hs = NA,
        ebr_hs = NA
      )
    })
  }, future.seed = TRUE)
}


### Sparse Models


pp_oos_sparse_parallel <- function(edgelists,
                                   m,
                                   parameters,
                                   covar,
                                   results_folder,
                                   statistics_folder,
                                   selected_folder,
                                   output_folder,
                                   num_cores,
                                   quantile = 0.95,
                                   new = 1000) {
  future::plan(future::multisession, workers = num_cores)
  
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        edgelists[[task$edgelist_index]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01a_estimates_dependent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          statistics_folder,
          "statistics_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          selected_folder,
          "selection_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      edgelist_oos <- edgelists[[task$edgelist_index]][(task$m_val + 1):(task$m_val + new), ]
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ]
      rm(statistics)
      
      
      mle_05_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output[["mle_coefs"]][, 1],!rownames(output[["mle_coefs"]]) %in% c("baseline", selected_vars$mle_05),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_95_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_95),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_mode_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_median_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_mean_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      if (task$m_val %in% c(100, 200)) {
        load(
          paste0(
            results_folder,
            "01b_estimates_ebr/estimates_brms_hs_m_",
            task$m_val,
            "_edgelist_",
            task$edgelist_index,
            ".RData"
          )
        )

        coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter))
        
        ebr_hs_95_pp <- predictive_performance_all(
          edgelist = edgelist_oos,
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_95),
            0
          ),
          reh = reh,
          statistics = statistics_oos,
          quantile = quantile
        )
        
        ebr_hs_mode_01_pp <- predictive_performance_all(
          edgelist = edgelist_oos,
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_mode_01),
            0
          ),
          reh = reh,
          statistics = statistics_oos,
          quantile = quantile
        )
        
        ebr_hs_median_01_pp <- predictive_performance_all(
          edgelist = edgelist_oos,
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_median_01),
            0
          ),
          reh = reh,
          statistics = statistics_oos,
          quantile = quantile
        )
        
        ebr_hs_mean_01_pp <- predictive_performance_all(
          edgelist = edgelist_oos,
          covar = covar,
          coefficients = replace(
            coefs_ebr_hs$MAP,!coefs_ebr_hs$Parameter %in% c("Intercept", selected_vars$ebr_hs_mean_01),
            0
          ),
          reh = reh,
          statistics = statistics_oos,
          quantile = quantile
        )
      } else {
        ebr_hs_95_pp <- NA
        ebr_hs_mode_01_pp <- NA
        ebr_hs_median_01_pp <- NA
        ebr_hs_mean_01_pp <- NA
      }
      
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle_05 = mle_05_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs_95 = abr_hs_95_pp,
        abr_hs_mode_01 = abr_hs_mode_01_pp,
        abr_hs_median_01 = abr_hs_median_01_pp,
        abr_hs_mean_01 = abr_hs_mean_01_pp,
        abr_ridge_95 = abr_ridge_95_pp,
        abr_ridge_mode_01 = abr_ridge_mode_01_pp,
        abr_ridge_median_01 = abr_ridge_median_01_pp,
        abr_ridge_mean_01 = abr_ridge_mean_01_pp,
        ebr_hs_95 = ebr_hs_95_pp,
        ebr_hs_mode_01 = ebr_hs_mode_01_pp,
        ebr_hs_median_01 = ebr_hs_median_01_pp,
        ebr_hs_mean_01 = ebr_hs_mean_01_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_sparse_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_ridge = NA,
        abr_hs_95 = NA,
        abr_hs_mode_01 = NA,
        abr_hs_median_01 = NA,
        abr_hs_mean_01 = NA,
        abr_ridge_95 = NA,
        abr_ridge_mode_01 = NA,
        abr_ridge_median_01 = NA,
        abr_ridge_mean_01 = NA,
        ebr_hs_95 = NA,
        ebr_hs_mode_01 = NA,
        ebr_hs_median_01 = NA,
        ebr_hs_mean_01 = NA
      )
    })
  }, future.seed = TRUE)
  
}


#### Independent


pp_oos_independent_parallel <- function(parameters,
                                        covar,
                                        results_folder,
                                        output_folder,
                                        num_cores,
                                        quantile = 0.95,
                                        task_list,
                                        new = 1000) {
  future::plan(future::multisession, workers = num_cores)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01c_estimates_independent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      edgelist_oos <- task$edgelist[(task$m_val + 1):(task$m_val + new), ]
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ]
      
      rm(statistics)
      
      
      mle_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      abr_hs_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_hs$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_independent_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_ridge = NA,
        abr_hs = NA
      )
    })
  }, future.seed = TRUE)
  
}


### Sparse Models


pp_oos_independent_sparse_parallel <- function(parameters,
                                               covar,
                                               results_folder,
                                               statistics_folder,
                                               selected_folder,
                                               output_folder,
                                               num_cores,
                                               quantile = 0.95,
                                               new = 1000,
                                               task_list) {
  future::plan(future::multisession, workers = num_cores)
  
  future_lapply(task_list, function(task) {
    tryCatch({
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      load(
        paste0(
          results_folder,
          "01c_estimates_independent/estimates_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      load(
        paste0(
          selected_folder,
          "selection_m_",
          task$m_val,
          "_edgelist_",
          task$edgelist_index,
          ".RData"
        )
      )
      
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      edgelist_oos <- task$edgelist[(task$m_val + 1):(task$m_val + new), ]
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ]
      
      rm(statistics)
      
      
      mle_05_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output[["mle_coefs"]][, 1],!rownames(output[["mle_coefs"]]) %in% c("baseline", selected_vars$mle_05),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$estimates$shrunk.mode,!rownames(output$shrink_hs$estimates) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_95_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_95),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_mode_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_median_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      abr_ridge_mean_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_ridge$estimates$shrunk.mode,!rownames(output$shrink_ridge$estimates) %in% c("baseline", selected_vars$abr_ridge_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle_05 = mle_05_pp,
        abr_hs_95 = abr_hs_95_pp,
        abr_hs_mode_01 = abr_hs_mode_01_pp,
        abr_hs_median_01 = abr_hs_median_01_pp,
        abr_hs_mean_01 = abr_hs_mean_01_pp,
        abr_ridge_95 = abr_ridge_95_pp,
        abr_ridge_mode_01 = abr_ridge_mode_01_pp,
        abr_ridge_median_01 = abr_ridge_median_01_pp,
        abr_ridge_mean_01 = abr_ridge_mean_01_pp
      )
      
      path <- paste0(
        output_folder ,
        "pp_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(pp, file = path)
      
    }, error = function(e) {
      err_path <- paste0(
        "../Results/06_errors/pp_sparse_error_m_",
        task$m_val,
        "_edgelist_",
        task$edgelist_index,
        ".RData"
      )
      save(e, file = err_path)
      list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = NA,
        abr_hs_95 = NA,
        abr_hs_mode_01 = NA,
        abr_hs_median_01 = NA,
        abr_hs_mean_01 = NA,
        abr_ridge_95 = NA,
        abr_ridge_mode_01 = NA,
        abr_ridge_median_01 = NA,
        abr_ridge_mean_01 = NA
      )
    })
  }, future.seed = TRUE)
  
}


## Plot Predictive Performance

### All Coefficients


plot_pp <- function(pp_folder, m, num_edgelists, title = "") {
  results <- list()
  
  for (i in m) {
    for (j in 1:num_edgelists) {
      file_path <- sprintf(paste0(pp_folder, "pp_m_%d_edgelist_%d.RData"), i, j)
      if (file.exists(file_path)) {
        load(file_path)
        if (exists("pp") && is.list(pp)) {
          temp_df <- data.frame(
            m         = if (!is.null(pp$m))
              pp$m
            else
              NA,
            iteration = if (!is.null(pp$iteration))
              pp$iteration
            else
              NA,
            mle       = if (!is.null(pp$mle))
              pp$mle
            else
              NA,
            abr_ridge     = if (!is.null(pp$abr_ridge))
              pp$abr_ridge
            else
              NA,
            abr_hs        = if (!is.null(pp$abr_hs))
              pp$abr_hs
            else
              NA,
            ebr_hs        = if (!is.null(pp$ebr_hs))
              pp$ebr_hs
            else
              NA
          )
          results[[length(results) + 1]] <- temp_df
        }
        rm(pp)
      }
    }
  }
  
  final_df <- do.call(rbind, results)
  
  pp_is <- final_df %>%
    group_by(m) %>%
    summarize(
      mean_mle   = mean(mle, na.rm = TRUE),
      mean_abr_ridge = mean(abr_ridge, na.rm = TRUE),
      mean_abr_hs    = mean(abr_hs, na.rm = TRUE),
      mean_ebr_hs    = mean(ebr_hs, na.rm = TRUE),
      se_mle     = sd(mle, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge   = sd(abr_ridge, na.rm = TRUE) / sqrt(n()),
      se_abr_hs      = sd(abr_hs, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs      = sd(ebr_hs, na.rm = TRUE) / sqrt(n()),
    ) %>%
    mutate(
      lower_mle   = mean_mle - qnorm(0.975) * se_mle,
      upper_mle   = mean_mle + qnorm(0.975) * se_mle,
      lower_abr_ridge = mean_abr_ridge - qnorm(0.975) * se_abr_ridge,
      upper_abr_ridge = mean_abr_ridge + qnorm(0.975) * se_abr_ridge,
      lower_abr_hs    = mean_abr_hs - qnorm(0.975) * se_abr_hs,
      upper_abr_hs    = mean_abr_hs + qnorm(0.975) * se_abr_hs,
      lower_ebr_hs    = mean_ebr_hs - qnorm(0.975) * se_ebr_hs,
      upper_ebr_hs    = mean_ebr_hs + qnorm(0.975) * se_ebr_hs
      
    )
  
  pp_long <- pp_is %>%
    pivot_longer(
      cols = matches("^(mean|lower|upper)_(mle|abr_ridge|abr_hs|ebr_hs)$"),
      names_to = c("stat", "method"),
      names_pattern = "(mean|lower|upper)_(mle|abr_ridge|abr_hs|ebr_hs)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::select(m, method, mean, lower, upper)
  
  legend_labels <- c(
    "mle" = "MLE",
    "abr_hs" = "ABR HS",
    "abr_ridge" = "ABR Ridge",
    "ebr_hs" = "EBR HS"
  )
  
  method_colors <- c(
    "mle" = "#0072B2",
    "abr_hs" = "#FFA500",
    "abr_ridge" = "#3CB371",
    "ebr_hs" = "#DF00FF"
  )
  
  method_linetypes <- c(
    "mle" = "solid",
    "abr_hs" = "dashed",
    "abr_ridge" = "dotdash",
    "ebr_hs" = "dotted"
  )
  
  p <- ggplot(pp_long, aes(x = m, y = mean, group = method)) +
    geom_ribbon(aes(
      ymin = lower,
      ymax = upper,
      fill = method
    ),
    alpha = 0.2,
    color = NA) +
    geom_line(aes(color = method, linetype = method), size = 1) +
    geom_point(aes(color = method), size = 1.5) +
    scale_x_log10(breaks = m) +
    scale_color_manual(values = method_colors, labels = legend_labels) +
    scale_fill_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    labs(
      title = title,
      x = "Events (M)",
      y = "Predictive Performance",
      color = "Model",
      fill = "Model",
      linetype = "Model"
    )
  
  print(p)
}




### Sparse


plot_pp_sparse <- function(pp_folder, m, num_edgelists, title = "") {
  results <- list()
  for (i in m) {
    for (j in 1:num_edgelists) {
      file_path <- sprintf(paste0(pp_folder, "pp_m_%d_edgelist_%d.RData"), i, j)
      if (file.exists(file_path)) {
        load(file_path)
        if (exists("pp") && is.list(pp)) {
          temp_df <- data.frame(
            m                 = if (!is.null(pp$m))
              pp$m
            else
              NA,
            iteration         = if (!is.null(pp$iteration))
              pp$iteration
            else
              NA,
            abr_ridge         = if (!is.null(pp$abr_ridge))
              pp$abr_ridge
            else
              NA,
            abr_hs_95         = if (!is.null(pp$abr_hs_95))
              pp$abr_hs_95
            else
              NA,
            abr_hs_mode_01    = if (!is.null(pp$abr_hs_mode_01))
              pp$abr_hs_mode_01
            else
              NA,
            abr_hs_mean_01    = if (!is.null(pp$abr_hs_mean_01))
              pp$abr_hs_mean_01
            else
              NA,
            abr_ridge_95      = if (!is.null(pp$abr_ridge_95))
              pp$abr_ridge_95
            else
              NA,
            abr_ridge_mode_01 = if (!is.null(pp$abr_ridge_mode_01))
              pp$abr_ridge_mode_01
            else
              NA,
            abr_ridge_mean_01 = if (!is.null(pp$abr_ridge_mean_01))
              pp$abr_ridge_mean_01
            else
              NA,
            ebr_hs_95         = if (!is.null(pp$ebr_hs_95))
              pp$ebr_hs_95
            else
              NA,
            ebr_hs_mode_01    = if (!is.null(pp$ebr_hs_mode_01))
              pp$ebr_hs_mode_01
            else
              NA,
            ebr_hs_median_01    = if (!is.null(pp$ebr_hs_median_01))
              pp$ebr_hs_median_01
            else
              NA,
            ebr_hs_mean_01    = if (!is.null(pp$ebr_hs_mean_01))
              pp$ebr_hs_mean_01
            else
              NA,
            mle_05            = if (!is.null(pp$mle_05))
              pp$mle_05
            else
              NA
          )
          results[[length(results) + 1]] <- temp_df
        }
        rm(pp)
      }
    }
  }
  final_df <- do.call(rbind, results)
  pp_is <- final_df %>%
    group_by(m) %>%
    summarize(
      mean_abr_ridge         = mean(abr_ridge, na.rm = TRUE),
      mean_abr_hs_95         = mean(abr_hs_95, na.rm = TRUE),
      mean_abr_hs_mode_01    = mean(abr_hs_mode_01, na.rm = TRUE),
      mean_abr_hs_mean_01    = mean(abr_hs_mean_01, na.rm = TRUE),
      mean_abr_ridge_95      = mean(abr_ridge_95, na.rm = TRUE),
      mean_abr_ridge_mode_01 = mean(abr_ridge_mode_01, na.rm = TRUE),
      mean_abr_ridge_mean_01 = mean(abr_ridge_mean_01, na.rm = TRUE),
      mean_ebr_hs_95         = mean(ebr_hs_95, na.rm = TRUE),
      mean_ebr_hs_mode_01    = mean(ebr_hs_mode_01, na.rm = TRUE),
      mean_ebr_hs_median_01    = mean(ebr_hs_median_01, na.rm = TRUE),
      mean_ebr_hs_mean_01    = mean(ebr_hs_mean_01, na.rm = TRUE),
      mean_mle_05            = mean(mle_05, na.rm = TRUE),
      se_abr_ridge           = sd(abr_ridge, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_95           = sd(abr_hs_95, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_mode_01      = sd(abr_hs_mode_01, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_mean_01      = sd(abr_hs_mean_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_95        = sd(abr_ridge_95, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_mode_01   = sd(abr_ridge_mode_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_mean_01   = sd(abr_ridge_mean_01, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_95           = sd(ebr_hs_95, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_mode_01      = sd(ebr_hs_mode_01, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_median_01      = sd(ebr_hs_median_01, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_mean_01      = sd(ebr_hs_mean_01, na.rm = TRUE) / sqrt(n()),
      se_mle_05              = sd(mle_05, na.rm = TRUE) / sqrt(n()),
    ) %>%
    mutate(
      lower_abr_ridge         = mean_abr_ridge - qnorm(0.975) * se_abr_ridge,
      upper_abr_ridge         = mean_abr_ridge + qnorm(0.975) * se_abr_ridge,
      lower_abr_hs_95         = mean_abr_hs_95 - qnorm(0.975) * se_abr_hs_95,
      upper_abr_hs_95         = mean_abr_hs_95 + qnorm(0.975) * se_abr_hs_95,
      lower_abr_hs_mode_01    = mean_abr_hs_mode_01 - qnorm(0.975) * se_abr_hs_mode_01,
      upper_abr_hs_mode_01    = mean_abr_hs_mode_01 + qnorm(0.975) * se_abr_hs_mode_01,
      lower_abr_hs_mean_01    = mean_abr_hs_mean_01 - qnorm(0.975) * se_abr_hs_mean_01,
      upper_abr_hs_mean_01    = mean_abr_hs_mean_01 + qnorm(0.975) * se_abr_hs_mean_01,
      lower_abr_ridge_95      = mean_abr_ridge_95 - qnorm(0.975) * se_abr_ridge_95,
      upper_abr_ridge_95      = mean_abr_ridge_95 + qnorm(0.975) * se_abr_ridge_95,
      lower_abr_ridge_mode_01 = mean_abr_ridge_mode_01 - qnorm(0.975) * se_abr_ridge_mode_01,
      upper_abr_ridge_mode_01 = mean_abr_ridge_mode_01 + qnorm(0.975) * se_abr_ridge_mode_01,
      lower_abr_ridge_mean_01 = mean_abr_ridge_mean_01 - qnorm(0.975) * se_abr_ridge_mean_01,
      upper_abr_ridge_mean_01 = mean_abr_ridge_mean_01 + qnorm(0.975) * se_abr_ridge_mean_01,
      lower_ebr_hs_95         = mean_ebr_hs_95 - qnorm(0.975) * se_ebr_hs_95,
      upper_ebr_hs_95         = mean_ebr_hs_95 + qnorm(0.975) * se_ebr_hs_95,
      lower_ebr_hs_mode_01    = mean_ebr_hs_mode_01 - qnorm(0.975) * se_ebr_hs_mode_01,
      upper_ebr_hs_mode_01    = mean_ebr_hs_mode_01 + qnorm(0.975) * se_ebr_hs_mode_01,
      lower_ebr_hs_median_01    = mean_ebr_hs_median_01 - qnorm(0.975) * se_ebr_hs_median_01,
      upper_ebr_hs_median_01    = mean_ebr_hs_median_01 + qnorm(0.975) * se_ebr_hs_median_01,
      lower_ebr_hs_mean_01    = mean_ebr_hs_mean_01 - qnorm(0.975) * se_ebr_hs_mean_01,
      upper_ebr_hs_mean_01    = mean_ebr_hs_mean_01 + qnorm(0.975) * se_ebr_hs_mean_01,
      lower_mle_05            = mean_mle_05 - qnorm(0.975) * se_mle_05,
      upper_mle_05            = mean_mle_05 + qnorm(0.975) * se_mle_05
    )
  pp_long <- pp_is %>%
    pivot_longer(
      cols = matches("^(mean|lower|upper)_"),
      names_to = c("stat", "method"),
      names_pattern = "(mean|lower|upper)_(.*)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::select(m, method, mean, lower, upper)
  
  legend_labels <- c(
    "abr_ridge"         = "ABR Ridge (Full)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95%-CI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "mle_05"              = "\u03B1 = 0.05 (MLE)"
  )
  
  method_colors <- c(
    "abr_ridge" = "grey20",
    
    # ABR HS methods
    "abr_hs_95"          = "#FFA500",
    "abr_hs_mode_01"     = "#FF4500",
    "abr_hs_median_01"   = "#FFA07A",
    "abr_hs_mean_01"     = "#FFD700",
    
    # ABR Ridge methods
    "abr_ridge_95"       = "#3CB371",
    "abr_ridge_mode_01"  = "#AFEEEE",
    "abr_ridge_median_01" = "#20B2AA",
    "abr_ridge_mean_01"  = "#40E0D0",
    
    # EBR HS
    "ebr_hs_95" = "#DF00FF",
    "ebr_hs_mean_01" = "#C54B8C",
    "ebr_hs_median_01" = "#5D3F6A",
    "ebr_hs_mode_01" = "#D88FD8",
    
    # MLE
    "mle_05"             = "#0072B2"
  )
  
  
  
  # different line types for each method
  method_linetypes <- c(
    "abr_ridge" = "solid",
    "abr_hs_95"          = "solid",
    "abr_hs_mode_01"     = "dotted",
    "abr_hs_median_01"   = "dashed",
    "abr_hs_mean_01"     = "dotdash",
    "abr_ridge_95"       = "solid",
    "abr_ridge_mode_01"  = "dotted",
    "abr_ridge_median_01" = "dashed",
    "abr_ridge_mean_01"  = "dotdash",
    "ebr_hs_95" = "solid",
    "ebr_hs_mean_01" = "dashed",
    "ebr_hs_median_01" = "dotdash",
    "ebr_hs_mode_01" = "dotted",
    "mle_05"             = "solid"
  )
  
  pp_long$method <- factor(
    pp_long$method,
    levels = c(
      "abr_ridge",
      "abr_hs_95",
      "abr_hs_mean_01",
      "abr_hs_median_01",
      "abr_hs_mode_01",
      "abr_ridge_95",
      "abr_ridge_mean_01",
      "abr_ridge_median_01",
      "abr_ridge_mode_01",
      "ebr_hs_95",
      "ebr_hs_mean_01",
      "ebr_hs_median_01",
      "ebr_hs_mode_01",
      "mle_05"
    )
  )
  
  p <- ggplot(pp_long, aes(x = m, y = mean, group = method)) +
    #geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.2, color = NA) +
    geom_line(aes(color = method, linetype = method), size = 1) +
    geom_point(aes(color = method), size = 1.5) +
    scale_x_log10(breaks = m) +
    scale_color_manual(values = method_colors, labels = legend_labels) +
    scale_fill_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    labs(
      title = title,
      x = "Events (M)",
      y = "Predictive Performance",
      color = "Model",
      fill = "Model",
      linetype = "Model"
    )
  print(p)
}




## Table Predictive Performance


table_pp <- function(pp_folder, subfolders, m, num_edgelists) {
  column_labels <- c(
    "mle" = "MLE",
    "abr_hs" = "ABR HS",
    "abr_ridge" = "ABR Ridge",
    "ebr_hs" = "EBR HS"
  )
  
  results <- list()
  
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # build file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive a simple label like "top5" from "01_top5"
          threshold <- sub(".*_", "", subdir)
          
          temp_df <- data.frame(
            m         = pp$m         %||% NA,
            iteration = pp$iteration %||% NA,
            mle       = pp$mle       %||% NA,
            abr_ridge = pp$abr_ridge %||% NA,
            abr_hs    = pp$abr_hs    %||% NA,
            ebr_hs    = pp$ebr_hs    %||% NA,
            threshold = threshold,
            stringsAsFactors = FALSE
          )
          results[[length(results) + 1]] <- temp_df
        }
        rm(pp)
      }
    }
  }
  
  final_df <- do.call(rbind, results)
  
  
  pp_is <- final_df %>%
    group_by(m, threshold) %>%
    summarize(
      mean_mle   = mean(mle, na.rm = TRUE),
      mean_abr_ridge = mean(abr_ridge, na.rm = TRUE),
      mean_abr_hs    = mean(abr_hs, na.rm = TRUE),
      mean_ebr_hs    = mean(ebr_hs, na.rm = TRUE),
      se_mle     = sd(mle, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge   = sd(abr_ridge, na.rm = TRUE) / sqrt(n()),
      se_abr_hs      = sd(abr_hs, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs      = sd(ebr_hs, na.rm = TRUE) / sqrt(n())
    )
  
  pp_long <- pp_is %>%
    pivot_longer(
      cols = matches("^(mean|se)_(mle|abr_ridge|abr_hs|ebr_hs)$"),
      names_to = c("stat", "method"),
      names_pattern = "(mean|se)_(mle|abr_ridge|abr_hs|ebr_hs)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    ungroup() %>%
    dplyr::select(m, method, mean, se, threshold)
  
  pp_wide <- pp_long %>%
    mutate(
      method = recode(method, !!!column_labels),
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>%
    pivot_wider(names_from   = method, values_from  = mean_se) %>%
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m)
  
  
  latex_table <- pp_wide %>%
    rowwise() %>%
    mutate(across(-c(threshold, m), ~ {
      mean_val <- as.numeric(str_extract(.x, "^[0-9\\.]+"))
      row_vals <- c_across(-c(threshold, m)) %>%
        str_extract("^[0-9\\.]+") %>%
        as.numeric()
      max_val <- max(row_vals, na.rm = TRUE)
      if (!is.na(mean_val) && mean_val == max_val) {
        cell_spec(.x, format = "latex", bold = TRUE)
      } else {
        .x
      }
    })) %>%
    ungroup() %>%
    mutate(across(-c(threshold, m), ~ if_else(is.na(.x), "-", .x))) %>%
    kable(format   = "latex",
          booktabs = TRUE,
          escape   = FALSE)
  
  return(list(latex = latex_table, data = pp_wide))
}


### Sparse


table_pp_sparse <- function(pp_folder, subfolders, m, num_edgelists) {
  column_labels <- c(
    "abr_hs_95"           = "95\\% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| $\\ge 0.1$ (ABR HS)",
    "abr_hs_median_01"    = "|Median| $\\ge 0.1$ (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| $\\ge 0.1$ (ABR HS)",
    "abr_ridge_95"        = "95\\%-CI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| $\\ge 0.1$ (ABR Ridge)",
    "ebr_hs_95"           = "95\\% HDI (EBR HS)",
    "ebr_hs_mode_01"      = "|Mode| $\\ge 0.1$ (EBR HS)",
    "ebr_hs_mean_01"      = "|Mean| $\\ge 0.1$ (EBR HS)",
    "ebr_hs_median_01"    = "|Median| $\\ge 0.1$ (EBR HS)",
    "mle_05"              = "$\\alpha = 0.05$ (MLE)"
  )
  
  results <- list()
  
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # build file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive a simple label like "top5" from "01_top5"
          threshold <- sub("^[^_]+_([^_]+)_.*$", "\\1", subdir)
          
          temp_df <- data.frame(
            m                 = if (!is.null(pp$m))
              pp$m
            else
              NA,
            iteration         = if (!is.null(pp$iteration))
              pp$iteration
            else
              NA,
            mle_05            = if (!is.null(pp$mle_05))
              pp$mle_05
            else
              NA,
            abr_hs_95         = if (!is.null(pp$abr_hs_95))
              pp$abr_hs_95
            else
              NA,
            abr_hs_mode_01    = if (!is.null(pp$abr_hs_mode_01))
              pp$abr_hs_mode_01
            else
              NA,
            abr_hs_median_01    = if (!is.null(pp$abr_hs_median_01))
              pp$abr_hs_median_01
            else
              NA,
            abr_hs_mean_01    = if (!is.null(pp$abr_hs_mean_01))
              pp$abr_hs_mean_01
            else
              NA,
            abr_ridge_95      = if (!is.null(pp$abr_ridge_95))
              pp$abr_ridge_95
            else
              NA,
            abr_ridge_mode_01 = if (!is.null(pp$abr_ridge_mode_01))
              pp$abr_ridge_mode_01
            else
              NA,
            abr_ridge_median_01 = if (!is.null(pp$abr_ridge_median_01))
              pp$abr_ridge_median_01
            else
              NA,
            abr_ridge_mean_01 = if (!is.null(pp$abr_ridge_mean_01))
              pp$abr_ridge_mean_01
            else
              NA,
            ebr_hs_95         = if (!is.null(pp$ebr_hs_95))
              pp$ebr_hs_95
            else
              NA,
            ebr_hs_mode_01    = if (!is.null(pp$ebr_hs_mode_01))
              pp$ebr_hs_mode_01
            else
              NA,
            ebr_hs_median_01    = if (!is.null(pp$ebr_hs_median_01))
              pp$ebr_hs_median_01
            else
              NA,
            ebr_hs_mean_01    = if (!is.null(pp$ebr_hs_mean_01))
              pp$ebr_hs_mean_01
            else
              NA,
            threshold = threshold,
            stringsAsFactors = FALSE
          )
          results[[length(results) + 1]] <- temp_df
        }
        rm(pp)
      }
    }
  }
  
  final_df <- do.call(rbind, results)
  
  
  pp <- final_df %>%
    group_by(m, threshold) %>%
    summarize(
      mean_mle_05            = mean(mle_05, na.rm = TRUE),
      mean_abr_hs_95         = mean(abr_hs_95, na.rm = TRUE),
      mean_abr_hs_mode_01    = mean(abr_hs_mode_01, na.rm = TRUE),
      mean_abr_hs_median_01    = mean(abr_hs_median_01, na.rm = TRUE),
      mean_abr_hs_mean_01    = mean(abr_hs_mean_01, na.rm = TRUE),
      mean_abr_ridge_95      = mean(abr_ridge_95, na.rm = TRUE),
      mean_abr_ridge_mode_01 = mean(abr_ridge_mode_01, na.rm = TRUE),
      mean_abr_ridge_median_01 = mean(abr_ridge_median_01, na.rm = TRUE),
      mean_abr_ridge_mean_01 = mean(abr_ridge_mean_01, na.rm = TRUE),
      mean_ebr_hs_95         = mean(ebr_hs_95, na.rm = TRUE),
      mean_ebr_hs_mode_01    = mean(ebr_hs_mode_01, na.rm = TRUE),
      mean_ebr_hs_median_01    = mean(ebr_hs_median_01, na.rm = TRUE),
      mean_ebr_hs_mean_01    = mean(ebr_hs_mean_01, na.rm = TRUE),
      se_mle_05              = sd(mle_05, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_95           = sd(abr_hs_95, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_mode_01      = sd(abr_hs_mode_01, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_median_01      = sd(abr_hs_median_01, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_mean_01      = sd(abr_hs_mean_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_95        = sd(abr_ridge_95, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_mode_01   = sd(abr_ridge_mode_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_median_01   = sd(abr_ridge_median_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_mean_01   = sd(abr_ridge_mean_01, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_95           = sd(ebr_hs_95, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_mode_01      = sd(ebr_hs_mode_01, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_median_01      = sd(ebr_hs_median_01, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs_mean_01      = sd(ebr_hs_mean_01, na.rm = TRUE) / sqrt(n())
    )
  
  pp_long <- pp %>%
    pivot_longer(
      cols = -c(m, threshold),
      names_to = c("stat", "method"),
      names_pattern = "^(mean|se)_(.+)$",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    ungroup() %>%
    dplyr::select(m, method, mean, se, threshold)
  
  pp_wide <- pp_long %>%
    mutate(
      method = recode(method, !!!column_labels),
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>%
    pivot_wider(names_from   = method, values_from  = mean_se) %>%
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m)
  
  latex_table <- pp_wide %>%
    rowwise() %>%
    mutate(across(-c(threshold, m), ~ {
      mean_val <- as.numeric(str_extract(.x, "^[0-9\\.]+"))
      row_vals <- c_across(-c(threshold, m)) %>%
        str_extract("^[0-9\\.]+") %>%
        as.numeric()
      max_val <- max(row_vals, na.rm = TRUE)
      if (!is.na(mean_val) && mean_val == max_val) {
        cell_spec(.x, format = "latex", bold = TRUE)
      } else {
        .x
      }
    })) %>%
    ungroup() %>%
    mutate(across(-c(threshold, m), ~ if_else(is.na(.x), "-", .x))) %>%
    kable(format   = "latex",
          booktabs = TRUE,
          escape   = FALSE)
  
  return(list(latex = latex_table, data = pp_wide))
}


## Independent


table_pp_independent <- function(pp_folder, subfolders, m, num_edgelists) {
  column_labels <- c("mle" = "MLE",
                     "abr_hs" = "ABR HS",
                     "abr_ridge" = "ABR Ridge")
  
  results <- list()
  
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # build file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive a simple label like "top5" from "01_top5"
          threshold <- sub(".*_", "", subdir)
          
          temp_df <- data.frame(
            m         = pp$m         %||% NA,
            iteration = pp$iteration %||% NA,
            mle       = pp$mle       %||% NA,
            abr_ridge = pp$abr_ridge %||% NA,
            abr_hs    = pp$abr_hs    %||% NA,
            threshold = threshold,
            stringsAsFactors = FALSE
          )
          results[[length(results) + 1]] <- temp_df
        }
        rm(pp)
      }
    }
  }
  
  final_df <- do.call(rbind, results)
  
  
  pp_is <- final_df %>%
    group_by(m, threshold) %>%
    summarize(
      mean_mle   = mean(mle, na.rm = TRUE),
      mean_abr_ridge = mean(abr_ridge, na.rm = TRUE),
      mean_abr_hs    = mean(abr_hs, na.rm = TRUE),
      se_mle     = sd(mle, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge   = sd(abr_ridge, na.rm = TRUE) / sqrt(n()),
      se_abr_hs      = sd(abr_hs, na.rm = TRUE) / sqrt(n())
    )
  
  pp_wide <- pp_is %>%
    pivot_longer(
      cols = matches("^(mean|se)_(mle|abr_ridge|abr_hs)$"),
      names_to = c("stat", "method"),
      names_pattern = "(mean|se)_(mle|abr_ridge|abr_hs)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    ungroup() %>%
    dplyr::select(m, method, mean, se, threshold)
  
  latex_table <- pp_wide %>%
    mutate(
      method = recode(method, !!!column_labels),
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>%
    pivot_wider(names_from   = method, values_from  = mean_se) %>%
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m) %>%
    rowwise() %>%
    mutate(across(-c(threshold, m), ~ {
      mean_val <- as.numeric(str_extract(.x, "^[0-9\\.]+"))
      row_vals <- c_across(-c(threshold, m)) %>%
        str_extract("^[0-9\\.]+") %>%
        as.numeric()
      max_val <- max(row_vals, na.rm = TRUE)
      if (!is.na(mean_val) && mean_val == max_val) {
        cell_spec(.x, format = "latex", bold = TRUE)
      } else {
        .x
      }
    })) %>%
    ungroup() %>%
    mutate(across(-c(threshold, m), ~ if_else(is.na(.x), "-", .x))) %>%
    kable(format   = "latex",
          booktabs = TRUE,
          escape   = FALSE)
  
  return(list(latex = latex_table, data = pp_wide))
}


### Sparse


table_pp_sparse_independent <- function(pp_folder, subfolders, m, num_edgelists) {
  column_labels <- c(
    "abr_hs_95"           = "95\\% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| $\\ge 0.1$ (ABR HS)",
    "abr_hs_median_01"    = "|Median| $\\ge 0.1$ (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| $\\ge 0.1$ (ABR HS)",
    "abr_ridge_95"        = "95\\%-CI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| $\\ge 0.1$ (ABR Ridge)",
    "mle_05"              = "$\\alpha = 0.05$ (MLE)"
  )
  
  results <- list()
  
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # build file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive a simple label like "top5" from "01_top5"
          threshold <- sub("^[^_]+_([^_]+)_.*$", "\\1", subdir)
          
          temp_df <- data.frame(
            m                 = if (!is.null(pp$m))
              pp$m
            else
              NA,
            iteration         = if (!is.null(pp$iteration))
              pp$iteration
            else
              NA,
            mle_05            = if (!is.null(pp$mle_05))
              pp$mle_05
            else
              NA,
            abr_hs_95         = if (!is.null(pp$abr_hs_95))
              pp$abr_hs_95
            else
              NA,
            abr_hs_mode_01    = if (!is.null(pp$abr_hs_mode_01))
              pp$abr_hs_mode_01
            else
              NA,
            abr_hs_median_01    = if (!is.null(pp$abr_hs_median_01))
              pp$abr_hs_median_01
            else
              NA,
            abr_hs_mean_01    = if (!is.null(pp$abr_hs_mean_01))
              pp$abr_hs_mean_01
            else
              NA,
            abr_ridge_95      = if (!is.null(pp$abr_ridge_95))
              pp$abr_ridge_95
            else
              NA,
            abr_ridge_mode_01 = if (!is.null(pp$abr_ridge_mode_01))
              pp$abr_ridge_mode_01
            else
              NA,
            abr_ridge_median_01 = if (!is.null(pp$abr_ridge_median_01))
              pp$abr_ridge_median_01
            else
              NA,
            abr_ridge_mean_01 = if (!is.null(pp$abr_ridge_mean_01))
              pp$abr_ridge_mean_01
            else
              NA,
            threshold = threshold,
            stringsAsFactors = FALSE
          )
          results[[length(results) + 1]] <- temp_df
        }
        rm(pp)
      }
    }
  }
  
  final_df <- do.call(rbind, results)
  
  
  pp <- final_df %>%
    group_by(m, threshold) %>%
    summarize(
      mean_mle_05            = mean(mle_05, na.rm = TRUE),
      mean_abr_hs_95         = mean(abr_hs_95, na.rm = TRUE),
      mean_abr_hs_mode_01    = mean(abr_hs_mode_01, na.rm = TRUE),
      mean_abr_hs_median_01    = mean(abr_hs_median_01, na.rm = TRUE),
      mean_abr_hs_mean_01    = mean(abr_hs_mean_01, na.rm = TRUE),
      mean_abr_ridge_95      = mean(abr_ridge_95, na.rm = TRUE),
      mean_abr_ridge_mode_01 = mean(abr_ridge_mode_01, na.rm = TRUE),
      mean_abr_ridge_median_01 = mean(abr_ridge_median_01, na.rm = TRUE),
      mean_abr_ridge_mean_01 = mean(abr_ridge_mean_01, na.rm = TRUE),
      se_mle_05              = sd(mle_05, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_95           = sd(abr_hs_95, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_mode_01      = sd(abr_hs_mode_01, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_median_01      = sd(abr_hs_median_01, na.rm = TRUE) / sqrt(n()),
      se_abr_hs_mean_01      = sd(abr_hs_mean_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_95        = sd(abr_ridge_95, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_mode_01   = sd(abr_ridge_mode_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_median_01   = sd(abr_ridge_median_01, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge_mean_01   = sd(abr_ridge_mean_01, na.rm = TRUE) / sqrt(n())
    )
  
  pp_wide <- pp %>%
    pivot_longer(
      cols = -c(m, threshold),
      names_to = c("stat", "method"),
      names_pattern = "^(mean|se)_(.+)$",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    ungroup() %>%
    dplyr::select(m, method, mean, se, threshold)
  
  latex_table <- pp_wide %>%
    mutate(
      method = recode(method, !!!column_labels),
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>%
    pivot_wider(names_from   = method, values_from  = mean_se) %>%
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m) %>%
    rowwise() %>%
    mutate(across(-c(threshold, m), ~ {
      mean_val <- as.numeric(str_extract(.x, "^[0-9\\.]+"))
      row_vals <- c_across(-c(threshold, m)) %>%
        str_extract("^[0-9\\.]+") %>%
        as.numeric()
      max_val <- max(row_vals, na.rm = TRUE)
      if (!is.na(mean_val) && mean_val == max_val) {
        cell_spec(.x, format = "latex", bold = TRUE)
      } else {
        .x
      }
    })) %>%
    ungroup() %>%
    mutate(across(-c(threshold, m), ~ if_else(is.na(.x), "-", .x))) %>%
    kable(format   = "latex",
          booktabs = TRUE,
          escape   = FALSE)
  
  return(list(latex = latex_table, data = pp_wide))
}




# Calculate Statistics


calculate_statistics_parallel <- function(edgelists,
                                          parameters,
                                          covar,
                                          results_folder,
                                          num_cores) {
  future::plan(future::multisession, workers = num_cores)
  
  
  future_lapply(seq_along(edgelists), function(idx) {
    tryCatch({
      reh <- remify(
        edgelists[[idx]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      
      file_path <- paste0(results_folder,
                          "06_statistics/statistics_edgelist_",
                          idx,
                          ".RData")
      save(statistics, file = file_path)
      list(status = "success",
           edgelist_index = idx,
           file = file_path)
    }, error = function(e) {
      error_file <- paste0(results_folder,
                           "06_errors/error_statistics_edgelist_",
                           idx,
                           ".RData")
      save(e, file = error_file)
      list(
        status = "error",
        edgelist_index = idx,
        error_file = error_file
      )
    })
  }, future.seed = TRUE)
}
