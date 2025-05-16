# functions.R

# This file contains all functions necessary to reproduce the analyses in the folder /Script/02_analyses/
# At the beginning of every Qmd file this script is sourced.

# Generate Edgelists

## Dependent

# Function to generate dependent edgelists in parallel

generate_edgelists_parallel <- function(effects, # effects from remulate
                                        effect_sizes, # effect sizes
                                        covar, # data frame with covariates
                                        num_events, # number of events to generate
                                        num_cores, # number of cores to use
                                        start, # start index for the edgelist
                                        end, # end index for the edgelist
                                        folder # folder to save the edgelists
                                        ) {
  
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  environment(effects) <- environment() # make sure the function is in the environment
  
  # Generate edgelists in parallel
  future.apply::future_lapply(start:end, function(i) {
    set.seed(i) # set seed for reproducibility
    
    # Generate edgelist with remulate
    edgelist <- remulate::remulateTie(effects, # effects from remulate
                                      actors = 1:50, # indices of actors
                                      events = num_events, # number of events to generate
                                      time   = 10000 # time frame
                                      )
    
    # Remove attributes from the edgelist
    edgelist <- as.data.frame(lapply(edgelist, function(x) {
      attributes(x) <- NULL
      x
    }))
    
    # save the edgelist to a file in the specified folder
    save(edgelist, file = paste0(folder, "edgelist_", i, ".RData"))
    
    # clean up memory
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

# Function to generate independent edgelists in parallel

generate_edgelists_independent_parallel <- function(effects, # effects from remulate
                                                    effect_sizes, # effect sizes
                                                    covar, # data frame with covariates
                                                    num_events, # number of events to generate
                                                    num_cores, # number of cores to use
                                                    start, # start index for the edgelist
                                                    end, # end index for the edgelist
                                                    folder # folder to save the edgelists
                                                    ) {
  # Set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  environment(effects) <- environment() # make sure the function is in the environment
  
  # Generate edgelists in parallel
  future.apply::future_lapply(start:end, function(i) {
    set.seed(num_events + i) # set seed for reproducibility
    edgelist <- remulate::remulateTie(effects, # effects from remulate
                                      actors = 1:50, # indices of actors
                                      events = num_events, # number of events to generate
                                      time   = 10000 # time frame
                                      )
    
    # Remove attributes from the edgelist
    edgelist <- as.data.frame(lapply(edgelist, function(x) {
      attributes(x) <- NULL
      x
    }))
    
    # save the edgelist to a file in the specified folder
    save(edgelist, file = paste0(folder, "edgelist_", i, ".RData"))
    
    # clean up memory
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



# Formula from Parameter List

# Function to generate formula from parameter list

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
          return("") # Handle or skip unrecognized components
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

# Function to subset edgelists for dependent analyses

subset_edgelist_single <- function(edgelist, m) {
  return(edgelist[1:m, ]) # Subset the first m rows
}


# Estimating MLE and ABR

# Function to estimate MLE and ABR in parallel


estimate_edgelists_parallel <- function(edgelists,
                                              parameters,
                                              covar,
                                              m,
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

# Function to transform data frame to a format which can be used for Poisson regression


poisson_df <- function(events, # data frame with events (edgelist)
                       tie_stats, # data frame with statistics
                       tie_reh, # remify object
                       t0 = 0 # start time
                       ) {
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


# Estimating brms Models (EBR HS)


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
        formula = glm_formula, # formula
        data = df_poisson, # data frame
        family = poisson(link = "log"), # family for Poisson regression
        prior = set_prior(
          horseshoe(
            df = 3, # 3 degrees of freedom
            scale_global = 1, # global scale
            df_global = 3, # 3 degrees of freedom for global scale
            scale_slab = 2, # slab scale
            df_slab = 4, # 4 degrees of freedom for slab scale
            par_ratio = NULL, # ratio of slab to global scale
            autoscale = TRUE # set to TRUE for automatic scaling
          ),
          class = "b" # only apply to coefficients
        ),
        backend = "cmdstanr", # use cmdstanr backend
        cores = 4 # run each chain on one core
      )
      
      coefs_ebr_hs <- bayestestR::describe_posterior(model, centrality = "all") # summary of posterior distribution
      
      
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
      
      fls <- list.files(tempdir(), full.names = T) # list temporary files
      file.remove(fls) # remove temporary files
      
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


# Selecting Variables

# Function to select variables from estimates in parallel


select_variables_parallel <- function(m, # values of m
                                      n_edgelists, # number of edgelists
                                      folder, # folder with estimates
                                      num_cores # number of cores to use
                                      ) {
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

        # Only get coefficients from parameters
        coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter))
        
        # Selected variables from EBR with Horseshoe
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
        "02a_selected_variables_dependent/selection_m_",
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



# Function to select variables from estimates for independent edgelists in parallel

select_variables_independent_parallel <- function(m, # values of m
                                                  n_edgelists, # number of edgelists
                                                  folder, # folder with estimates
                                                  num_cores # number of cores to use
                                                  ) {
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


# Extracting Estimates from Files

# A function that turns all files with estimates into a data frame

extract_estimates <- function(m, # values of m
                              n_edgelists, # number of edgelists
                              results_folder, # folder with estimates
                              num_cores # number of cores to use
                              ) {
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
      
      # Load the estimates
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
      
      ncoefs <- length(output[["mle_coefs"]][, 1]) # number of coefficients
      
      # Create a data frame to store estimates
      df <- data.frame(matrix(nrow = ncoefs, ncol = 7))
      
      # Set column names
      colnames(df) <- c("m",
                        "iteration",
                        "coef",
                        "mle",
                        "abr_hs",
                        "abr_ridge",
                        "ebr_hs")
      
      # Fill in the data frame with estimates
      df$m <- task$m_val # value of m
      df$iteration <- task$edgelist_index # edgelist index
      
      df$coef <- names(output[["mle_coefs"]][, 1]) # coefficient names
      df$mle <- output[["mle_coefs"]][, 1] # MLE estimates
      df$abr_hs <- output[["shrink_hs"]]$shrunk.mode # ABR HS estimates
      df$abr_ridge <- output[["shrink_ridge"]]$shrunk.mode # ABR Ridge estimates
      
      # Load EBR HS estimates only for m = 100 and 200
      if (task$m_val %in% c(100, 200)) {
        
        # Load the EBR HS estimates
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
        
        # Fill coefficients into the data frame
        df$ebr_hs <- coefs_ebr_hs$MAP
      } else {
        df$ebr_hs <- NA # No EBR HS estimates for other m values
      }
      
      
      return(df)
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  df_estimates <- do.call(rbind, estimates) # Combine all data frames into one
  
  return(df_estimates)
}


# Function to extract estimates from independent edgelists

extract_estimates_independent <- function(m, # values of m
                                          n_edgelists, # number of edgelists
                                          results_folder, # folder with estimates
                                          num_cores # number of cores to use
                                          ) {
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
      # Load the estimates
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
      
      ncoefs <- length(output[["mle_coefs"]][, 1]) # number of coefficients
      df <- data.frame(matrix(nrow = ncoefs, ncol = 6)) # initialize data frame
      colnames(df) <- c("m", "iteration", "coef", "mle", "abr_hs", "abr_ridge") # set column names
      
      # Fill in the data frame with estimates
      df$m <- task$m_val # value of m
      df$iteration <- task$edgelist_index # edgelist index
      
      df$coef <- names(output[["mle_coefs"]][, 1]) # coefficient names
      df$mle <- output[["mle_coefs"]][, 1] # MLE estimates
      df$abr_hs <- output[["shrink_hs"]]$shrunk.mode # ABR HS estimates
      df$abr_ridge <- output[["shrink_ridge"]]$shrunk.mode # ABR Ridge estimates
      
      
      return(df)
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  
  df_estimates <- do.call(rbind, estimates) # Combine all data frames into one
  
  return(df_estimates)
}




# Plot Estimates


plot_estimates <- function(df_estimates, # data frame with estimates
                           parameters, # true parameters
                           m_val, # value of m
                           estimate = c("mle", "abr_hs", "abr_ridge"), # estimate to plot
                           uncertainty = c("ci", "iqr") # uncertainty measure (either confidence interval or interquartile range)
                           ) { 
  
  param_vec <- unlist(parameters) # extract parameters
  
  # Create a data frame with true parameters
  df_true <- data.frame(coef = names(param_vec), true_param = param_vec) %>%
    filter(coef != "baseline") # remove baseline
  
  if (uncertainty == "ci") {
    # Calculate confidence intervals
    
    df_summarized <- df_estimates %>%
      filter(m == m_val, coef != "baseline") %>% # filter for the current m value and remove baseline
      group_by(coef = factor(coef, levels = names(param_vec))) %>% # group by coefficient
      summarise(
        coef_mean = mean(get(estimate)), # mean estimate
        coef_sd = sd(get(estimate)), # standard deviation
        n = n(), # number of estimates
        .groups = "drop"
      ) %>%
      mutate(
        se = coef_sd / sqrt(n), # standard error
        lower_ci = coef_mean - qnorm(0.975) * se, # lower confidence interval limit
        upper_ci = coef_mean + qnorm(0.975) * se # upper confidence interval limit
      )
  } else if (uncertainty == "iqr") {
    # Calculate interquartile range
    
    df_summarized <- df_estimates %>%
      filter(m == m_val, coef != "baseline") %>% # filter for the current m value and remove baseline
      group_by(coef = factor(coef, levels = names(param_vec))) %>% # group by coefficient
      summarise(
        coef_mean = mean(get(estimate)), # mean estimate
        lower_ci = quantile(get(estimate), 0.25), # lower iqr limit
        upper_ci = quantile(get(estimate), 0.75), # upper iqr limit
        .groups = "drop"
      )
  }
  
  df_plot <- left_join(df_summarized, df_true, by = "coef") # join with true parameters
  df_plot$coef <- factor(df_plot$coef, levels = rev(names(param_vec)[names(param_vec) != "baseline"])) # reverse order of coefficients for plot
  
  # Create the plot
  ggplot(df_plot, aes(y = coef)) +
    # Points of true parameters and estimates
    geom_point(aes(x = true_param, color = "True"),
               alpha = 0.5,
               size = 1.5) +
    geom_point(aes(x = coef_mean, color = "Estimate"),
               alpha = 0.5,
               size = 1.5) +
    
    # Error bars for confidence intervals or interquartile range
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci),
                   height = 0.2,
                   color = "black") +
    
    # Dashed lines to emphasize thresholds
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
    
    # limits and labels
    xlim(c(-6, 6)) +
    labs(title = paste0("M=", m_val),
         x = "Value",
         y = "") +
    
    # Theme
    theme_linedraw() +
    
    # Colors and font sizes
    scale_color_manual(values = c("Estimate" = "black", "True" = "red")) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.title = element_blank(),
      title = element_text(size = 10)
    )
}





# Assess Bias


bias_estimates <- function(parameters, # true parameters
                           df_estimates, # data frame with estimates
                           method = c("mle", "abr_hs", "abr_ridge", "ebr_hs") # methods to assess bias
                           ) {
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

# Functions to calculate discovery rates for dependent and independent edgelists


discovery_rate <- function(folder, # folder with estimates
                           parameters, # true parameters
                           m, # values of m
                           num_edgelists, # number of edgelists
                           endo = TRUE, # if TRUE, get endogenous predictors
                           tdr = TRUE, # if TRUE, calculate true discovery rate
                           effect_size = NULL, # if not NULL, filter parameters by effect size
                           num_cores # number of cores to use
                           ) {
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

# Function to calculate discovery rates for independent edgelists


discovery_rate_independent <- function(folder, # folder with estimates
                                       parameters, # true parameters
                                       m, # values of m
                                       num_edgelists, # number of edgelists
                                       endo = TRUE, # if TRUE, get endogenous predictors
                                       tdr = TRUE, # if TRUE, calculate true discovery rate
                                       num_cores # number of cores to use
                                       ) {
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
          "02a_selected_variables_dependent/selection_m_",
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

# Function to plot discovery rates for dependent and independent edgelists

plot_discovery_rate <- function(data, criteria, title = "") {
  
  # Define legend labels
  legend_labels <- c(
    "mle_05"              = "\u03B1 = 0.05 (MLE)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95% HDI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)"
  )
  
  # Define distinct colors for each method
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

# Function to plot distance metric for dependent and independent edgelists


plot_distance_metric <- function(data, criteria, title = "") {
  
  # Define legend labels
  legend_labels <- c(
    "mle_05"              = "\u03B1 = 0.05 (MLE)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95% HDI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)"
  )
  
  # Define distinct colors for each method
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

# Function to calculate Matthews' correlation coefficient (MCC) for dependent and independent edgelists


mcc <- function(data, # data frame with selection results
                parameters, # true parameters
                endo = TRUE, # if TRUE, get endogenous predictors
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
  
  selection_cols <- grep("^selected", names(data), value = TRUE) # detect selection columns
  
  # Calculate Matthews' correlation coefficient (MCC)
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
      
      # Compute the numerator and denominator for MCC
      numerator <- tp * tn - fp * fn
      denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      
      # Check if denominator is NA or zero
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
    # Reformat data to have method and stat as separate columns
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(lower = mean - qnorm(0.975) * se,
           upper = mean + qnorm(0.975) * se) %>%
    # Select and order the desired columns: m, method, mean, sd, lower, upper.
    dplyr::select(m, method, mean, se, lower, upper)
  
  return(results)
}


## Plot MCC


plot_mcc <- function(data, criteria, title = "") {
  
  # Define legend labels
  legend_labels <- c(
    "mle_05"              = "\u03B1 = 0.05 (MLE)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95% HDI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)"
  )
  
  # Define distinct colors for each method
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
    dplyr::filter(method %in% criteria) %>% # filter data by criteria
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

# Function to calculate predictive performance for dependent edgelists
predictive_performance_all <- function(edgelist, # list of edgelists
                                       covar, # data frame with covariates
                                       coefficients, # coefficients
                                       reh, # remify object
                                       statistics, # array with statistics
                                       quantile = 0.95, # quantile to evaluate
                                       warnings = TRUE # display warnings
                                       ) {
  top5 <- rep(FALSE, nrow(edgelist)) # Initialize vector to store results
  
  if (!warnings) {
    suppressWarnings({
      for (i in 1:nrow(edgelist)) {
        # Calculate event rates for every actor in riskset
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


# In-Sample Predictive Performance

## All Coefficients

# Function to calculate in-sample predictive performance for dependent edgelists

pp_is_parallel <- function(edgelists, # list of edgelists
                           m, # vector of m values
                           parameters, # true parameters
                           covar, # data frame with covariates
                           results_folder, # folder with results
                           statistics_folder, # folder with statistics
                           output_folder, # folder to save output
                           num_cores, # number of cores to use
                           quantile = 0.95 # quantile to evaluate
                           ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  
  # Set up task list
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE) # flatten the list
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      
      # remify object
      reh <- remify(
        edgelists[[task$edgelist_index]], # edgelist
        directed = TRUE, # directed network
        model = "tie", # tie-oriented model
        riskset = "full" # full riskset
      )
      
      # load estimates and statistics
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
      
      # Predictive performance for MLE
      mle_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      # Predictive performance for ABR HS
      abr_hs_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_hs$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      
      # If m is 100 or 200, load EBR HS estimates
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
        
        # Predictive performance for EBR HS
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
        ebr_hs_pp <- NA # if m is not 100 or 200
      }
      
      # combine results into a list
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp,
        ebr_hs = ebr_hs_pp
      )
      
      # save results to file
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

# Function to calculate in-sample predictive performance for dependent edgelists with sparse models after variable selection

pp_is_sparse_parallel <- function(edgelists, # list of edgelists
                                  m, # vector of m values
                                  parameters, # true parameters
                                  covar, # data frame with covariates
                                  results_folder, # folder with results
                                  statistics_folder, # folder with statistics
                                  selected_folder, # folder with selected variables
                                  output_folder, # folder to save output
                                  num_cores, # number of cores to use
                                  quantile = 0.95 # quantile to evaluate
                                  ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Set up task list
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE) # flatten the list
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      
      # remify object
      reh <- remify(
        edgelists[[task$edgelist_index]], # edgelist
        directed = TRUE, # directed network
        model = "tie", # tie-oriented model
        riskset = "full" # full riskset
      )
      
      # load estimates, statistics and selected variables
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
      
      statistics_is <- statistics[1:task$m_val, , ] # statistics for in-sample
      
      rm(statistics) # remove statistics to save memory
      
      # Predictive performance for MLE
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
      
      # Predictive performance for ABR Ridge
      
      # 95% HDI
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # abs Mode >= 0.1
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # abs Median >= 0.1
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # abs Mean >= 0.1
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      
      # 95% HDI
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
      
      # abs Mode >= 0.1
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
      
      # abs Median >= 0.1
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
      
      # abs Mean >= 0.1
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
      
      # Predictive performance for EBR HS
      # only if m is 100 or 200
      if (task$m_val %in% c(100, 200)) {
        # load estimates
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

        coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter)) # only keep coefficients
        
        # 95% HDI
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
        
        # abs Mode >= 0.1
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
        
        # abs Median >= 0.1
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
        
        # abs Mean >= 0.1
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
      
      # Full ABR Ridge
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelists[[task$edgelist_index]][1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_is,
        quantile = quantile,
        warnings = FALSE
      )
      
      # Combine results into a list
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
      
      # save results to file
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

# Function to calculate in-sample predictive performance for independent edgelists

pp_is_independent_parallel <- function(parameters, # true parameters
                                       covar, # data frame with covariates
                                       results_folder, # folder with results
                                       output_folder, # folder to save output
                                       num_cores, # number of cores to use
                                       quantile = 0.95, # quantile to evaluate
                                       task_list # list of tasks with edgelists
                                       ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      
      # remify object
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      # load estimates
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
      
      # calculate statistics
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      
      # Predictive performance for MLE
      mle_pp <- predictive_performance_all(
        edgelist = task$edgelist[1:task$m_val, ],
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      abr_ridge_pp <- predictive_performance_all(
        edgelist = task$edgelist[1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      # Predictive performance for ABR HS
      abr_hs_pp <- predictive_performance_all(
        edgelist = task$edgelist[1:task$m_val, ],
        covar = covar,
        coefficients = output$shrink_hs$shrunk.mode,
        reh = reh,
        statistics = statistics[1:task$m_val, , ],
        quantile = quantile,
        warnings = FALSE
      )
      
      
      # List to store results
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp
      )
      
      # save results to file
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

# Function to calculate in-sample predictive performance for independent edgelists with sparse models after variable selection

pp_is_independent_sparse_parallel <- function(parameters, # true parameters
                                              covar, # data frame with covariates
                                              results_folder, # folder with results
                                              selected_folder, # folder with selected variables
                                              output_folder, # folder to save output
                                              num_cores, # number of cores to use
                                              quantile = 0.95, # quantile to evaluate
                                              task_list # list of tasks with edgelists
                                              ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      
      # remify object
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      # calculate statistics
      statistics_is <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      
      # load estimates
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
      
      # load selected variables
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
      
      # Predictive performance for MLE
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
      
      # Predictive performance for ABR Ridge
      
      # 95% HDI
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # abs Mode >= 0.1
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # abs Median >= 0.1
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # abs Mean >= 0.1
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = task$edgelist,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_is,
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      
      # 95% HDI
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
      
      # abs Mode >= 0.1
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
      
      # abs Median >= 0.1
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
      
      # abs Mean >= 0.1
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
      
      # Combine results into a list
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
      
      # save results to file
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



# Out-Of-Sample Predictive Performance

## Full Models

# Function to calculate out-of-sample predictive performance for dependent edgelists

pp_oos_parallel <- function(edgelists, # list of edgelists
                            m, # vector of m values
                            parameters, # true parameters
                            covar, # data frame with covariates
                            results_folder, # folder with results
                            statistics_folder, # folder with statistics
                            output_folder, # folder to save output
                            num_cores, # number of cores to use
                            quantile = 0.95, # quantile to evaluate
                            new = 1000 # number of out-of-sample observations to evaluate
                            ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Set up task list
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE) # flatten the list
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      
      # remify object
      reh <- remify(
        edgelists[[task$edgelist_index]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      # load estimates and statistics
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
      
      edgelist_oos <- edgelists[[task$edgelist_index]][(task$m_val + 1):(task$m_val + new), ] # out-of-sample edgelist
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ] # out-of-sample statistics
      
      # Predictive performance for MLE
      mle_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      # Predictive performance for ABR Ridge
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      # Predictive performance for ABR HS
      abr_hs_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_hs$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      # Predictive performance for EBR HS
      # only if m is 100 or 200
      if (task$m_val %in% c(100, 200)) {
        ä# load estimates
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
        ebr_hs_pp <- NA # if m is not 100 or 200, set to NA
      }
      
      # combine results into a list
      pp <- list(
        m = task$m_val,
        iteration = task$edgelist_index,
        mle = mle_pp,
        abr_ridge = abr_ridge_pp,
        abr_hs = abr_hs_pp,
        ebr_hs = ebr_hs_pp
      )
      
      # save results to file
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


## Sparse Models

# Function to calculate out-of-sample predictive performance for dependent edgelists with sparse models after variable selection


pp_oos_sparse_parallel <- function(edgelists, # list of edgelists
                                   m, # vector of m values
                                   parameters, # true parameters
                                   covar, # data frame with covariates
                                   results_folder, # folder with results
                                   statistics_folder, # folder with statistics
                                   selected_folder, # folder with selected variables
                                   output_folder, # folder to save output
                                   num_cores, # number of cores to use
                                   quantile = 0.95, # quantile to evaluate
                                   new = 1000 # number of out-of-sample observations to evaluate
                                   ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Set up task list
  task_list <- lapply(m, function(m_val) {
    lapply(seq_along(edgelists), function(edgelist_index) {
      list(m_val = m_val, edgelist_index = edgelist_index)
    })
  })
  
  task_list <- unlist(task_list, recursive = FALSE) # flatten the list
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      # remify object
      reh <- remify(
        edgelists[[task$edgelist_index]],
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      # load estimates, statistics and selected variables
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
      
      edgelist_oos <- edgelists[[task$edgelist_index]][(task$m_val + 1):(task$m_val + new), ] # out-of-sample edgelist
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ] # out-of-sample statistics
      rm(statistics)
      
      # Predictive performance for MLE
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
      
      # Predictive performance for ABR HS
      
      # 95% HDI
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # abs Mode >= 0.1
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # abs Median >= 0.1
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # abs Mean >= 0.1
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      
      # 95% HDI
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
      
      # abs Mode >= 0.1
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
      
      # abs Median >= 0.1
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
      
      # abs Mean >= 0.1
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
      
      # Predictive performance for EBR HS
      # only if m is 100 or 200
      if (task$m_val %in% c(100, 200)) {
        # load estimates
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

        coefs_ebr_hs$Parameter <- sub("1$", "", sub("^b_", "", coefs_ebr_hs$Parameter)) # only coefficients
        
        # 95% HDI
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
        
        # abs Mode >= 0.1
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
        
        # abs Median >= 0.1
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
        
        # abs Mean >= 0.1
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
      
      # Full ABR Ridge
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      
      # Combine results into a list
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
      
      # save results to file
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


## Independent

# Function to calculate out-of-sample predictive performance for independent edgelists

pp_oos_independent_parallel <- function(parameters, # true parameters
                                        covar, # data frame with covariates
                                        results_folder, # folder with results
                                        output_folder, # folder to save output
                                        num_cores, # number of cores to use
                                        quantile = 0.95, # quantile to evaluate
                                        task_list, # list of tasks with edgelists
                                        new = 1000 # number of out-of-sample observations to evaluate
                                        ) {
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      
      # remify object
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      # load estimates
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
      
      # calculate statistics
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      edgelist_oos <- task$edgelist[(task$m_val + 1):(task$m_val + new), ] # out-of-sample edgelist
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ] # out-of-sample statistics
      
      rm(statistics) # remove statistics to save memory
      
      # Predictive performance for MLE
      mle_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output[["mle_coefs"]][, 1],
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      abr_ridge_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_ridge$estimates$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      # Predictive performance for ABR HS
      abr_hs_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = output$shrink_hs$shrunk.mode,
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile,
        warnings = FALSE
      )
      
      
      # Combine results into a list
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


pp_oos_independent_sparse_parallel <- function(parameters, # true parameters
                                               covar, # data frame with covariates
                                               results_folder, # folder with results
                                               statistics_folder, # folder with statistics
                                               selected_folder, # folder with selected variables
                                               output_folder, # folder to save output
                                               num_cores, # number of cores to use
                                               quantile = 0.95, # quantile to evaluate
                                               new = 1000, # number of out-of-sample observations to evaluate
                                               task_list # list of tasks with edgelists
                                               ) {
  
  # set up parallel backend
  future::plan(future::multisession, workers = num_cores)
  
  # Run tasks in parallel
  future_lapply(task_list, function(task) {
    tryCatch({
      # remify object
      reh <- remify(
        task$edgelist,
        directed = TRUE,
        model = "tie",
        riskset = "full"
      )
      
      # load estimates and selected variables
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
      
      # calculate statistics
      statistics <- remstats(
        reh = reh,
        tie_effects = generate_formula(parameters),
        attr_actors = covar
      )
      edgelist_oos <- task$edgelist[(task$m_val + 1):(task$m_val + new), ] # out-of-sample edgelist
      statistics_oos <- statistics[(task$m_val + 1):(task$m_val + new), , ] # out-of-sample statistics
      
      rm(statistics) # remove statistics to save memory
      
      # Predictive performance for MLE
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
      
      # Predictive performance for ABR HS
      
      # 95% HDI
      abr_hs_95_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_95),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # abs Mode >= 0.1
      abr_hs_mode_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mode_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # abs Median >= 0.1
      abr_hs_median_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_median_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # abs Mean >= 0.1
      abr_hs_mean_01_pp <- predictive_performance_all(
        edgelist = edgelist_oos,
        covar = covar,
        coefficients = replace(
          output$shrink_hs$shrunk.mode,!rownames(output$shrink_hs) %in% c("baseline", selected_vars$abr_hs_mean_01),
          0
        ),
        reh = reh,
        statistics = statistics_oos,
        quantile = quantile
      )
      
      # Predictive performance for ABR Ridge
      
      # 95% HDI
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
      
      # abs Mode >= 0.1
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
      
      # abs Median >= 0.1
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
      
      # abs Mean >= 0.1
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
      
      # combine results into a list
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
      
      # save results to file
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


# Plot Predictive Performance

## All Coefficients

# Function to plot predictive performance for full models

plot_pp <- function(pp_folder, # folder with predictive performance results
                    m, # vector of m values
                    num_edgelists, # number of edgelists
                    title = "" # title of the plot
                    ) {
  # Create empty list to store results
  results <- list()
  
  # Fill list with results
  for (i in m) {
    for (j in 1:num_edgelists) {
      file_path <- sprintf(paste0(pp_folder, "pp_m_%d_edgelist_%d.RData"), i, j) # file path
      if (file.exists(file_path)) {
        load(file_path) # load file
        
        # check if pp is a list and exists
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
  
  final_df <- do.call(rbind, results) # combine list to data frame
  
  # Summarize results
  pp_is <- final_df %>%
    group_by(m) %>% # group by m
    summarize(
      # calculate mean and standard error for each method
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
  
  # Reshape data for plotting
  pp_long <- pp_is %>%
    pivot_longer(
      cols = matches("^(mean|lower|upper)_(mle|abr_ridge|abr_hs|ebr_hs)$"),
      names_to = c("stat", "method"),
      names_pattern = "(mean|lower|upper)_(mle|abr_ridge|abr_hs|ebr_hs)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::select(m, method, mean, lower, upper)
  
  # labels for legend on the plot
  legend_labels <- c(
    "mle" = "MLE",
    "abr_hs" = "ABR HS",
    "abr_ridge" = "ABR Ridge",
    "ebr_hs" = "EBR HS"
  )
  
  # colors and linetypes for each method
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
  
  # create plot
  p <- ggplot(pp_long, aes(x = m, y = mean, group = method)) +
    # confidence intervals
    geom_ribbon(aes(
      ymin = lower,
      ymax = upper,
      fill = method
    ),
    alpha = 0.2,
    color = NA) +
    # lines and points
    geom_line(aes(color = method, linetype = method), size = 1) +
    geom_point(aes(color = method), size = 1.5) +
    
    # log the x axis
    scale_x_log10(breaks = m) +
    
    # customize colors and linetypes
    scale_color_manual(values = method_colors, labels = legend_labels) + 
    scale_fill_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    
    # customize theme
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    
    # customize labels
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
  
  # Create empty list to store results
  results <- list()
  
  # Fill list with results
  for (i in m) {
    for (j in 1:num_edgelists) {
      file_path <- sprintf(paste0(pp_folder, "pp_m_%d_edgelist_%d.RData"), i, j) # file path
      if (file.exists(file_path)) {
        load(file_path) # load file
        if (exists("pp") && is.list(pp)) { # check if pp is a list and exists
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
  final_df <- do.call(rbind, results) # combine list to data frame
  
  # Summarize results
  pp_is <- final_df %>%
    group_by(m) %>% # group by m
    summarize(
      
      # calculate mean and standard error for each method
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
      # calculate confidence intervals
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
  
  # Reshape data for plotting
  pp_long <- pp_is %>%
    pivot_longer(
      cols = matches("^(mean|lower|upper)_"),
      names_to = c("stat", "method"),
      names_pattern = "(mean|lower|upper)_(.*)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::select(m, method, mean, lower, upper)
  
  # labels for legend on the plot
  legend_labels <- c(
    "abr_ridge"         = "ABR Ridge (Full)",
    "abr_hs_95"           = "95% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| \u2265 0.1 (ABR HS)",
    "abr_hs_median_01"    = "|Median| \u2265 0.1 (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| \u2265 0.1 (ABR HS)",
    "abr_ridge_95"        = "95% HDI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| \u2265 0.1 (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| \u2265 0.1 (ABR Ridge)",
    "ebr_hs_95" = "95% HDI (EBR HS)",
    "ebr_hs_mean_01" = "|Mean| \u2265 0.1 (EBR HS)",
    "ebr_hs_median_01" = "|Median| \u2265 0.1 (EBR HS)",
    "ebr_hs_mode_01" = "|Mode| \u2265 0.1 (EBR HS)",
    "mle_05"              = "\u03B1 = 0.05 (MLE)"
  )
  
  # colors for each method
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
  
  # reorder method factor levels
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
  
  # create plot
  p <- ggplot(pp_long, aes(x = m, y = mean, group = method)) +
    # line and point for each method
    geom_line(aes(color = method, linetype = method), size = 1) +
    geom_point(aes(color = method), size = 1.5) +
    
    # log the x axis
    scale_x_log10(breaks = m) +
    
    # customize colors and linetypes
    scale_color_manual(values = method_colors, labels = legend_labels) +
    scale_fill_manual(values = method_colors, labels = legend_labels) +
    scale_linetype_manual(values = method_linetypes, labels = legend_labels) +
    
    # customize theme
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    
    # customize labels
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




# Table Predictive Performance

# Function to create a table of predictive performance results


table_pp <- function(pp_folder, # folder with predictive performance results
                     subfolders, # subfolders with different thresholds
                     m, # vector of m values
                     num_edgelists # number of edgelists
                     ) {
  
  # labels for the columns
  column_labels <- c(
    "mle" = "MLE",
    "abr_hs" = "ABR HS",
    "abr_ridge" = "ABR Ridge",
    "ebr_hs" = "EBR HS"
  )
  
  # create empty list to store results
  results <- list()
  
  # fill list with results
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive threshold from subdir name
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
  
  final_df <- do.call(rbind, results) # combine list to data frame
  
  # Summarize results
  pp_is <- final_df %>%
    group_by(m, threshold) %>% # group by m and threshold
    summarize(
      # calculate mean and standard error for each method
      mean_mle   = mean(mle, na.rm = TRUE),
      mean_abr_ridge = mean(abr_ridge, na.rm = TRUE),
      mean_abr_hs    = mean(abr_hs, na.rm = TRUE),
      mean_ebr_hs    = mean(ebr_hs, na.rm = TRUE),
      se_mle     = sd(mle, na.rm = TRUE) / sqrt(n()),
      se_abr_ridge   = sd(abr_ridge, na.rm = TRUE) / sqrt(n()),
      se_abr_hs      = sd(abr_hs, na.rm = TRUE) / sqrt(n()),
      se_ebr_hs      = sd(ebr_hs, na.rm = TRUE) / sqrt(n())
    )
  
  # Reshape data for table
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
  
  # Prepare data for table
  pp_wide <- pp_long %>%
    mutate(
      method = recode(method, !!!column_labels),
      # combine mean and se into a single string
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>% # select relevant columns
    pivot_wider(names_from   = method, values_from  = mean_se) %>% # turn to wide format
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m) # sort by threshold and m
  
  # Create LaTeX table
  latex_table <- pp_wide %>%
    rowwise() %>% # apply row-wise operations
    
    # apply bold formatting to the maximum value in each row
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


table_pp_sparse <- function(pp_folder, # folder with predictive performance results
                            subfolders, # subfolders with different thresholds
                            m, # vector of m values
                            num_edgelists # number of edgelists
                            ) {
  # labels for the columns
  column_labels <- c(
    "abr_hs_95"           = "95\\% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| $\\ge 0.1$ (ABR HS)",
    "abr_hs_median_01"    = "|Median| $\\ge 0.1$ (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| $\\ge 0.1$ (ABR HS)",
    "abr_ridge_95"        = "95\\% HDI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| $\\ge 0.1$ (ABR Ridge)",
    "ebr_hs_95"           = "95\\% HDI (EBR HS)",
    "ebr_hs_mode_01"      = "|Mode| $\\ge 0.1$ (EBR HS)",
    "ebr_hs_mean_01"      = "|Mean| $\\ge 0.1$ (EBR HS)",
    "ebr_hs_median_01"    = "|Median| $\\ge 0.1$ (EBR HS)",
    "mle_05"              = "$\\alpha = 0.05$ (MLE)"
  )
  
  # create empty list to store results
  results <- list()
  
  # fill list with results
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive threshold from subdir name
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
  
  final_df <- do.call(rbind, results) # combine list to data frame
  
  
  pp <- final_df %>%
    group_by(m, threshold) %>% # group by m and threshold
    summarize(
      # calculate mean and standard error for each method
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
  
  # Summarize results
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
      # combine mean and se into a single string
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>% # select relevant columns
    pivot_wider(names_from   = method, values_from  = mean_se) %>% # turn to wide format
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m) # sort by threshold and m
  
  # Create LaTeX table
  latex_table <- pp_wide %>%
    rowwise() %>%
    # apply bold formatting to the maximum value in each row
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

# Function to create a table of predictive performance results for independent data

table_pp_independent <- function(pp_folder, # folder with predictive performance results
                                 subfolders, # subfolders with different thresholds
                                 m, # vector of m values
                                 num_edgelists # number of edgelists
                                 ) {
  
  # labels for the columns
  column_labels <- c("mle" = "MLE",
                     "abr_hs" = "ABR HS",
                     "abr_ridge" = "ABR Ridge")
  
  # create empty list to store results
  results <- list()
  
  # fill list with results
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path)) # load the data
        if (is.list(pp)) {
          # derive threshold from subdir name
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
  
  final_df <- do.call(rbind, results) # combine list to data frame
  
  # Summarize results
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
  
  # Reshape data for table
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
  
  # Create LaTeX table
  latex_table <- pp_wide %>%
    mutate(
      method = recode(method, !!!column_labels),
      # combine mean and se into a single string
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>% # select relevant columns
    pivot_wider(names_from   = method, values_from  = mean_se) %>% # turn to wide format
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m) %>% # sort by threshold and m
    rowwise() %>%
    
    # apply bold formatting to the maximum value in each row
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


table_pp_sparse_independent <- function(pp_folder, # folder with predictive performance results
                                        subfolders, # subfolders with different thresholds
                                        m, # vector of m values
                                        num_edgelists # number of edgelists
                                        ) {
  
  # labels for the columns
  column_labels <- c(
    "abr_hs_95"           = "95\\% HDI (ABR HS)",
    "abr_hs_mode_01"      = "|Mode| $\\ge 0.1$ (ABR HS)",
    "abr_hs_median_01"    = "|Median| $\\ge 0.1$ (ABR HS)",
    "abr_hs_mean_01"      = "|Mean| $\\ge 0.1$ (ABR HS)",
    "abr_ridge_95"        = "95\\% HDI (ABR Ridge)",
    "abr_ridge_mode_01"   = "|Mode| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_median_01" = "|Median| $\\ge 0.1$ (ABR Ridge)",
    "abr_ridge_mean_01"   = "|Mean| $\\ge 0.1$ (ABR Ridge)",
    "mle_05"              = "$\\alpha = 0.05$ (MLE)"
  )
  
  # create empty list to store results
  results <- list() 
  
  # fill list with results
  for (i in m) {
    for (j in seq_len(num_edgelists)) {
      for (subdir in subfolders) {
        # file path
        file_path <- sprintf("%s%s/pp_m_%d_edgelist_%d.RData",
                             pp_folder,
                             subdir,
                             i,
                             j)
        if (!file.exists(file_path))
          next
        
        pp <- get(load(file_path))
        if (is.list(pp)) {
          # derive threshold from subdir name
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
  
  final_df <- do.call(rbind, results) # combine list to data frame
  
  # Summarize results
  pp <- final_df %>%
    group_by(m, threshold) %>% # group by m and threshold
    summarize(
      # calculate mean and standard error for each method
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
  
  # Reshape data for table
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
  
  # Create LaTeX table
  latex_table <- pp_wide %>%
    mutate(
      method = recode(method, !!!column_labels),
      # combine mean and se into a single string
      mean_se = if_else(
        is.na(mean) | is.na(se),
        NA_character_,
        sprintf("%.3f (%.3f)", mean, se)
      )
    ) %>%
    dplyr::select(threshold, m, method, mean_se) %>% # select relevant columns
    pivot_wider(names_from   = method, values_from  = mean_se) %>% # turn to wide format
    arrange(factor(threshold, levels = c("top20", "top10", "top5")), m) %>% # sort by threshold and m
    rowwise() %>%
    # apply bold formatting to the maximum value in each row
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

# Function to calculate statistics for each edgelist in parallel


calculate_statistics_parallel <- function(edgelists, # list of edgelists
                                          parameters, # true parameters
                                          covar, # data frame with covariates
                                          results_folder, # folder to save results
                                          num_cores # number of cores to use
                                          ) {
  
  # set up parallel processing
  future::plan(future::multisession, workers = num_cores)
  
  # run the remify and remstats functions in parallel
  future_lapply(seq_along(edgelists), function(idx) {
    tryCatch({
      
      # remify the edgelist and calculate statistics
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
      
      # save the statistics
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
