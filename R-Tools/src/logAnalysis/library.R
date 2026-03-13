library(tidyverse)
library(nptest)
library(knitr)
library(TOSTER)
library(data.table)

#
#
#  H E L P E R S  
#
#


boot.lower <- function(x,y,Iter){
  npbs <- np.boot(x = x, statistic = y, R = Iter)
  return(ifelse(is.nan(npbs$bca[2,1]),0,npbs$bca[2,1]))
}

boot.upper <- function(x,y,Iter){
  npbs <- np.boot(x = x, statistic = y, R = Iter)
  return(ifelse(is.nan(npbs$bca[2,2]),0,npbs$bca[2,2]))
  return(npbs$bca[2,2])
}

format_simtime <- function(simtime_ms) {
  total_seconds <- simtime_ms / 1000
  hours   <- floor(total_seconds / 3600)
  minutes <- floor((total_seconds %% 3600) / 60)
  seconds <- floor(total_seconds %% 60)
  millis  <- round((total_seconds - floor(total_seconds)) * 1000)
  
  sprintf("%02d:%02d:%02d.%03d", hours, minutes, seconds, millis)
}

auto_time_to_ms <- function(x) {
  
  # If it's already numeric → treat as milliseconds
  if (is.numeric(x)) {
    return(as.integer(x))
  }
  
  # Remove accidental spaces
  x <- trimws(x)
  
  # Is it a pure integer string?
  if (grepl("^[0-9]+$", x)) {
    return(as.integer(x))
  }
  
  # Otherwise parse as time string
  parts <- unlist(strsplit(x, ":"))
  
  # Case 1: SS.mmm
  if (length(parts) == 1) {
    sec <- as.numeric(parts[1])
    return(round(sec * 1000))
  }
  
  # Case 2: MM:SS.mmm
  if (length(parts) == 2) {
    min <- as.numeric(parts[1])
    sec <- as.numeric(parts[2])
    return(round((min * 60 + sec) * 1000))
  }
  
  # Case 3: HH:MM:SS.mmm
  if (length(parts) == 3) {
    hour <- as.numeric(parts[1])
    min  <- as.numeric(parts[2])
    sec  <- as.numeric(parts[3])
    return(round((hour * 3600 + min * 60 + sec) * 1000))
  }
  
  stop(paste("Unrecognized time format:", x))
}


#
#
# ATTACK THEORETICAL PROBABILITY
#
#

attacker_success_probability <- function(q, z) {
  if (q <= 0) return(0)
  if (q >= 1) return(1)
  
  p <- 1 - q
  lambda <- z * (q / p)
  
  sum <- 1.0
  
  for (k in 0:z) {
    poisson <- dpois(k, lambda)   # numerically stable Poisson probability
    sum <- sum - poisson * (1 - (q / p)^(z - k + 1))
    #sum <- sum - poisson * (1 - (q / p)^(z - k))
  }
  
  return(sum)
}




equivalence <- function(trials, successes, theoretical, delta = 0.03, alpha = 0.05) {
  n  <- trials
  x  <- successes
  phat <- x/n
  p0 <- theoretical

  low  <- p0 - delta
  high <- p0 + delta
  
  # Standard error under observed proportion
  se <- sqrt(phat * (1 - phat) / n)
  
  # Lower test: H0: p <= low
  z1 <- (phat - low) / se
  p1 <- 1 - pnorm(z1)
  
  # Upper test: H0: p >= high
  z2 <- (phat - high) / se
  p2 <- pnorm(z2)
  
  
  if ((p1 < alpha) && (p2 < alpha)){
    result <- "pass"
  } else {
    result <- "fail"
  }
  
  ci90 <- prop.test(x, n,
                    conf.level = 1 - 2*alpha,
                    correct = FALSE)$conf.int
  
  if (result == "pass") {
    cat(paste0(
      "PASS. An equivalence test was conducted to evaluate whether the observed proportion ",
      "phat differed from the theoretical value p0 by more than delta. ",
      "The observed proportion was ", phat,
      " (90% CI [", round(ci90[1],5), ", ", round(ci90[2],5), "]). ",
      "The entire confidence interval was contained within the equivalence bounds [",
      round(low,5), ", ", round(high,5), 
      "], supporting statistical equivalence at significance level alpha."
    ),"\n")
  } else {
    cat(paste0(
      "FAIL. An equivalence test was conducted to evaluate whether the observed proportion ",
      "phat differed from the theoretical value p0 by more than delta. ",
      "The observed proportion was ", phat,
      " (90% CI [", round(ci90[1],5), ", ", round(ci90[2],5), "]). ",
      "The confidence interval overlapped the equivalence bounds [",
      round(low,5), ", ", round(high,5), "], ",
      "so statistical equivalence at significance level alpha could not be established."
    ),"\n")
  }
}
  

equivalence_vec <- function(trials, successes, theoretical, delta = 0.03, alpha = 0.05) {
  n  <- trials
  x  <- successes
  phat <- x / n
  p0 <- theoretical
  
  
  
  low  <- max(0,p0 - delta)
  high <- min(1,p0 + delta)
  
  se <- sqrt(phat * (1 - phat) / n)
  se1 <- sqrt(low * (1 - low) / n)
  se2 <- sqrt(high * (1 - high) / n)
  
  z1 <- (phat - low) / se1
  p1 <- 1 - pnorm(z1)
  
  z2 <- (phat - high) / se2
  p2 <- pnorm(z2)
  
  result <- ifelse((p1 < alpha) & (p2 < alpha), "PASS", "FAIL")
  
  ci90 <- prop.test(x, n,
                    conf.level = 1 - 2*alpha,
                    correct = FALSE)$conf.int
  
  tibble::tibble(
    lower_bound = low,
    ci_low = ci90[1],
    observed = phat,
    theoretical = p0,
    ci_high = ci90[2],
    upper_bound = high,
    p_lower = p1,
    p_upper = p2,
    result = result,
    delta = delta
  )
}






#
#
# PACE ANALYSIS
#
#


producePaceData <- function(dataFolder,experiment){
  pacedata <- read_csv(paste0(dataFolder,experiment,"/EventLog - ",experiment,".csv")) %>% 
    filter(EventType == "Event_ContainerValidation") %>%
    arrange(SimID, SimTime) %>%
    group_by(SimID) %>%
    mutate(SimDiff = SimTime - lag(SimTime, default = 0)) %>%  # use 0 for first lag
    ungroup() %>%
    select(SimID, SimDiff)
  
  write_csv(pacedata,paste0(dataFolder,experiment,"/PaceData - ",experiment,".csv"))
}


paceInfo <- function(dataFolder,experiment) {
  
  # Ensure experiment is a character vector
  experiment <- unlist(experiment)
  
  
  # Read and merge all CSV files
  pacedata <- purrr::map_dfr(experiment, function(exp) {
    readr::read_csv(
      paste0(dataFolder, exp, "/PaceData - ", exp, ".csv")
    )
  })

  pace = pacedata %>%
  arrange(SimID, SimDiff) %>%             # ensure ordered
  summarise(
    `Block Time (mins) - Mean....:` = mean(SimDiff)/60000,
    `Block Time (mins) - St. Dev.:`   = sd(SimDiff)/60000,
    Samples = n()
  )
  
  return(pace)
}



pace_theoretical <- function(n, difficulty,power) {
  p <- 1 / difficulty
  u <- runif(n)
  trials = floor(log1p(-u) / log1p(-p))
  return (trials / power)
}

pace_theoretical2 <- function(n, difficulty,power) {
  trials = rgeom(n, p = 1/difficulty)
  return (trials / power)
}



#
#
# SIMULATION PERFORMANCE
#
#

getRunTime <- function(dataFolder,experiment) {
  runTimes <- read_csv(paste0(dataFolder,experiment,"/EventLog - ",experiment,".csv")) %>% 
    group_by(SimID) %>% summarise(
      sysTime = as.integer(max(SysTime) - min(SysTime)),
      sysTime_formated = format_simtime(sysTime))
  return(runTimes)
}



#
#
#  M A I N   F U N C T I O N S 
#
#

getTxArrivalTimes <- function(data,txVector) {
  
  if ("EventType" %in% names(data)) { # Event Log is used
    return(data %>% 
             filter(EventType == "Event_NewTransactionArrival", ObjectID %in% txVector) %>% 
             select(Simulation = SimID,Transaction = ObjectID, Time = SimTime))
  } else { # Input data is used
    return(data %>% 
             filter(TxID %in% txVector) %>% 
             select(Simulation = SimID,Transaction = TxID, Time = ArrivalTime))
    
  }
}

timeAlignBeliefData <- function(beliefData, arrivalTimes, offset = 0) {
  
  if (missing(beliefData)){
    rlang::abort("Argument `beliefData` is required and was not supplied.")
  }
  
  if (missing(arrivalTimes)){
    rlang::abort("Argument `arrivalTimes` is required and was not supplied.")
  }
  
  
  augData <- beliefData %>% 
    inner_join(arrivalTimes, by = c("Simulation", "Transaction"),
               suffix = c("", "_arrival"))  %>%
    filter(Time > Time_arrival)
  
  return (augData %>%
    group_by(Simulation, Transaction) %>%
    mutate(Time = (Time - min(Time) + offset)) %>%
    ungroup())
}



prepareGraphData <- function(beliefData, VaR = FALSE, 
                             boot = FALSE, R_Boot = 10000, 
                             alignTimes = FALSE, 
                             timeAlignmentOffset = 0,
                             arrivalTimes) {
  
  if (nrow(beliefData) == 0) {
    rlang::abort("beliefData does not contain rows. Are txList transactions included in the set?")
  }
  
  if (alignTimes) {
    if (missing(arrivalTimes)) {
      rlang::abort("For timeAlign = TRUE, argument `arrivalTimes` is required and was not supplied.")
    } else {
      beliefData <- timeAlignBeliefData(beliefData, arrivalTimes, timeAlignmentOffset) 
    }
  }
  
  
  # Bounding dataset 
  print(paste0("[",as.character(sys.call(0)[[1]]),"]: Calculating time span."))
  endTime = min(beliefData %>% group_by(Simulation) %>% summarise(maxTime = max(Time),.groups = "drop") %>% select(maxTime))
  net = beliefData %>% filter(Time<=endTime)

  print(paste0("[",as.character(sys.call(0)[[1]]),"]: Grouping by Simulation ID."))  
  if (boot) {
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Bootstrapping enabled."))  
    if (VaR) {
      confs = net %>% group_by(Time,Transaction) %>%
        summarise(avgConf = mean(Belief), 
                  sdConf = sd(Belief), 
                  medConf = median(Belief),
                  .groups = "drop") %>% 
        mutate(
          lwr = boot.lower(avgConf,mean,R_Boot),upr = boot.upper(avgConf,mean,R_Boot),
          VaR = quantile(avgConf,0.05, type = 7))
    } else {
      confs = net %>% group_by(Time,Transaction) %>% 
        summarise(avgConf = mean(Belief), 
                  sdConf = sd(Belief), 
                  medConf = median(Belief),
                  .groups = "drop") %>%
        mutate(
          lwr = boot.lower(avgConf,mean,R_Boot),upr = boot.upper(avgConf,mean,R_Boot))
    }
    
  } else {
    
    if (VaR) {
      part1 = net %>% group_by(Time,Transaction) %>%
        summarise(avgConf = mean(Belief), 
                  sdConf = sd(Belief), 
                  medConf = median(Belief),
                  VaR = quantile(Belief, 0.05),
                  .groups = "drop")  %>%
        mutate(Transaction = as.factor(Transaction)) 
      
      part2 = part1 %>%
        group_by(Time) %>%
        summarise(
          Transaction = "VaR",
          avgConf     = min(VaR),
          sdConf      = NA_real_,      # optional, can leave NA
          medConf     = NA_real_,
          .groups = "drop"
        ) %>%
        mutate(Transaction = as.factor(Transaction)) 
      
      confs <- bind_rows(part1, part2) %>%
        arrange(Time, Transaction) %>% select(-VaR)
      
    } else {
      confs = net %>% group_by(Time,Transaction) %>%
        summarise(avgConf = mean(Belief), 
                  sdConf = sd(Belief), 
                  medConf = median(Belief),
                  .groups = "drop")
    }

  }
  
  print(paste0("[",as.character(sys.call(0)[[1]]),"]: Preparing dataset."))
  confs2plot <- confs %>% 
    mutate(Transaction = as.factor(Transaction)) %>%
    mutate(Time = Time)
  
  return(confs2plot) 
  
}

getBeliefGraph <- function(graphData,threshold, faceted = FALSE, VaR = "include",
                           xlims,timeUnit = "min",txVector) {
  
  # Parameter Processing and Validation
  if (!(timeUnit %in% c("min","sec","ms"))) {
    rlang::abort(paste0("timeUnit '",timeUnit, "' not recognized."))
  }
  
  timeDivider <- switch(timeUnit, min = 60000, sec = 1000, ms = 1)
  
  graphData <- graphData %>% mutate(Time = Time/timeDivider)
  
  if (missing(xlims)){
    xlim_min = 0
    xlim_max = ceiling(max(graphData$Time))
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: x axis bounds set to (",xlim_min,",",xlim_max,") ", timeUnit))
  } else {
    xlim_min = auto_time_to_ms(xlims[1])/timeDivider
    xlim_max = auto_time_to_ms(xlims[2])/timeDivider
  }
  
  if (!missing(txVector)) {
    graphData <- graphData %>% filter(Transaction %in% txVector)
  }
  
  if ((nrow(graphData) == 0) || missing(graphData)) {
    rlang::abort("graphData does not contain rows or is missing.")
  }
  

  if (!missing("VaR")) {
    
    if (!(VaR %in% c("include","exclude","only"))) {
      rlang::abort("Uknown value for paramter VaR. Choose one of 'include','exclude','only'")  
    }
    
    if (VaR == "exclude") {
      graphData <- graphData %>% filter(Transaction != "VaR")
    }
    if (VaR == "only"){
      graphData <- graphData %>% filter(Transaction == "VaR")
    }
  }
  
    
  # Get existing levels

  print(paste0("[",as.character(sys.call(0)[[1]]),"]: Producing graph."))
  
  if (faceted) {
    p1 <- ggplot(data = graphData, aes(x = Time, y=avgConf, label = avgConf)) + 
      geom_line() + 
      facet_wrap(~ Transaction, scales = "free_y") +
      ylab("Confidence") + 
      xlab(paste0("Time (",timeUnit,")")) +
      xlim(xlim_min,xlim_max) +
      ggtitle(label = "Confidence in f over time")
  } else {
    p1 <- ggplot(data = graphData, aes(x = Time, y=avgConf, label = avgConf)) + 
      geom_line(aes(color = Transaction)) + 
      ylab("Confidence") + xlab(paste0("Time (",timeUnit,")")) +
      xlim(xlim_min,xlim_max) +
      ggtitle(label = "Confidence in f over time")
  }
  
  
  if (!missing(threshold)) {
    p1 = p1 + geom_hline(yintercept = threshold) +   
      annotate(
      "text",
      x = xlim_min,
      y = threshold,
      label = "Belief threshold",
      vjust = -0.5,
      hjust = 0,                  # flush left
      color = "darkred"
    )
  }
  
  if (rlang::has_name(graphData,"VaR")) {
    p1 = p1 + geom_line(aes(y = VaR, color = "VaR"), linetype="dashed")
  }
  
  if (rlang::has_name(graphData,"lwr")) {
    p2 <- ggplot(data = graphData, aes(x = Time, y=avgConf, label = avgConf)) + 
      geom_ribbon(aes(ymin=lwr,ymax=upr,fill = "Confidence Interval"),alpha=0.2) + 
      ylab("Confidence") + xlab(paste0("Time (",timeUnit,")")) +
      xlim(xlim_min,xlim_max) +
      ggtitle(label = "Confidence in f over time") +
      scale_color_manual(name = "Legend", values = c("VaR" = "red")) +
      scale_fill_manual(name = "Legend", values = c("Confidence interval" = "blue")) + 
    if (!missing(threshold)) {
      p2 = p2 + geom_hline(yintercept = threshold)
    }
    if (rlang::has_name(graphData,"VaR")) {
      p2 = p2 +  geom_line(aes(y = VaR, color = "VaR"), linewidth = 1)
    }
    return(list(txGraph = p1, bootGraph = p2)) 
  } else {
    return(graph = p1) 
  }
}



getFinality <- function(
    beliefData,
    horizon,
    alignTimes = TRUE,
    arrivalTimes,
    timeAlignmentOffset = 0,
    threshold = 0.9,
    thresTest_p = 0.95,
    alpha = 0.05,
    txVector) {

  
  if (alignTimes) {
    if (missing(arrivalTimes)) {
      rlang::abort("For timeAlign = TRUE, argument `arrivalTimes` is required and was not supplied.")
    } else {
      beliefData <- timeAlignBeliefData(beliefData, arrivalTimes, timeAlignmentOffset) 
    }
  }
  
  
  h = max(beliefData$Time)
  if (missing(horizon) || (auto_time_to_ms(horizon) > h)) {
    horizon = h
  } else {
    horizon = auto_time_to_ms(horizon)
  }
  

  
  
  result <- beliefData %>%
    # Step 1: keep only times <= horizon
    filter(Time <= horizon) %>%
    
    # Step 2: for each Transaction, find the maximum Time
    group_by(Simulation,Transaction) %>%
    filter(Time == max(Time)) %>%
    group_by(Transaction) %>%
    # Step 3: calculate mean Belief over simulations for that Transaction × Time
    summarise(
      Time_ms = first(Time),                    # the max time
      sampleMean = mean(Belief, na.rm = TRUE),
      sampleSD = sd(Belief, na.rm = TRUE),
      n = n(),
      n_above = sum(Belief > threshold, na.rm = TRUE),
      #thres_test = list(binom.test(n_above,n,p=thresTest_p,alternative = "greater")),
      #test = list(binom.test(n_above,n)),
      .groups = "drop"
    ) %>% 
    mutate(
      thres_test = map2(n_above, n,
                        ~ binom.test(.x, .y, p = thresTest_p, alternative = "greater")),
      test       = map2(n_above, n,
                        ~ binom.test(.x, .y))
    ) %>% 
    mutate(Time = format_simtime(Time_ms), 
           threstest_p_val      = map_dbl(thres_test, ~ .x$p.value),
           threstest_confIntLow = map_dbl(thres_test, ~ .x$conf.int[1]),
           threstest_confIntHigh= map_dbl(thres_test, ~ .x$conf.int[2]),
           threstest_Final      = threstest_p_val <= alpha,
           point_estimate       = map_dbl(test, ~ .x$estimate),
           ci_low               = map_dbl(test, ~ .x$conf.int[1]),
           ci_high              = map_dbl(test, ~ .x$conf.int[2])
           ) %>% 
    select(Transaction, Time, Time_ms,
           sampleMean, sampleSD,
           trials = n,
           successes = n_above,
           threstest_p_val, threstest_confIntLow, threstest_confIntHigh, threstest_Final,
           point_estimate,ci_low,ci_high
           )
  
  if (missing(txVector)) {
    return(result)
  } else {
    return(result %>% filter(Transaction %in% txVector))
  }
  
}


getTimeToFinality <- function(graphData,threshold) {
  result <- graphData %>%
    filter(!is.na(avgConf)) %>%                # ignore NA rows
    filter(avgConf > threshold) %>%                   # exceed threshold
    group_by(Transaction) %>%
    summarise(
      first_time = format_simtime(min(Time)),
      .groups = "drop"
    )
  return(result)
}

getCulpritSims <- function(
    beliefData,
    horizon,
    alignTimes = TRUE,
    arrivalTimes,
    timeAlignmentOffset = 0,
    threshold = 0.9,
    txVector) {
  
  
  if (alignTimes) {
    if (missing(arrivalTimes)) {
      rlang::abort("For timeAlign = TRUE, argument `arrivalTimes` is required and was not supplied.")
    } else {
      beliefData <- timeAlignBeliefData(beliefData, arrivalTimes, timeAlignmentOffset) 
    }
  }
  
  
  h = max(beliefData$Time)
  if (missing(horizon) || (auto_time_to_ms(horizon) > h)) {
    horizon = h
  } else {
    horizon = auto_time_to_ms(horizon)
  }
  
  result <- beliefData %>%
    # Step 1: keep only times <= horizon
    filter(Time <= horizon) %>%
    
    # Step 2: for each Transaction, find the maximum Time
    group_by(Simulation,Transaction) %>%
    filter(Time == max(Time)) %>%
    filter(Belief <= threshold) %>% 
    select(Transaction, Simulation, Time, Belief) %>%
    arrange(Transaction)

    
  if (missing(txVector)) {
    return(result)
  } else {
    return(result %>% filter(Transaction %in% txVector))
  }
}


getTransactionHistory <- function(
    eventData,
    blockData,
    txVector,
    simVector) {
  
  
  # Function to map EventType to text
  event_text <- function(type) {
    if(type == "Event_NewTransactionArrival") {
      "(*arrival*) arrives at node"
    } else if(type == "Event_TransactionPropagation") {
      "propagates to node"
    } else if(type == "Node Completes Validation") {
      "(*validation*) appears in a block that was just validated by Node "
    } else if(type == "Node Receives Propagated Block") {
      "is received in a validated block by Node "
    } else if(type == "Appended On Chain (w/ parent)") {
      "(*belief*) is appended on local chain by receiving Node "
    } else if(type == "Appended On Chain (parentless)") {
      "(*belief*) is appended on local chain by validating Node "
    } else {
      type  # fallback
    }
  }
  
  # History according to EventLog
  ed_out <- tibble()
  if (!missing(eventData)) {
    
    eD = eventData %>% filter(ObjectID %in% txVector) %>% arrange(SimTime)
    
    eD_out <- eD %>%
      arrange(SimID, SimTime, EventID) %>%
      rowwise() %>%
      mutate(
        Time = SimTime,
        Transaction = ObjectID,
        Event = paste0(event_text(EventType), " ", NodeID)
      ) %>%
      ungroup() %>%
      select(SimID, Time, Transaction, Event)  # keep only the columns you want
  }
  
  # History according to BlockLog
  bD_out <- tibble() 
  if (!missing(blockData)) {
    
    for (tx in txVector) {
      bD <- blockData %>%
        rowwise() %>%
        filter(length(intersect(as.integer(str_split(str_remove_all(BlockContent, "[{}]"), ";")[[1]]), tx)) > 0) %>%
        ungroup()
      
      bD_out <- bD_out %>%
        bind_rows(
          bD %>%
            arrange(SimID, SimTime) %>%
            rowwise() %>%
            mutate(
              Time = SimTime,
              Transaction = tx,
              Event = paste0(event_text(EventType)," ", NodeID)
            ) %>%
            ungroup() %>%
            select(SimID, Time, Transaction, Event)  # keep only the columns you want
        )
    }
    
  }
  
  all_out <- bind_rows(bD_out, eD_out)
  
  if (!missing(simVector)) {
    all_out <- all_out %>% filter(SimID %in% simVector)  
  }
  
  
  for (sim in unique(all_out$SimID)) {
    
    cat("Simulation", sim, ":\n")
    cat(paste0("    HH:MM:SS.mmm\n"))
    all_out %>%
      filter(SimID == sim) %>%
      arrange(Time, Transaction) %>%  # optional: order events
      rowwise() %>%
      mutate(
        line = paste0("    ", format_simtime(Time), ": Transaction ", Transaction, " ", Event)
      ) %>%
      pull(line) %>%
      paste(collapse = "\n") %>%
      cat("\n")
    
    cat("\n")  # extra line between simulations
  }
}




# concatenate_results.R
# Concatenates CSV/TXT files with matching prefixes across multiple result folders.
#
# Usage:
#   source("concatenate_results.R")
#   concatenate_results(
#     dir_names = c("folder1", "folder2", ...),
#     output_folder_name = "folderALL",
#     base_dir = getwd()
#   )

# concatenate_results.R
# Concatenates CSV/TXT files with matching prefixes across multiple result folders.
#
# Usage:
#   source("concatenate_results.R")
#   concatenate_results(
#     dir_names = c("folder1", "folder2", ...),
#     output_folder_name = "folderALL",
#     input_dir = "path/to/results"
#   )
concatenate_results <- function(dir_names, output_folder_name, input_dir) {
  output_dir <- file.path(input_dir, output_folder_name)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Discover files from the first directory
  first_dir <- file.path(input_dir, dir_names[1])
  all_files <- list.files(first_dir)
  
  # Extract prefix (everything before the first " - ")
  get_prefix <- function(filename) str_replace(filename, " - .*$", "")
  prefixes <- unique(map_chr(all_files, get_prefix))
  
  cat("Found prefixes:", paste(prefixes, collapse = ", "), "\n")
  cat("Output folder:", output_dir, "\n\n")
  
  for (prefix in prefixes) {
    cat("Processing:", prefix, "... ")
    
    all_lines <- character(0)
    is_first <- TRUE
    
    for (dir_name in dir_names) {
      dir_path <- file.path(input_dir, dir_name)
      files_in_dir <- list.files(dir_path, full.names = TRUE)
      match <- files_in_dir[str_starts(basename(files_in_dir), fixed(str_c(prefix, " - ")))]
      
      if (length(match) == 0) {
        warning(str_glue("No file with prefix {prefix} found in {dir_name}"))
        next
      }
      
      lines <- read_lines(match[1])
      
      if (is_first) {
        all_lines <- lines
        is_first <- FALSE
      } else if (length(lines) > 1) {
        # Skip header line, append data lines only
        all_lines <- c(all_lines, lines[-1])
      }
    }
    
    if (length(all_lines) == 0) {
      cat("SKIPPED (no files found)\n")
      next
    }
    
    # Preserve original file extension
    original_file <- list.files(first_dir, pattern = str_c("^", fixed(prefix), " - "))
    ext <- str_extract(original_file[1], "\\.[^.]+$")
    output_file <- file.path(output_dir, str_c(prefix, " - ", output_folder_name, ext))
    
    write_lines(all_lines, output_file)
    cat(length(all_lines) - 1, "data lines written\n")
  }
  
  cat("\nDone.\n")
}
# ── Example usage ──────────────────────────────────────────────────────────────
# concatenate_results(
#   dir_names = c(
#     "confirmation.attack.0.40.4.1.200 - 2026.03.12 17.59.17",
#     "confirmation.attack.0.40.4.201.400 - 2026.03.12 18.09.19",
#     "confirmation.attack.0.40.4.401.600 - 2026.03.12 18.19.51",
#     "confirmation.attack.0.40.4.601.800 - 2026.03.12 18.30.11",
#     "confirmation.attack.0.40.4.801.1000 - 2026.03.12 18.40.36"
#   ),
#   output_folder_name = "confirmation.attack.0.40.4.ALL",
#   input_dir = "F:/4. Software/cnsim-bitcoin/examples/results"
# )
