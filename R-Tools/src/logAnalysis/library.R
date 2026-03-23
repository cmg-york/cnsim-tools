#' Load and install required libraries
#'
#' Checks if each library is available; if not, installs it from CRAN
#' and then loads it.
#'
#' @param libs Character vector of package names to load.
#' @return Called for side effects (loading packages). Returns invisible NULL.
#' @export
load_libraries <- function(libs) {
  for (lib in libs) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib)
      library(lib, character.only = TRUE)
    }
  }
}

load_libraries(c("tidyverse", "nptest", "knitr", "xtable","TOSTER", "data.table"))


#
#
#  H E L P E R S  
#
#


#' Bootstrap lower confidence bound
#'
#' Computes the lower bound of a BCa bootstrap confidence interval
#' using non-parametric bootstrap.
#'
#' @param x Numeric vector of data.
#' @param y A function (statistic) to apply to the data.
#' @param Iter Integer number of bootstrap replicates.
#' @return Numeric lower bound of the BCa confidence interval; returns 0 if NaN.
#' @export
boot.lower <- function(x,y,Iter){
  npbs <- np.boot(x = x, statistic = y, R = Iter)
  return(ifelse(is.nan(npbs$bca[2,1]),0,npbs$bca[2,1]))
}

#' Bootstrap upper confidence bound
#'
#' Computes the upper bound of a BCa bootstrap confidence interval
#' using non-parametric bootstrap.
#'
#' @param x Numeric vector of data.
#' @param y A function (statistic) to apply to the data.
#' @param Iter Integer number of bootstrap replicates.
#' @return Numeric upper bound of the BCa confidence interval; returns 0 if NaN.
#' @export
boot.upper <- function(x,y,Iter){
  npbs <- np.boot(x = x, statistic = y, R = Iter)
  return(ifelse(is.nan(npbs$bca[2,2]),0,npbs$bca[2,2]))
  return(npbs$bca[2,2])
}

#' Format simulation time as HH:MM:SS.mmm
#'
#' Converts a time value in milliseconds to a human-readable string
#' in the format \code{HH:MM:SS.mmm}.
#'
#' @param simtime_ms Numeric time in milliseconds.
#' @return Character string formatted as \code{"HH:MM:SS.mmm"}.
#' @export
format_simtime <- function(simtime_ms) {
  total_seconds <- simtime_ms / 1000
  hours   <- floor(total_seconds / 3600)
  minutes <- floor((total_seconds %% 3600) / 60)
  seconds <- floor(total_seconds %% 60)
  millis  <- round((total_seconds - floor(total_seconds)) * 1000)
  
  sprintf("%02d:%02d:%02d.%03d", hours, minutes, seconds, millis)
}

#' Auto-detect time format and convert to milliseconds
#'
#' Accepts numeric values (treated as milliseconds), pure integer strings,
#' or time strings in \code{SS.mmm}, \code{MM:SS.mmm}, or \code{HH:MM:SS.mmm}
#' format and converts them to integer milliseconds.
#'
#' @param x A numeric value or character string representing a time.
#' @return Integer time in milliseconds.
#' @export
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

#' Nakamoto attacker success probability
#'
#' Calculates the probability that an attacker with hash-power fraction \code{q}
#' can catch up and double-spend after \code{z} confirmations, using the
#' original Nakamoto formula.
#'
#' @param q Numeric attacker's fraction of total network hash power (0 to 1).
#' @param z Integer number of confirmations.
#' @return Numeric probability of a successful attack.
#' @export
attacker_success_probability <- function(q, z) {
  if (q <= 0) return(0)
  if (q >= 1) return(1)
  
  p <- 1 - q
  lambda <- z * (q / p)
  
  sum <- 1.0
  
  for (k in 0:z) {
    poisson <- dpois(k, lambda)   # numerically stable Poisson probability
    sum <- sum - poisson * (1 - (q / p)^(z - k + 1))
  }
  
  return(sum)
}

#' Modified Nakamoto catch-up probability
#'
#' Calculates the probability that an attacker can surpass the honest chain
#' after z confirmations. Differs from the original Nakamoto algorithm in that
#' the attacker starts at a deficit of z+1 (rather than z) and must strictly
#' surpass the honest chain (rather than merely equal it).
#'
#' @param q Attacker's fraction of total network hash power (0 < q < 0.5)
#' @param z Number of confirmations
#' @return Probability that the attacker successfully surpasses the honest chain
nakamoto_catchup_confirmations <- function(q, z) {
  # Honest network's fraction of hash power
  p <- 1.0 - q
  
  # Expected number of blocks the attacker mines during z confirmations
  lambda <- z * q / p
  
  # Step probability ratio for the random-walk catch-up model
  ratio <- q / p
  
  # Poisson PMF terms P(X = k) for k = 0, 1, ..., z+1
  # where X is the number of blocks the attacker has mined
  poisson <- dpois(0:(z + 1), lambda)
  
  # Component 1 — P_sim: probability the attacker has already mined more than
  # z+1 blocks, placing them strictly ahead without further catch-up needed.
  # This is the upper tail of the Poisson distribution: P(X > z+1)
  p_sim <- 1.0 - sum(poisson)
  
  # Component 2 — P_tail: for each k <= z+1, the attacker has a remaining
  # deficit of (z - k + 2) blocks to overcome. The probability of closing
  # this gap via a random walk is (q/p)^(z - k + 2), weighted by the
  # Poisson probability of having mined exactly k blocks.
  k <- 0:(z + 1)
  p_tail <- sum(poisson * ratio^(z - k + 2))
  
  # Total success probability is the sum of both components
  p_total <- p_sim + p_tail
  return(p_total)
}



#' TOST equivalence test for proportions (console output)
#'
#' Performs a two one-sided tests (TOST) procedure to assess whether an
#' observed proportion is statistically equivalent to a theoretical value
#' within a margin of \code{delta}. Prints PASS/FAIL result to the console.
#'
#' @param trials Integer total number of trials.
#' @param successes Integer number of successes observed.
#' @param theoretical Numeric theoretical proportion to test against.
#' @param delta Numeric equivalence margin (default 0.03).
#' @param alpha Numeric significance level (default 0.05).
#' @return Called for side effects (prints result). Returns invisible NULL.
#' @export
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
  

#' TOST equivalence test for proportions (tibble output)
#'
#' Vectorized version of the TOST equivalence test. Returns a tibble with
#' bounds, confidence intervals, p-values, and PASS/FAIL result instead
#' of printing to the console.
#'
#' @param trials Integer total number of trials.
#' @param successes Integer number of successes observed.
#' @param theoretical Numeric theoretical proportion to test against.
#' @param delta Numeric equivalence margin (default 0.03).
#' @param alpha Numeric significance level (default 0.05).
#' @return A \code{\link[tibble]{tibble}} with columns: \code{lower_bound},
#'   \code{ci_low}, \code{observed}, \code{theoretical}, \code{ci_high},
#'   \code{upper_bound}, \code{p_lower}, \code{p_upper}, \code{result}, \code{delta}.
#' @export
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


#' Produce pace data from event log
#'
#' Reads the event log for a given experiment, filters for container validation
#' events, computes inter-event time differences per simulation, and writes
#' the result to a CSV file.
#'
#' @param dataFolder Character path to the base data directory.
#' @param experiment Character name of the experiment.
#' @return Called for side effects (writes CSV). Returns invisible NULL.
#' @export
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


#' Summarize pace statistics
#'
#' Reads pace data CSVs for one or more experiments and computes the mean
#' and standard deviation of block time (in minutes), along with the
#' total sample count.
#'
#' @param dataFolder Character path to the base data directory.
#' @param experiment Character vector of one or more experiment names.
#' @return A tibble with mean block time, standard deviation, and sample count.
#' @export
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



#' Simulate theoretical block times (inverse-CDF method)
#'
#' Generates \code{n} random block times from a geometric distribution
#' using the inverse-CDF (log) transform of uniform random variates,
#' scaled by mining power.
#'
#' @param n Integer number of samples to generate.
#' @param difficulty Numeric mining difficulty (reciprocal of success probability).
#' @param power Numeric mining hash rate (trials per time unit).
#' @return Numeric vector of simulated block times.
#' @export
pace_theoretical <- function(n, difficulty,power) {
  p <- 1 / difficulty
  u <- runif(n)
  trials = floor(log1p(-u) / log1p(-p))
  return (trials / power)
}

#' Simulate theoretical block times (rgeom method)
#'
#' Generates \code{n} random block times from a geometric distribution
#' using R's built-in \code{\link[stats]{rgeom}}, scaled by mining power.
#'
#' @param n Integer number of samples to generate.
#' @param difficulty Numeric mining difficulty (reciprocal of success probability).
#' @param power Numeric mining hash rate (trials per time unit).
#' @return Numeric vector of simulated block times.
#' @export
pace_theoretical2 <- function(n, difficulty,power) {
  trials = rgeom(n, p = 1/difficulty)
  return (trials / power)
}



#
#
# SIMULATION PERFORMANCE
#
#

#' Get simulation run times
#'
#' Reads the event log and computes the wall-clock run time for each
#' simulation as the difference between the maximum and minimum system
#' timestamps.
#'
#' @param dataFolder Character path to the base data directory.
#' @param experiment Character name of the experiment.
#' @return A tibble with columns \code{SimID}, \code{sysTime} (ms),
#'   and \code{sysTime_formated} (HH:MM:SS.mmm string).
#' @export
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

#' Get transaction arrival times
#'
#' Extracts arrival times for specified transactions from either an event log
#' or input data frame, depending on which columns are present.
#'
#' @param data A data frame — either an event log (with \code{EventType} column)
#'   or input data (with \code{TxID} column).
#' @param txVector Character or integer vector of transaction IDs to filter.
#' @return A tibble with columns \code{Simulation}, \code{Transaction}, and \code{Time}.
#' @export
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

#' Time-align belief data to transaction arrival times
#'
#' Joins belief data with arrival times, filters to keep only observations
#' after each transaction's arrival, and re-zeros the time axis so that
#' \code{Time = 0} corresponds to each transaction's arrival (plus an
#' optional offset).
#'
#' @param beliefData A data frame with columns \code{Simulation}, \code{Transaction},
#'   \code{Time}, and \code{Belief}.
#' @param arrivalTimes A data frame with columns \code{Simulation}, \code{Transaction},
#'   and \code{Time} (as returned by \code{\link{getTxArrivalTimes}}).
#' @param offset Numeric offset in milliseconds added after re-zeroing (default 0).
#' @return A tibble with the same structure as \code{beliefData}, time-aligned.
#' @export
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



#' Prepare aggregated graph data from belief data
#'
#' Aggregates raw belief data across simulations by computing mean, standard
#' deviation, and median confidence at each time point per transaction.
#' Optionally adds bootstrap confidence intervals, Value-at-Risk (5th
#' percentile), and time alignment to arrival times.
#'
#' @param beliefData A data frame with columns \code{Simulation}, \code{Transaction},
#'   \code{Time}, and \code{Belief}.
#' @param VaR Logical; if \code{TRUE}, compute the 5th-percentile VaR of belief
#'   and append a synthetic \code{"VaR"} transaction (default \code{FALSE}).
#' @param boot Logical; if \code{TRUE}, compute BCa bootstrap confidence
#'   intervals for the mean (default \code{FALSE}).
#' @param R_Boot Integer number of bootstrap replicates (default 10000).
#' @param alignTimes Logical; if \code{TRUE}, time-align belief data to
#'   transaction arrival times before aggregation (default \code{FALSE}).
#' @param timeAlignmentOffset Numeric offset in ms added after re-zeroing
#'   (default 0).
#' @param arrivalTimes A data frame of arrival times (required when
#'   \code{alignTimes = TRUE}).
#' @return A tibble with columns \code{Time}, \code{Transaction} (factor),
#'   \code{avgConf}, \code{sdConf}, \code{medConf}, and optionally
#'   \code{lwr}, \code{upr}, \code{VaR}.
#' @export
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

#' Produce belief confidence graph(s)
#'
#' Generates one or two \pkg{ggplot2} plots of average confidence over time,
#' with optional faceting by transaction, threshold line, VaR overlay, and
#' bootstrap ribbon. Returns a single plot or a named list of plots when
#' bootstrap data is present.
#'
#' @param graphData A data frame as returned by \code{\link{prepareGraphData}}.
#' @param threshold Numeric belief threshold; if supplied, a horizontal line
#'   and label are added to the plot.
#' @param faceted Logical; if \code{TRUE}, facet the plot by transaction
#'   (default \code{FALSE}).
#' @param VaR Character controlling VaR display: \code{"include"} (default),
#'   \code{"exclude"}, or \code{"only"}.
#' @param xlims A two-element vector of time limits (in any format accepted
#'   by \code{\link{auto_time_to_ms}}). If omitted, limits are derived from data.
#' @param timeUnit Character time unit for the x-axis: \code{"min"} (default),
#'   \code{"sec"}, or \code{"ms"}.
#' @param txVector Optional vector of transaction IDs to include in the plot.
#' @return A \code{ggplot} object, or a named list with elements \code{txGraph}
#'   and \code{bootGraph} when bootstrap data is present.
#' @export
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



#' Compute finality statistics per transaction
#'
#' For each transaction, evaluates belief at the end of the observation
#' horizon across all simulations. Computes the sample mean, standard
#' deviation, proportion of simulations exceeding the belief threshold,
#' and a one-sided binomial test for whether that proportion meets
#' \code{thresTest_p}.
#'
#' @param beliefData A data frame with columns \code{Simulation},
#'   \code{Transaction}, \code{Time}, and \code{Belief}.
#' @param horizon Time horizon (any format accepted by
#'   \code{\link{auto_time_to_ms}}). If omitted or larger than available
#'   data, defaults to the maximum time in the dataset.
#' @param alignTimes Logical; if \code{TRUE}, time-align belief data to
#'   arrival times (default \code{TRUE}).
#' @param arrivalTimes A data frame of arrival times (required when
#'   \code{alignTimes = TRUE}).
#' @param timeAlignmentOffset Numeric offset in ms (default 0).
#' @param threshold Numeric belief threshold for counting successes
#'   (default 0.9).
#' @param thresTest_p Numeric null-hypothesis proportion for the binomial
#'   test (default 0.95).
#' @param alpha Numeric significance level (default 0.05).
#' @param txVector Optional vector of transaction IDs to filter results.
#' @return A tibble with per-transaction finality statistics including
#'   sample mean, SD, trial/success counts, binomial test p-value,
#'   confidence intervals, and a PASS/FAIL indicator.
#' @export
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


#' Get time to finality per transaction
#'
#' Finds the earliest time at which the average confidence exceeds the
#' given threshold for each transaction.
#'
#' @param graphData A data frame as returned by \code{\link{prepareGraphData}},
#'   with columns \code{Time}, \code{Transaction}, and \code{avgConf}.
#' @param threshold Numeric belief threshold.
#' @return A tibble with columns \code{Transaction} and \code{first_time}
#'   (formatted as HH:MM:SS.mmm).
#' @export
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

#' Find simulations where belief stays below the threshold
#'
#' Identifies individual simulations in which the final belief for a
#' transaction does not exceed the specified threshold by the end of
#' the observation horizon.
#'
#' @param beliefData A data frame with columns \code{Simulation},
#'   \code{Transaction}, \code{Time}, and \code{Belief}.
#' @param horizon Time horizon (any format accepted by
#'   \code{\link{auto_time_to_ms}}).
#' @param alignTimes Logical; if \code{TRUE}, time-align to arrival times
#'   (default \code{TRUE}).
#' @param arrivalTimes A data frame of arrival times (required when
#'   \code{alignTimes = TRUE}).
#' @param timeAlignmentOffset Numeric offset in ms (default 0).
#' @param threshold Numeric belief threshold (default 0.9).
#' @param txVector Optional vector of transaction IDs to filter results.
#' @return A tibble with columns \code{Transaction}, \code{Simulation},
#'   \code{Time}, and \code{Belief} for the culprit simulation/transaction
#'   pairs.
#' @export
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


#' Print human-readable transaction event history
#'
#' Combines event log and block log data to produce a chronological,
#' human-readable narrative of each transaction's lifecycle (arrival,
#' propagation, validation, chain append) for selected simulations.
#' Output is printed to the console.
#'
#' @param eventData A data frame of event log data (with \code{EventType},
#'   \code{ObjectID}, \code{SimTime}, \code{NodeID}, \code{SimID} columns).
#' @param blockData A data frame of block log data (with \code{BlockContent},
#'   \code{EventType}, \code{SimTime}, \code{NodeID}, \code{SimID} columns).
#' @param txVector Integer or character vector of transaction IDs to trace.
#' @param simVector Optional vector of simulation IDs to filter.
#' @return Called for side effects (prints to console). Returns invisible NULL.
#' @export
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
#' Concatenate CSV/TXT result files across experiment folders
#'
#' Discovers files by prefix in the first directory, then concatenates
#' matching files from all listed directories (skipping duplicate headers).
#' Writes combined output to a new folder.
#'
#' @param dir_names Character vector of folder names containing result files.
#' @param output_folder_name Character name for the output folder.
#' @param input_dir Character path to the parent directory containing
#'   \code{dir_names}.
#' @return Called for side effects (writes concatenated files). Returns
#'   invisible NULL.
#' @export
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


#' Load all experiment data into a named list
#'
#' Reads the block log, event log, input data, and belief log for a given
#' experiment and returns them bundled together with metadata in a list.
#'
#' @param outputFolder Character path to the base results directory.
#' @param experiment Character name of the experiment.
#' @param txVector Vector of transaction IDs of interest.
#' @return A named list with elements: \code{outputFolder}, \code{experiment},
#'   \code{blockData}, \code{eventData}, \code{inputData}, \code{beliefData},
#'   and \code{txVector}.
#' @export
setVars <- function(outputFolder, experiment, txVector) {
  # Read the CSV files
  blockData <- safe_read_csv(paste0(outputFolder, experiment, "/BlockLog - ", experiment, ".csv"))
  eventData <- read_csv(paste0(outputFolder, experiment, "/EventLog - ", experiment, ".csv"))
  inputData <- read_csv(paste0(outputFolder, experiment, "/Input - ", experiment, ".csv"))
  beliefData <- read_csv(paste0(outputFolder, experiment, "/BeliefLogShort - ", experiment, ".csv")) %>%
    rename(Simulation = SimID,
           Time = `Time (ms from start)`,
           Transaction = `Transaction ID`)
  
  # Build a list indexed by code
  result <- list(
    outputFolder = outputFolder,
    experiment = experiment,
    blockData = blockData,
    eventData = eventData,
    inputData = inputData,
    beliefData = beliefData,
    txVector = txVector
  )
  
  return(result)
}




#' Safely read a CSV file
#'
#' Wraps \code{\link[readr]{read_csv}} in a \code{tryCatch} so that file-read
#' errors return \code{NULL} instead of stopping execution.
#'
#' @param path Character file path to the CSV.
#' @return A tibble if the file is read successfully, or \code{NULL} on error.
#' @export
safe_read_csv <- function(path) {
  tryCatch(
    readr::read_csv(path),
    error = function(e) NULL
  )
}

