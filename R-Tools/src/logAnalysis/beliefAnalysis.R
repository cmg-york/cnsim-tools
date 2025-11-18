library(tidyverse)
library(nptest)



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
#  M A I N   F U N C T I O N S 
#
#

getTxArrivalTimes <- function(eventData,txVector) {
  return(eventData %>% 
    filter(EventType == "Event_NewTransactionArrival", Object %in% txVector) %>% 
    select(Simulation = SimID,Transaction = Object, Time = SimTime))
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
      confs = net %>% group_by(Time,Transaction) %>%
        summarise(avgConf = mean(Belief), 
                  sdConf = sd(Belief), 
                  medConf = median(Belief),
                  .groups = "drop") %>% 
        mutate(
          VaR = quantile(avgConf,0.05))
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

getBeliefGraph <- function(graphData,threshold, faceted = FALSE,
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
    p1 = p1 + geom_hline(yintercept = threshold)
  }
  
  if (rlang::has_name(graphData,"VaR")) {
    p1 = p1 + geom_line(aes(y = VaR, color = "VaR"))
  }
  
  if (rlang::has_name(graphData,"lwr")) {
    p2 <- ggplot(data = graphData, aes(x = Time, y=avgConf, label = avgConf)) + 
      geom_ribbon(aes(ymin=lwr,ymax=upr,fill = "Confidence Interval"),alpha=0.2) + 
      ylab("Confidence") + xlab(paste0("Time (",timeUnit,")")) +
      xlim(xlim_min,xlim_max) +
      ggtitle(label = "Confidence in f over time") +
      scale_color_manual(name = "Legend", values = c("VaR" = "red")) +
      scale_fill_manual(name = "Legend", values = c("Confidence interval" = "blue"))
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
      n_above = sum(Belief >= threshold),
      thres_test = list(binom.test(n_above,n,p=thresTest_p,alternative = "greater")),
      test = list(binom.test(n_above,n)),
      .groups = "drop"
    ) %>% 
    mutate(Time = format_simtime(Time_ms), 
           threstest_p_val = thres_test[[1]]$p.value,
           threstest_confIntLow = thres_test[[1]]$conf.int[1],
           threstest_confIntHigh = thres_test[[1]]$conf.int[2],
           threstest_Final = (threstest_p_val <= alpha),
           point_estimate = test[[1]]$estimate,
           ci_low = test[[1]]$conf.int[1],
           ci_high = test[[1]]$conf.int[2]
           ) %>% 
    select(Transaction, Time, Time_ms,
           sampleMean, sampleSD,
           threstest_p_val, threstest_confIntLow, threstest_confIntHigh, threstest_Final,
           point_estimate,ci_low,ci_high
           )
  
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
    
    eD = eventData %>% filter(Object %in% txVector) %>% arrange(SimTime)
    
    eD_out <- eD %>%
      arrange(SimID, SimTime, EventID) %>%
      rowwise() %>%
      mutate(
        Time = SimTime,
        Transaction = Object,
        Event = paste0(event_text(EventType), " ", Node)
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
              Event = paste0(event_text(EvtType)," ", NodeID)
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






getBeliefGraph(gr,0.9,timeUnit = "min",xlims = c(0,"20:90.000"))$txGraph


#
#
# D E P R E C A T E D
#
#

getBeliefGraphOLD <- function(beliefData, 
                           VaR = FALSE,
                           boot = FALSE, R_Boot = 10000,
                           threshold,
                           xlims,
                           timeUnit = "min",
                           alignTimes = FALSE,
                           timeAlignmentOffset = 0,
                           arrivalTimes,
                           txVector
) {
  
  
  # Parameter Processing and Validation
  
  if (!(timeUnit %in% c("min","sec","ms"))) {
    rlang::abort(paste0("timeUnit '",timeUnit, "' not recognized."))
  }
  
  timeDivider <- switch(timeUnit, min = 60000, sec = 1000, ms = 1)
  
  if (missing(xlims)){
    xlim_min = 0
    xlim_max = ceiling(endTime/timeDivider)
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: x axis bounds set to (",xlim_min,",",xlim_max,") ", timeUnit))
  } else {
    xlim_min = xlims[1]
    xlim_max = xlims[2]
  }
  
  if (!missing(txVector)) {
    beliefData <- beliefData %>% filter(Transaction %in% txVector)
  }
  
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
  endTime = min(beliefData %>% group_by(Simulation) %>% summarise(maxTime = max(Time)) %>% select(maxTime))
  net = beliefData %>% filter(Time<=endTime)
  
  if (boot) {
    
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Grouping by Simulation ID. Bootstraping for confidence band."))
    confs = net %>% group_by(Time,Transaction) %>% 
      summarise(avgConf = mean(Belief), sd = sd(Belief), medConf = median(Belief)) %>%
      mutate(
        lwr = boot.lower(avgConf,mean,R_Boot),upr = boot.upper(avgConf,mean,R_Boot),
        VaR = quantile(avgConf,0.05))
    
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Preparing dataset."))
    confs2plot <- confs %>% mutate(Transaction = as.factor(Transaction)) %>%
      mutate(Time = Time / timeDivider)
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Producing graph."))
    p1 <- ggplot(data = confs2plot, aes(x = Time, y=avgConf, label = avgConf)) + 
      geom_line(aes(color = Transaction)) + 
      ylab("Confidence") + xlab(paste0("Time (",timeUnit,")")) +
      xlim(xlim_min,xlim_max) +
      ggtitle(label = "Confidence in f over time")
    
    if (!missing(threshold)) {
      p1 = p1 + geom_hline(yintercept = threshold)
    }
    
    p2 <- ggplot(data = confs2plot, aes(x = Time, y=avgConf, label = avgConf)) + 
      geom_line(aes(y = VaR, color = "VaR")) + 
      geom_ribbon(aes(ymin=lwr,ymax=upr,fill = "Confidence Interval"),alpha=0.1) + 
      ylab("Confidence") + xlab(paste0("Time (",timeUnit,")")) +
      xlim(xlim_min,xlim_max) +
      ggtitle(label = "Confidence in f over time") +
      scale_color_manual(name = "Legend", values = c("VaR" = "red")) +
      scale_fill_manual(name = "Legend", values = c("Confidence interval" = "lightblue"))
    
    
    
    p2 <- ggplot(data = confs2plot, aes(x = Time)) +
      # Ribbon for the confidence interval
      geom_ribbon(
        aes(ymin = lwr, ymax = upr, fill = "Confidence interval"),
        alpha = 0.2
      ) +
      
      # Line for VaR
      geom_line(
        aes(y = VaR, color = "VaR"),
        linewidth = 1
      ) +
      
      # Axis labels and title
      ylab("Confidence") +
      xlab(paste0("Time (", timeUnit, ")")) +
      xlim(xlim_min, xlim_max) +
      ggtitle("Confidence in f over time") +
      
      # Manual color and fill scales — *names must match* the strings in aes()
      scale_color_manual(
        name   = "Legend",
        values = c("VaR" = "red")
      ) +
      scale_fill_manual(
        name   = "Legend",
        values = c("Confidence interval" = "blue")
      )
    
    
    if (!missing(threshold)) {
      p2 = p2 + geom_hline(yintercept = threshold)
    }
    
    return(list(txGraph = p1, VaRGraph = p2, data = confs2plot)) 
    
  } else {
    
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Grouping by Simulation ID"))
    
    if (VaR) {
      confs = net %>% group_by(Time,Transaction) %>%
        summarise(avgConf = mean(Belief), sd = sd(Belief), medConf = median(Belief)) %>% 
        mutate(
        VaR = quantile(avgConf,0.05))
    } else {
      confs = net %>% group_by(Time,Transaction) %>%
        summarise(avgConf = mean(Belief), sd = sd(Belief), medConf = median(Belief))
    }
    
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Preparing dataset."))
    confs2plot <- confs %>% mutate(Transaction = as.factor(Transaction)) %>%
      mutate(Time = Time / timeDivider)
    
    print(paste0("[",as.character(sys.call(0)[[1]]),"]: Producing graph."))
    p <- ggplot(data = confs2plot, aes(x = Time, y=avgConf)) + 
      geom_line(aes(color = Transaction)) + 
      ylab("Confidence") + xlab(paste0("Time (",timeUnit,")")) +
      ggtitle(label = "Confidence in f over time") + 
      xlim(xlim_min,xlim_max)
    
    if (VaR) {
      p <- p + geom_line(
        aes(y = VaR, color = "VaR"),
        linewidth = 1
      ) + 
        # Legend titles
        scale_color_discrete(name = "Transaction") +
        scale_linetype_manual(
          name = "Reference",
          values = c("VaR" = "solid")
        )
    }

    
    
    if (!missing(threshold)) {
      p = p + geom_hline(yintercept = threshold)
    }
    return(list(graph = p, data = confs2plot %>% mutate(Time = Time*timeDivider))) 
  }
}



redrawGrapOLD <- function(data, 
                        faceted = FALSE,
                        timeUnit,
                        xmin, xmax) {
  
  if (missing(xmin)) {
    xmin = min(data$Time)
  }
  
  if (missing(xmax)) {
    xmax = max(data$Time)
  }
  
  if (!(timeUnit %in% c("min","sec","ms"))) {
    rlang::abort(paste0("timeUnit '",timeUnit, "' not recognized."))
  }

  timeDivider <- switch(timeUnit, min = 60000, sec = 1000, ms = 1)
  
  if (faceted) {
    ggplot(data, aes(x = Time/timeDivider, y = avgConf)) +
      geom_line() +
      facet_wrap(~ Transaction, scales = "free_y") +
      ylab("Confidence") +
      xlab(paste0("Time (", timeUnit, ")")) +
      ggtitle("Confidence in f over time") +
      xlim(xmin,xmax)
  } else {
    ggplot(data = data, aes(x = Time/timeDivider, y=avgConf)) + 
      geom_line(aes(color = Transaction)) + 
      ylab("Confidence") + 
      xlab(paste0("Time (",timeUnit,")")) + 
      ggtitle(label = "Confidence in f over time")+
      xlim(xmin,xmax)
  }
}




#
# A D V A N C E D
#


#Aggregate over simulations
# confs = net %>% group_by(Time,Transaction) %>% 
#   summarise(avgConf = mean(conf), sd = sd(conf), medConf = median(conf), 
#             lwr = boot.lower(conf,mean,Bootstrap_R),upr = boot.upper(conf,mean,Bootstrap_R),
#             VaR = quantile(conf,0.05))
# 
# confs2plot <- confs %>% mutate(Condition = as.factor(Transaction))
# 
# 
# ggplot(data = confs2plot, aes(x = Time, y=avgConf)) + 
#   geom_line(aes(color = Condition)) + 
#   geom_ribbon(aes(ymin=lwr,ymax=upr,fill = Condition),alpha=0.1) + 
#   ylab("Confidence") +
#   ggtitle(label = "Confidence in f over time (mean and 95% CI)") + 
#   geom_hline(yintercept = 0.8)
# 
# 
# View(confs2plot)
