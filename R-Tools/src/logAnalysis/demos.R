source("beliefAnalysis.R")


#
#
# LOW RATE 1
#
#

outputRelFolder = "../../../../cnsim-bitcoin/examples/results/"
experiment = "lowrate.5node.30sim - 2025.11.13 20.06.40"
blockData = read_csv(paste0(outputRelFolder,experiment,"/BlockLog - ",experiment,".csv"))
eventData = read_csv(paste0(outputRelFolder,experiment,"/EventLog - ",experiment,".csv"))
beliefData = read_csv(paste0(outputRelFolder,experiment,"/BeliefLogShort - ",experiment,".csv")) %>% 
  rename(Simulation = SimID, Time = `Time (ms from start)`, Transaction = `Transaction ID`)
txVector = c(200,210,220,230,240,250,260,270.280,290)


#
# Investigate transaction history for the txVector
#
getTransactionHistory(eventData,blockData,txVector = c(210), simVector = 5)


#
# Generate belief graph
#
lowrate1.plain = getBeliefGraph(beliefData, 
                                timeUnit = "min",
                                threshold = 0.8,
                                alignTimes = TRUE,
                                arrivalTimes = getTxArrivalTimes(eventData,txVector),
                                txVector = txVector)

lowrate1.plain$graph
redrawGraph(lowrate1.plain$data, 
        faceted = TRUE, 
        timeUnit = "min",
        0, 40)


getFinality(beliefData, alignTimes = TRUE, arrivalTimes = getTxArrivalTimes(eventData,txVector), 
            horizon = "19:00.000")
lowrate1.plain = getBeliefGraph(beliefData, VaR=TRUE, timeUnit = "min",threshold = 0.8,txVector = txVector)


#lowrate1.boot = getBeliefGraph(beliefData, VaR=TRUE, boot=TRUE, R_Boot = 1000, timeUnit = "min",threshold = 0.8,txVector = txVector)



#
#
# LOW RATE 2 - Fee differences
#
#

outputRelFolder = "../../examples/results/"
experiment = "lowrate.5node.30sim.txvalue - 2025.11.14 12.45.15"
blockData = read_csv(paste0(outputRelFolder,experiment,"/BlockLog - ",experiment,".csv"))
eventData = read_csv(paste0(outputRelFolder,experiment,"/EventLog - ",experiment,".csv"))
beliefData = read_csv(paste0(outputRelFolder,experiment,"/BeliefLogShort - ",experiment,".csv")) %>% 
  rename(Simulation = SimID, Time = `Time (ms from start)`, Transaction = `Transaction ID`)
txVector = c(200,210,220,230,240,250,260,270.280,290)

getTransactionHistory(eventData,blockData,txVector = c(290))
lowrate2.plain = getBeliefGraph(beliefData, timeUnit = "min",threshold = 0.8,txVector = txVector)
redrawGraph(lowrate2.plain$data,xmax = 20,faceted = TRUE)
getFinality(beliefData,txVector = txVector)



#
#
# LOW RATE 3 - Network Slowness
#
#

outputRelFolder = "../../../../cnsim-bitcoin/examples/results/"
experiment = "lowrate.5node.30sim.slownet - 2025.11.14 14.57.57"
blockData = read_csv(paste0(outputRelFolder,experiment,"/BlockLog - ",experiment,".csv"))
eventData = read_csv(paste0(outputRelFolder,experiment,"/EventLog - ",experiment,".csv"))
beliefData = read_csv(paste0(outputRelFolder,experiment,"/BeliefLogShort - ",experiment,".csv")) %>% 
  rename(Simulation = SimID, Time = `Time (ms from start)`, Transaction = `Transaction ID`)
txVector = seq(1000,9000,by=1000)

lowrate2.plain = getBeliefGraph(beliefData, timeUnit = "min",threshold = 0.8,txVector = txVector)
lowrate2.plain$graph
redrawGraph(lowrate2.plain$data,faceted = TRUE)



# Redesign getFinality to take into account the entry of the Transaction
getFinality(beliefData,txVector = txVector,horizon = "17:12")

# Clean-Up Events that are not transaction related
getTransactionHistory(eventData,blockData,txVector = c(2000))

# Complete PACING package

# Caclualte Branching factor from StructureLog

# Visualize from StructureLog


#
#
# FAITHFUL 1
#
#



outputRelFolder = "../../examples/results/"
experiment_part1 = "faithful.10node.30sim - 2025.11.05 14.08.34"
experiment_part2 = "faithful.10node.30sim - 2025.11.05 14.09.34"
experiment_part3 = "faithful.10node.30sim - 2025.11.05 14.10.12"

beliefData <- 
  rbind(
    read_csv(paste0(outputRelFolder,experiment_part1,"/BeliefLogShort - ",experiment_part1,".csv")) %>% rename(Simulation = SimID, Transaction = `Transaction ID`,Time = `Time (ms from start)`),
    read_csv(paste0(outputRelFolder,experiment_part2,"/BeliefLogShort - ",experiment_part2,".csv")) %>% rename(Simulation = SimID, Transaction = `Transaction ID`,Time = `Time (ms from start)`),
    read_csv(paste0(outputRelFolder,experiment_part3,"/BeliefLogShort - ",experiment_part3,".csv")) %>% rename(Simulation = SimID, Transaction = `Transaction ID`,Time = `Time (ms from start)`)
  )


eventData <- 
  rbind(
    read_csv(paste0(outputRelFolder,experiment_part1,"/EventLog - ",experiment_part1,".csv")), 
    read_csv(paste0(outputRelFolder,experiment_part2,"/EventLog - ",experiment_part2,".csv")),
    read_csv(paste0(outputRelFolder,experiment_part3,"/EventLog - ",experiment_part3,".csv")) 
  )


txVector = c(200,210,220,230,240,250,260,270.280,290)

getTransactionHistory(eventData,txVector = c(210))



faithful1.plain = getBeliefGraph(beliefData, boot=FALSE, R_Boot = 10, timeUnit = "min",threshold = 0.8,txVector = txVector)


faithful1.plain$graph

faithful1.band = getBeliefGraph(beliefData, boot=TRUE, R_Boot = 100, xlims = c(0,180))


