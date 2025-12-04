library(tidyverse)
library(knitr)
library(TOSTER)

source("beliefAnalysis.R")

#
# Produce Pace Data
#

folder.10 = "../../../../cnsim-bitcoin/examples/results/faithful/"  

# 5 Node
#df_5node = "faithful.5node.30sim - 2025.11.04 12.07.10"
#producePaceData(folder, df_5node)
#Sys.sleep(1)

# 10 Node
df_10Node_1 = "faithful.10node.30sim.1-5 - 2025.12.03 17.58.53"
producePaceData(folder.10, df_10Node_1)
df_10Node_2 = "faithful.10node.30sim.6-10 - 2025.12.03 18.03.06"
producePaceData(folder.10, df_10Node_2)
df_10Node_3 = "faithful.10node.30sim.11-15 - 2025.12.03 18.03.32"
producePaceData(folder.10, df_10Node_3)
df_10Node_4 = "faithful.10node.30sim.16-20 - 2025.12.03 18.03.56"
producePaceData(folder.10, df_10Node_4)
df_10Node_5 = "faithful.10node.30sim.21-25 - 2025.12.03 18.04.34"
producePaceData(folder.10, df_10Node_5)
df_10Node_6 = "faithful.10node.30sim.26-30 - 2025.12.03 18.12.24"
producePaceData(folder.10, df_10Node_6)


runTimes <- rbind(
  getRunTime(folder.10, df_10Node_1),
  getRunTime(folder.10, df_10Node_2),
  getRunTime(folder.10, df_10Node_3),
  getRunTime(folder.10, df_10Node_4),
  getRunTime(folder.10, df_10Node_5),
  getRunTime(folder.10, df_10Node_6)
  )

runTimes %>% summarise(`Run Time (mean)` = format_simtime(as.integer(mean(sysTime))), 
                       `Run Time (sd)` =  format_simtime(as.integer(sd(sysTime))))

# 30 Node
folder.30 = "../../../../cnsim-bitcoin/examples/results/lowrate/"  

df_30Node_1 = "lowrate.30node.100sim.1-20 - 2025.12.03 19.27.10"
producePaceData(folder.30, df_30Node_1)
df_30Node_2 = "lowrate.30node.100sim.21-40 - 2025.12.03 19.27.35"
producePaceData(folder.30, df_30Node_2)
df_30Node_3 = "lowrate.30node.100sim.41-60 - 2025.12.03 19.27.57"
producePaceData(folder.30, df_30Node_3)
df_30Node_4 = "lowrate.30node.100sim.61-80 - 2025.12.03 19.29.59"
producePaceData(folder.30, df_30Node_4)
df_30Node_5 = "lowrate.30node.100sim.81-100 - 2025.12.03 19.32.21"
producePaceData(folder.30, df_30Node_5)


#paste0(folder,df_10Node_1,"/EventLog - ",df_10Node_1,".csv")

#df_10Node_2 = "faithful.10node.30sim - 2025.11.05 14.09.34"
#df_10Node_3 = "faithful.10node.30sim - 2025.11.05 14.10.12"
#Sys.sleep(1)
#producePaceData(folder, df_10Node_2)
#Sys.sleep(1)
#producePaceData(folder, df_10Node_3)
#Sys.sleep(1)

#
# Compare with authoritative
#

auth = 8.521
low_eqbound = -1.5
high_eqbound = 1.5

# 5 Node

# pacedata.5 = read_csv(paste0(folder,df_5node,"/PaceData - ",df_5node,".csv"))
# 
# pace.5 = pacedata.5 %>%
#   arrange(SimID, SimDiff) %>%             # ensure ordered
#   summarise(
#     `Block Time (mins) - Mean....:` = mean(SimDiff)/60000,
#     `Block Time (mins) - St. Dev.:`   = sd(SimDiff)/60000
#   )
# t(pace.5)
# 
# 
# distances.5 = pacedata.5$SimDiff/60000
# 
# 
# TOSTone.raw(
#         m = mean(distances.5),
#         sd = sd(distances.5),
#         n = length(distances.5),
#         mu = auth,
#         low_eqbound = low_eqbound,
#         high_eqbound = high_eqbound,
#         alpha = 0.05)


# 10 Node

pacedata.10 = rbind(read_csv(paste0(folder,df_10Node_1,"/PaceData - ",df_10Node_1,".csv")),
read_csv(paste0(folder,df_10Node_2,"/PaceData - ",df_10Node_2,".csv")),
read_csv(paste0(folder,df_10Node_3,"/PaceData - ",df_10Node_3,".csv")),
read_csv(paste0(folder,df_10Node_4,"/PaceData - ",df_10Node_4,".csv")),
read_csv(paste0(folder,df_10Node_5,"/PaceData - ",df_10Node_5,".csv")),
read_csv(paste0(folder,df_10Node_6,"/PaceData - ",df_10Node_6,".csv")))

pace.10 = pacedata.10 %>%
  arrange(SimID, SimDiff) %>%             # ensure ordered
  summarise(
    `Block Time (mins) - Mean....:` = mean(SimDiff)/60000,
    `Block Time (mins) - St. Dev.:`   = sd(SimDiff)/60000
  )
t(pace.10)



distances.10 = pacedata.10$SimDiff/60000


TOSTone.raw(
  m = mean(distances.10),
  sd = sd(distances.10),
  n = length(distances.10),
  mu = auth,
  low_eqbound = low_eqbound,
  high_eqbound = high_eqbound,
  alpha = 0.05)



# 30 Node

pacedata.30 = rbind(read_csv(paste0(folder.30,df_30Node_1,"/PaceData - ",df_30Node_1,".csv")),
                    read_csv(paste0(folder.30,df_30Node_2,"/PaceData - ",df_30Node_2,".csv")),
                    read_csv(paste0(folder.30,df_30Node_3,"/PaceData - ",df_30Node_3,".csv")),
                    read_csv(paste0(folder.30,df_30Node_4,"/PaceData - ",df_30Node_4,".csv")),
                    read_csv(paste0(folder.30,df_30Node_5,"/PaceData - ",df_30Node_5,".csv")))

pace.30 = pacedata.30 %>%
  arrange(SimID, SimDiff) %>%             # ensure ordered
  summarise(
    `Block Time (mins) - Mean....:` = mean(SimDiff)/60000,
    `Block Time (mins) - St. Dev.:`   = sd(SimDiff)/60000
  )
t(pace.30)



distances.30 = pacedata.30$SimDiff/60000


TOSTone.raw(
  m = mean(distances.30),
  sd = sd(distances.30),
  n = length(distances.30),
  mu = auth,
  low_eqbound = low_eqbound,
  high_eqbound = high_eqbound,
  alpha = 0.05)


