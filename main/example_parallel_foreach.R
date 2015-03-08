# Simple example for working loops in parallel on cluster or on systems with more than 1 CPU
# install.packages('doParallel')
library(doParallel)

# create R cluster; Please register number of CPU you want to use here.
cl<-makeCluster(2)
registerDoParallel(cl)

# Work task in foreach loop in parallel: %dopar% 
# For working in sequential mode: %do%

cluster_test <- foreach(i=1:1E6) %dopar% {sqrt(i)}

# Information about employed CPU workers and closing the cluster
workers <- getDoParWorkers()
print('CPU workers employed: ')
print(workers)
stopCluster(cl)
