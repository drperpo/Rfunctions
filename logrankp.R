# this is a function to compute the log-rank test of the differencies in survival per pair
# when you have more than two groups
# 22-2-2006

log.rankp <- function(time,status,x){
# Runs a log-rank test for survival differencies per 
# pairs 
# Original version written by Aris Perperoglou, for
# LUMC, Leiden University, The Netherlands
#22 Feb 2006.
i <- 1:(length(table(x))-1)
j <- 2:length(table(x))
xn <-substitute(x)
if(length(table(x))==1) stop("There is only one group, cannot compare")
if(length(table(x))==2) print(survdiff(Surv(time,status)~x)) else
for(i in i){
    for(j in j){
        cat("Testing",paste(xn,i,sep="")," with",paste(xn,j,sep=""),"\n")
       print( survdiff(Surv(time[x==i|x==j],status[x==i|x==j])~x[x==i|x==j]))
        cat("============================","\n")}}}


#log.rankp(time,death,grade)
