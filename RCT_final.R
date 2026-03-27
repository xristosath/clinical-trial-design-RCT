set.seed(123)
OC_3plus3= function(p, doses=NULL){
  if (is.null(doses)) {
    doses= paste("Dose",seq_along(p))
  }
  Pi= dbinom(0,size = 3, prob = p)+dbinom(1,size=3,prob=p)*dbinom(0,size = 3, prob = p)
  Qi=cumprod(Pi)
  OC=1-Qi
  data.frame(
    dose= doses,
    p=p,
    Pi=Pi,
    Qi= Qi,
    OC= OC
  )
}
res_ScenA=OC_3plus3(c(0.08,0.18,0.30,0.45))
res_ScenB=OC_3plus3(c(0.12,0.22,0.28,0.40))
res_ScenA
res_ScenB

#MTD Calc
#we want to calculate the difference OC[i]-OC[i-1] 
#aka the probability of stopping the study exactly at dose i

#Scenario A
stop_probA= diff(c(0,res_ScenA$OC))
stop_probA
MTD_probs_A= c(
  None= stop_probA[1], #stop at dose 1 -> no MTD
  Dose1=stop_probA[2], #stop at dose 2 -> MTD Dose= 1
  Dose2=stop_probA[3], # stop at Dose 3 -> MTD Dose= 2
  Dose3= stop_probA[4], # stop at Dose 4 -> MTD Dose= 3
  Dose4= 1- res_ScenA$OC[4] #No stopping at all -> MTD= Dose 4
)
MTD_probs_A

# Scenario B
stop_probs_B = diff(c(0, res_ScenB$OC))
MTD_probs_B = c(
  None  = stop_probs_B[1],
  Dose1 = stop_probs_B[2],
  Dose2 = stop_probs_B[3],
  Dose3 = stop_probs_B[4],
  Dose4 = 1 - res_ScenB$OC[4]
)

#results
round(MTD_probs_A, 3)
round(MTD_probs_B, 3)

##### Askisi 2 ####
library(clinfun)
#Scenario A
single_stage= function(p0,p1,alpha,power) {
  n=1
  while(TRUE) {
    #find the smallest k : P(X>k|p0)<=alpha
    kvec=0:n
    ptail0= pbinom(kvec,size = n, prob = p0, lower.tail = F)
    idx= which(ptail0<=alpha)[1]
    
    if (!is.na(idx)) {
      k_crit=kvec[idx] #if X>k_crit we reject H0
      #we calculate the power P(X>k_crit|p1)
      ptail1= pbinom(k_crit,size = n,prob = p1,lower.tail = F)
      
      if (ptail1>=power) {
        return(c(n=n,k=k_crit,alpha_exact=ptail0[idx],power_exact=ptail1))
      }
    }
    n=n+1
  }
}

res_single_A <-single_stage(0.15, 0.40, 0.10, 0.80)
round(res_single_A,3)
res_single_B <- single_stage(0.20, 0.40, 0.05, 0.90)
round(res_single_B,3)

results_table <- rbind(res_single_A, res_single_B)


### Two Stage Design ###
library(clinfun)
#Scenario A
simon_A= ph2simon(pu=0.15,pa=0.40,ep1 = 0.10,ep2 = 0.20,nmax = 100)
#Scenario B
simon_B=ph2simon(pu=0.20,pa=0.40,ep1 = 0.05,ep2 = 0.10,nmax = 100)

#Results
simon_A
simon_B



##find power if p=0.3
p_true=0.3
# Scenario A
# Parameters: n1=7, r1=1, Total_N=18, Total_r=4
# Rule: Stop if Stage 1 responses (x1) <= 1. Continue if x1 > 1.
# Logic: We sum the probabilities ONLY for cases where the trial continues (x1 from 2 to 7).
n1=7
r1=1
N=18
r=4
n2=N-n1
x1_val= (r1+1):n1
power_A=sum(dbinom(x1_val,n1,p_true) # probability of getting enough responses in stage 2 to reach total r>4
            *pbinom(r-x1_val,n2,p_true,lower.tail =F )) #we need x2< (4-x1) from the remaining 18-7 patients
round(power_A,3)

#Scenario B 
# Parameters: n1=24, r1=5, Total_N=45, Total_r=13
# We stop if Stage 1 responses (x1) <= 5. Continue if x1 > 5
# Logic: The loop runs for x1 = 6 to 24 
n1 = 24
r1 = 5
N  = 45
r  = 13
n2 = N - n1
x1_val= (r1+1):n1
power_B <- sum(
  # Probability of observing x1 responses in Stage 1 (x1 must be > 5 to continue)
  dbinom(x1_val, n1,p_true) * # Probability of success in Stage 2 given x1
    # We need total responses > 13, so x2 > (13 - x1)
    pbinom(r-x1_val, n2, p_true, lower.tail=F)
)
round(power_B,3)


### Askisi 3 ####
set.seed(123) #for reproducibility
N=80
prob_exact=dbinom(x=N/2,size = N,p=0.5)
round(prob_exact,4)

# Simulation 
#Simulation Parameters
N=80
K=1000

results= data.frame(
  replication = 1:K,
  NumberA= NA_integer_
)

#loop
for (i in 1:K) {
  u= runif(N)
  #We directly calculate the number in Arm A.
  results$NumberA[i]=sum(u<0.5)
}
ci_95= quantile(results$NumberA/N, probs = c(0.025,0.975))
round(ci_95,4)

#Scenario B Block Randomization
#Block size=4
library(blockrand)
block_rand=blockrand(
  n=N,
  num.levels = 2,
  levels = c("A","B"),
  block.sizes = 2 # num.levels*block.sizes=4
)
first_20 =head(block_rand,20)
# Check Balance
table(block_rand$treatment)

wrong_size_test <- blockrand(
  n = N, 
  num.levels = 2, 
  levels = c("A", "B"), 
  block.sizes = 3 
)

#check block.size number
head(wrong_size_test, 10)


#### Askisi ####
library(gsDesign)
#Scenario A O-F
scenA= gsDesign(
  k=3,
  test.type = 2, # two sided symmetric
  alpha = 0.025,
  beta = 0.2,
  timing = c(0.3,0.6,1),
  sfu = "OF"
)
scenA
gsBoundSummary(scenA)

#Scenario B Pockock
scenB=gsDesign(
  k=3,
  test.type = 2,
  alpha = 0.025,
  beta = 0.2,
  timing = c(0.3,0.6,1),
  sfu = "Pocock"
)
scenB
gsBoundSummary(scenB)
