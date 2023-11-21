
#------------------------------------------------------------------------------#
# Bloomsbury Policy Lab 2023, last updated: 8th of Sept 2023
# Cost-effectiveness of mSTR versus comparator for MDR/RR TB in EECA
#------------------------------------------------------------------------------#

library(truncnorm)
library(ggExtra)
library(ggplot2)
library(markovchain)
library(tidyverse)
library(shape)
library(diagram)
library(Hmisc)
library(yarrr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)

#------------------------------------------------------------------------------#
#General parameters
#------------------------------------------------------------------------------#
n_t <- 240 # number of years * by months, 20y * 12
n_states <- 8 #number of health states
#n_c <- 1000 #number of patients
Pregular_mortalityR <- 0.0001 #regular mortality rate
discount <- 0.03 #discount rate
discount_m <- discount/12 #discount rate per month
N_psa <- 1000 #number of repetitions for probability sensitivity analyses
#------------------------------------------------------------------------------#
##I. DETERMINISTIC ANALYSIS                                                   #
#------------------------------------------------------------------------------#
##CREATING MODEL FUNCTION; we can change the "params" and run it again per country.
model <- function(.params) { 
  with(.params, {
    n_t <- 240
    n_s <- 8
    v_state_names <- c("mSTR","SLtreatment","LostFollowUp","Unresolved","TreatmentCompleted","AdverseEffects","Cured", "Dead")
    m_p <- matrix(0, nrow= 8, ncol=8, dimnames= list(from= v_state_names, to = v_state_names))
    
    #prob matrix   
    m_p["mSTR", "mSTR"] <- 1  -p_mSTR_AdverseEffects -p_mSTR_SLtreatment -p_mSTR_LostFollowUp - p_mSTR_TreatmentCompleted -p_mSTR_Dead
    m_p["mSTR", "SLtreatment"] <- p_mSTR_SLtreatment
    m_p["mSTR", "LostFollowUp"] <- p_mSTR_LostFollowUp
    m_p["mSTR", "Unresolved"] <- p_mSTR_Unresolved
    m_p["mSTR", "TreatmentCompleted"] <- p_mSTR_TreatmentCompleted 
    m_p["mSTR", "AdverseEffects"] <- p_mSTR_AdverseEffects
    m_p["mSTR", "Cured"] <- p_mSTR_Cured
    m_p["mSTR", "Dead"] <-  p_mSTR_Dead
    
    m_p["SLtreatment", "mSTR"] <- p_SLtreatment_mSTR
    m_p["SLtreatment", "SLtreatment"] <- 1 -p_SLtreatment_Unresolved -p_SLtreatment_TreatmentCompleted -( p_SLtreatment_Dead )
    m_p["SLtreatment", "LostFollowUp"] <- p_SLtreatment_LostFollowUp
    m_p["SLtreatment", "Unresolved"] <- p_SLtreatment_Unresolved
    m_p["SLtreatment", "TreatmentCompleted"] <- p_SLtreatment_TreatmentCompleted
    m_p["SLtreatment", "AdverseEffects"] <- p_SLtreatment_AdverseEffects
    m_p["SLtreatment", "Cured"] <- p_SLtreatment_Cured
    m_p["SLtreatment", "Dead"] <-  p_SLtreatment_Dead 
    
    m_p["LostFollowUp", "mSTR"] <- p_LostFollowUp_mSTR
    m_p["LostFollowUp", "SLtreatment"] <- 0
    m_p["LostFollowUp", "LostFollowUp"] <- 1 - p_LostFollowUp_mSTR -( p_LostFollowUp_Dead)
    m_p["LostFollowUp", "Unresolved"] <- 0
    m_p["LostFollowUp", "TreatmentCompleted"] <-0
    m_p["LostFollowUp", "AdverseEffects"] <-  0
    m_p["LostFollowUp", "Cured"] <- 0
    m_p["LostFollowUp", "Dead"] <-  p_LostFollowUp_Dead
    
    m_p["Unresolved", "mSTR"] <- 0
    m_p["Unresolved", "SLtreatment"] <- 0
    m_p["Unresolved", "LostFollowUp"] <- 0  
    m_p["Unresolved", "Unresolved"] <- 1 - (p_Unresolved_Dead)
    m_p["Unresolved", "TreatmentCompleted"] <- 0
    m_p["Unresolved", "AdverseEffects"] <- 0
    m_p["Unresolved", "Cured"] <- 0
    m_p["Unresolved", "Dead"] <-  (0.75)*p_Unresolved_Dead
    
    m_p["TreatmentCompleted", "mSTR"] <- 0
    m_p["TreatmentCompleted", "SLtreatment"] <- p_TreatmentCompleted_SLtreatment
    m_p["TreatmentCompleted", "LostFollowUp"] <- 0
    m_p["TreatmentCompleted", "Unresolved"] <-  0
    m_p["TreatmentCompleted", "TreatmentCompleted"] <-  1- p_TreatmentCompleted_SLtreatment -p_TreatmentCompleted_Cured -(p_TreatmentCompleted_Dead)
    m_p["TreatmentCompleted", "AdverseEffects"] <-  0
    m_p["TreatmentCompleted", "Cured"] <-  p_TreatmentCompleted_Cured
    m_p["TreatmentCompleted", "Dead"] <-   p_TreatmentCompleted_Dead
    
    m_p["AdverseEffects", "mSTR"] <- 0
    m_p["AdverseEffects", "SLtreatment"] <- p_AdverseEffects_SLtreatment
    m_p["AdverseEffects", "LostFollowUp"] <- 0
    m_p["AdverseEffects", "Unresolved"] <-  p_AdverseEffects_Unresolved
    m_p["AdverseEffects", "TreatmentCompleted"] <-  0
    m_p["AdverseEffects", "AdverseEffects"] <-  1- p_AdverseEffects_Unresolved - p_AdverseEffects_SLtreatment -(p_AdverseEffects_Dead)
    m_p["AdverseEffects", "Cured"] <-  0
    m_p["AdverseEffects", "Dead"] <-   p_AdverseEffects_Dead
    
    m_p["Cured", "mSTR"] <- 0
    m_p["Cured", "SLtreatment"] <- 0
    m_p["Cured", "LostFollowUp"] <- 0
    m_p["Cured", "Unresolved"] <- 0
    m_p["Cured", "TreatmentCompleted"] <- 0
    m_p["Cured", "AdverseEffects"] <- 0
    m_p["Cured", "Cured"] <-  1 - (p_Cured_Dead)
    m_p["Cured", "Dead"] <-   p_Cured_Dead
    
    m_p["Dead", "mSTR"] <- 0
    m_p["Dead", "SLtreatment"] <- 0
    m_p["Dead", "LostFollowUp"] <- 0
    m_p["Dead", "Unresolved"] <- 0
    m_p["Dead", "TreatmentCompleted"] <- 0
    m_p["Dead", "AdverseEffects"] <- 0
    m_p["Dead", "Cured"] <- 0
    m_p["Dead", "Dead"] <-  1
    
    #Start drafting the markov matrix
    state_membership <- array(NA_real_, dim= c(n_t, n_s), dimnames= list(cycle = 1:n_t, state = v_state_names))
    #Start filling in first row; only with treatment box.
    state_membership[1,]= c(n_c, 0, 0, 0, 0, 0, 0, 0)
    for (i in 2:n_t) {state_membership[i, ] <- state_membership[i-1, ] %*% m_p
    }
    #Insert costs per healthy, diseased and death, let's wait for Tom
    m_payoffs <- matrix(0, nrow = 8, ncol = 2, dimnames= list(state= v_state_names, 
                                                              payoff= c("Costs", "QALYs")))
    
    m_payoffsCost<- matrix(0, nrow=8, ncol=1,
                           dimnames= list(state= v_state_names, 
                                          payoff= c("Costs")))
    m_payoffsQALY <- matrix(0, nrow=8, ncol=1,
                            dimnames= list(state= v_state_names, 
                                           payoff= c("QALYs")))
    
    m_payoffs["mSTR",  "Costs"] <-  Cost_mSTR
    m_payoffs["SLtreatment", "Costs"] <-  Cost_SLtreatment
    m_payoffs["LostFollowUp", "Costs"] <- Cost_LostFollowUp
    m_payoffs["Unresolved", "Costs"] <- Cost_Unresolved
    m_payoffs["TreatmentCompleted", "Costs"] <- Cost_TreatmentCompleted
    m_payoffs["AdverseEffects", "Costs"] <-  Cost_AdverseEffects
    m_payoffs["Cured", "Costs"] <- Cost_Cured
    m_payoffs["Dead", "Costs"] <- Cost_Dead
    
    m_payoffsCost["mSTR",  "Costs"] <-  Cost_mSTR
    m_payoffsCost["SLtreatment", "Costs"] <-  Cost_SLtreatment
    m_payoffsCost["LostFollowUp", "Costs"] <- Cost_LostFollowUp
    m_payoffsCost["Unresolved", "Costs"] <- Cost_Unresolved
    m_payoffsCost["TreatmentCompleted", "Costs"] <- Cost_TreatmentCompleted
    m_payoffsCost["AdverseEffects", "Costs"] <-  Cost_AdverseEffects
    m_payoffsCost["Cured", "Costs"] <- Cost_Cured
    m_payoffsCost["Dead", "Costs"] <- Cost_Dead
    
    m_payoffs["mSTR",  "QALYs"] <-  qaly_mSTR/12
    m_payoffs["SLtreatment", "QALYs"] <-  qaly_SLtreatment/12
    m_payoffs["LostFollowUp", "QALYs"] <- qaly_LostFollowUp/12
    m_payoffs["Unresolved", "QALYs"] <- qaly_Unresolved/12
    m_payoffs["TreatmentCompleted", "QALYs"] <- qaly_TreatmentCompleted/12
    m_payoffs["AdverseEffects", "QALYs"] <-  qaly_AdverseEffects/12
    m_payoffs["Cured", "QALYs"] <- qaly_Cured/12
    m_payoffs["Dead", "QALYs"] <- qaly_Dead/12
    
    m_payoffsQALY["mSTR",  "QALYs"] <-  qaly_mSTR/12
    m_payoffsQALY["SLtreatment", "QALYs"] <-  qaly_SLtreatment/12
    m_payoffsQALY["LostFollowUp", "QALYs"] <- qaly_LostFollowUp/12
    m_payoffsQALY["Unresolved", "QALYs"] <- qaly_Unresolved/12
    m_payoffsQALY["TreatmentCompleted", "QALYs"] <- qaly_TreatmentCompleted/12
    m_payoffsQALY["AdverseEffects", "QALYs"] <-  qaly_AdverseEffects/12
    m_payoffsQALY["Cured", "QALYs"] <- qaly_Cured/12
    m_payoffsQALY["Dead", "QALYs"] <- qaly_Dead/12
    
    #Check probability matrix by typing : m_payoffs
    payoff_trace <- state_membership %*% m_payoffs
    payoff_trace_d <- matrix(0, nrow = 240, ncol = 2, dimnames= list(state= 1:240, payoff= c("Costs", "QALYs")))
    for (i in 1:n_t) {payoff_trace_d[i, ] <- payoff_trace[i, ] * (1/((1+discount_m)^(i)))}
    
    payoff_traceCost <- state_membership %*% m_payoffsCost
    payoff_traceCost_d <- matrix(0, nrow = 240, ncol = 1, dimnames= list(state= 1:240, payoff= c("Costs")))
    for (i in 1:n_t) {payoff_traceCost_d[i, ] <- payoff_traceCost[i, ] * (1/((1+discount_m)^(i)))}
    
    payoff_traceQaly <- state_membership %*% m_payoffsQALY
    payoff_traceQaly_d <- matrix(0, nrow = 240, ncol = 1, dimnames= list(state= 1:240, payoff= c("QALYs")))
    for (i in 1:n_t) {payoff_traceQaly_d[i, ] <- payoff_traceQaly[i, ] * (1/((1+discount_m)^(i)))}
    
    #payoff_traceSUCCt <- state_membership$TreatmentCompleted %*% m_payoffs
    #payoff_traceSUCCt_d <- matrix(0, nrow = 240, ncol = 2, dimnames= list(state= 1:240, payoff= c("Costs","QALYs")))
    #for (i in 1:n_t) {payoff_traceSUCCt_d[i, ] <- payoff_traceSUCCt[i, ] * (1/((1+discount_m)^(i)))}
    #check results by typing: payoff_trace, payoff_traceCost, payoff_traceQaly
    #matplot(1:n_t, state_membership, type= 'l')
    #legend("topright", colnames(state_membership),col=seq_len(8),cex=0.6,fill=seq_len(8), lwd=1)
    
    #Outputs, discount rate applied.
    list(
      summary_resultsAll = colSums(payoff_trace_d) / n_c,
      summary_resultsperStateCost = colSums(payoff_traceCost_d) / n_c,
      summary_resultsperStateQaly = colSums(payoff_traceQaly_d) / n_c,
      summary_resultsSUCC_t = colSums(payoff_trace_d) / state_membership[240, 7],
      payoff_trace_perCycle = payoff_trace_d,
      payoff_trace_perCycle_cost = payoff_traceCost_d,
      payoff_trace_perCycle_qaly = payoff_traceQaly_d,
      individual_shifts = view(state_membership)
    )
  })
}
#COUNTRY-SPECIFIC parameters
#[[[[[ Notes: There were data available only for 4 countries]]]

# Parameters Kazakhstan--------------------------------------------------------------
#PARAMETERS LIST BELOW (change costs accordingly, QALYS are the same for everycountry)
n_c<- 3755
paramsmSTR_kaz <- list(
  #Costs & QALYS
  Cost_mSTR = 635.14,
  Cost_SLtreatment = 658.4,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  635.14*(1-0.26),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0338*(1-0.0412)),
  p_mSTR_LostFollowUp = (1/9)*(0.0068*(1-0.0412)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.5203+0.4189)*(1-0.0412)),
  p_mSTR_AdverseEffects = (1/9)*(0.0412),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.0203)*(1-0.0412)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.0840+0.0336),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.7267,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.1252,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.28/12))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.625),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.375), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

paramsCON_kaz <- list(
  #Costs & QALYS
  Cost_mSTR = 658.40,
  Cost_SLtreatment = 658.4,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  658.4*(1-0.37),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/20)*(0.033*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/20)*(0.1035*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/20)*((0.7614)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/20)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/20)*((0.1018)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.0840+0.0336),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.7267,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.1252,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)



# Parameters Moldova--------------------------------------------------------------
paramsmSTR_mol <- list(
  #Costs & QALYS
  Cost_mSTR = 693.52,
  Cost_SLtreatment = 451.92,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  693.52*(1-0.16),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0388*(1-0.0991)),
  p_mSTR_LostFollowUp = (1/9)*(0.0583*(1-0.0991)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.0388+0.8544)*(1-0.0991)),
  p_mSTR_AdverseEffects = (1/9)*(0.0991),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.0097)*(1-0.0991)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/18)*(0.1892+0.1351),
  p_SLtreatment_TreatmentCompleted = (1/18)*0.5676,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/18)*0.1081,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.28/12))**0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.6154+0.1538),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.2308), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)
# now conventional treatment, mSTR labels below means Conventional.
paramsCON_mol <- list(
  #Costs & QALYS
  Cost_mSTR = 451.92,
  Cost_SLtreatment = 451.92,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  451.92*(1-0.2),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/18)*(0.0626*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/18)*(0.1091*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/18)*((0.6923)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/18)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/18)*((0.1360)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/18)*(0.1892+0.1351),
  p_SLtreatment_TreatmentCompleted = (1/18)*0.5676,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/18)*0.1081,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

# Parameters Georgia--------------------------------------------------------------
paramsmSTR_geo <- list(
  #Costs & QALYS
  Cost_mSTR = 786.92,
  Cost_SLtreatment = 643.49,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  643.49*(1-0.24),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0101*(1-0.0187)),
  p_mSTR_LostFollowUp = (1/9)*(0.1515*(1-0.0187)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.7980+0.0303)*(1-0.0187)),
  p_mSTR_AdverseEffects = (1/9)*(0.0187),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.0101)*(1-0.0187)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/17)*(0.0943+0.2642),
  p_SLtreatment_TreatmentCompleted = (1/17)*0.5660,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/17)*0.0755,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.28/12))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.5),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.5), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

paramsCON_geo <- list(
  #Costs & QALYS
  Cost_mSTR = 643.49,
  Cost_SLtreatment = 643.49,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  643.49*(1-0.32),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/17)*(0.0183*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/17)*(0.1376*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/17)*((0.7844)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/17)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/17)*((0.0596)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/17)*(0.0943+0.2642),
  p_SLtreatment_TreatmentCompleted = (1/17)*0.5660,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/17)*0.0755,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
)
# Parameters Belarus--------------------------------------------------------------
n_c<- 801
paramsmSTR_bel <- list(
  #Costs & QALYS
  Cost_mSTR = 844.07,
  Cost_SLtreatment = 933.53,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects = 844.07*(1-0.16),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0261*(1-0.1963)),
  p_mSTR_LostFollowUp = (1/9)*(0.0279*(1-0.1963)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.0112+0.8957)*(1-0.1963)),
  p_mSTR_AdverseEffects = (1/9)*(0.1963),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.0391)*(1-0.1963)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/16)*(0.0662+0.1324),
  p_SLtreatment_TreatmentCompleted = (1/16)*0.7178,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/16)*0.0836,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.28/12))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.6471),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.0074+0.0294),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.3162), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

paramsCON_bel <- list(
  #Costs & QALYS
  Cost_mSTR = 933.533,
  Cost_SLtreatment = 933.533,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  933.533*(1-0.32),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/16)*(0.0168*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/16)*(0.0712*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/16)*((0.8290)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/16)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/16)*((0.0829)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/16)*(0.0662+0.1324),
  p_SLtreatment_TreatmentCompleted = (1/16)*0.7178,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/16)*0.0836,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.109), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
)






#Rest countries with no data below#############

# Parameters Armenia---
paramsmSTR_arm <- list(
  
)

paramsCON_arm <- list(
  
)

#Parameters Azerbaijan ---
paramsmSTR_aze <- list(
  
)

paramsCON_aze <- list(
  
)
# Parameters Ukraine---
paramsmSTR_ukr <- list(
  
)

paramsCON_ukr <- list(
  
)
# Parameters Uzbekistan---
paramsmSTR_uzb <- list(
  
)

paramsCON_uzb <- list(
  
)
# Parameters Tajikistan---
paramsmSTR_taj <- list(
  
)

paramsCON_taj <- list(
  
)
# Parameters Turkmenistan---
paramsmSTR_tur <- list(
  
)

paramsCON_tur <- list(
  
)
# Parameters Azerbaijan---
paramsmSTR_aze <- list(
  
)

paramsCON_aze <- list(
  
)
# Parameters Kyrgyzstan---
paramsmSTR_kyr <- list(
  
)

paramsCON_kyr <- list(
  
)

################################################################################
#------------------------------------------------------------------------------#
#I.  DETERMINISTIC MODEL and FUNCTION                                          #
#------------------------------------------------------------------------------#
# Function 
#CREATING MODEL FUNCTION; we can change the "params" and run it again per country.
model <- function(.params) { 
  with(.params, {
    n_t <- 240
    n_s <- 8
    v_state_names <- c("mSTR","SLtreatment","LostFollowUp","Unresolved","TreatmentCompleted","AdverseEffects","Cured", "Dead")
    m_p <- matrix(0, nrow= 8, ncol=8, dimnames= list(from= v_state_names, to = v_state_names))
   
    #prob matrix   
    m_p["mSTR", "mSTR"] <- 1  -p_mSTR_AdverseEffects -p_mSTR_SLtreatment -p_mSTR_LostFollowUp - p_mSTR_TreatmentCompleted -p_mSTR_Dead
    m_p["mSTR", "SLtreatment"] <- p_mSTR_SLtreatment
    m_p["mSTR", "LostFollowUp"] <- p_mSTR_LostFollowUp
    m_p["mSTR", "Unresolved"] <- p_mSTR_Unresolved
    m_p["mSTR", "TreatmentCompleted"] <- p_mSTR_TreatmentCompleted 
    m_p["mSTR", "AdverseEffects"] <- p_mSTR_AdverseEffects
    m_p["mSTR", "Cured"] <- p_mSTR_Cured
    m_p["mSTR", "Dead"] <-  p_mSTR_Dead
    
    m_p["SLtreatment", "mSTR"] <- p_SLtreatment_mSTR
    m_p["SLtreatment", "SLtreatment"] <- 1 -p_SLtreatment_Unresolved -p_SLtreatment_TreatmentCompleted -( p_SLtreatment_Dead )
    m_p["SLtreatment", "LostFollowUp"] <- p_SLtreatment_LostFollowUp
    m_p["SLtreatment", "Unresolved"] <- p_SLtreatment_Unresolved
    m_p["SLtreatment", "TreatmentCompleted"] <- p_SLtreatment_TreatmentCompleted
    m_p["SLtreatment", "AdverseEffects"] <- p_SLtreatment_AdverseEffects
    m_p["SLtreatment", "Cured"] <- p_SLtreatment_Cured
    m_p["SLtreatment", "Dead"] <-  p_SLtreatment_Dead 
    
    m_p["LostFollowUp", "mSTR"] <- p_LostFollowUp_mSTR
    m_p["LostFollowUp", "SLtreatment"] <- 0
    m_p["LostFollowUp", "LostFollowUp"] <- 1 - p_LostFollowUp_mSTR -( p_LostFollowUp_Dead)
    m_p["LostFollowUp", "Unresolved"] <- 0
    m_p["LostFollowUp", "TreatmentCompleted"] <-0
    m_p["LostFollowUp", "AdverseEffects"] <-  0
    m_p["LostFollowUp", "Cured"] <- 0
    m_p["LostFollowUp", "Dead"] <-  p_LostFollowUp_Dead
    
    m_p["Unresolved", "mSTR"] <- 0
    m_p["Unresolved", "SLtreatment"] <- 0
    m_p["Unresolved", "LostFollowUp"] <- 0  
    m_p["Unresolved", "Unresolved"] <- 1 - (p_Unresolved_Dead)
    m_p["Unresolved", "TreatmentCompleted"] <- 0
    m_p["Unresolved", "AdverseEffects"] <- 0
    m_p["Unresolved", "Cured"] <- 0
    m_p["Unresolved", "Dead"] <-  (0.75)*p_Unresolved_Dead
    
    m_p["TreatmentCompleted", "mSTR"] <- 0
    m_p["TreatmentCompleted", "SLtreatment"] <- p_TreatmentCompleted_SLtreatment
    m_p["TreatmentCompleted", "LostFollowUp"] <- 0
    m_p["TreatmentCompleted", "Unresolved"] <-  0
    m_p["TreatmentCompleted", "TreatmentCompleted"] <-  1- p_TreatmentCompleted_SLtreatment -p_TreatmentCompleted_Cured -(p_TreatmentCompleted_Dead)
    m_p["TreatmentCompleted", "AdverseEffects"] <-  0
    m_p["TreatmentCompleted", "Cured"] <-  p_TreatmentCompleted_Cured
    m_p["TreatmentCompleted", "Dead"] <-   p_TreatmentCompleted_Dead
    
    m_p["AdverseEffects", "mSTR"] <- 0
    m_p["AdverseEffects", "SLtreatment"] <- p_AdverseEffects_SLtreatment
    m_p["AdverseEffects", "LostFollowUp"] <- 0
    m_p["AdverseEffects", "Unresolved"] <-  p_AdverseEffects_Unresolved
    m_p["AdverseEffects", "TreatmentCompleted"] <-  0
    m_p["AdverseEffects", "AdverseEffects"] <-  1- p_AdverseEffects_Unresolved - p_AdverseEffects_SLtreatment -(p_AdverseEffects_Dead)
    m_p["AdverseEffects", "Cured"] <-  0
    m_p["AdverseEffects", "Dead"] <-   p_AdverseEffects_Dead
    
    m_p["Cured", "mSTR"] <- 0
    m_p["Cured", "SLtreatment"] <- 0
    m_p["Cured", "LostFollowUp"] <- 0
    m_p["Cured", "Unresolved"] <- 0
    m_p["Cured", "TreatmentCompleted"] <- 0
    m_p["Cured", "AdverseEffects"] <- 0
    m_p["Cured", "Cured"] <-  1 - (p_Cured_Dead)
    m_p["Cured", "Dead"] <-   p_Cured_Dead
    
    m_p["Dead", "mSTR"] <- 0
    m_p["Dead", "SLtreatment"] <- 0
    m_p["Dead", "LostFollowUp"] <- 0
    m_p["Dead", "Unresolved"] <- 0
    m_p["Dead", "TreatmentCompleted"] <- 0
    m_p["Dead", "AdverseEffects"] <- 0
    m_p["Dead", "Cured"] <- 0
    m_p["Dead", "Dead"] <-  1
    
    #Start drafting the markov matrix
    state_membership <- array(NA_real_, dim= c(n_t, n_s), dimnames= list(cycle = 1:n_t, state = v_state_names))
    #Start filling in first row; only with treatment box.
    state_membership[1,]= c(n_c, 0, 0, 0, 0, 0, 0, 0)
    for (i in 2:n_t) {state_membership[i, ] <- state_membership[i-1, ] %*% m_p
    }
    #Insert costs per healthy, diseased and death, let's wait for Tom
    m_payoffs <- matrix(0, nrow = 8, ncol = 2, dimnames= list(state= v_state_names, 
                payoff= c("Costs", "QALYs")))
    
    m_payoffsCost<- matrix(0, nrow=8, ncol=1,
                        dimnames= list(state= v_state_names, 
                        payoff= c("Costs")))
    m_payoffsQALY <- matrix(0, nrow=8, ncol=1,
                        dimnames= list(state= v_state_names, 
                        payoff= c("QALYs")))

    m_payoffs["mSTR",  "Costs"] <-  Cost_mSTR
    m_payoffs["SLtreatment", "Costs"] <-  Cost_SLtreatment
    m_payoffs["LostFollowUp", "Costs"] <- Cost_LostFollowUp
    m_payoffs["Unresolved", "Costs"] <- Cost_Unresolved
    m_payoffs["TreatmentCompleted", "Costs"] <- Cost_TreatmentCompleted
    m_payoffs["AdverseEffects", "Costs"] <-  Cost_AdverseEffects
    m_payoffs["Cured", "Costs"] <- Cost_Cured
    m_payoffs["Dead", "Costs"] <- Cost_Dead
    
    m_payoffsCost["mSTR",  "Costs"] <-  Cost_mSTR
    m_payoffsCost["SLtreatment", "Costs"] <-  Cost_SLtreatment
    m_payoffsCost["LostFollowUp", "Costs"] <- Cost_LostFollowUp
    m_payoffsCost["Unresolved", "Costs"] <- Cost_Unresolved
    m_payoffsCost["TreatmentCompleted", "Costs"] <- Cost_TreatmentCompleted
    m_payoffsCost["AdverseEffects", "Costs"] <-  Cost_AdverseEffects
    m_payoffsCost["Cured", "Costs"] <- Cost_Cured
    m_payoffsCost["Dead", "Costs"] <- Cost_Dead
    
    m_payoffs["mSTR",  "QALYs"] <-  qaly_mSTR/12
    m_payoffs["SLtreatment", "QALYs"] <-  qaly_SLtreatment/12
    m_payoffs["LostFollowUp", "QALYs"] <- qaly_LostFollowUp/12
    m_payoffs["Unresolved", "QALYs"] <- qaly_Unresolved/12
    m_payoffs["TreatmentCompleted", "QALYs"] <- qaly_TreatmentCompleted/12
    m_payoffs["AdverseEffects", "QALYs"] <-  qaly_AdverseEffects/12
    m_payoffs["Cured", "QALYs"] <- qaly_Cured/12
    m_payoffs["Dead", "QALYs"] <- qaly_Dead/12
    
    m_payoffsQALY["mSTR",  "QALYs"] <-  qaly_mSTR/12
    m_payoffsQALY["SLtreatment", "QALYs"] <-  qaly_SLtreatment/12
    m_payoffsQALY["LostFollowUp", "QALYs"] <- qaly_LostFollowUp/12
    m_payoffsQALY["Unresolved", "QALYs"] <- qaly_Unresolved/12
    m_payoffsQALY["TreatmentCompleted", "QALYs"] <- qaly_TreatmentCompleted/12
    m_payoffsQALY["AdverseEffects", "QALYs"] <-  qaly_AdverseEffects/12
    m_payoffsQALY["Cured", "QALYs"] <- qaly_Cured/12
    m_payoffsQALY["Dead", "QALYs"] <- qaly_Dead/12
    
    #Check probability matrix by typing : m_payoffs
    payoff_trace <- state_membership %*% m_payoffs
    payoff_trace_d <- matrix(0, nrow = 240, ncol = 2, dimnames= list(state= 1:240, payoff= c("Costs", "QALYs")))
    for (i in 1:n_t) {payoff_trace_d[i, ] <- payoff_trace[i, ] * (1/((1+discount_m)^(i)))}
    
    payoff_traceCost <- state_membership %*% m_payoffsCost
    payoff_traceCost_d <- matrix(0, nrow = 240, ncol = 1, dimnames= list(state= 1:240, payoff= c("Costs")))
    for (i in 1:n_t) {payoff_traceCost_d[i, ] <- payoff_traceCost[i, ] * (1/((1+discount_m)^(i)))}
    
    payoff_traceQaly <- state_membership %*% m_payoffsQALY
    payoff_traceQaly_d <- matrix(0, nrow = 240, ncol = 1, dimnames= list(state= 1:240, payoff= c("QALYs")))
    for (i in 1:n_t) {payoff_traceQaly_d[i, ] <- payoff_traceQaly[i, ] * (1/((1+discount_m)^(i)))}
    
    #payoff_traceSUCCt <- state_membership$TreatmentCompleted %*% m_payoffs
    #payoff_traceSUCCt_d <- matrix(0, nrow = 240, ncol = 2, dimnames= list(state= 1:240, payoff= c("Costs","QALYs")))
    #for (i in 1:n_t) {payoff_traceSUCCt_d[i, ] <- payoff_traceSUCCt[i, ] * (1/((1+discount_m)^(i)))}
    #check results by typing: payoff_trace, payoff_traceCost, payoff_traceQaly
    #matplot(1:n_t, state_membership, type= 'l')
    #legend("topright", colnames(state_membership),col=seq_len(8),cex=0.6,fill=seq_len(8), lwd=1)
    
    #Outputs, discount rate applied.
    list(
      summary_resultsAll = colSums(payoff_trace_d) / n_c,
      summary_resultsperStateCost = colSums(payoff_traceCost_d) / n_c,
      summary_resultsperStateQaly = colSums(payoff_traceQaly_d) / n_c,
      summary_resultsSUCC_t = colSums(payoff_trace_d) / state_membership[240, 7],
      payoff_trace_perCycle = payoff_trace_d,
      payoff_trace_perCycle_cost = payoff_traceCost_d,
      payoff_trace_perCycle_qaly = payoff_traceQaly_d,
      individual_shifts = view(state_membership)
    )
  })
}

# Exploratory analyses by country: Graphs and diagrams.
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/Bloomsbury_PolicyLab/0_Article/Figures")

#------------------------------------####
#PLOTING DISTRIBUTIONS per country  #
#------------------------------------#

#Kazakhstan###################
#Calling outputs
n_c<- 3755
resultsKAZ_mSTR <- model(paramsmSTR_kaz)
matrix1 <- resultsKAZ_mSTR$individual_shifts
resultsKAZ_CON <- model(paramsCON_kaz)
matrix2 <- resultsKAZ_CON$individual_shifts
#Plot Dynamics
theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="Kazahstan_dynamics.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(1:n_t, matrix1$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 4000), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250), main="mSTR treatment, Kazakhstan")
lines(1:n_t,matrix1$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix1$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix1$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix1$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix1$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix1$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix1$Dead, lwd=3, col="sienna1", type="l")
legend(x = 80, y = 3600, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), pt.cex=0.6, lwd=c(3,3),cex=0.75, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
plot(1:n_t, matrix2$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 4000), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250),  main="Conventional treatment, Kazakhstan")
lines(1:n_t,matrix2$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix2$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix2$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix2$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix2$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix2$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix2$Dead, lwd=3, col="sienna1", type="l")
#legend(x = 80, y = 3800, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), lwd=c(3,3), cex=0.5,  pch = 1, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
dev.off()

#Moldova##########
#Calling outputs
n_c<- 593
resultsmol_mSTR <- model(paramsmSTR_mol)
matrix1 <- resultsmol_mSTR$individual_shifts
resultsmol_CON <- model(paramsCON_mol)
matrix2 <- resultsmol_CON$individual_shifts
#Plot Dynamics
theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="Moldova_dynamics.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(1:n_t, matrix1$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 600), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250), main="mSTR treatment, Moldova")
lines(1:n_t,matrix1$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix1$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix1$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix1$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix1$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix1$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix1$Dead, lwd=3, col="sienna1", type="l")
legend(x = 80, y = 500, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), pt.cex=0.6, lwd=c(3,3),cex=0.75, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
plot(1:n_t, matrix2$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 600), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250),  main="Conventional treatment, Moldova")
lines(1:n_t,matrix2$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix2$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix2$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix2$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix2$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix2$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix2$Dead, lwd=3, col="sienna1", type="l")
#legend(x = 80, y = 3800, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), lwd=c(3,3), cex=0.5,  pch = 1, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
dev.off()


#Belarus####
#Calling outputs
n_c<- 801
resultsbel_mSTR <- model(paramsmSTR_bel)
matrix1 <- resultsbel_mSTR$individual_shifts
resultsbel_CON <- model(paramsCON_bel)
matrix2 <- resultsbel_CON$individual_shifts
#Plot Dynamics
theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="Belarus_dynamics.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(1:n_t, matrix1$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 800), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250), main="mSTR treatment, Belarus")
lines(1:n_t,matrix1$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix1$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix1$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix1$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix1$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix1$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix1$Dead, lwd=3, col="sienna1", type="l")
legend(x = 80, y = 600, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), pt.cex=0.6, lwd=c(3,3),cex=0.75, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
plot(1:n_t, matrix2$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 800), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250),  main="Conventional treatment, Belarus")
lines(1:n_t,matrix2$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix2$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix2$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix2$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix2$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix2$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix2$Dead, lwd=3, col="sienna1", type="l")
#legend(x = 80, y = 3800, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), lwd=c(3,3), cex=0.5,  pch = 1, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
dev.off()


#Georgia####
n_c<- 187
resultsGEO_mSTR <- model(paramsmSTR_geo)
matrix1 <- resultsGEO_mSTR$individual_shifts
resultsGEO_CON <- model(paramsCON_geo)
matrix2 <- resultsGEO_CON$individual_shifts
#Plot Dynamics
theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="Georgia_dynamics.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(1:n_t, matrix1$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 200), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250), main="mSTR treatment, Georgia")
lines(1:n_t,matrix1$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix1$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix1$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix1$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix1$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix1$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix1$Dead, lwd=3, col="sienna1", type="l")
legend(x = 80, y = 140, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), pt.cex=0.6, lwd=c(3,3),cex=0.75, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
plot(1:n_t, matrix2$mSTR, type="l", lwd=3, col="blue", ylim=c(0, 200), xaxs="i", yaxs="i", xlab= "Time in months", ylab="Incidence of MDR/RR TB (n)", las = 1, xlim=c(0, 250),  main="Conventional treatment, Georgia")
lines(1:n_t,matrix2$SLtreatment, lwd=3, col="red", type="l")
lines(1:n_t,matrix2$LostFollowUp, lwd=3, col="navy", type="l")
lines(1:n_t,matrix2$Unresolved, lwd=3, col="springgreen4", type="l")
lines(1:n_t,matrix2$TreatmentCompleted, lwd=3, col="darkred", type="l")
lines(1:n_t,matrix2$AdverseEffects, lwd=3, col="gold", type="l")
lines(1:n_t,matrix2$Cured, lwd=3, col="orchid3", type="l")
lines(1:n_t,matrix2$Dead, lwd=3, col="sienna1", type="l")
#legend(x = 80, y = 3800, legend=c("mSTR","SLtreatment", "LTFU", "Unresolved", "Treatment completed", "AE", "Cured", "Dead"), lwd=c(3,3), cex=0.5,  pch = 1, bty = "n", col=c("blue","red", "navy", "springgreen4", "darkred", "gold", "orchid3", "sienna1"))
dev.off()





################################################################################
#RESULTS GENERAL, deterministic, by country ##################################
#Results GENERAL########################:
#KAZAHSTAN!
n_c<- 3755
resultsKAZ_mSTR <- model(paramsmSTR_kaz)
resultsKAZ_CON <- model(paramsCON_kaz)
ICER_kaz <- (resultsKAZ_mSTR$summary_resultsperStateCost - resultsKAZ_CON$summary_resultsperStateCost)/( resultsKAZ_mSTR$summary_resultsperStateQaly -resultsKAZ_CON$summary_resultsperStateQaly)
resultsKAZ_mSTR$summary_resultsperStateCost
resultsKAZ_CON$summary_resultsperStateCost
resultsKAZ_mSTR$summary_resultsperStateQaly
resultsKAZ_CON$summary_resultsperStateQaly
excess_mort <- resultsKAZ_CON$individual_shifts$Dead[240]-resultsKAZ_mSTR$individual_shifts$Dead[240]
resultsKAZ_mSTR$summary_resultsSUCC_t
resultsKAZ_CON$summary_resultsSUCC_t
resultsKAZ_mSTR$individual_shifts$Dead[240]/n_c
resultsKAZ_CON$individual_shifts$Dead[240]/n_c
resultsKAZ_mSTR$individual_shifts$Cured[240]/n_c
resultsKAZ_CON$individual_shifts$Cured[240]/n_c


#Moldova!
n_c<- 593
resultsmol_mSTR <- model(paramsmSTR_mol)
resultsmol_CON <- model(paramsCON_mol)
ICER_mol <- (resultsmol_mSTR$summary_resultsperStateCost - resultsmol_CON$summary_resultsperStateCost)/( resultsmol_mSTR$summary_resultsperStateQaly -resultsmol_CON$summary_resultsperStateQaly)
resultsmol_mSTR$summary_resultsperStateCost
resultsmol_CON$summary_resultsperStateCost
resultsmol_mSTR$summary_resultsperStateQaly
resultsmol_CON$summary_resultsperStateQaly
excess_mort <- resultsmol_CON$individual_shifts$Dead[240]-resultsmol_mSTR$individual_shifts$Dead[240]
resultsmol_mSTR$summary_resultsSUCC_t
resultsmol_CON$summary_resultsSUCC_t
resultsmol_mSTR$individual_shifts$Dead[240]/n_c
resultsmol_CON$individual_shifts$Dead[240]/n_c
resultsmol_mSTR$individual_shifts$Cured[240]/n_c
resultsmol_CON$individual_shifts$Cured[240]/n_c

#Belarus
n_c<- 801
resultsbel_mSTR <- model(paramsmSTR_bel)
resultsbel_CON <- model(paramsCON_bel)
ICER_bel<- (resultsbel_mSTR$summary_resultsperStateCost - resultsbel_CON$summary_resultsperStateCost)/( resultsbel_mSTR$summary_resultsperStateQaly -resultsbel_CON$summary_resultsperStateQaly)
resultsbel_mSTR$summary_resultsperStateCost
resultsbel_CON$summary_resultsperStateCost
resultsbel_mSTR$summary_resultsperStateQaly
resultsbel_CON$summary_resultsperStateQaly
excess_mort <- resultsbel_CON$individual_shifts$Dead[240]-resultsbel_mSTR$individual_shifts$Dead[240]
resultsbel_mSTR$summary_resultsSUCC_t
resultsbel_CON$summary_resultsSUCC_t
resultsbel_mSTR$individual_shifts$Dead[240]/n_c
resultsbel_CON$individual_shifts$Dead[240]/n_c
resultsbel_mSTR$individual_shifts$Cured[240]/n_c
resultsbel_CON$individual_shifts$Cured[240]/n_c

#Georgia
n_c<- 187
resultsgeo_mSTR <- model(paramsmSTR_geo)
resultsgeo_CON <- model(paramsCON_geo)
ICER_geo<- (resultsgeo_mSTR$summary_resultsperStateCost - resultsgeo_CON$summary_resultsperStateCost)/( resultsgeo_mSTR$summary_resultsperStateQaly -resultsgeo_CON$summary_resultsperStateQaly)
resultsgeo_mSTR$summary_resultsperStateCost
resultsgeo_CON$summary_resultsperStateCost
resultsgeo_mSTR$summary_resultsperStateQaly
resultsgeo_CON$summary_resultsperStateQaly
excess_mort <- resultsgeo_CON$individual_shifts$Dead[240]-resultsgeo_mSTR$individual_shifts$Dead[240]
resultsgeo_mSTR$summary_resultsSUCC_t
resultsgeo_CON$summary_resultsSUCC_t

resultsgeo_mSTR$individual_shifts$Dead[240]/n_c
resultsgeo_CON$individual_shifts$Dead[240]/n_c
resultsgeo_mSTR$individual_shifts$Cured[240]/n_c
resultsgeo_CON$individual_shifts$Cured[240]/n_c

#####Rest of the Countries_with no Data######
#Armenia!
n_c<- 64
resultsarm_mSTR <- model(paramsmSTR_arm)
resultsarm_CON <- model(paramsCON_arm)
ICER_arm <- (resultsarm_mSTR$summary_resultsperStateCost - resultsarm_CON$summary_resultsperStateCost)/( resultsarm_mSTR$summary_resultsperStateQaly -resultsarm_CON$summary_resultsperStateQaly)
resultsarm_mSTR$summary_resultsperStateCost
resultsarm_CON$summary_resultsperStateCost
resultsarm_mSTR$summary_resultsperStateQaly
resultsarm_CON$summary_resultsperStateQaly
excess_mort <- resultsarm_CON$individual_shifts$Dead[240]-resultsarm_mSTR$individual_shifts$Dead[240]
resultsarm_mSTR$summary_resultsSUCC_t
resultsarm_CON$summary_resultsSUCC_t
resultsarm_mSTR$individual_shifts$Dead[240]/n_c
resultsarm_CON$individual_shifts$Dead[240]/n_c
resultsarm_mSTR$individual_shifts$Cured[240]/n_c
resultsarm_CON$individual_shifts$Cured[240]/n_c

#Azerbaijan!
n_c<- 1040
resultsaze_mSTR <- model(paramsmSTR_aze)
resultsaze_CON <- model(paramsCON_aze)
ICER_aze<- (resultsaze_mSTR$summary_resultsperStateCost - resultsaze_CON$summary_resultsperStateCost)/( resultsaze_mSTR$summary_resultsperStateQaly -resultsaze_CON$summary_resultsperStateQaly)
resultsaze_mSTR$summary_resultsperStateCost
resultsaze_CON$summary_resultsperStateCost
resultsaze_mSTR$summary_resultsperStateQaly
resultsaze_CON$summary_resultsperStateQaly
excess_mort <- resultsaze_CON$individual_shifts$Dead[240]-resultsaze_mSTR$individual_shifts$Dead[240]
resultsaze_mSTR$summary_resultsSUCC_t
resultsaze_CON$summary_resultsSUCC_t
resultsaze_mSTR$individual_shifts$Dead[240]/n_c
resultsaze_CON$individual_shifts$Dead[240]/n_c
resultsaze_mSTR$individual_shifts$Cured[240]/n_c
resultsaze_CON$individual_shifts$Cured[240]/n_c



#Ukraine
n_c<- 4046
resultsukr_mSTR <- model(paramsmSTR_ukr)
resultsukr_CON <- model(paramsCON_ukr)
ICER_ukr<- (resultsukr_mSTR$summary_resultsperStateCost - resultsukr_CON$summary_resultsperStateCost)/( resultsukr_mSTR$summary_resultsperStateQaly -resultsukr_CON$summary_resultsperStateQaly)
resultsukr_mSTR$summary_resultsperStateCost
resultsukr_CON$summary_resultsperStateCost
resultsukr_mSTR$summary_resultsperStateQaly
resultsukr_CON$summary_resultsperStateQaly
excess_mort <- resultsukr_CON$individual_shifts$Dead[240]-resultsukr_mSTR$individual_shifts$Dead[240]
resultsukr_mSTR$summary_resultsSUCC_t
resultsukr_CON$summary_resultsSUCC_t

resultsukr_mSTR$individual_shifts$Dead[240]/n_c
resultsukr_CON$individual_shifts$Dead[240]/n_c
resultsukr_mSTR$individual_shifts$Cured[240]/n_c
resultsukr_CON$individual_shifts$Cured[240]/n_c

#Uzbekistan
n_c<- 2147
resultsuzb_mSTR <- model(paramsmSTR_uzb)
resultsuzb_CON <- model(paramsCON_uzb)
ICER_uzb<- (resultsuzb_mSTR$summary_resultsperStateCost - resultsuzb_CON$summary_resultsperStateCost)/( resultsuzb_mSTR$summary_resultsperStateQaly -resultsuzb_CON$summary_resultsperStateQaly)
resultsuzb_mSTR$summary_resultsperStateCost
resultsuzb_CON$summary_resultsperStateCost
resultsuzb_mSTR$summary_resultsperStateQaly
resultsuzb_CON$summary_resultsperStateQaly
excess_mort <- resultsuzb_CON$individual_shifts$Dead[240]-resultsuzb_mSTR$individual_shifts$Dead[240]
resultsuzb_mSTR$summary_resultsSUCC_t
resultsuzb_CON$summary_resultsSUCC_t
resultsuzb_mSTR$individual_shifts$Dead[240]/n_c
resultsuzb_CON$individual_shifts$Dead[240]/n_c
resultsuzb_mSTR$individual_shifts$Cured[240]/n_c
resultsuzb_CON$individual_shifts$Cured[240]/n_c

#Tajikistan
n_c<- 606
resultstaj_mSTR <- model(paramsmSTR_taj)
resultstaj_CON <- model(paramsCON_taj)
ICER_taj<- (resultstaj_mSTR$summary_resultsperStateCost - resultstaj_CON$summary_resultsperStateCost)/( resultstaj_mSTR$summary_resultsperStateQaly -resultstaj_CON$summary_resultsperStateQaly)
resultstaj_mSTR$summary_resultsperStateCost
resultstaj_CON$summary_resultsperStateCost
resultstaj_mSTR$summary_resultsperStateQaly
resultstaj_CON$summary_resultsperStateQaly
excess_mort <- resultstaj_CON$individual_shifts$Dead[240]-resultstaj_mSTR$individual_shifts$Dead[240]
resultstaj_mSTR$summary_resultsSUCC_t
resultstaj_CON$summary_resultsSUCC_t
resultstaj_mSTR$individual_shifts$Dead[240]/n_c
resultstaj_CON$individual_shifts$Dead[240]/n_c
resultstaj_mSTR$individual_shifts$Cured[240]/n_c
resultstaj_CON$individual_shifts$Cured[240]/n_c

#Turkmenistan
n_c<- 808
resultstur_mSTR <- model(paramsmSTR_tur)
resultstur_CON <- model(paramsCON_tur)
ICER_tur<- (resultstur_mSTR$summary_resultsperStateCost - resultstur_CON$summary_resultsperStateCost)/( resultstur_mSTR$summary_resultsperStateQaly -resultstur_CON$summary_resultsperStateQaly)
resultstur_mSTR$summary_resultsperStateCost
resultstur_CON$summary_resultsperStateCost
resultstur_mSTR$summary_resultsperStateQaly
resultstur_CON$summary_resultsperStateQaly
excess_mort <- resultstur_CON$individual_shifts$Dead[240]-resultstur_mSTR$individual_shifts$Dead[240]
resultstur_mSTR$summary_resultsSUCC_t
resultstur_CON$summary_resultsSUCC_t

resultstur_mSTR$individual_shifts$Dead[240]/n_c
resultstur_CON$individual_shifts$Dead[240]/n_c
resultstur_mSTR$individual_shifts$Cured[240]/n_c
resultstur_CON$individual_shifts$Cured[240]/n_c



#Kyrgistan
n_c<- 918
resultskyr_mSTR <- model(paramsmSTR_kyr)
resultskyr_CON <- model(paramsCON_kyr)
ICER_kyr<- (resultskyr_mSTR$summary_resultsperStateCost - resultskyr_CON$summary_resultsperStateCost)/( resultskyr_mSTR$summary_resultsperStateQaly -resultskyr_CON$summary_resultsperStateQaly)
resultskyr_mSTR$summary_resultsperStateCost
resultskyr_CON$summary_resultsperStateCost
resultskyr_mSTR$summary_resultsperStateQaly
resultskyr_CON$summary_resultsperStateQaly
excess_mort <- resultskyr_CON$individual_shifts$Dead[240]-resultskyr_mSTR$individual_shifts$Dead[240]
resultskyr_mSTR$summary_resultsSUCC_t
resultskyr_CON$summary_resultsSUCC_t

resultskyr_mSTR$individual_shifts$Dead[240]/n_c
resultskyr_CON$individual_shifts$Dead[240]/n_c
resultskyr_mSTR$individual_shifts$Cured[240]/n_c
resultskyr_CON$individual_shifts$Cured[240]/n_c




################################################################################

#------------------------------------------------------------------------------#
#II UNIVARIATE SENSITIVITY ANLYSES###############################################
#------------------------------------------------------------------------------#
#Costs of main treatments separately, Costs of SL treatment .25, .5, discount rates 0 and 6%, return to care from LTFU.
#Kazakhstan###################

#(1) Treatment costs
n_c<- 3755
  paramsmSTR_kaz$Cost_mSTR<- (paramsmSTR_kaz$Cost_mSTR*0.26)+((1-0.26)*paramsmSTR_kaz$Cost_mSTR*0.75)
  resultsKAZ_mSTR_LB <- model(paramsmSTR_kaz)
  resultsKAZ_CON <- model(paramsCON_kaz)
  qaly_kaz_mSTR_LB <- (resultsKAZ_mSTR_LB$summary_resultsperStateQaly -resultsKAZ_CON$summary_resultsperStateQaly)
  cost_kaz_mSTR_LB <- (resultsKAZ_mSTR_LB$summary_resultsperStateCost - resultsKAZ_CON$summary_resultsperStateCost)
  paramsmSTR_kaz$Cost_mSTR<-  635.14
  
  paramsCON_kaz$Cost_mSTR<- (0.37*paramsCON_kaz$Cost_mSTR)+((1-0.37)*paramsCON_kaz$Cost_mSTR*0.75)
  resultsKAZ_CON_LB <- model(paramsCON_kaz)
  resultsKAZ_mSTR <- model(paramsmSTR_kaz)
  qaly_kaz_CON_LB <- ( resultsKAZ_mSTR$summary_resultsperStateQaly -resultsKAZ_CON_LB$summary_resultsperStateQaly)
  cost_kaz_CON_LB <- (resultsKAZ_mSTR$summary_resultsperStateCost - resultsKAZ_CON_LB$summary_resultsperStateCost)
  paramsCON_kaz$Cost_mSTR<- 658.40
  
  paramsmSTR_kaz$Cost_mSTR<-  (paramsmSTR_kaz$Cost_mSTR*0.26)+((1-0.26)*paramsmSTR_kaz$Cost_mSTR*1.25)
  resultsKAZ_mSTR_UB <- model(paramsmSTR_kaz)
  resultsKAZ_CON <- model(paramsCON_kaz)
  qaly_kaz_mSTR_UB <- (resultsKAZ_mSTR_UB$summary_resultsperStateQaly -resultsKAZ_CON$summary_resultsperStateQaly)
  cost_kaz_mSTR_UB <- (resultsKAZ_mSTR_UB$summary_resultsperStateCost - resultsKAZ_CON$summary_resultsperStateCost)
  paramsmSTR_kaz$Cost_mSTR<-  635.14
  
  paramsCON_kaz$Cost_mSTR<- (0.37*paramsCON_kaz$Cost_mSTR)+((1-0.37)*paramsCON_kaz$Cost_mSTR*1.25)
  resultsKAZ_CON_UB <- model(paramsCON_kaz)
  resultsKAZ_mSTR <- model(paramsmSTR_kaz)
  qaly_kaz_CON_UB <- (resultsKAZ_mSTR$summary_resultsperStateQaly -resultsKAZ_CON_UB$summary_resultsperStateQaly)
  cost_kaz_CON_UB <- (resultsKAZ_mSTR$summary_resultsperStateCost - resultsKAZ_CON_UB$summary_resultsperStateCost)
  paramsCON_kaz$Cost_mSTR<- 658.40

qaly_kaz_CON_LB
qaly_kaz_CON_UB
qaly_kaz_mSTR_LB
qaly_kaz_mSTR_UB
cost_kaz_CON_LB
cost_kaz_CON_UB
cost_kaz_mSTR_LB
cost_kaz_mSTR_UB

#(2) Costs for SL treatment
paramsmSTR_kaz$Cost_SLtreatment<- paramsmSTR_kaz$Cost_SLtreatment*0.75
paramsCON_kaz$Cost_SLtreatment<- paramsCON_kaz$Cost_SLtreatment*0.75
resultsKAZ_mSTR_SL075 <- model(paramsmSTR_kaz)
resultsKAZ_CON_SL075 <- model(paramsCON_kaz)
cost_kaz_mSTR_SL075 <- (resultsKAZ_mSTR_SL075$summary_resultsperStateCost - resultsKAZ_CON_SL075$summary_resultsperStateCost)
qaly_kaz_mSTR_SL075 <- (resultsKAZ_mSTR_SL075$summary_resultsperStateQaly -resultsKAZ_CON_SL075$summary_resultsperStateQaly)
paramsmSTR_kaz$Cost_SLtreatment<- paramsmSTR_kaz$Cost_SLtreatment/0.75
paramsCON_kaz$Cost_SLtreatment<- paramsCON_kaz$Cost_SLtreatment/0.75

paramsmSTR_kaz$Cost_SLtreatment<- paramsmSTR_kaz$Cost_SLtreatment*1.25
paramsCON_kaz$Cost_SLtreatment<- paramsCON_kaz$Cost_SLtreatment*1.25
resultsKAZ_mSTR_SL125 <- model(paramsmSTR_kaz)
resultsKAZ_CON_SL125 <- model(paramsCON_kaz)
cost_kaz_mSTR_SL125 <- (resultsKAZ_mSTR_SL125$summary_resultsperStateCost - resultsKAZ_CON_SL125$summary_resultsperStateCost)
qaly_kaz_mSTR_SL125 <- (resultsKAZ_mSTR_SL125$summary_resultsperStateQaly -resultsKAZ_CON_SL125$summary_resultsperStateQaly)
paramsmSTR_kaz$Cost_SLtreatment<- paramsmSTR_kaz$Cost_SLtreatment/1.25
paramsCON_kaz$Cost_SLtreatment<- paramsCON_kaz$Cost_SLtreatment/1.25

qaly_kaz_mSTR_SL075
cost_kaz_mSTR_SL075
qaly_kaz_mSTR_SL125
cost_kaz_mSTR_SL125


#(3) discount rate
discount <- 0.0 
discount_m <- discount/12
resultsKAZ_mSTR_d0 <- model(paramsmSTR_kaz)
resultsKAZ_CON_d0 <- model(paramsCON_kaz)
qaly_kaz_mSTR_d0 <- (resultsKAZ_mSTR_d0$summary_resultsperStateQaly -resultsKAZ_CON_d0$summary_resultsperStateQaly)
cost_kaz_mSTR_d0 <- (resultsKAZ_mSTR_d0$summary_resultsperStateCost - resultsKAZ_CON_d0$summary_resultsperStateCost)
discount <- 0.06 
discount_m <- discount/12
resultsKAZ_mSTR_d06 <- model(paramsmSTR_kaz)
resultsKAZ_CON_d06 <- model(paramsCON_kaz)
qaly_kaz_mSTR_d06 <- ( resultsKAZ_mSTR_d06$summary_resultsperStateQaly -resultsKAZ_CON_d06$summary_resultsperStateQaly)
cost_kaz_mSTR_d06 <- (resultsKAZ_mSTR_d06$summary_resultsperStateCost - resultsKAZ_CON_d06$summary_resultsperStateCost)
discount <-0.03
discount_m <- discount/12

qaly_kaz_mSTR_d0
cost_kaz_mSTR_d0
qaly_kaz_mSTR_d06
cost_kaz_mSTR_d06

#(4) Health state utility after cured TB or and LTFU/others
paramsmSTR_kaz$qaly_AdverseEffects <-paramsmSTR_kaz$qaly_AdverseEffects + 0.14
paramsCON_kaz$qaly_AdverseEffects <- paramsCON_kaz$qaly_AdverseEffects + 0.14
paramsmSTR_kaz$qaly_Unresolved <-paramsmSTR_kaz$qaly_Unresolved + 0.14
paramsCON_kaz$qaly_Unresolved <- paramsCON_kaz$qaly_Unresolved + 0.14
paramsmSTR_kaz$qaly_LostFollowUp <-paramsmSTR_kaz$qaly_LostFollowUp + 0.14
paramsCON_kaz$qaly_LostFollowUp <- paramsCON_kaz$qaly_LostFollowUp +  0.14
resultsKAZ_mSTR_Ud <- model(paramsmSTR_kaz)
resultsKAZ_CON_Ud <- model(paramsCON_kaz)
qaly_kaz_mSTR_Ud <- (resultsKAZ_mSTR_Ud$summary_resultsperStateQaly -resultsKAZ_CON_Ud$summary_resultsperStateQaly)
cost_kaz_mSTR_Ud <- (resultsKAZ_mSTR_Ud$summary_resultsperStateCost - resultsKAZ_CON_Ud$summary_resultsperStateCost)

paramsmSTR_kaz$qaly_AdverseEffects <-paramsmSTR_kaz$qaly_AdverseEffects - (0.28)
paramsCON_kaz$qaly_AdverseEffects <- paramsCON_kaz$qaly_AdverseEffects - (0.28)
paramsmSTR_kaz$qaly_Unresolved <-paramsmSTR_kaz$qaly_Unresolved - (0.28)
paramsCON_kaz$qaly_Unresolved <- paramsCON_kaz$qaly_Unresolved - (0.28)
paramsmSTR_kaz$qaly_LostFollowUp <-paramsmSTR_kaz$qaly_LostFollowUp - (0.28)
paramsCON_kaz$qaly_LostFollowUp <- paramsCON_kaz$qaly_LostFollowUp - (0.28)
resultsKAZ_mSTR_Dd <- model(paramsmSTR_kaz)
resultsKAZ_CON_Dd <- model(paramsCON_kaz)
qaly_kaz_mSTR_Dd <- ( resultsKAZ_mSTR_Dd$summary_resultsperStateQaly -resultsKAZ_CON_Dd$summary_resultsperStateQaly)
cost_kaz_mSTR_Dd <- (resultsKAZ_mSTR_Dd$summary_resultsperStateCost - resultsKAZ_CON_Dd$summary_resultsperStateCost)

paramsmSTR_kaz$qaly_AdverseEffects <-paramsmSTR_kaz$qaly_AdverseEffects + (0.28)
paramsCON_kaz$qaly_AdverseEffects <- paramsCON_kaz$qaly_AdverseEffects + (0.28)
paramsmSTR_kaz$qaly_Unresolved <-paramsmSTR_kaz$qaly_Unresolved + (0.28)
paramsCON_kaz$qaly_Unresolved <- paramsCON_kaz$qaly_Unresolved + (0.28)
paramsmSTR_kaz$qaly_LostFollowUp <-paramsmSTR_kaz$qaly_LostFollowUp + (0.28)
paramsCON_kaz$qaly_LostFollowUp <- paramsCON_kaz$qaly_LostFollowUp + (0.28)


paramsmSTR_kaz$qaly_Cured <-paramsmSTR_kaz$qaly_Cured + 0.04
paramsCON_kaz$qaly_Cured <- paramsCON_kaz$qaly_Cured + 0.04
paramsmSTR_kaz$qaly_Cured <-paramsmSTR_kaz$qaly_TreatmentCompleted + 0.04
paramsCON_kaz$qaly_Cured <- paramsCON_kaz$qaly_TreatmentCompleted + 0.04
paramsmSTR_kaz$qaly_mSTR<-paramsmSTR_kaz$qaly_mSTR + 0.04
paramsCON_kaz$qaly_mSTR <- paramsCON_kaz$qaly_mSTR + 0.04
paramsmSTR_kaz$qaly_SLtreatment<-paramsmSTR_kaz$qaly_SLtreatment + 0.04
paramsCON_kaz$qaly_SLtreatment <- paramsCON_kaz$qaly_SLtreatment + 0.04
resultsKAZ_mSTR_Ut <- model(paramsmSTR_kaz)
resultsKAZ_CON_Ut <- model(paramsCON_kaz)
qaly_kaz_mSTR_Ut <- ( resultsKAZ_mSTR_Ut$summary_resultsperStateQaly -resultsKAZ_CON_Ut$summary_resultsperStateQaly)
cost_kaz_mSTR_Ut <- (resultsKAZ_mSTR_Ut$summary_resultsperStateCost - resultsKAZ_CON_Ut$summary_resultsperStateCost)

paramsmSTR_kaz$qaly_Cured <-paramsmSTR_kaz$qaly_Cured - (0.08)
paramsCON_kaz$qaly_Cured <- paramsCON_kaz$qaly_Cured - (0.08)
paramsmSTR_kaz$qaly_Cured <-paramsmSTR_kaz$qaly_TreatmentCompleted - (0.08)
paramsCON_kaz$qaly_Cured <- paramsCON_kaz$qaly_TreatmentCompleted - (0.08)
paramsmSTR_kaz$qaly_mSTR<-paramsmSTR_kaz$qaly_mSTR - (0.08)
paramsCON_kaz$qaly_mSTR <- paramsCON_kaz$qaly_mSTR - (0.08)
paramsmSTR_kaz$qaly_SLtreatment<-paramsmSTR_kaz$qaly_SLtreatment - (0.08)
paramsCON_kaz$qaly_SLtreatment <- paramsCON_kaz$qaly_SLtreatment - (0.08)
resultsKAZ_mSTR_Dt <- model(paramsmSTR_kaz)
resultsKAZ_CON_Dt <- model(paramsCON_kaz)
cost_kaz_mSTR_Dt <- (resultsKAZ_mSTR_Dt$summary_resultsperStateCost - resultsKAZ_CON_Dt$summary_resultsperStateCost)
qaly_kaz_mSTR_Dt <- (resultsKAZ_mSTR_Dt$summary_resultsperStateQaly -resultsKAZ_CON_Dt$summary_resultsperStateQaly)

paramsmSTR_kaz$qaly_Cured <-paramsmSTR_kaz$qaly_Cured + (0.08)
paramsCON_kaz$qaly_Cured <- paramsCON_kaz$qaly_Cured + (0.08)
paramsmSTR_kaz$qaly_Cured <-paramsmSTR_kaz$qaly_TreatmentCompleted + (0.08)
paramsCON_kaz$qaly_Cured <- paramsCON_kaz$qaly_TreatmentCompleted + (0.08)
paramsmSTR_kaz$qaly_mSTR<-paramsmSTR_kaz$qaly_mSTR + (0.08)
paramsCON_kaz$qaly_mSTR <- paramsCON_kaz$qaly_mSTR + (0.08)
paramsmSTR_kaz$qaly_SLtreatment<-paramsmSTR_kaz$qaly_SLtreatment + (0.08)
paramsCON_kaz$qaly_SLtreatment <- paramsCON_kaz$qaly_SLtreatment + (0.08)


qaly_kaz_mSTR_Ud
qaly_kaz_mSTR_Dd
qaly_kaz_mSTR_Ut 
qaly_kaz_mSTR_Dt 
cost_kaz_mSTR_Ud
cost_kaz_mSTR_Dd
cost_kaz_mSTR_Ut 
cost_kaz_mSTR_Dt 



#time horizon
n_t <- 120
resultsKAZ_mSTR_120 <- model(paramsmSTR_kaz)
resultsKAZ_CON_120 <- model(paramsCON_kaz)
qaly_kaz_mSTR_120 <- (resultsKAZ_mSTR_120$summary_resultsperStateQaly -resultsKAZ_CON_120$summary_resultsperStateQaly)
cost_kaz_mSTR_120 <- (resultsKAZ_mSTR_120$summary_resultsperStateCost - resultsKAZ_CON_120$summary_resultsperStateCost)

n_t <- 360
resultsKAZ_mSTR_360 <- model(paramsmSTR_kaz)
resultsKAZ_CON_360 <- model(paramsCON_kaz)
qaly_kaz_mSTR_360 <- (resultsKAZ_mSTR_360$summary_resultsperStateQaly -resultsKAZ_CON_360$summary_resultsperStateQaly)
cost_kaz_mSTR_360 <- (resultsKAZ_mSTR_360$summary_resultsperStateCost - resultsKAZ_CON_360$summary_resultsperStateCost)

n_t <- 480
resultsKAZ_mSTR_480 <- model(paramsmSTR_kaz)
resultsKAZ_CON_480 <- model(paramsCON_kaz)
qaly_kaz_mSTR_480 <- (resultsKAZ_mSTR_480$summary_resultsperStateQaly -resultsKAZ_CON_480$summary_resultsperStateQaly)
cost_kaz_mSTR_480 <- (resultsKAZ_mSTR_480$summary_resultsperStateCost - resultsKAZ_CON_480$summary_resultsperStateCost)
n_t <- 240

qaly_kaz_mSTR_120
cost_kaz_mSTR_120
qaly_kaz_mSTR_360
cost_kaz_mSTR_360
qaly_kaz_mSTR_480
cost_kaz_mSTR_480




#Moldova###################
n_c<- 593
paramsmSTR_mol$Cost_mSTR<- (paramsmSTR_mol$Cost_mSTR*0.16)+((1-0.16)*paramsmSTR_mol$Cost_mSTR*0.75)
resultsMOL_mSTR_LB <- model(paramsmSTR_mol)
resultsMOL_CON <- model(paramsCON_mol)
qaly_mol_mSTR_LB <- (resultsMOL_mSTR_LB$summary_resultsperStateQaly -resultsMOL_CON$summary_resultsperStateQaly)
cost_mol_mSTR_LB <- (resultsMOL_mSTR_LB$summary_resultsperStateCost - resultsMOL_CON$summary_resultsperStateCost)
paramsmSTR_mol$Cost_mSTR<-  693.52

paramsCON_mol$Cost_mSTR<- (paramsCON_mol$Cost_mSTR*0.2)+((1-0.2)*paramsCON_mol$Cost_mSTR*0.75)
resultsMOL_CON_LB <- model(paramsCON_mol)
resultsMOL_mSTR <- model(paramsmSTR_mol)
qaly_mol_CON_LB <- ( resultsMOL_mSTR$summary_resultsperStateQaly -resultsMOL_CON_LB$summary_resultsperStateQaly)
cost_mol_CON_LB <- (resultsMOL_mSTR$summary_resultsperStateCost - resultsMOL_CON_LB$summary_resultsperStateCost)
paramsCON_mol$Cost_mSTR<- 451.92

paramsmSTR_mol$Cost_mSTR<- (paramsmSTR_mol$Cost_mSTR*0.16)+((1-0.16)*paramsmSTR_mol$Cost_mSTR*1.25)
resultsMOL_mSTR_UB <- model(paramsmSTR_mol)
resultsMOL_CON <- model(paramsCON_mol)
qaly_mol_mSTR_UB <- (resultsMOL_mSTR_UB$summary_resultsperStateQaly -resultsMOL_CON$summary_resultsperStateQaly)
cost_mol_mSTR_UB <- (resultsMOL_mSTR_UB$summary_resultsperStateCost - resultsMOL_CON$summary_resultsperStateCost)
paramsmSTR_mol$Cost_mSTR<- 693.52

paramsCON_mol$Cost_mSTR<-(paramsCON_mol$Cost_mSTR*0.2)+((1-0.2)*paramsCON_mol$Cost_mSTR*1.25)
resultsMOL_CON_UB <- model(paramsCON_mol)
resultsMOL_mSTR <- model(paramsmSTR_mol)
qaly_mol_CON_UB <- (resultsMOL_mSTR$summary_resultsperStateQaly -resultsMOL_CON_UB$summary_resultsperStateQaly)
cost_mol_CON_UB <- (resultsMOL_mSTR$summary_resultsperStateCost - resultsMOL_CON_UB$summary_resultsperStateCost)
paramsCON_mol$Cost_mSTR<- 451.92

qaly_mol_CON_LB
qaly_mol_CON_UB
qaly_mol_mSTR_LB
qaly_mol_mSTR_UB
cost_mol_CON_LB
cost_mol_CON_UB
cost_mol_mSTR_LB
cost_mol_mSTR_UB

#(2) Costs for SL treatment
paramsmSTR_mol$Cost_SLtreatment<- paramsmSTR_mol$Cost_SLtreatment*0.75
paramsCON_mol$Cost_SLtreatment<- paramsCON_mol$Cost_SLtreatment*0.75
resultsMOL_mSTR_SL075 <- model(paramsmSTR_mol)
resultsMOL_CON_SL075 <- model(paramsCON_mol)
cost_mol_mSTR_SL075 <- (resultsMOL_mSTR_SL075$summary_resultsperStateCost - resultsMOL_CON_SL075$summary_resultsperStateCost)
qaly_mol_mSTR_SL075 <- (resultsMOL_mSTR_SL075$summary_resultsperStateQaly -resultsMOL_CON_SL075$summary_resultsperStateQaly)
paramsmSTR_mol$Cost_SLtreatment<- paramsmSTR_mol$Cost_SLtreatment/0.75
paramsCON_mol$Cost_SLtreatment<- paramsCON_mol$Cost_SLtreatment/0.75

paramsmSTR_mol$Cost_SLtreatment<- paramsmSTR_mol$Cost_SLtreatment*1.25
paramsCON_mol$Cost_SLtreatment<- paramsCON_mol$Cost_SLtreatment*1.25
resultsMOL_mSTR_SL125 <- model(paramsmSTR_mol)
resultsMOL_CON_SL125 <- model(paramsCON_mol)
cost_mol_mSTR_SL125 <- (resultsMOL_mSTR_SL125$summary_resultsperStateCost - resultsMOL_CON_SL125$summary_resultsperStateCost)
qaly_mol_mSTR_SL125 <- (resultsMOL_mSTR_SL125$summary_resultsperStateQaly -resultsMOL_CON_SL125$summary_resultsperStateQaly)
paramsmSTR_mol$Cost_SLtreatment<- paramsmSTR_mol$Cost_SLtreatment/1.25
paramsCON_mol$Cost_SLtreatment<- paramsCON_mol$Cost_SLtreatment/1.25

qaly_mol_mSTR_SL075
cost_mol_mSTR_SL075
qaly_mol_mSTR_SL125
cost_mol_mSTR_SL125


#(3) discount rate
discount <- 0.0 
discount_m <- discount/12
resultsMOL_mSTR_d0 <- model(paramsmSTR_mol)
resultsMOL_CON_d0 <- model(paramsCON_mol)
qaly_mol_mSTR_d0 <- (resultsMOL_mSTR_d0$summary_resultsperStateQaly -resultsMOL_CON_d0$summary_resultsperStateQaly)
cost_mol_mSTR_d0 <- (resultsMOL_mSTR_d0$summary_resultsperStateCost - resultsMOL_CON_d0$summary_resultsperStateCost)
discount <- 0.06 
discount_m <- discount/12
resultsMOL_mSTR_d06 <- model(paramsmSTR_mol)
resultsMOL_CON_d06 <- model(paramsCON_mol)
qaly_mol_mSTR_d06 <- ( resultsMOL_mSTR_d06$summary_resultsperStateQaly -resultsMOL_CON_d06$summary_resultsperStateQaly)
cost_mol_mSTR_d06 <- (resultsMOL_mSTR_d06$summary_resultsperStateCost - resultsMOL_CON_d06$summary_resultsperStateCost)
discount <-0.03
discount_m <- discount/12

qaly_mol_mSTR_d0
cost_mol_mSTR_d0
qaly_mol_mSTR_d06
cost_mol_mSTR_d06

#(4) Health state utility after cured TB or and LTFU/others
paramsmSTR_mol$qaly_AdverseEffects <-paramsmSTR_mol$qaly_AdverseEffects + 0.14
paramsCON_mol$qaly_AdverseEffects <- paramsCON_mol$qaly_AdverseEffects + 0.14
paramsmSTR_mol$qaly_Unresolved <-paramsmSTR_mol$qaly_Unresolved + 0.14
paramsCON_mol$qaly_Unresolved <- paramsCON_mol$qaly_Unresolved + 0.14
paramsmSTR_mol$qaly_LostFollowUp <-paramsmSTR_mol$qaly_LostFollowUp + 0.14
paramsCON_mol$qaly_LostFollowUp <- paramsCON_mol$qaly_LostFollowUp +  0.14
resultsMOL_mSTR_Ud <- model(paramsmSTR_mol)
resultsMOL_CON_Ud <- model(paramsCON_mol)
qaly_mol_mSTR_Ud <- (resultsMOL_mSTR_Ud$summary_resultsperStateQaly -resultsMOL_CON_Ud$summary_resultsperStateQaly)
cost_mol_mSTR_Ud <- (resultsMOL_mSTR_Ud$summary_resultsperStateCost - resultsMOL_CON_Ud$summary_resultsperStateCost)

paramsmSTR_mol$qaly_AdverseEffects <-paramsmSTR_mol$qaly_AdverseEffects - (0.28)
paramsCON_mol$qaly_AdverseEffects <- paramsCON_mol$qaly_AdverseEffects - (0.28)
paramsmSTR_mol$qaly_Unresolved <-paramsmSTR_mol$qaly_Unresolved - (0.28)
paramsCON_mol$qaly_Unresolved <- paramsCON_mol$qaly_Unresolved - (0.28)
paramsmSTR_mol$qaly_LostFollowUp <-paramsmSTR_mol$qaly_LostFollowUp - (0.28)
paramsCON_mol$qaly_LostFollowUp <- paramsCON_mol$qaly_LostFollowUp - (0.28)
resultsMOL_mSTR_Dd <- model(paramsmSTR_mol)
resultsMOL_CON_Dd <- model(paramsCON_mol)
qaly_mol_mSTR_Dd <- ( resultsMOL_mSTR_Dd$summary_resultsperStateQaly -resultsMOL_CON_Dd$summary_resultsperStateQaly)
cost_mol_mSTR_Dd <- (resultsMOL_mSTR_Dd$summary_resultsperStateCost - resultsMOL_CON_Dd$summary_resultsperStateCost)

paramsmSTR_mol$qaly_AdverseEffects <-paramsmSTR_mol$qaly_AdverseEffects + (0.28)
paramsCON_mol$qaly_AdverseEffects <- paramsCON_mol$qaly_AdverseEffects + (0.28)
paramsmSTR_mol$qaly_Unresolved <-paramsmSTR_mol$qaly_Unresolved + (0.28)
paramsCON_mol$qaly_Unresolved <- paramsCON_mol$qaly_Unresolved + (0.28)
paramsmSTR_mol$qaly_LostFollowUp <-paramsmSTR_mol$qaly_LostFollowUp + (0.28)
paramsCON_mol$qaly_LostFollowUp <- paramsCON_mol$qaly_LostFollowUp + (0.28)


paramsmSTR_mol$qaly_Cured <-paramsmSTR_mol$qaly_Cured + 0.04
paramsCON_mol$qaly_Cured <- paramsCON_mol$qaly_Cured + 0.04
paramsmSTR_mol$qaly_Cured <-paramsmSTR_mol$qaly_TreatmentCompleted + 0.04
paramsCON_mol$qaly_Cured <- paramsCON_mol$qaly_TreatmentCompleted + 0.04
paramsmSTR_mol$qaly_mSTR<-paramsmSTR_mol$qaly_mSTR + 0.04
paramsCON_mol$qaly_mSTR <- paramsCON_mol$qaly_mSTR + 0.04
paramsmSTR_mol$qaly_SLtreatment<-paramsmSTR_mol$qaly_SLtreatment + 0.04
paramsCON_mol$qaly_SLtreatment <- paramsCON_mol$qaly_SLtreatment + 0.04
resultsMOL_mSTR_Ut <- model(paramsmSTR_mol)
resultsMOL_CON_Ut <- model(paramsCON_mol)
qaly_mol_mSTR_Ut <- ( resultsMOL_mSTR_Ut$summary_resultsperStateQaly -resultsMOL_CON_Ut$summary_resultsperStateQaly)
cost_mol_mSTR_Ut <- (resultsMOL_mSTR_Ut$summary_resultsperStateCost - resultsMOL_CON_Ut$summary_resultsperStateCost)

paramsmSTR_mol$qaly_Cured <-paramsmSTR_mol$qaly_Cured - (0.08)
paramsCON_mol$qaly_Cured <- paramsCON_mol$qaly_Cured - (0.08)
paramsmSTR_mol$qaly_Cured <-paramsmSTR_mol$qaly_TreatmentCompleted - (0.08)
paramsCON_mol$qaly_Cured <- paramsCON_mol$qaly_TreatmentCompleted - (0.08)
paramsmSTR_mol$qaly_mSTR<-paramsmSTR_mol$qaly_mSTR - (0.08)
paramsCON_mol$qaly_mSTR <- paramsCON_mol$qaly_mSTR - (0.08)
paramsmSTR_mol$qaly_SLtreatment<-paramsmSTR_mol$qaly_SLtreatment - (0.08)
paramsCON_mol$qaly_SLtreatment <- paramsCON_mol$qaly_SLtreatment - (0.08)
resultsMOL_mSTR_Dt <- model(paramsmSTR_mol)
resultsMOL_CON_Dt <- model(paramsCON_mol)
cost_mol_mSTR_Dt <- (resultsMOL_mSTR_Dt$summary_resultsperStateCost - resultsMOL_CON_Dt$summary_resultsperStateCost)
qaly_mol_mSTR_Dt <- (resultsMOL_mSTR_Dt$summary_resultsperStateQaly -resultsMOL_CON_Dt$summary_resultsperStateQaly)

paramsmSTR_mol$qaly_Cured <-paramsmSTR_mol$qaly_Cured + (0.08)
paramsCON_mol$qaly_Cured <- paramsCON_mol$qaly_Cured + (0.08)
paramsmSTR_mol$qaly_Cured <-paramsmSTR_mol$qaly_TreatmentCompleted + (0.08)
paramsCON_mol$qaly_Cured <- paramsCON_mol$qaly_TreatmentCompleted + (0.08)
paramsmSTR_mol$qaly_mSTR<-paramsmSTR_mol$qaly_mSTR + (0.08)
paramsCON_mol$qaly_mSTR <- paramsCON_mol$qaly_mSTR + (0.08)
paramsmSTR_mol$qaly_SLtreatment<-paramsmSTR_mol$qaly_SLtreatment + (0.08)
paramsCON_mol$qaly_SLtreatment <- paramsCON_mol$qaly_SLtreatment + (0.08)


qaly_mol_mSTR_Ud
qaly_mol_mSTR_Dd
qaly_mol_mSTR_Ut 
qaly_mol_mSTR_Dt 
cost_mol_mSTR_Ud
cost_mol_mSTR_Dd
cost_mol_mSTR_Ut 
cost_mol_mSTR_Dt


#Belarus###################
n_c<- 801
paramsmSTR_bel$Cost_mSTR<-(paramsmSTR_bel$Cost_mSTR*0.16)+((1-0.16)*paramsmSTR_bel$Cost_mSTR*0.75)
resultsBEL_mSTR_LB <- model(paramsmSTR_bel)
resultsBEL_CON <- model(paramsCON_bel)
qaly_bel_mSTR_LB <- (resultsBEL_mSTR_LB$summary_resultsperStateQaly -resultsBEL_CON$summary_resultsperStateQaly)
cost_bel_mSTR_LB <- (resultsBEL_mSTR_LB$summary_resultsperStateCost - resultsBEL_CON$summary_resultsperStateCost)
paramsmSTR_bel$Cost_mSTR<- 844.07

paramsCON_bel$Cost_mSTR<- (paramsCON_bel$Cost_mSTR*0.32)+((1-0.32)*paramsCON_bel$Cost_mSTR*0.75)
resultsBEL_CON_LB <- model(paramsCON_bel)
resultsBEL_mSTR <- model(paramsmSTR_bel)
qaly_bel_CON_LB <- ( resultsBEL_mSTR$summary_resultsperStateQaly -resultsBEL_CON_LB$summary_resultsperStateQaly)
cost_bel_CON_LB <- (resultsBEL_mSTR$summary_resultsperStateCost - resultsBEL_CON_LB$summary_resultsperStateCost)
paramsCON_bel$Cost_mSTR<- 933.53

paramsmSTR_bel$Cost_mSTR<- (paramsmSTR_bel$Cost_mSTR*0.16)+((1-0.16)*paramsmSTR_bel$Cost_mSTR*1.25)
resultsBEL_mSTR_UB <- model(paramsmSTR_bel)
resultsBEL_CON <- model(paramsCON_bel)
qaly_bel_mSTR_UB <- (resultsBEL_mSTR_UB$summary_resultsperStateQaly -resultsBEL_CON$summary_resultsperStateQaly)
cost_bel_mSTR_UB <- (resultsBEL_mSTR_UB$summary_resultsperStateCost - resultsBEL_CON$summary_resultsperStateCost)
paramsmSTR_bel$Cost_mSTR<- 844.07

paramsCON_bel$Cost_mSTR<- (paramsCON_bel$Cost_mSTR*0.32)+((1-0.32)*paramsCON_bel$Cost_mSTR*1.25)
resultsBEL_CON_UB <- model(paramsCON_bel)
resultsBEL_mSTR <- model(paramsmSTR_bel)
qaly_bel_CON_UB <- (resultsBEL_mSTR$summary_resultsperStateQaly -resultsBEL_CON_UB$summary_resultsperStateQaly)
cost_bel_CON_UB <- (resultsBEL_mSTR$summary_resultsperStateCost - resultsBEL_CON_UB$summary_resultsperStateCost)
paramsCON_bel$Cost_mSTR<- 933.53

qaly_bel_CON_LB
qaly_bel_CON_UB
qaly_bel_mSTR_LB
qaly_bel_mSTR_UB
cost_bel_CON_LB
cost_bel_CON_UB
cost_bel_mSTR_LB
cost_bel_mSTR_UB

#(2) Costs for SL treatment
paramsmSTR_bel$Cost_SLtreatment<- paramsmSTR_bel$Cost_SLtreatment*0.75
paramsCON_bel$Cost_SLtreatment<- paramsCON_bel$Cost_SLtreatment*0.75
resultsBEL_mSTR_SL075 <- model(paramsmSTR_bel)
resultsBEL_CON_SL075 <- model(paramsCON_bel)
cost_bel_mSTR_SL075 <- (resultsBEL_mSTR_SL075$summary_resultsperStateCost - resultsBEL_CON_SL075$summary_resultsperStateCost)
qaly_bel_mSTR_SL075 <- (resultsBEL_mSTR_SL075$summary_resultsperStateQaly -resultsBEL_CON_SL075$summary_resultsperStateQaly)
paramsmSTR_bel$Cost_SLtreatment<- paramsmSTR_bel$Cost_SLtreatment/0.75
paramsCON_bel$Cost_SLtreatment<- paramsCON_bel$Cost_SLtreatment/0.75

paramsmSTR_bel$Cost_SLtreatment<- paramsmSTR_bel$Cost_SLtreatment*1.25
paramsCON_bel$Cost_SLtreatment<- paramsCON_bel$Cost_SLtreatment*1.25
resultsBEL_mSTR_SL125 <- model(paramsmSTR_bel)
resultsBEL_CON_SL125 <- model(paramsCON_bel)
cost_bel_mSTR_SL125 <- (resultsBEL_mSTR_SL125$summary_resultsperStateCost - resultsBEL_CON_SL125$summary_resultsperStateCost)
qaly_bel_mSTR_SL125 <- (resultsBEL_mSTR_SL125$summary_resultsperStateQaly -resultsBEL_CON_SL125$summary_resultsperStateQaly)
paramsmSTR_bel$Cost_SLtreatment<- paramsmSTR_bel$Cost_SLtreatment/1.25
paramsCON_bel$Cost_SLtreatment<- paramsCON_bel$Cost_SLtreatment/1.25

qaly_bel_mSTR_SL075
cost_bel_mSTR_SL075
qaly_bel_mSTR_SL125
cost_bel_mSTR_SL125


#(3) discount rate
discount <- 0.0 
discount_m <- discount/12
resultsBEL_mSTR_d0 <- model(paramsmSTR_bel)
resultsBEL_CON_d0 <- model(paramsCON_bel)
qaly_bel_mSTR_d0 <- (resultsBEL_mSTR_d0$summary_resultsperStateQaly -resultsBEL_CON_d0$summary_resultsperStateQaly)
cost_bel_mSTR_d0 <- (resultsBEL_mSTR_d0$summary_resultsperStateCost - resultsBEL_CON_d0$summary_resultsperStateCost)
discount <- 0.06 
discount_m <- discount/12
resultsBEL_mSTR_d06 <- model(paramsmSTR_bel)
resultsBEL_CON_d06 <- model(paramsCON_bel)
qaly_bel_mSTR_d06 <- ( resultsBEL_mSTR_d06$summary_resultsperStateQaly -resultsBEL_CON_d06$summary_resultsperStateQaly)
cost_bel_mSTR_d06 <- (resultsBEL_mSTR_d06$summary_resultsperStateCost - resultsBEL_CON_d06$summary_resultsperStateCost)
discount <-0.03
discount_m <- discount/12

qaly_bel_mSTR_d0
cost_bel_mSTR_d0
qaly_bel_mSTR_d06
cost_bel_mSTR_d06

#(4) Health state utility after cured TB or and LTFU/others
paramsmSTR_bel$qaly_AdverseEffects <-paramsmSTR_bel$qaly_AdverseEffects + 0.14
paramsCON_bel$qaly_AdverseEffects <- paramsCON_bel$qaly_AdverseEffects + 0.14
paramsmSTR_bel$qaly_Unresolved <-paramsmSTR_bel$qaly_Unresolved + 0.14
paramsCON_bel$qaly_Unresolved <- paramsCON_bel$qaly_Unresolved + 0.14
paramsmSTR_bel$qaly_LostFollowUp <-paramsmSTR_bel$qaly_LostFollowUp + 0.14
paramsCON_bel$qaly_LostFollowUp <- paramsCON_bel$qaly_LostFollowUp +  0.14
resultsBEL_mSTR_Ud <- model(paramsmSTR_bel)
resultsBEL_CON_Ud <- model(paramsCON_bel)
qaly_bel_mSTR_Ud <- (resultsBEL_mSTR_Ud$summary_resultsperStateQaly -resultsBEL_CON_Ud$summary_resultsperStateQaly)
cost_bel_mSTR_Ud <- (resultsBEL_mSTR_Ud$summary_resultsperStateCost - resultsBEL_CON_Ud$summary_resultsperStateCost)

paramsmSTR_bel$qaly_AdverseEffects <-paramsmSTR_bel$qaly_AdverseEffects - (0.28)
paramsCON_bel$qaly_AdverseEffects <- paramsCON_bel$qaly_AdverseEffects - (0.28)
paramsmSTR_bel$qaly_Unresolved <-paramsmSTR_bel$qaly_Unresolved - (0.28)
paramsCON_bel$qaly_Unresolved <- paramsCON_bel$qaly_Unresolved - (0.28)
paramsmSTR_bel$qaly_LostFollowUp <-paramsmSTR_bel$qaly_LostFollowUp - (0.28)
paramsCON_bel$qaly_LostFollowUp <- paramsCON_bel$qaly_LostFollowUp - (0.28)
resultsBEL_mSTR_Dd <- model(paramsmSTR_bel)
resultsBEL_CON_Dd <- model(paramsCON_bel)
qaly_bel_mSTR_Dd <- ( resultsBEL_mSTR_Dd$summary_resultsperStateQaly -resultsBEL_CON_Dd$summary_resultsperStateQaly)
cost_bel_mSTR_Dd <- (resultsBEL_mSTR_Dd$summary_resultsperStateCost - resultsBEL_CON_Dd$summary_resultsperStateCost)

paramsmSTR_bel$qaly_AdverseEffects <-paramsmSTR_bel$qaly_AdverseEffects + (0.28)
paramsCON_bel$qaly_AdverseEffects <- paramsCON_bel$qaly_AdverseEffects + (0.28)
paramsmSTR_bel$qaly_Unresolved <-paramsmSTR_bel$qaly_Unresolved + (0.28)
paramsCON_bel$qaly_Unresolved <- paramsCON_bel$qaly_Unresolved + (0.28)
paramsmSTR_bel$qaly_LostFollowUp <-paramsmSTR_bel$qaly_LostFollowUp + (0.28)
paramsCON_bel$qaly_LostFollowUp <- paramsCON_bel$qaly_LostFollowUp + (0.28)


paramsmSTR_bel$qaly_Cured <-paramsmSTR_bel$qaly_Cured + 0.04
paramsCON_bel$qaly_Cured <- paramsCON_bel$qaly_Cured + 0.04
paramsmSTR_bel$qaly_Cured <-paramsmSTR_bel$qaly_TreatmentCompleted + 0.04
paramsCON_bel$qaly_Cured <- paramsCON_bel$qaly_TreatmentCompleted + 0.04
paramsmSTR_bel$qaly_mSTR<-paramsmSTR_bel$qaly_mSTR + 0.04
paramsCON_bel$qaly_mSTR <- paramsCON_bel$qaly_mSTR + 0.04
paramsmSTR_bel$qaly_SLtreatment<-paramsmSTR_bel$qaly_SLtreatment + 0.04
paramsCON_bel$qaly_SLtreatment <- paramsCON_bel$qaly_SLtreatment + 0.04
resultsBEL_mSTR_Ut <- model(paramsmSTR_bel)
resultsBEL_CON_Ut <- model(paramsCON_bel)
qaly_bel_mSTR_Ut <- ( resultsBEL_mSTR_Ut$summary_resultsperStateQaly -resultsBEL_CON_Ut$summary_resultsperStateQaly)
cost_bel_mSTR_Ut <- (resultsBEL_mSTR_Ut$summary_resultsperStateCost - resultsBEL_CON_Ut$summary_resultsperStateCost)

paramsmSTR_bel$qaly_Cured <-paramsmSTR_bel$qaly_Cured - (0.08)
paramsCON_bel$qaly_Cured <- paramsCON_bel$qaly_Cured - (0.08)
paramsmSTR_bel$qaly_Cured <-paramsmSTR_bel$qaly_TreatmentCompleted - (0.08)
paramsCON_bel$qaly_Cured <- paramsCON_bel$qaly_TreatmentCompleted - (0.08)
paramsmSTR_bel$qaly_mSTR<-paramsmSTR_bel$qaly_mSTR - (0.08)
paramsCON_bel$qaly_mSTR <- paramsCON_bel$qaly_mSTR - (0.08)
paramsmSTR_bel$qaly_SLtreatment<-paramsmSTR_bel$qaly_SLtreatment - (0.08)
paramsCON_bel$qaly_SLtreatment <- paramsCON_bel$qaly_SLtreatment - (0.08)
resultsBEL_mSTR_Dt <- model(paramsmSTR_bel)
resultsBEL_CON_Dt <- model(paramsCON_bel)
cost_bel_mSTR_Dt <- (resultsBEL_mSTR_Dt$summary_resultsperStateCost - resultsBEL_CON_Dt$summary_resultsperStateCost)
qaly_bel_mSTR_Dt <- (resultsBEL_mSTR_Dt$summary_resultsperStateQaly -resultsBEL_CON_Dt$summary_resultsperStateQaly)

paramsmSTR_bel$qaly_Cured <-paramsmSTR_bel$qaly_Cured + (0.08)
paramsCON_bel$qaly_Cured <- paramsCON_bel$qaly_Cured + (0.08)
paramsmSTR_bel$qaly_Cured <-paramsmSTR_bel$qaly_TreatmentCompleted + (0.08)
paramsCON_bel$qaly_Cured <- paramsCON_bel$qaly_TreatmentCompleted + (0.08)
paramsmSTR_bel$qaly_mSTR<-paramsmSTR_bel$qaly_mSTR + (0.08)
paramsCON_bel$qaly_mSTR <- paramsCON_bel$qaly_mSTR + (0.08)
paramsmSTR_bel$qaly_SLtreatment<-paramsmSTR_bel$qaly_SLtreatment + (0.08)
paramsCON_bel$qaly_SLtreatment <- paramsCON_bel$qaly_SLtreatment + (0.08)


qaly_bel_mSTR_Ud
qaly_bel_mSTR_Dd
qaly_bel_mSTR_Ut 
qaly_bel_mSTR_Dt 
cost_bel_mSTR_Ud
cost_bel_mSTR_Dd
cost_bel_mSTR_Ut 
cost_bel_mSTR_Dt










#Georgia###################
n_c<- 187
paramsmSTR_geo$Cost_mSTR<-(paramsmSTR_geo$Cost_mSTR*0.16)+((1-0.16)*paramsmSTR_geo$Cost_mSTR*0.75)
resultsGEO_mSTR_LB <- model(paramsmSTR_geo)
resultsGEO_CON <- model(paramsCON_geo)
qaly_geo_mSTR_LB <- (resultsGEO_mSTR_LB$summary_resultsperStateQaly -resultsGEO_CON$summary_resultsperStateQaly)
cost_geo_mSTR_LB <- (resultsGEO_mSTR_LB$summary_resultsperStateCost - resultsGEO_CON$summary_resultsperStateCost)
paramsmSTR_geo$Cost_mSTR<- 786.92

paramsCON_geo$Cost_mSTR<- (paramsCON_geo$Cost_mSTR*0.32)+((1-0.32)*paramsCON_geo$Cost_mSTR*0.75)
resultsGEO_CON_LB <- model(paramsCON_geo)
resultsGEO_mSTR <- model(paramsmSTR_geo)
qaly_geo_CON_LB <- ( resultsGEO_mSTR$summary_resultsperStateQaly -resultsGEO_CON_LB$summary_resultsperStateQaly)
cost_geo_CON_LB <- (resultsGEO_mSTR$summary_resultsperStateCost - resultsGEO_CON_LB$summary_resultsperStateCost)
paramsCON_geo$Cost_mSTR<- 643.49

paramsmSTR_geo$Cost_mSTR<- (paramsmSTR_geo$Cost_mSTR*0.16)+((1-0.16)*paramsmSTR_geo$Cost_mSTR*1.25)
resultsGEO_mSTR_UB <- model(paramsmSTR_geo)
resultsGEO_CON <- model(paramsCON_geo)
qaly_geo_mSTR_UB <- (resultsGEO_mSTR_UB$summary_resultsperStateQaly -resultsGEO_CON$summary_resultsperStateQaly)
cost_geo_mSTR_UB <- (resultsGEO_mSTR_UB$summary_resultsperStateCost - resultsGEO_CON$summary_resultsperStateCost)
paramsmSTR_geo$Cost_mSTR<- 786.92

paramsCON_geo$Cost_mSTR<- (paramsCON_geo$Cost_mSTR*0.32)+((1-0.32)*paramsCON_geo$Cost_mSTR*1.25)
resultsGEO_CON_UB <- model(paramsCON_geo)
resultsGEO_mSTR <- model(paramsmSTR_geo)
qaly_geo_CON_UB <- (resultsGEO_mSTR$summary_resultsperStateQaly -resultsGEO_CON_UB$summary_resultsperStateQaly)
cost_geo_CON_UB <- (resultsGEO_mSTR$summary_resultsperStateCost - resultsGEO_CON_UB$summary_resultsperStateCost)
paramsCON_geo$Cost_mSTR<- 643.49

qaly_geo_CON_LB
qaly_geo_CON_UB
qaly_geo_mSTR_LB
qaly_geo_mSTR_UB
cost_geo_CON_LB
cost_geo_CON_UB
cost_geo_mSTR_LB
cost_geo_mSTR_UB

#(2) Costs for SL treatment
paramsmSTR_geo$Cost_SLtreatment<- paramsmSTR_geo$Cost_SLtreatment*0.75
paramsCON_geo$Cost_SLtreatment<- paramsCON_geo$Cost_SLtreatment*0.75
resultsGEO_mSTR_SL075 <- model(paramsmSTR_geo)
resultsGEO_CON_SL075 <- model(paramsCON_geo)
cost_geo_mSTR_SL075 <- (resultsGEO_mSTR_SL075$summary_resultsperStateCost - resultsGEO_CON_SL075$summary_resultsperStateCost)
qaly_geo_mSTR_SL075 <- (resultsGEO_mSTR_SL075$summary_resultsperStateQaly -resultsGEO_CON_SL075$summary_resultsperStateQaly)
paramsmSTR_geo$Cost_SLtreatment<- paramsmSTR_geo$Cost_SLtreatment/0.75
paramsCON_geo$Cost_SLtreatment<- paramsCON_geo$Cost_SLtreatment/0.75

paramsmSTR_geo$Cost_SLtreatment<- paramsmSTR_geo$Cost_SLtreatment*1.25
paramsCON_geo$Cost_SLtreatment<- paramsCON_geo$Cost_SLtreatment*1.25
resultsGEO_mSTR_SL125 <- model(paramsmSTR_geo)
resultsGEO_CON_SL125 <- model(paramsCON_geo)
cost_geo_mSTR_SL125 <- (resultsGEO_mSTR_SL125$summary_resultsperStateCost - resultsGEO_CON_SL125$summary_resultsperStateCost)
qaly_geo_mSTR_SL125 <- (resultsGEO_mSTR_SL125$summary_resultsperStateQaly -resultsGEO_CON_SL125$summary_resultsperStateQaly)
paramsmSTR_geo$Cost_SLtreatment<- paramsmSTR_geo$Cost_SLtreatment/1.25
paramsCON_geo$Cost_SLtreatment<- paramsCON_geo$Cost_SLtreatment/1.25

qaly_geo_mSTR_SL075
cost_geo_mSTR_SL075
qaly_geo_mSTR_SL125
cost_geo_mSTR_SL125


#(3) discount rate
discount <- 0.0 
discount_m <- discount/12
resultsGEO_mSTR_d0 <- model(paramsmSTR_geo)
resultsGEO_CON_d0 <- model(paramsCON_geo)
qaly_geo_mSTR_d0 <- (resultsGEO_mSTR_d0$summary_resultsperStateQaly -resultsGEO_CON_d0$summary_resultsperStateQaly)
cost_geo_mSTR_d0 <- (resultsGEO_mSTR_d0$summary_resultsperStateCost - resultsGEO_CON_d0$summary_resultsperStateCost)
discount <- 0.06 
discount_m <- discount/12
resultsGEO_mSTR_d06 <- model(paramsmSTR_geo)
resultsGEO_CON_d06 <- model(paramsCON_geo)
qaly_geo_mSTR_d06 <- ( resultsGEO_mSTR_d06$summary_resultsperStateQaly -resultsGEO_CON_d06$summary_resultsperStateQaly)
cost_geo_mSTR_d06 <- (resultsGEO_mSTR_d06$summary_resultsperStateCost - resultsGEO_CON_d06$summary_resultsperStateCost)
discount <-0.03
discount_m <- discount/12

qaly_geo_mSTR_d0
cost_geo_mSTR_d0
qaly_geo_mSTR_d06
cost_geo_mSTR_d06

#(4) Health state utility after cured TB or and LTFU/others
paramsmSTR_geo$qaly_AdverseEffects <-paramsmSTR_geo$qaly_AdverseEffects + 0.14
paramsCON_geo$qaly_AdverseEffects <- paramsCON_geo$qaly_AdverseEffects + 0.14
paramsmSTR_geo$qaly_Unresolved <-paramsmSTR_geo$qaly_Unresolved + 0.14
paramsCON_geo$qaly_Unresolved <- paramsCON_geo$qaly_Unresolved + 0.14
paramsmSTR_geo$qaly_LostFollowUp <-paramsmSTR_geo$qaly_LostFollowUp + 0.14
paramsCON_geo$qaly_LostFollowUp <- paramsCON_geo$qaly_LostFollowUp +  0.14
resultsGEO_mSTR_Ud <- model(paramsmSTR_geo)
resultsGEO_CON_Ud <- model(paramsCON_geo)
qaly_geo_mSTR_Ud <- (resultsGEO_mSTR_Ud$summary_resultsperStateQaly -resultsGEO_CON_Ud$summary_resultsperStateQaly)
cost_geo_mSTR_Ud <- (resultsGEO_mSTR_Ud$summary_resultsperStateCost - resultsGEO_CON_Ud$summary_resultsperStateCost)

paramsmSTR_geo$qaly_AdverseEffects <-paramsmSTR_geo$qaly_AdverseEffects - (0.28)
paramsCON_geo$qaly_AdverseEffects <- paramsCON_geo$qaly_AdverseEffects - (0.28)
paramsmSTR_geo$qaly_Unresolved <-paramsmSTR_geo$qaly_Unresolved - (0.28)
paramsCON_geo$qaly_Unresolved <- paramsCON_geo$qaly_Unresolved - (0.28)
paramsmSTR_geo$qaly_LostFollowUp <-paramsmSTR_geo$qaly_LostFollowUp - (0.28)
paramsCON_geo$qaly_LostFollowUp <- paramsCON_geo$qaly_LostFollowUp - (0.28)
resultsGEO_mSTR_Dd <- model(paramsmSTR_geo)
resultsGEO_CON_Dd <- model(paramsCON_geo)
qaly_geo_mSTR_Dd <- ( resultsGEO_mSTR_Dd$summary_resultsperStateQaly -resultsGEO_CON_Dd$summary_resultsperStateQaly)
cost_geo_mSTR_Dd <- (resultsGEO_mSTR_Dd$summary_resultsperStateCost - resultsGEO_CON_Dd$summary_resultsperStateCost)

paramsmSTR_geo$qaly_AdverseEffects <-paramsmSTR_geo$qaly_AdverseEffects + (0.28)
paramsCON_geo$qaly_AdverseEffects <- paramsCON_geo$qaly_AdverseEffects + (0.28)
paramsmSTR_geo$qaly_Unresolved <-paramsmSTR_geo$qaly_Unresolved + (0.28)
paramsCON_geo$qaly_Unresolved <- paramsCON_geo$qaly_Unresolved + (0.28)
paramsmSTR_geo$qaly_LostFollowUp <-paramsmSTR_geo$qaly_LostFollowUp + (0.28)
paramsCON_geo$qaly_LostFollowUp <- paramsCON_geo$qaly_LostFollowUp + (0.28)


paramsmSTR_geo$qaly_Cured <-paramsmSTR_geo$qaly_Cured + 0.04
paramsCON_geo$qaly_Cured <- paramsCON_geo$qaly_Cured + 0.04
paramsmSTR_geo$qaly_Cured <-paramsmSTR_geo$qaly_TreatmentCompleted + 0.04
paramsCON_geo$qaly_Cured <- paramsCON_geo$qaly_TreatmentCompleted + 0.04
paramsmSTR_geo$qaly_mSTR<-paramsmSTR_geo$qaly_mSTR + 0.04
paramsCON_geo$qaly_mSTR <- paramsCON_geo$qaly_mSTR + 0.04
paramsmSTR_geo$qaly_SLtreatment<-paramsmSTR_geo$qaly_SLtreatment + 0.04
paramsCON_geo$qaly_SLtreatment <- paramsCON_geo$qaly_SLtreatment + 0.04
resultsGEO_mSTR_Ut <- model(paramsmSTR_geo)
resultsGEO_CON_Ut <- model(paramsCON_geo)
qaly_geo_mSTR_Ut <- ( resultsGEO_mSTR_Ut$summary_resultsperStateQaly -resultsGEO_CON_Ut$summary_resultsperStateQaly)
cost_geo_mSTR_Ut <- (resultsGEO_mSTR_Ut$summary_resultsperStateCost - resultsGEO_CON_Ut$summary_resultsperStateCost)

paramsmSTR_geo$qaly_Cured <-paramsmSTR_geo$qaly_Cured - (0.08)
paramsCON_geo$qaly_Cured <- paramsCON_geo$qaly_Cured - (0.08)
paramsmSTR_geo$qaly_Cured <-paramsmSTR_geo$qaly_TreatmentCompleted - (0.08)
paramsCON_geo$qaly_Cured <- paramsCON_geo$qaly_TreatmentCompleted - (0.08)
paramsmSTR_geo$qaly_mSTR<-paramsmSTR_geo$qaly_mSTR - (0.08)
paramsCON_geo$qaly_mSTR <- paramsCON_geo$qaly_mSTR - (0.08)
paramsmSTR_geo$qaly_SLtreatment<-paramsmSTR_geo$qaly_SLtreatment - (0.08)
paramsCON_geo$qaly_SLtreatment <- paramsCON_geo$qaly_SLtreatment - (0.08)
resultsGEO_mSTR_Dt <- model(paramsmSTR_geo)
resultsGEO_CON_Dt <- model(paramsCON_geo)
cost_geo_mSTR_Dt <- (resultsGEO_mSTR_Dt$summary_resultsperStateCost - resultsGEO_CON_Dt$summary_resultsperStateCost)
qaly_geo_mSTR_Dt <- (resultsGEO_mSTR_Dt$summary_resultsperStateQaly -resultsGEO_CON_Dt$summary_resultsperStateQaly)

paramsmSTR_geo$qaly_Cured <-paramsmSTR_geo$qaly_Cured + (0.08)
paramsCON_geo$qaly_Cured <- paramsCON_geo$qaly_Cured + (0.08)
paramsmSTR_geo$qaly_Cured <-paramsmSTR_geo$qaly_TreatmentCompleted + (0.08)
paramsCON_geo$qaly_Cured <- paramsCON_geo$qaly_TreatmentCompleted + (0.08)
paramsmSTR_geo$qaly_mSTR<-paramsmSTR_geo$qaly_mSTR + (0.08)
paramsCON_geo$qaly_mSTR <- paramsCON_geo$qaly_mSTR + (0.08)
paramsmSTR_geo$qaly_SLtreatment<-paramsmSTR_geo$qaly_SLtreatment + (0.08)
paramsCON_geo$qaly_SLtreatment <- paramsCON_geo$qaly_SLtreatment + (0.08)


qaly_geo_mSTR_Ud
qaly_geo_mSTR_Dd
qaly_geo_mSTR_Ut 
qaly_geo_mSTR_Dt 
cost_geo_mSTR_Ud
cost_geo_mSTR_Dd
cost_geo_mSTR_Ut 
cost_geo_mSTR_Dt



#####-------------------------------######

#DATASET FOR UNIVARIATE ANALYSES####
#Creating a matrix with all univariate sensitivity analyses
labe_names <- c("Country","caracteristic","value1q","Value2q","Value1c","Value2c")
univ_res <- matrix(NA_real_, nrow= 24, ncol=6, dimnames= list(from= 1:24, to = labe_names))
univ_res[1:6, "Country"]<- 'Kazakhstan'
univ_res[7:12, "Country"]<- 'Moldova'
univ_res[13:18, "Country"]<- 'Belarus'
univ_res[19:24, "Country"]<- 'Georgia'

labe2_names <- c("Total non-drug costs mSTR treatment*","Total non-drug costs conventional treatment*","Total costs for SL treatment","Discount rate","Utility MDR/RR TB unresolved","Utility MDR/RR TB cured")
univ_res[1:6,2] <-labe2_names
for (x in 1:3) {
  univ_res[(6*x+1):(6*x+6),2] <-labe2_names
}
#kazakhstan
univ_res[1,3]<- round(qaly_kaz_mSTR_LB,1)
univ_res[1,4]<- round(qaly_kaz_mSTR_UB,1)
univ_res[2,3]<- round(qaly_kaz_CON_LB,1)
univ_res[2,4]<- round(qaly_kaz_CON_UB,1)
univ_res[3,3]<- round(qaly_kaz_mSTR_SL075,1)
univ_res[3,4]<- round(qaly_kaz_mSTR_SL125,1)
univ_res[4,3]<- round(qaly_kaz_mSTR_d0,1)
univ_res[4,4]<- round(qaly_kaz_mSTR_d06,1)
univ_res[5,3]<- round(qaly_kaz_mSTR_Dd,1)
univ_res[5,4]<- round(qaly_kaz_mSTR_Ud,1)
univ_res[6,3]<- round(qaly_kaz_mSTR_Dt,1)
univ_res[6,4]<- round(qaly_kaz_mSTR_Ut,1)

univ_res[1,5]<- round(cost_kaz_mSTR_LB,0)
univ_res[1,6]<- round(cost_kaz_mSTR_UB,0)
univ_res[2,5]<- round(cost_kaz_CON_LB,0)
univ_res[2,6]<- round(cost_kaz_CON_UB,0)
univ_res[3,5]<- round(cost_kaz_mSTR_SL075,0)
univ_res[3,6]<- round(cost_kaz_mSTR_SL125,0)
univ_res[4,5]<- round(cost_kaz_mSTR_d0,0)
univ_res[4,6]<- round(cost_kaz_mSTR_d06,0)
univ_res[5,5]<- round(cost_kaz_mSTR_Dd,0)
univ_res[5,6]<- round(cost_kaz_mSTR_Ud,0)
univ_res[6,5]<- round(cost_kaz_mSTR_Dt,0)
univ_res[6,6]<- round(cost_kaz_mSTR_Ut,0)
#moldova
univ_res[7,3]<- round(qaly_mol_mSTR_LB,1)
univ_res[7,4]<- round(qaly_mol_mSTR_UB,1)
univ_res[8,3]<- round(qaly_mol_CON_LB,1)
univ_res[8,4]<- round(qaly_mol_CON_UB,1)
univ_res[9,3]<- round(qaly_mol_mSTR_SL075,1)
univ_res[9,4]<- round(qaly_mol_mSTR_SL125,1)
univ_res[10,3]<- round(qaly_mol_mSTR_d0,1)
univ_res[10,4]<- round(qaly_mol_mSTR_d06,1)
univ_res[11,3]<- round(qaly_mol_mSTR_Dd,1)
univ_res[11,4]<- round(qaly_mol_mSTR_Ud,1)
univ_res[12,3]<- round(qaly_mol_mSTR_Dt,1)
univ_res[12,4]<- round(qaly_mol_mSTR_Ut,1)

univ_res[7,5]<- round(cost_mol_mSTR_LB,0)
univ_res[7,6]<- round(cost_mol_mSTR_UB,0)
univ_res[8,5]<- round(cost_mol_CON_LB,0)
univ_res[8,6]<- round(cost_mol_CON_UB,0)
univ_res[9,5]<- round(cost_mol_mSTR_SL075,0)
univ_res[9,6]<- round(cost_mol_mSTR_SL125,0)
univ_res[10,5]<- round(cost_mol_mSTR_d0,0)
univ_res[10,6]<- round(cost_mol_mSTR_d06,0)
univ_res[11,5]<- round(cost_mol_mSTR_Dd,0)
univ_res[11,6]<- round(cost_mol_mSTR_Ud,0)
univ_res[12,5]<- round(cost_mol_mSTR_Dt,0)
univ_res[12,6]<- round(cost_mol_mSTR_Ut,0)
#belarus
univ_res[13,3]<- round(qaly_bel_mSTR_LB,1)
univ_res[13,4]<- round(qaly_bel_mSTR_UB,1)
univ_res[14,3]<- round(qaly_bel_CON_LB,1)
univ_res[14,4]<- round(qaly_bel_CON_UB,1)
univ_res[15,3]<- round(qaly_bel_mSTR_SL075,1)
univ_res[15,4]<- round(qaly_bel_mSTR_SL125,1)
univ_res[16,3]<- round(qaly_bel_mSTR_d0,1)
univ_res[16,4]<- round(qaly_bel_mSTR_d06,1)
univ_res[17,3]<- round(qaly_bel_mSTR_Dd,1)
univ_res[17,4]<- round(qaly_bel_mSTR_Ud,1)
univ_res[18,3]<- round(qaly_bel_mSTR_Dt,1)
univ_res[18,4]<- round(qaly_bel_mSTR_Ut,1)

univ_res[13,5]<- round(cost_bel_mSTR_LB,0)
univ_res[13,6]<- round(cost_bel_mSTR_UB,0)
univ_res[14,5]<- round(cost_bel_CON_LB,0)
univ_res[14,6]<- round(cost_bel_CON_UB,0)
univ_res[15,5]<- round(cost_bel_mSTR_SL075,0)
univ_res[15,6]<- round(cost_bel_mSTR_SL125,0)
univ_res[16,5]<- round(cost_bel_mSTR_d0,0)
univ_res[16,6]<- round(cost_bel_mSTR_d06,0)
univ_res[17,5]<- round(cost_bel_mSTR_Dd,0)
univ_res[17,6]<- round(cost_bel_mSTR_Ud,0)
univ_res[18,5]<- round(cost_bel_mSTR_Dt,0)
univ_res[18,6]<- round(cost_bel_mSTR_Ut,0)

#Georgia
univ_res[19,3]<- round(qaly_geo_mSTR_LB,1)
univ_res[19,4]<- round(qaly_geo_mSTR_UB,1)
univ_res[20,3]<- round(qaly_geo_CON_LB,1)
univ_res[20,4]<- round(qaly_geo_CON_UB,1)
univ_res[21,3]<- round(qaly_geo_mSTR_SL075,1)
univ_res[21,4]<- round(qaly_geo_mSTR_SL125,1)
univ_res[22,3]<- round(qaly_geo_mSTR_d0,1)
univ_res[22,4]<- round(qaly_geo_mSTR_d06,1)
univ_res[23,3]<- round(qaly_geo_mSTR_Dd,1)
univ_res[23,4]<- round(qaly_geo_mSTR_Ud,1)
univ_res[24,3]<- round(qaly_geo_mSTR_Dt,1)
univ_res[24,4]<- round(qaly_geo_mSTR_Ut,1)

univ_res[19,5]<- round(cost_geo_mSTR_LB,0)
univ_res[19,6]<- round(cost_geo_mSTR_UB,0)
univ_res[20,5]<- round(cost_geo_CON_LB,0)
univ_res[20,6]<- round(cost_geo_CON_UB,0)
univ_res[21,5]<- round(cost_geo_mSTR_SL075,0)
univ_res[21,6]<- round(cost_geo_mSTR_SL125,0)
univ_res[22,5]<- round(cost_geo_mSTR_d0,0)
univ_res[22,6]<- round(cost_geo_mSTR_d06,0)
univ_res[23,5]<- round(cost_geo_mSTR_Dd,0)
univ_res[23,6]<- round(cost_geo_mSTR_Ud,0)
univ_res[24,5]<- round(cost_geo_mSTR_Dt,0)
univ_res[24,6]<- round(cost_geo_mSTR_Ut,0)

view(univ_res)
univ_resDf <- as.data.frame(univ_res)
#Graph Univariate sensitivity analyses #####
# import data
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
data_bel <- as.data.frame(filter(univ_resDf, Country=="Belarus"))
data_mol <- as.data.frame(filter(univ_resDf, Country=="Moldova"))
data_kaz <- as.data.frame(filter(univ_resDf, Country=="Kazakhstan"))
data_geo <- as.data.frame(filter(univ_resDf, Country=="Georgia"))


#costs & qalys separated
data_bel_cost <- as.data.frame(filter(data_bel, caracteristic == "Total non-drug costs mSTR treatment*" |  caracteristic =="Total non-drug costs conventional treatment*" | caracteristic == "Total costs for SL treatment" | caracteristic =="Discount rate"))
data_bel_qaly <- as.data.frame(filter(data_bel, caracteristic == "Utility MDR/RR TB cured" | caracteristic =="Utility MDR/RR TB unresolved" | caracteristic =="Discount rate"))
data_mol_cost <- as.data.frame(filter(data_mol, caracteristic == "Total non-drug costs mSTR treatment*" | caracteristic =="Total non-drug costs conventional treatment*"| caracteristic == "Total costs for SL treatment"  | caracteristic =="Discount rate"))
data_mol_qaly <- as.data.frame(filter(data_mol, caracteristic == "Utility MDR/RR TB cured" | caracteristic =="Utility MDR/RR TB unresolved" | caracteristic =="Discount rate"))
data_kaz_cost <- as.data.frame(filter(data_kaz, caracteristic == "Total non-drug costs mSTR treatment*" | caracteristic =="Total non-drug costs conventional treatment*"| caracteristic == "Total costs for SL treatment"  | caracteristic =="Discount rate"))
data_kaz_qaly <- as.data.frame(filter(data_kaz, caracteristic == "Utility MDR/RR TB cured" | caracteristic =="Utility MDR/RR TB unresolved" | caracteristic =="Discount rate"))
data_geo_cost <- as.data.frame(filter(data_geo, caracteristic == "Total non-drug costs mSTR treatment*" | caracteristic =="Total non-drug costs conventional treatment*"| caracteristic == "Total costs for SL treatment"  | caracteristic =="Discount rate"))
data_geo_qaly <- as.data.frame(filter(data_geo, caracteristic == "Utility MDR/RR TB cured" | caracteristic =="Utility MDR/RR TB unresolved" | caracteristic =="Discount rate"))

#Costs incremental
a1<- ggplot(data_bel_cost) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=Value1c, yend=Value2c), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=Value1c), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2c), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Belarus") +
  xlab("") +
  ylab("Incremental economic costs")+
  theme(plot.title = element_text(face = "bold"))

a2<- ggplot(data_kaz_cost) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=Value1c, yend=Value2c), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=Value1c), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2c), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Kazakhstan") +
  xlab("") +
  ylab("Incremental economic costs")+
  theme(plot.title = element_text(face = "bold"))

a3<- ggplot(data_mol_cost) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=Value2c, yend=Value1c), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=Value1c), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2c), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Moldova") +
  xlab("") +
  ylab("Incremental economic costs")+
  theme(plot.title = element_text(face = "bold"))

a4<- ggplot(data_geo_cost) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=Value2c, yend=Value1c), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=Value1c), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2c), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Georgia") +
  xlab("") +
  ylab("Incremental economic costs")+
  theme(plot.title = element_text(face = "bold"))

setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/Bloomsbury_PolicyLab/0_Article/Figures")
tiff(file="bel_mol_kaz_geo_univSA_incCost.tiff",
     width=14, height=9.5, units="in", res=1000)
ggarrange(a1, a2, a3, a4,
          labels = c("(A)", "(B)", "(C)", "(D)"),
          ncol = 2, nrow = 2,common.legend = TRUE)
dev.off()
#Qalys incremental
b1<- ggplot(data_bel_qaly) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=value1q, yend=Value2q), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=value1q), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2q), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Belarus") +
  xlab("") +
  ylab("Incremental QALYs")+
  theme(plot.title = element_text(face = "bold"))

b2<- ggplot(data_kaz_qaly) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=value1q, yend=Value2q), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=value1q), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2q), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Kazakhstan") +
  xlab("") +
  ylab("Incremental QALYs")+
  theme(plot.title = element_text(face = "bold"))

b3<- ggplot(data_mol_qaly) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=value1q, yend=Value2q), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=value1q), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2q), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Moldova") +
  xlab("") +
  ylab("Incremental QALYs")+
  theme(plot.title = element_text(face = "bold"))

b4<- ggplot(data_geo_qaly) +
  geom_segment( aes(x=caracteristic, xend=caracteristic, y=value1q, yend=Value2q), color="cornsilk3", linewidth=2) +
  geom_point( aes(x=caracteristic, y=value1q), color="salmon2", size=4 ) +
  geom_point( aes(x=caracteristic, y=Value2q), color="darkblue", size=4 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  labs(title="Georgia") +
  xlab("") +
  ylab("Incremental QALYs")+
  theme(plot.title = element_text(face = "bold"))

setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/Bloomsbury_PolicyLab/0_Article/Figures")
tiff(file="bel_mol_kaz_geo_univSA_incQaly.tiff",
     width=14, height=9, units="in", res=800)
ggarrange(b1, b2, b3, b4,
          labels = c("(A)", "(B)", "(C)", "(D)"),
          ncol = 2, nrow = 2)
dev.off()



################################################################################

#------------------------------------------------------------------------------#
# III. SCENARIO ANALYSES#########################################################
#------------------------------------------------------------------------------#
# Parameters Kazakhstan--------------------------------------------------------------
#PARAMETERS LIST BELOW (change costs accordingly, QALYS are the same for everycountry)
n_c<- 3755
paramsmSTR_kaz <- list(
  #Costs & QALYS
  Cost_mSTR = 635.14,
  Cost_SLtreatment = 658.4,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  635.14,
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.033*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/9)*(0.1035*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.7614)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/9)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.1018)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.0840+0.0336),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.7267,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.1252,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.28/12))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

paramsCON_kaz <- list(
  #Costs & QALYS
  Cost_mSTR = 658.40,
  Cost_SLtreatment = 658.4,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  658.4,
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/20)*(0.033*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/20)*(0.1035*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/20)*((0.7614)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/20)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/20)*((0.1018)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.0840+0.0336),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.7267,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.1252,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)




# Parameters Moldova--------------------------------------------------------------
paramsmSTR_mol <- list(
  #Costs & QALYS
  Cost_mSTR = 693.52,
  Cost_SLtreatment = 451.92,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  693.52,
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0626*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/9)*(0.1091*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.6923)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/9)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.1360)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.1892+0.1351),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.5676,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.1081,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.28/12))**0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)
# now conventional treatment, mSTR labels below means Conventional.
paramsCON_mol <- list(
  #Costs & QALYS
  Cost_mSTR = 451.92,
  Cost_SLtreatment = 451.92,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  451.92,
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/18)*(0.0626*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/18)*(0.1091*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/18)*((0.6923)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/18)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/18)*((0.1360)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.1892+0.1351),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.5676,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.1081,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)


# Parameters Belarus--------------------------------------------------------------
n_c<- 801
paramsmSTR_bel <- list(
  #Costs & QALYS
  Cost_mSTR = 844.07,
  Cost_SLtreatment = 933.53,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects = 933.53,
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0168*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/9)*(0.0712*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.8290)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/9)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.0829)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.0662+0.1324),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.7178,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.0836,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.28/12,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

paramsCON_bel <- list(
  #Costs & QALYS
  Cost_mSTR = 933.533,
  Cost_SLtreatment = 933.533,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  933.533,
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/15.57)*(0.0168*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/15.57)*(0.0712*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/15.57)*((0.8290)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/15.57)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/15.57)*((0.0829)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/20)*(0.0662+0.1324),
  p_SLtreatment_TreatmentCompleted = (1/20)*0.7178,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/20)*0.0836,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
)


# Parameters Georgia--------------------------------------------------------------
n_c<- 187
paramsmSTR_geo <- list(
  #Costs & QALYS
  Cost_mSTR = 786.92,
  Cost_SLtreatment = 643.49,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  643.49*(1-0.24),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/9)*(0.0183*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/9)*(0.1376*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/9)*((0.7844)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/9)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/9)*((0.0596)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/17)*(0.0943+0.2642),
  p_SLtreatment_TreatmentCompleted = (1/17)*0.5660,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/17)*0.0755,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
  
)

paramsCON_geo <- list(
  #Costs & QALYS
  Cost_mSTR = 643.49,
  Cost_SLtreatment = 643.49,
  Cost_LostFollowUp =  1,
  Cost_Unresolved =  0,
  Cost_TreatmentCompleted =  0,
  Cost_AdverseEffects =  643.49*(1-0.32),
  Cost_Cured =  0,
  Cost_Dead =  0,
  qaly_mSTR = 0.81,
  qaly_SLtreatment = 0.79,
  qaly_LostFollowUp = 0.68,
  qaly_Unresolved = 0.68,
  qaly_TreatmentCompleted = 0.81,
  qaly_AdverseEffects = 0.68,
  qaly_Cured = 0.81,
  qaly_Dead = 0,
  #Matrix of transition probabilities {Define here whether we incorporate monthly probs or not}
  p_mSTR_SLtreatment =  (1/17)*(0.0183*(1-0.2351)),
  p_mSTR_LostFollowUp = (1/17)*(0.1376*(1-0.2351)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = (1/17)*((0.7844)*(1-0.2351)),
  p_mSTR_AdverseEffects = (1/17)*(0.2351),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  (1/17)*((0.0596)*(1-0.2351)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = (1/17)*(0.0943+0.2642),
  p_SLtreatment_TreatmentCompleted = (1/17)*0.5660,
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = (1/17)*0.0755,  
  
  p_LostFollowUp_mSTR = 0.28/12,
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = 0.0686,
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  0.0686,
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = 0.023,
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = (1-(0.023))*0.08333, 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = 0.5*(0.697),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0.5*(0.195),
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = 0.5*(0.108), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =0  
)

#----------------------------------------------------------#
##Results scenario analyses#################################
#----------------------------------------------------------#
#KAZAHSTAN!
n_c<- 3755
resultsKAZ_mSTR <- model(paramsmSTR_kaz)
resultsKAZ_CON <- model(paramsCON_kaz)
ICER_kaz <- (resultsKAZ_mSTR$summary_resultsperStateCost - resultsKAZ_CON$summary_resultsperStateCost)/( resultsKAZ_mSTR$summary_resultsperStateQaly -resultsKAZ_CON$summary_resultsperStateQaly)
resultsKAZ_mSTR$summary_resultsperStateCost
resultsKAZ_CON$summary_resultsperStateCost
resultsKAZ_mSTR$summary_resultsperStateQaly
resultsKAZ_CON$summary_resultsperStateQaly
excess_mort <- resultsKAZ_CON$individual_shifts$Dead[240]-resultsKAZ_mSTR$individual_shifts$Dead[240]
resultsKAZ_mSTR$summary_resultsSUCC_t
resultsKAZ_CON$summary_resultsSUCC_t
resultsKAZ_mSTR$individual_shifts$Dead[240]/n_c
resultsKAZ_CON$individual_shifts$Dead[240]/n_c
resultsKAZ_mSTR$individual_shifts$Cured[240]/n_c
resultsKAZ_CON$individual_shifts$Cured[240]/n_c

#Moldova!
n_c<- 593
resultsmol_mSTR <- model(paramsmSTR_mol)
resultsmol_CON <- model(paramsCON_mol)
ICER_mol <- (resultsmol_mSTR$summary_resultsperStateCost - resultsmol_CON$summary_resultsperStateCost)/( resultsmol_mSTR$summary_resultsperStateQaly -resultsmol_CON$summary_resultsperStateQaly)
resultsmol_mSTR$summary_resultsperStateCost
resultsmol_CON$summary_resultsperStateCost
resultsmol_mSTR$summary_resultsperStateQaly
resultsmol_CON$summary_resultsperStateQaly
excess_mort <- resultsmol_CON$individual_shifts$Dead[240]-resultsmol_mSTR$individual_shifts$Dead[240]
resultsmol_mSTR$summary_resultsSUCC_t
resultsmol_CON$summary_resultsSUCC_t
resultsmol_mSTR$individual_shifts$Dead[240]/n_c
resultsmol_CON$individual_shifts$Dead[240]/n_c
resultsmol_mSTR$individual_shifts$Cured[240]/n_c
resultsmol_CON$individual_shifts$Cured[240]/n_c

#Belarus
n_c<- 801
resultsbel_mSTR <- model(paramsmSTR_bel)
resultsbel_CON <- model(paramsCON_bel)
ICER_bel<- (resultsbel_mSTR$summary_resultsperStateCost - resultsbel_CON$summary_resultsperStateCost)/( resultsbel_mSTR$summary_resultsperStateQaly -resultsbel_CON$summary_resultsperStateQaly)
resultsbel_mSTR$summary_resultsperStateCost
resultsbel_CON$summary_resultsperStateCost
resultsbel_mSTR$summary_resultsperStateQaly
resultsbel_CON$summary_resultsperStateQaly
excess_mort <- resultsbel_CON$individual_shifts$Dead[240]-resultsbel_mSTR$individual_shifts$Dead[240]
resultsbel_mSTR$summary_resultsSUCC_t
resultsbel_CON$summary_resultsSUCC_t
resultsbel_mSTR$individual_shifts$Dead[240]/n_c
resultsbel_CON$individual_shifts$Dead[240]/n_c
resultsbel_mSTR$individual_shifts$Cured[240]/n_c
resultsbel_CON$individual_shifts$Cured[240]/n_c

#Georgia
n_c<- 187
resultsgeo_mSTR <- model(paramsmSTR_geo)
resultsgeo_CON <- model(paramsCON_geo)
ICER_geo<- (resultsgeo_mSTR$summary_resultsperStateCost - resultsgeo_CON$summary_resultsperStateCost)/( resultsgeo_mSTR$summary_resultsperStateQaly -resultsgeo_CON$summary_resultsperStateQaly)
resultsgeo_mSTR$summary_resultsperStateCost
resultsgeo_CON$summary_resultsperStateCost
resultsgeo_mSTR$summary_resultsperStateQaly
resultsgeo_CON$summary_resultsperStateQaly
excess_mort <- resultsgeo_CON$individual_shifts$Dead[240]-resultsgeo_mSTR$individual_shifts$Dead[240]
resultsgeo_mSTR$summary_resultsSUCC_t
resultsgeo_CON$summary_resultsSUCC_t

resultsgeo_mSTR$individual_shifts$Dead[240]/n_c
resultsgeo_CON$individual_shifts$Dead[240]/n_c
resultsgeo_mSTR$individual_shifts$Cured[240]/n_c
resultsgeo_CON$individual_shifts$Cured[240]/n_c

##Other countries with no-data####
#Armenia!
n_c<- 64
resultsarm_mSTR <- model(paramsmSTR_arm)
resultsarm_CON <- model(paramsCON_arm)
ICER_arm <- (resultsarm_mSTR$summary_resultsperStateCost - resultsarm_CON$summary_resultsperStateCost)/( resultsarm_mSTR$summary_resultsperStateQaly -resultsarm_CON$summary_resultsperStateQaly)
resultsarm_mSTR$summary_resultsperStateCost
resultsarm_CON$summary_resultsperStateCost
resultsarm_mSTR$summary_resultsperStateQaly
resultsarm_CON$summary_resultsperStateQaly
excess_mort <- resultsarm_CON$individual_shifts$Dead[240]-resultsarm_mSTR$individual_shifts$Dead[240]
resultsarm_mSTR$summary_resultsSUCC_t
resultsarm_CON$summary_resultsSUCC_t
resultsarm_mSTR$individual_shifts$Dead[240]/n_c
resultsarm_CON$individual_shifts$Dead[240]/n_c
resultsarm_mSTR$individual_shifts$Cured[240]/n_c
resultsarm_CON$individual_shifts$Cured[240]/n_c


#Azerbaijan!
n_c<- 1040
resultsaze_mSTR <- model(paramsmSTR_aze)
resultsaze_CON <- model(paramsCON_aze)
ICER_aze<- (resultsaze_mSTR$summary_resultsperStateCost - resultsaze_CON$summary_resultsperStateCost)/( resultsaze_mSTR$summary_resultsperStateQaly -resultsaze_CON$summary_resultsperStateQaly)
resultsaze_mSTR$summary_resultsperStateCost
resultsaze_CON$summary_resultsperStateCost
resultsaze_mSTR$summary_resultsperStateQaly
resultsaze_CON$summary_resultsperStateQaly
excess_mort <- resultsaze_CON$individual_shifts$Dead[240]-resultsaze_mSTR$individual_shifts$Dead[240]
resultsaze_mSTR$summary_resultsSUCC_t
resultsaze_CON$summary_resultsSUCC_t
resultsaze_mSTR$individual_shifts$Dead[240]/n_c
resultsaze_CON$individual_shifts$Dead[240]/n_c
resultsaze_mSTR$individual_shifts$Cured[240]/n_c
resultsaze_CON$individual_shifts$Cured[240]/n_c


#Ukraine
n_c<- 4046
resultsukr_mSTR <- model(paramsmSTR_ukr)
resultsukr_CON <- model(paramsCON_ukr)
ICER_ukr<- (resultsukr_mSTR$summary_resultsperStateCost - resultsukr_CON$summary_resultsperStateCost)/( resultsukr_mSTR$summary_resultsperStateQaly -resultsukr_CON$summary_resultsperStateQaly)
resultsukr_mSTR$summary_resultsperStateCost
resultsukr_CON$summary_resultsperStateCost
resultsukr_mSTR$summary_resultsperStateQaly
resultsukr_CON$summary_resultsperStateQaly
excess_mort <- resultsukr_CON$individual_shifts$Dead[240]-resultsukr_mSTR$individual_shifts$Dead[240]
resultsukr_mSTR$summary_resultsSUCC_t
resultsukr_CON$summary_resultsSUCC_t

resultsukr_mSTR$individual_shifts$Dead[240]/n_c
resultsukr_CON$individual_shifts$Dead[240]/n_c
resultsukr_mSTR$individual_shifts$Cured[240]/n_c
resultsukr_CON$individual_shifts$Cured[240]/n_c

#Uzbekistan
n_c<- 2147
resultsuzb_mSTR <- model(paramsmSTR_uzb)
resultsuzb_CON <- model(paramsCON_uzb)
ICER_uzb<- (resultsuzb_mSTR$summary_resultsperStateCost - resultsuzb_CON$summary_resultsperStateCost)/( resultsuzb_mSTR$summary_resultsperStateQaly -resultsuzb_CON$summary_resultsperStateQaly)
resultsuzb_mSTR$summary_resultsperStateCost
resultsuzb_CON$summary_resultsperStateCost
resultsuzb_mSTR$summary_resultsperStateQaly
resultsuzb_CON$summary_resultsperStateQaly
excess_mort <- resultsuzb_CON$individual_shifts$Dead[240]-resultsuzb_mSTR$individual_shifts$Dead[240]
resultsuzb_mSTR$summary_resultsSUCC_t
resultsuzb_CON$summary_resultsSUCC_t
resultsuzb_mSTR$individual_shifts$Dead[240]/n_c
resultsuzb_CON$individual_shifts$Dead[240]/n_c
resultsuzb_mSTR$individual_shifts$Cured[240]/n_c
resultsuzb_CON$individual_shifts$Cured[240]/n_c

#Tajikistan
n_c<- 606
resultstaj_mSTR <- model(paramsmSTR_taj)
resultstaj_CON <- model(paramsCON_taj)
ICER_taj<- (resultstaj_mSTR$summary_resultsperStateCost - resultstaj_CON$summary_resultsperStateCost)/( resultstaj_mSTR$summary_resultsperStateQaly -resultstaj_CON$summary_resultsperStateQaly)
resultstaj_mSTR$summary_resultsperStateCost
resultstaj_CON$summary_resultsperStateCost
resultstaj_mSTR$summary_resultsperStateQaly
resultstaj_CON$summary_resultsperStateQaly
excess_mort <- resultstaj_CON$individual_shifts$Dead[240]-resultstaj_mSTR$individual_shifts$Dead[240]
resultstaj_mSTR$summary_resultsSUCC_t
resultstaj_CON$summary_resultsSUCC_t
resultstaj_mSTR$individual_shifts$Dead[240]/n_c
resultstaj_CON$individual_shifts$Dead[240]/n_c
resultstaj_mSTR$individual_shifts$Cured[240]/n_c
resultstaj_CON$individual_shifts$Cured[240]/n_c

#Turkmenistan
n_c<- 808
resultstur_mSTR <- model(paramsmSTR_tur)
resultstur_CON <- model(paramsCON_tur)
ICER_tur<- (resultstur_mSTR$summary_resultsperStateCost - resultstur_CON$summary_resultsperStateCost)/( resultstur_mSTR$summary_resultsperStateQaly -resultstur_CON$summary_resultsperStateQaly)
resultstur_mSTR$summary_resultsperStateCost
resultstur_CON$summary_resultsperStateCost
resultstur_mSTR$summary_resultsperStateQaly
resultstur_CON$summary_resultsperStateQaly
excess_mort <- resultstur_CON$individual_shifts$Dead[240]-resultstur_mSTR$individual_shifts$Dead[240]
resultstur_mSTR$summary_resultsSUCC_t
resultstur_CON$summary_resultsSUCC_t

resultstur_mSTR$individual_shifts$Dead[240]/n_c
resultstur_CON$individual_shifts$Dead[240]/n_c
resultstur_mSTR$individual_shifts$Cured[240]/n_c
resultstur_CON$individual_shifts$Cured[240]/n_c


#Kyrgistan
n_c<- 918
resultskyr_mSTR <- model(paramsmSTR_kyr)
resultskyr_CON <- model(paramsCON_kyr)
ICER_kyr<- (resultskyr_mSTR$summary_resultsperStateCost - resultskyr_CON$summary_resultsperStateCost)/( resultskyr_mSTR$summary_resultsperStateQaly -resultskyr_CON$summary_resultsperStateQaly)
resultskyr_mSTR$summary_resultsperStateCost
resultskyr_CON$summary_resultsperStateCost
resultskyr_mSTR$summary_resultsperStateQaly
resultskyr_CON$summary_resultsperStateQaly
excess_mort <- resultskyr_CON$individual_shifts$Dead[240]-resultskyr_mSTR$individual_shifts$Dead[240]
resultskyr_mSTR$summary_resultsSUCC_t
resultskyr_CON$summary_resultsSUCC_t

resultskyr_mSTR$individual_shifts$Dead[240]/n_c
resultskyr_CON$individual_shifts$Dead[240]/n_c
resultskyr_mSTR$individual_shifts$Cured[240]/n_c
resultskyr_CON$individual_shifts$Cured[240]/n_c

############################################################
#----------------------------------------------------------#

#------------------------------------------------------------------------------#
#IV. PROBABILITY SENSITIVITY ANALYSIS.                                         #
#------------------------------------------------------------------------------#
#Definition of parameters and distributions, by country
set.seed(123)
#Kazakhstan_param_mSTR##########

params_psa_mSTR_kaz <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = 635.14 /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=658.4/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=635.14/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_kaz$p_mSTR_SLtreatment*100, 100-(paramsmSTR_kaz$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_kaz$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_kaz$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_kaz$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_kaz$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_kaz$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_kaz$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_kaz$p_mSTR_Dead*100), 100-(paramsmSTR_kaz$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_kaz$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_kaz$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_kaz$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_kaz$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_kaz$p_SLtreatment_Dead*100), 100-(paramsmSTR_kaz$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_kaz$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_kaz$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_kaz$p_LostFollowUp_Dead*100), 100-(paramsmSTR_kaz$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_kaz$p_Unresolved_Dead*100), 100-(paramsmSTR_kaz$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_kaz$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_kaz$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_kaz$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_kaz$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_kaz$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_kaz$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_kaz$p_AdverseEffects_Dead*100), 100-(paramsmSTR_kaz$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_kaz <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = 658.40 /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=658.4/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=658.4/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_kaz$qaly_mSTR, 
  qaly_SLtreatment = params_psa_mSTR_kaz$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_kaz$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_kaz$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_kaz$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_kaz$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_kaz$qaly_Cured,
  qaly_Dead = params_psa_mSTR_kaz$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_kaz$p_mSTR_SLtreatment*100, 100-(paramsCON_kaz$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_kaz$p_mSTR_LostFollowUp*100), 100-(paramsCON_kaz$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_kaz$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_kaz$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_kaz$p_mSTR_AdverseEffects*100), 100-(paramsCON_kaz$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_kaz$p_mSTR_Dead*100), 100-(paramsmSTR_kaz$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_kaz$p_SLtreatment_Unresolved*100), 100-(paramsCON_kaz$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_kaz$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_kaz$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_kaz$p_SLtreatment_Dead*100), 100-(paramsCON_kaz$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_kaz$p_LostFollowUp_mSTR*100), 100-(paramsCON_kaz$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_kaz$p_LostFollowUp_Dead*100), 100-(paramsCON_kaz$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_kaz$p_Unresolved_Dead*100), 100-(paramsCON_kaz$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_kaz$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_kaz$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_kaz$p_TreatmentCompleted_Cured*100), 100-(paramsCON_kaz$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_kaz$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_kaz$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_kaz$p_AdverseEffects_Dead*100), 100-(paramsCON_kaz$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)




#summary(params)


#Belarus_param_mSTR#######
params_psa_mSTR_bel <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_bel$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_bel$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_bel$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_bel$p_mSTR_SLtreatment*100, 100-(paramsmSTR_bel$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_bel$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_bel$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_bel$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_bel$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_bel$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_bel$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_bel$p_mSTR_Dead*100), 100-(paramsmSTR_bel$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_bel$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_bel$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_bel$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_bel$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_bel$p_SLtreatment_Dead*100), 100-(paramsmSTR_bel$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_bel$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_bel$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_bel$p_LostFollowUp_Dead*100), 100-(paramsmSTR_bel$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_bel$p_Unresolved_Dead*100), 100-(paramsmSTR_bel$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_bel$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_bel$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_bel$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_bel$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_bel$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_bel$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_bel$p_AdverseEffects_Dead*100), 100-(paramsmSTR_bel$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_bel <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_bel$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_bel$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_bel$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_bel$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_bel$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_bel$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_bel$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_bel$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_bel$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_bel$qaly_Cured,
  qaly_Dead = params_psa_mSTR_bel$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_bel$p_mSTR_SLtreatment*100, 100-(paramsCON_bel$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_bel$p_mSTR_LostFollowUp*100), 100-(paramsCON_bel$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_bel$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_bel$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_bel$p_mSTR_AdverseEffects*100), 100-(paramsCON_bel$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_bel$p_mSTR_Dead*100), 100-(paramsmSTR_bel$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_bel$p_SLtreatment_Unresolved*100), 100-(paramsCON_bel$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_bel$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_bel$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_bel$p_SLtreatment_Dead*100), 100-(paramsCON_bel$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_bel$p_LostFollowUp_mSTR*100), 100-(paramsCON_bel$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_bel$p_LostFollowUp_Dead*100), 100-(paramsCON_bel$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_bel$p_Unresolved_Dead*100), 100-(paramsCON_bel$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_bel$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_bel$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_bel$p_TreatmentCompleted_Cured*100), 100-(paramsCON_bel$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_bel$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_bel$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_bel$p_AdverseEffects_Dead*100), 100-(paramsCON_bel$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

#Moldova_param_mSTR#######
params_psa_mSTR_mol <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_mol$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_mol$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_mol$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_mol$p_mSTR_SLtreatment*100, 100-(paramsmSTR_mol$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_mol$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_mol$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_mol$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_mol$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_mol$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_mol$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_mol$p_mSTR_Dead*100), 100-(paramsmSTR_mol$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_mol$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_mol$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_mol$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_mol$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_mol$p_SLtreatment_Dead*100), 100-(paramsmSTR_mol$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_mol$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_mol$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_mol$p_LostFollowUp_Dead*100), 100-(paramsmSTR_mol$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_mol$p_Unresolved_Dead*100), 100-(paramsmSTR_mol$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_mol$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_mol$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_mol$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_mol$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_mol$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_mol$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_mol$p_AdverseEffects_Dead*100), 100-(paramsmSTR_mol$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_mol <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_mol$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_mol$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_mol$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_mol$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_mol$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_mol$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_mol$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_mol$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_mol$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_mol$qaly_Cured,
  qaly_Dead = params_psa_mSTR_mol$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_mol$p_mSTR_SLtreatment*100, 100-(paramsCON_mol$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_mol$p_mSTR_LostFollowUp*100), 100-(paramsCON_mol$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_mol$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_mol$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_mol$p_mSTR_AdverseEffects*100), 100-(paramsCON_mol$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_mol$p_mSTR_Dead*100), 100-(paramsmSTR_mol$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_mol$p_SLtreatment_Unresolved*100), 100-(paramsCON_mol$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_mol$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_mol$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_mol$p_SLtreatment_Dead*100), 100-(paramsCON_mol$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_mol$p_LostFollowUp_mSTR*100), 100-(paramsCON_mol$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_mol$p_LostFollowUp_Dead*100), 100-(paramsCON_mol$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_mol$p_Unresolved_Dead*100), 100-(paramsCON_mol$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_mol$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_mol$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_mol$p_TreatmentCompleted_Cured*100), 100-(paramsCON_mol$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_mol$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_mol$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_mol$p_AdverseEffects_Dead*100), 100-(paramsCON_mol$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)



#Georgia_param_mSTR#######
params_psa_mSTR_geo <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_geo$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_geo$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_geo$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_geo$p_mSTR_SLtreatment*100, 100-(paramsmSTR_geo$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_geo$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_geo$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_geo$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_geo$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_geo$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_geo$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_geo$p_mSTR_Dead*100), 100-(paramsmSTR_geo$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_geo$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_geo$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_geo$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_geo$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_geo$p_SLtreatment_Dead*100), 100-(paramsmSTR_geo$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_geo$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_geo$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_geo$p_LostFollowUp_Dead*100), 100-(paramsmSTR_geo$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_geo$p_Unresolved_Dead*100), 100-(paramsmSTR_geo$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_geo$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_geo$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_geo$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_geo$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_geo$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_geo$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_geo$p_AdverseEffects_Dead*100), 100-(paramsmSTR_geo$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_geo <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_geo$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_geo$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_geo$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_geo$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_geo$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_geo$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_geo$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_geo$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_geo$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_geo$qaly_Cured,
  qaly_Dead = params_psa_mSTR_geo$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_geo$p_mSTR_SLtreatment*100, 100-(paramsCON_geo$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_geo$p_mSTR_LostFollowUp*100), 100-(paramsCON_geo$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_geo$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_geo$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_geo$p_mSTR_AdverseEffects*100), 100-(paramsCON_geo$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_geo$p_mSTR_Dead*100), 100-(paramsmSTR_geo$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_geo$p_SLtreatment_Unresolved*100), 100-(paramsCON_geo$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_geo$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_geo$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_geo$p_SLtreatment_Dead*100), 100-(paramsCON_geo$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_geo$p_LostFollowUp_mSTR*100), 100-(paramsCON_geo$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_geo$p_LostFollowUp_Dead*100), 100-(paramsCON_geo$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_geo$p_Unresolved_Dead*100), 100-(paramsCON_geo$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_geo$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_geo$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_geo$p_TreatmentCompleted_Cured*100), 100-(paramsCON_geo$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_geo$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_geo$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_geo$p_AdverseEffects_Dead*100), 100-(paramsCON_geo$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)


#Other countries with NO data#####
#Azerbaijan_param_mSTR#######
params_psa_mSTR_aze <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_aze$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_aze$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_aze$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_aze$p_mSTR_SLtreatment*100, 100-(paramsmSTR_aze$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_aze$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_aze$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_aze$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_aze$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_aze$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_aze$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_aze$p_mSTR_Dead*100), 100-(paramsmSTR_aze$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_aze$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_aze$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_aze$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_aze$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_aze$p_SLtreatment_Dead*100), 100-(paramsmSTR_aze$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_aze$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_aze$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_aze$p_LostFollowUp_Dead*100), 100-(paramsmSTR_aze$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_aze$p_Unresolved_Dead*100), 100-(paramsmSTR_aze$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_aze$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_aze$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_aze$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_aze$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_aze$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_aze$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_aze$p_AdverseEffects_Dead*100), 100-(paramsmSTR_aze$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_aze <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_aze$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_aze$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_aze$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_aze$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_aze$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_aze$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_aze$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_aze$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_aze$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_aze$qaly_Cured,
  qaly_Dead = params_psa_mSTR_aze$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_aze$p_mSTR_SLtreatment*100, 100-(paramsCON_aze$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_aze$p_mSTR_LostFollowUp*100), 100-(paramsCON_aze$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_aze$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_aze$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_aze$p_mSTR_AdverseEffects*100), 100-(paramsCON_aze$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_aze$p_mSTR_Dead*100), 100-(paramsmSTR_aze$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_aze$p_SLtreatment_Unresolved*100), 100-(paramsCON_aze$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_aze$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_aze$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_aze$p_SLtreatment_Dead*100), 100-(paramsCON_aze$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_aze$p_LostFollowUp_mSTR*100), 100-(paramsCON_aze$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_aze$p_LostFollowUp_Dead*100), 100-(paramsCON_aze$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_aze$p_Unresolved_Dead*100), 100-(paramsCON_aze$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_aze$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_aze$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_aze$p_TreatmentCompleted_Cured*100), 100-(paramsCON_aze$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_aze$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_aze$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_aze$p_AdverseEffects_Dead*100), 100-(paramsCON_aze$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)


#Ukraine_param_mSTR#######
params_psa_mSTR_ukr <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_ukr$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_ukr$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_ukr$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_ukr$p_mSTR_SLtreatment*100, 100-(paramsmSTR_ukr$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_ukr$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_ukr$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_ukr$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_ukr$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_ukr$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_ukr$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_ukr$p_mSTR_Dead*100), 100-(paramsmSTR_ukr$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_ukr$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_ukr$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_ukr$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_ukr$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_ukr$p_SLtreatment_Dead*100), 100-(paramsmSTR_ukr$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_ukr$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_ukr$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_ukr$p_LostFollowUp_Dead*100), 100-(paramsmSTR_ukr$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_ukr$p_Unresolved_Dead*100), 100-(paramsmSTR_ukr$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_ukr$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_ukr$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_ukr$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_ukr$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_ukr$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_ukr$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_ukr$p_AdverseEffects_Dead*100), 100-(paramsmSTR_ukr$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_ukr <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_ukr$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_ukr$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_ukr$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_ukr$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_ukr$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_ukr$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_ukr$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_ukr$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_ukr$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_ukr$qaly_Cured,
  qaly_Dead = params_psa_mSTR_ukr$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_ukr$p_mSTR_SLtreatment*100, 100-(paramsCON_ukr$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_ukr$p_mSTR_LostFollowUp*100), 100-(paramsCON_ukr$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_ukr$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_ukr$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_ukr$p_mSTR_AdverseEffects*100), 100-(paramsCON_ukr$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_ukr$p_mSTR_Dead*100), 100-(paramsmSTR_ukr$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_ukr$p_SLtreatment_Unresolved*100), 100-(paramsCON_ukr$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_ukr$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_ukr$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_ukr$p_SLtreatment_Dead*100), 100-(paramsCON_ukr$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_ukr$p_LostFollowUp_mSTR*100), 100-(paramsCON_ukr$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_ukr$p_LostFollowUp_Dead*100), 100-(paramsCON_ukr$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_ukr$p_Unresolved_Dead*100), 100-(paramsCON_ukr$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_ukr$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_ukr$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_ukr$p_TreatmentCompleted_Cured*100), 100-(paramsCON_ukr$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_ukr$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_ukr$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_ukr$p_AdverseEffects_Dead*100), 100-(paramsCON_ukr$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

#Uzbekistan_param_mSTR#######
params_psa_mSTR_uzb <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_uzb$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_uzb$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_uzb$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_uzb$p_mSTR_SLtreatment*100, 100-(paramsmSTR_uzb$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_uzb$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_uzb$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_uzb$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_uzb$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_uzb$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_uzb$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_uzb$p_mSTR_Dead*100), 100-(paramsmSTR_uzb$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_uzb$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_uzb$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_uzb$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_uzb$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_uzb$p_SLtreatment_Dead*100), 100-(paramsmSTR_uzb$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_uzb$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_uzb$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_uzb$p_LostFollowUp_Dead*100), 100-(paramsmSTR_uzb$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_uzb$p_Unresolved_Dead*100), 100-(paramsmSTR_uzb$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_uzb$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_uzb$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_uzb$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_uzb$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_uzb$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_uzb$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_uzb$p_AdverseEffects_Dead*100), 100-(paramsmSTR_uzb$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_uzb <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_uzb$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_uzb$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_uzb$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_uzb$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_uzb$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_uzb$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_uzb$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_uzb$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_uzb$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_uzb$qaly_Cured,
  qaly_Dead = params_psa_mSTR_uzb$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_uzb$p_mSTR_SLtreatment*100, 100-(paramsCON_uzb$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_uzb$p_mSTR_LostFollowUp*100), 100-(paramsCON_uzb$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_uzb$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_uzb$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_uzb$p_mSTR_AdverseEffects*100), 100-(paramsCON_uzb$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_uzb$p_mSTR_Dead*100), 100-(paramsmSTR_uzb$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_uzb$p_SLtreatment_Unresolved*100), 100-(paramsCON_uzb$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_uzb$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_uzb$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_uzb$p_SLtreatment_Dead*100), 100-(paramsCON_uzb$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_uzb$p_LostFollowUp_mSTR*100), 100-(paramsCON_uzb$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_uzb$p_LostFollowUp_Dead*100), 100-(paramsCON_uzb$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_uzb$p_Unresolved_Dead*100), 100-(paramsCON_uzb$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_uzb$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_uzb$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_uzb$p_TreatmentCompleted_Cured*100), 100-(paramsCON_uzb$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_uzb$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_uzb$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_uzb$p_AdverseEffects_Dead*100), 100-(paramsCON_uzb$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

#Tajikistan_param_mSTR#######
params_psa_mSTR_taj <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_taj$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_taj$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_taj$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_taj$p_mSTR_SLtreatment*100, 100-(paramsmSTR_taj$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_taj$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_taj$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_taj$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_taj$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_taj$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_taj$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_taj$p_mSTR_Dead*100), 100-(paramsmSTR_taj$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_taj$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_taj$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_taj$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_taj$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_taj$p_SLtreatment_Dead*100), 100-(paramsmSTR_taj$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_taj$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_taj$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_taj$p_LostFollowUp_Dead*100), 100-(paramsmSTR_taj$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_taj$p_Unresolved_Dead*100), 100-(paramsmSTR_taj$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_taj$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_taj$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_taj$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_taj$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_taj$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_taj$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_taj$p_AdverseEffects_Dead*100), 100-(paramsmSTR_taj$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_taj <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_taj$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_taj$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_taj$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_taj$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_taj$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_taj$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_taj$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_taj$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_taj$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_taj$qaly_Cured,
  qaly_Dead = params_psa_mSTR_taj$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_taj$p_mSTR_SLtreatment*100, 100-(paramsCON_taj$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_taj$p_mSTR_LostFollowUp*100), 100-(paramsCON_taj$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_taj$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_taj$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_taj$p_mSTR_AdverseEffects*100), 100-(paramsCON_taj$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_taj$p_mSTR_Dead*100), 100-(paramsmSTR_taj$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_taj$p_SLtreatment_Unresolved*100), 100-(paramsCON_taj$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_taj$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_taj$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_taj$p_SLtreatment_Dead*100), 100-(paramsCON_taj$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_taj$p_LostFollowUp_mSTR*100), 100-(paramsCON_taj$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_taj$p_LostFollowUp_Dead*100), 100-(paramsCON_taj$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_taj$p_Unresolved_Dead*100), 100-(paramsCON_taj$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_taj$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_taj$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_taj$p_TreatmentCompleted_Cured*100), 100-(paramsCON_taj$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_taj$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_taj$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_taj$p_AdverseEffects_Dead*100), 100-(paramsCON_taj$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

#Turkmenistan_param_mSTR#######
params_psa_mSTR_tur <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_tur$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_tur$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_tur$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_tur$p_mSTR_SLtreatment*100, 100-(paramsmSTR_tur$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_tur$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_tur$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_tur$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_tur$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_tur$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_tur$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_tur$p_mSTR_Dead*100), 100-(paramsmSTR_tur$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_tur$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_tur$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_tur$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_tur$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_tur$p_SLtreatment_Dead*100), 100-(paramsmSTR_tur$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_tur$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_tur$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_tur$p_LostFollowUp_Dead*100), 100-(paramsmSTR_tur$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_tur$p_Unresolved_Dead*100), 100-(paramsmSTR_tur$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_tur$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_tur$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_tur$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_tur$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_tur$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_tur$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_tur$p_AdverseEffects_Dead*100), 100-(paramsmSTR_tur$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_tur <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_tur$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_tur$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_tur$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_tur$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_tur$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_tur$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_tur$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_tur$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_tur$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_tur$qaly_Cured,
  qaly_Dead = params_psa_mSTR_tur$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_tur$p_mSTR_SLtreatment*100, 100-(paramsCON_tur$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_tur$p_mSTR_LostFollowUp*100), 100-(paramsCON_tur$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_tur$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_tur$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_tur$p_mSTR_AdverseEffects*100), 100-(paramsCON_tur$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_tur$p_mSTR_Dead*100), 100-(paramsmSTR_tur$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_tur$p_SLtreatment_Unresolved*100), 100-(paramsCON_tur$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_tur$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_tur$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_tur$p_SLtreatment_Dead*100), 100-(paramsCON_tur$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_tur$p_LostFollowUp_mSTR*100), 100-(paramsCON_tur$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_tur$p_LostFollowUp_Dead*100), 100-(paramsCON_tur$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_tur$p_Unresolved_Dead*100), 100-(paramsCON_tur$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_tur$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_tur$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_tur$p_TreatmentCompleted_Cured*100), 100-(paramsCON_tur$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_tur$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_tur$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_tur$p_AdverseEffects_Dead*100), 100-(paramsCON_tur$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

#Armenia_param_mSTR#######

params_psa_mSTR_arm <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_arm$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_arm$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_arm$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_arm$p_mSTR_SLtreatment*100, 100-(paramsmSTR_arm$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_arm$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_arm$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_arm$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_arm$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_arm$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_arm$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_arm$p_mSTR_Dead*100), 100-(paramsmSTR_arm$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_arm$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_arm$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_arm$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_arm$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_arm$p_SLtreatment_Dead*100), 100-(paramsmSTR_arm$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_arm$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_arm$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_arm$p_LostFollowUp_Dead*100), 100-(paramsmSTR_arm$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_arm$p_Unresolved_Dead*100), 100-(paramsmSTR_arm$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_arm$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_arm$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_arm$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_arm$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_arm$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_arm$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_arm$p_AdverseEffects_Dead*100), 100-(paramsmSTR_arm$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_arm <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_arm$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_arm$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_arm$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_arm$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_arm$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_arm$qaly_LostFollowup,
  qaly_Unresolved = params_psa_mSTR_arm$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_arm$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_arm$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_arm$qaly_Cured,
  qaly_Dead = params_psa_mSTR_arm$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_arm$p_mSTR_SLtreatment*100, 100-(paramsCON_arm$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_arm$p_mSTR_LostFollowUp*100), 100-(paramsCON_arm$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_arm$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_arm$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_arm$p_mSTR_AdverseEffects*100), 100-(paramsCON_arm$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_arm$p_mSTR_Dead*100), 100-(paramsmSTR_arm$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_arm$p_SLtreatment_Unresolved*100), 100-(paramsCON_arm$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_arm$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_arm$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_arm$p_SLtreatment_Dead*100), 100-(paramsCON_arm$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_arm$p_LostFollowUp_mSTR*100), 100-(paramsCON_arm$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_arm$p_LostFollowUp_Dead*100), 100-(paramsCON_arm$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_arm$p_Unresolved_Dead*100), 100-(paramsCON_arm$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_arm$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_arm$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_arm$p_TreatmentCompleted_Cured*100), 100-(paramsCON_arm$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_arm$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_arm$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_arm$p_AdverseEffects_Dead*100), 100-(paramsCON_arm$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)



#Kyrgystan_param_mSTR#######
params_psa_mSTR_kyr <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsmSTR_kyr$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale= paramsmSTR_kyr$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale= paramsmSTR_kyr$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = rbeta(N_psa, 81, 19),
  qaly_SLtreatment = rbeta(N_psa, 79, 21),
  qaly_LostFollowUp = rbeta(N_psa, 68, 32),
  qaly_Unresolved = rbeta(N_psa, 68, 32),
  qaly_TreatmentCompleted = rbeta(N_psa, 81, 19),
  qaly_AdverseEffects = rbeta(N_psa, 68, 32),
  qaly_Cured = rbeta(N_psa, 81, 19),
  qaly_Dead = rbeta(N_psa, 0, 100),
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsmSTR_kyr$p_mSTR_SLtreatment*100, 100-(paramsmSTR_kyr$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsmSTR_kyr$p_mSTR_LostFollowUp*100), 100-(paramsmSTR_kyr$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_kyr$p_mSTR_TreatmentCompleted*100), 100-(paramsmSTR_kyr$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsmSTR_kyr$p_mSTR_AdverseEffects*100), 100-(paramsmSTR_kyr$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsmSTR_kyr$p_mSTR_Dead*100), 100-(paramsmSTR_kyr$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsmSTR_kyr$p_SLtreatment_Unresolved*100), 100-(paramsmSTR_kyr$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsmSTR_kyr$p_SLtreatment_TreatmentCompleted*100), 100-(paramsmSTR_kyr$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsmSTR_kyr$p_SLtreatment_Dead*100), 100-(paramsmSTR_kyr$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsmSTR_kyr$p_LostFollowUp_mSTR*100), 100-(paramsmSTR_kyr$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsmSTR_kyr$p_LostFollowUp_Dead*100), 100-(paramsmSTR_kyr$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsmSTR_kyr$p_Unresolved_Dead*100), 100-(paramsmSTR_kyr$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsmSTR_kyr$p_TreatmentCompleted_SLtreatment*100), 100-(paramsmSTR_kyr$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsmSTR_kyr$p_TreatmentCompleted_Cured*100), 100-(paramsmSTR_kyr$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsmSTR_kyr$p_AdverseEffects_SLtreatment*100), 100-(paramsmSTR_kyr$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsmSTR_kyr$p_AdverseEffects_Dead*100), 100-(paramsmSTR_kyr$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

params_psa_CON_kyr <- data.frame(
  #Costs
  Cost_mSTR = rgamma(N_psa, shape = 25, scale = paramsCON_kyr$Cost_mSTR /25),
  Cost_SLtreatment = rgamma(N_psa, shape=25, scale=paramsCON_kyr$Cost_SLtreatment/25),
  Cost_LostFollowUp =  rgamma(N_psa, shape=25, scale=1/25),
  Cost_Unresolved =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_TreatmentCompleted =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_AdverseEffects =  rgamma(N_psa, shape=25, scale=paramsCON_kyr$Cost_AdverseEffects/25),
  Cost_Cured =  rgamma(N_psa, shape=25, scale=0/25),
  Cost_Dead =  rgamma(N_psa, shape=25, scale=0/25),
  qaly_mSTR = params_psa_mSTR_kyr$qaly_mSTR,
  qaly_SLtreatment = params_psa_mSTR_kyr$qaly_SLtreatment,
  qaly_LostFollowUp = params_psa_mSTR_kyr$qaly_LostFollowUp,
  qaly_Unresolved = params_psa_mSTR_kyr$qaly_Unresolved,
  qaly_TreatmentCompleted = params_psa_mSTR_kyr$qaly_TreatmentCompleted,
  qaly_AdverseEffects = params_psa_mSTR_kyr$qaly_AdverseEffects,
  qaly_Cured = params_psa_mSTR_kyr$qaly_Cured,
  qaly_Dead = params_psa_mSTR_kyr$qaly_Dead,
  
  #Matrix of transition probabilities (add the corresponding probabilities from above but using the decimal number and then substracting it 100)
  p_mSTR_SLtreatment = rbeta(N_psa,paramsCON_kyr$p_mSTR_SLtreatment*100, 100-(paramsCON_kyr$p_mSTR_SLtreatment*100)),
  p_mSTR_LostFollowUp = rbeta(N_psa,(paramsCON_kyr$p_mSTR_LostFollowUp*100), 100-(paramsCON_kyr$p_mSTR_LostFollowUp*100)),
  p_mSTR_Unresolved = 0,
  p_mSTR_TreatmentCompleted = rbeta(N_psa,(paramsCON_kyr$p_mSTR_TreatmentCompleted*100), 100-(paramsCON_kyr$p_mSTR_TreatmentCompleted*100)),
  p_mSTR_AdverseEffects = rbeta(N_psa,(paramsCON_kyr$p_mSTR_AdverseEffects*100), 100-(paramsCON_kyr$p_mSTR_AdverseEffects*100)),
  p_mSTR_Cured = 0,
  p_mSTR_Dead =  rbeta(N_psa,(paramsCON_kyr$p_mSTR_Dead*100), 100-(paramsmSTR_kyr$p_mSTR_Dead*100)),
  
  p_SLtreatment_mSTR = 0,
  p_SLtreatment_LostFollowUp = 0,
  p_SLtreatment_Unresolved = rbeta(N_psa,(paramsCON_kyr$p_SLtreatment_Unresolved*100), 100-(paramsCON_kyr$p_SLtreatment_Unresolved*100)),
  p_SLtreatment_TreatmentCompleted = rbeta(N_psa,(paramsCON_kyr$p_SLtreatment_TreatmentCompleted*100), 100-(paramsCON_kyr$p_SLtreatment_TreatmentCompleted*100)),
  p_SLtreatment_AdverseEffects = 0,
  p_SLtreatment_Cured = 0,
  p_SLtreatment_Dead = rbeta(N_psa,(paramsCON_kyr$p_SLtreatment_Dead*100), 100-(paramsCON_kyr$p_SLtreatment_Dead*100)),  
  
  p_LostFollowUp_mSTR = rbeta(N_psa,(paramsCON_kyr$p_LostFollowUp_mSTR*100), 100-(paramsCON_kyr$p_LostFollowUp_mSTR*100)),
  p_LostFollowUp_SLtreatment = 0,
  p_LostFollowUp_Unresolved = 0,
  p_LostFollowUp_TreatmentCompleted = 0,
  p_LostFollowUp_AdverseEffects = 0,
  p_LostFollowUp_Cured = 0,
  p_LostFollowUp_Dead = rbeta(N_psa,(paramsCON_kyr$p_LostFollowUp_Dead*100), 100-(paramsCON_kyr$p_LostFollowUp_Dead*100)),
  
  p_Unresolved_mSTR = 0,
  p_Unresolved_SLtreatment = 0,
  p_Unresolved_LostFollowUp = 0,
  p_Unresolved_TreatmentCompleted = 0,
  p_Unresolved_AdverseEffects = 0,
  p_Unresolved_Cured = 0,
  p_Unresolved_Dead =  rbeta(N_psa,(paramsCON_kyr$p_Unresolved_Dead*100), 100-(paramsCON_kyr$p_Unresolved_Dead*100)),
  
  p_TreatmentCompleted_mSTR = 0,
  p_TreatmentCompleted_SLtreatment = rbeta(N_psa,(paramsCON_kyr$p_TreatmentCompleted_SLtreatment*100), 100-(paramsCON_kyr$p_TreatmentCompleted_SLtreatment*100)),
  p_TreatmentCompleted_LostFollowUp = 0,
  p_TreatmentCompleted_Unresolved =  0,
  p_TreatmentCompleted_AdverseEffects = 0, 
  p_TreatmentCompleted_Cured = rbeta(N_psa,(paramsCON_kyr$p_TreatmentCompleted_Cured*100), 100-(paramsCON_kyr$p_TreatmentCompleted_Cured*100)), 
  p_TreatmentCompleted_Dead = 0,  
  
  p_AdverseEffects_mSTR = 0,
  p_AdverseEffects_SLtreatment = rbeta(N_psa,(paramsCON_kyr$p_AdverseEffects_SLtreatment*100), 100-(paramsCON_kyr$p_AdverseEffects_SLtreatment*100)),
  p_AdverseEffects_LostFollowUp = 0,
  p_AdverseEffects_Unresolved =  0,
  p_AdverseEffects_TreatmentCompleted = 0,
  p_AdverseEffects_Cured = 0,
  p_AdverseEffects_Dead = rbeta(N_psa,(paramsCON_kyr$p_AdverseEffects_Dead*100), 100-(paramsCON_kyr$p_AdverseEffects_Dead*100)), 
  
  p_Cured_mSTR = 0,
  p_Cured_SLtreatment = 0,
  p_Cured_LostFollowUp = 0,
  p_Cured_Unresolved = 0 ,
  p_Cured_TreatmentCompleted = 0,
  p_Cured_AdverseEffects= 0,
  p_Cured_Dead =  0
  
)

##########

#------------------------------------------------------------------------------#
#IV.I  MODEL FUNCTION for probability sensitivity analyses (PSA)
modelpsa <- function(.params) { 
  with(.params, {
    n_t <- 240
    n_s <- 8
    #n_c <- 1000
    v_state_names <- c("mSTR","SLtreatment","LostFollowUp","Unresolved","TreatmentCompleted","AdverseEffects","Cured", "Dead")
    m_p <- matrix(0, nrow= n_states, ncol=n_states, byrow= TRUE, dimnames= list(from= v_state_names, to = v_state_names))
    
    #prob matrix   
    m_p["mSTR", "mSTR"] <- 1  -p_mSTR_AdverseEffects -p_mSTR_SLtreatment -p_mSTR_LostFollowUp - p_mSTR_TreatmentCompleted -p_mSTR_Dead
    m_p["mSTR", "SLtreatment"] <- p_mSTR_SLtreatment
    m_p["mSTR", "LostFollowUp"] <- p_mSTR_LostFollowUp
    m_p["mSTR", "Unresolved"] <- p_mSTR_Unresolved
    m_p["mSTR", "TreatmentCompleted"] <- p_mSTR_TreatmentCompleted 
    m_p["mSTR", "AdverseEffects"] <- p_mSTR_AdverseEffects
    m_p["mSTR", "Cured"] <- p_mSTR_Cured
    m_p["mSTR", "Dead"] <-  p_mSTR_Dead
    
    m_p["SLtreatment", "mSTR"] <- p_SLtreatment_mSTR
    m_p["SLtreatment", "SLtreatment"] <- 1 -p_SLtreatment_Unresolved -p_SLtreatment_TreatmentCompleted -( p_SLtreatment_Dead )
    m_p["SLtreatment", "LostFollowUp"] <- p_SLtreatment_LostFollowUp
    m_p["SLtreatment", "Unresolved"] <- p_SLtreatment_Unresolved
    m_p["SLtreatment", "TreatmentCompleted"] <- p_SLtreatment_TreatmentCompleted
    m_p["SLtreatment", "AdverseEffects"] <- p_SLtreatment_AdverseEffects
    m_p["SLtreatment", "Cured"] <- p_SLtreatment_Cured
    m_p["SLtreatment", "Dead"] <-  p_SLtreatment_Dead 
    
    m_p["LostFollowUp", "mSTR"] <- p_LostFollowUp_mSTR
    m_p["LostFollowUp", "SLtreatment"] <- 0
    m_p["LostFollowUp", "LostFollowUp"] <- 1 - p_LostFollowUp_mSTR -( p_LostFollowUp_Dead)
    m_p["LostFollowUp", "Unresolved"] <- 0
    m_p["LostFollowUp", "TreatmentCompleted"] <-0
    m_p["LostFollowUp", "AdverseEffects"] <-  0
    m_p["LostFollowUp", "Cured"] <- 0
    m_p["LostFollowUp", "Dead"] <-  p_LostFollowUp_Dead
    
    m_p["Unresolved", "mSTR"] <- 0
    m_p["Unresolved", "SLtreatment"] <- 0
    m_p["Unresolved", "LostFollowUp"] <- 0  
    m_p["Unresolved", "Unresolved"] <- 1 - (p_Unresolved_Dead)
    m_p["Unresolved", "TreatmentCompleted"] <- 0
    m_p["Unresolved", "AdverseEffects"] <- 0
    m_p["Unresolved", "Cured"] <- 0
    m_p["Unresolved", "Dead"] <-  (0.75)*p_Unresolved_Dead
    
    m_p["TreatmentCompleted", "mSTR"] <- 0
    m_p["TreatmentCompleted", "SLtreatment"] <- p_TreatmentCompleted_SLtreatment
    m_p["TreatmentCompleted", "LostFollowUp"] <- 0
    m_p["TreatmentCompleted", "Unresolved"] <-  0
    m_p["TreatmentCompleted", "TreatmentCompleted"] <-  1- p_TreatmentCompleted_SLtreatment -p_TreatmentCompleted_Cured -(p_TreatmentCompleted_Dead)
    m_p["TreatmentCompleted", "AdverseEffects"] <-  0
    m_p["TreatmentCompleted", "Cured"] <-  p_TreatmentCompleted_Cured
    m_p["TreatmentCompleted", "Dead"] <-   p_TreatmentCompleted_Dead
    
    m_p["AdverseEffects", "mSTR"] <- 0
    m_p["AdverseEffects", "SLtreatment"] <- p_AdverseEffects_SLtreatment
    m_p["AdverseEffects", "LostFollowUp"] <- 0
    m_p["AdverseEffects", "Unresolved"] <-  p_AdverseEffects_Unresolved
    m_p["AdverseEffects", "TreatmentCompleted"] <-  0
    m_p["AdverseEffects", "AdverseEffects"] <-  1- p_AdverseEffects_Unresolved - p_AdverseEffects_SLtreatment -(p_AdverseEffects_Dead)
    m_p["AdverseEffects", "Cured"] <-  0
    m_p["AdverseEffects", "Dead"] <-   p_AdverseEffects_Dead
    
    m_p["Cured", "mSTR"] <- 0
    m_p["Cured", "SLtreatment"] <- 0
    m_p["Cured", "LostFollowUp"] <- 0
    m_p["Cured", "Unresolved"] <- 0
    m_p["Cured", "TreatmentCompleted"] <- 0
    m_p["Cured", "AdverseEffects"] <- 0
    m_p["Cured", "Cured"] <-  1 - (p_Cured_Dead)
    m_p["Cured", "Dead"] <-   p_Cured_Dead
    
    m_p["Dead", "mSTR"] <- 0
    m_p["Dead", "SLtreatment"] <- 0
    m_p["Dead", "LostFollowUp"] <- 0
    m_p["Dead", "Unresolved"] <- 0
    m_p["Dead", "TreatmentCompleted"] <- 0
    m_p["Dead", "AdverseEffects"] <- 0
    m_p["Dead", "Cured"] <- 0
    m_p["Dead", "Dead"] <-  1
    
    
    #Start drafting the markov matrix
    state_membership <- array(NA_real_, dim= c(n_t, n_s), dimnames= list(cycle = 1:n_t, state = v_state_names))
    view(state_membership)
    #Start filling in first row; only with treatment box.
    state_membership[1,]= c(n_c, 0, 0, 0, 0, 0, 0, 0)
    for (i in 2:n_t) {state_membership[i, ] <- state_membership[i-1, ] %*% m_p}
    
    #Insert costs per healthy, diseased and death, let's wait for Tom
    m_payoffs <- matrix(0, nrow = 8, ncol = 2, dimnames= list(state= v_state_names, 
                                                              payoff= c("Costs", "QALYs")))
    
    m_payoffsCost<- matrix(0, nrow=8, ncol=1,
                           dimnames= list(state= v_state_names, 
                                          payoff= c("Costs")))
    m_payoffsQALY <- matrix(0, nrow=8, ncol=1,
                            dimnames= list(state= v_state_names, 
                                           payoff= c("QALYs")))
    
    m_payoffs["mSTR",  "Costs"] <-  Cost_mSTR
    m_payoffs["SLtreatment", "Costs"] <-  Cost_SLtreatment
    m_payoffs["LostFollowUp", "Costs"] <- Cost_LostFollowUp
    m_payoffs["Unresolved", "Costs"] <- Cost_Unresolved
    m_payoffs["TreatmentCompleted", "Costs"] <- Cost_TreatmentCompleted
    m_payoffs["AdverseEffects", "Costs"] <-  Cost_AdverseEffects
    m_payoffs["Cured", "Costs"] <- Cost_Cured
    m_payoffs["Dead", "Costs"] <- Cost_Dead
    
    m_payoffsCost["mSTR",  "Costs"] <-  Cost_mSTR
    m_payoffsCost["SLtreatment", "Costs"] <-  Cost_SLtreatment
    m_payoffsCost["LostFollowUp", "Costs"] <- Cost_LostFollowUp
    m_payoffsCost["Unresolved", "Costs"] <- Cost_Unresolved
    m_payoffsCost["TreatmentCompleted", "Costs"] <- Cost_TreatmentCompleted
    m_payoffsCost["AdverseEffects", "Costs"] <-  Cost_AdverseEffects
    m_payoffsCost["Cured", "Costs"] <- Cost_Cured
    m_payoffsCost["Dead", "Costs"] <- Cost_Dead
    
    m_payoffs["mSTR",  "QALYs"] <-  qaly_mSTR/12
    m_payoffs["SLtreatment", "QALYs"] <-  qaly_SLtreatment/12
    m_payoffs["LostFollowUp", "QALYs"] <- qaly_LostFollowUp/12
    m_payoffs["Unresolved", "QALYs"] <- qaly_Unresolved/12
    m_payoffs["TreatmentCompleted", "QALYs"] <- qaly_TreatmentCompleted/12
    m_payoffs["AdverseEffects", "QALYs"] <-  qaly_AdverseEffects/12
    m_payoffs["Cured", "QALYs"] <- qaly_Cured/12
    m_payoffs["Dead", "QALYs"] <- qaly_Dead/12
    
    m_payoffsQALY["mSTR",  "QALYs"] <-  qaly_mSTR/12
    m_payoffsQALY["SLtreatment", "QALYs"] <-  qaly_SLtreatment/12
    m_payoffsQALY["LostFollowUp", "QALYs"] <- qaly_LostFollowUp/12
    m_payoffsQALY["Unresolved", "QALYs"] <- qaly_Unresolved/12
    m_payoffsQALY["TreatmentCompleted", "QALYs"] <- qaly_TreatmentCompleted/12
    m_payoffsQALY["AdverseEffects", "QALYs"] <-  qaly_AdverseEffects/12
    m_payoffsQALY["Cured", "QALYs"] <- qaly_Cured/12
    m_payoffsQALY["Dead", "QALYs"] <- qaly_Dead/12
    
    #Check probability matrix by typing : m_payoffs
    payoff_trace <- state_membership %*% m_payoffs
    payoff_trace_d <- matrix(0, nrow = 240, ncol = 2, dimnames= list(state= 1:240, payoff= c("Costs", "QALYs")))
    for (i in 1:n_t) {payoff_trace_d[i, ] <- payoff_trace[i, ] * (1/((1+discount_m)^(i)))}
    
    payoff_traceCost <- state_membership %*% m_payoffsCost
    payoff_traceCost_d <- matrix(0, nrow = 240, ncol = 1, dimnames= list(state= 1:240, payoff= c("Costs")))
    for (i in 1:n_t) {payoff_traceCost_d[i, ] <- payoff_traceCost[i, ] * (1/((1+discount_m)^(i)))}
    
    payoff_traceQaly <- state_membership %*% m_payoffsQALY
    payoff_traceQaly_d <- matrix(0, nrow = 240, ncol = 1, dimnames= list(state= 1:240, payoff= c("QALYs")))
    for (i in 1:n_t) {payoff_traceQaly_d[i, ] <- payoff_traceQaly[i, ] * (1/((1+discount_m)^(i)))}
  
    #Outputs, discount rate applied.

    #colSums(payoff_trace) / n_c
    
    summary_resultsAll_psaD = colSums(payoff_trace_d) / n_c
    #summary_resultsAll_psaDcost = colSums(payoff_traceCost_d) / n_c
    #summary_resultsAll_psaDqaly = colSums(payoff_traceQaly_d) / n_c
    
  })
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#IV.II RESULTS after applying PSA                                                #
#------------------------------------------------------------------------------#
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/Bloomsbury_PolicyLab/0_Article/Figures")
#(1):KAZAHSTAN##########
n_c<- 3755
psa_results_mSTR_kaz <- t(sapply(
  X = split(params_psa_mSTR_kaz, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
psa_results_CON_kaz <- t(sapply(
  X = split(params_psa_CON_kaz, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
#Calculations

theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="kazahstan_psa.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(psa_results_mSTR_kaz[,2], psa_results_mSTR_kaz[,1], type='p', cex.lab = 1,  cex.axis = 0.8, main ="mSTR treatment, Kazakhstan", 
  xlab = "QALYs gained",
  ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(6, 14), ylim=c(0, 20000),
  col = transparent("salmon", .5),
  #bg = "deepskyblue",   # Fill color
  #col = "blue", # Border color
  cex = 1,      # Symbol size
  lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_mSTR_kaz[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_mSTR_kaz[,1])),col="dark green", lty=8)
points(mean(psa_results_mSTR_kaz[,2]), mean(psa_results_mSTR_kaz[,1]), col="dark green", cex=2)

plot(psa_results_CON_kaz[,2], psa_results_CON_kaz[,1], type='p', cex.lab = 1,  cex.axis = 0.8,  main ="Conventional treatment, Kazakhstan",
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(3, 14), ylim=c(0, 60000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_CON_kaz[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_CON_kaz[,1])),col="dark green", lty=8)
points(mean(psa_results_CON_kaz[,2]), mean(psa_results_CON_kaz[,1]), col="dark green", cex=2)
dev.off()

#GRAPH ICERS
ICER_psa_Kaz <- (psa_results_mSTR_kaz[,1] - psa_results_CON_kaz[,1])/( psa_results_mSTR_kaz[,2] -psa_results_CON_kaz[,2])
ICER_psa_kaz1 <- psa_results_mSTR_kaz[,1] - psa_results_CON_kaz[,1]
ICER_psa_kaz2 <- psa_results_mSTR_kaz[,2] -psa_results_CON_kaz[,2]
#view(ICER_psa_Kaz)

tiff(file="kazahstan_psaICER.tiff",
     width=10, height=7, units="in", res=400)
par(mfrow=c(1,1))
a31 <- plot(ICER_psa_kaz1, ICER_psa_kaz2, pch = 21, main= "Incremental outcomes/costs between mSTR and conevntional treatment for MDR/RR TB, Kazahstan",
     bg = "bisque",   # Fill color
     col = "darkgoldenrod1", # Border color
     cex = 2,      # Symbol size
     lwd = 3, 
     type='p', cex.lab = 1,  cex.axis = 0.8,
     ylab = "Incremental QALYs gained",
     xlab = "Incremental economic costs in 2022 USD", las = 1, xlim=c(-50000, 20000), ylim=c(-5, 10))      
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_kaz1),mean(ICER_psa_kaz2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-36000, -40, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
dev.off()


#Graph2;effectiveness as percentage depending on scenario(4 scenarios)
name <- list("+QALYs, - Costs", "-QALYs, -Costs", "+QALYs, +Costs", "-QALYs, +Costs")
value <- list(sum(ICER_psa_kaz1<0.01 & ICER_psa_kaz2>0.01), sum(ICER_psa_kaz1<0.01 & ICER_psa_kaz2<0.01), sum(ICER_psa_kaz1>0.01 & ICER_psa_kaz2>0.01), sum(ICER_psa_kaz1>0.01 & ICER_psa_kaz2<0.01))
datakaznodf <- cbind( name, value)
datacekaz <- data.frame(datakaznodf)
datacekaz[1,3]= 866/1000*100
datacekaz[2,3]= 107/1000*100
datacekaz[3,3]= 22/1000*100
datacekaz[4,3]= 5/1000*100
# Uniform color GRAPH!
tiff(file="kazahstan_psaICER_varA.tiff",
     width=10, height=9, units="in", res=400)
par(mfrow=c(1,1))
barplot(height=datacekaz$V3, names=datacekaz$name, 
        col= "bisque1", cex.axis=1.2, cex.names=1.4,
        horiz=T, las=1, yaxt="n", xlab = "Percentage", ylab = "", main= "Kazahstan", xlim=c(0, 100))
text(60,1:4,name,cex=1.5)
mtext("Variation in Costs and QALYs values between mSTR and conventional treatment", side=2, line=2.2, cex=1.3)
dev.off()
#Graph3; effectiveness using WTP-thresholds
effective_ICER_kaz <- matrix(0, nrow= 1000, ncol=1)
wtp_kaz <- matrix(0, nrow= 1000, ncol=1)
wtp_kaz[1,1] <-0
effective_ICER_kaz[1,1] <- (sum(ICER_psa_Kaz < 0, na.rm=TRUE))/1000
for (i in 2:1000) {
effective_ICER_kaz[i,1] <- (sum(ICER_psa_Kaz < (i-1)*10, na.rm=TRUE))/1000
wtp_kaz[i,1] <- (i-1)*10
}
# Create DataFrame for Plotting
DF <- data.frame(wtp_kaz, effective_ICER_kaz)
tiff(file="kazahstan_psaICER_varB.tiff",
     width=10, height=8, units="in", res=400)
par(mfrow=c(1,1))
ggplot(DF, aes(wtp_kaz, effective_ICER_kaz)) +                                     
  geom_line(color = "dodgerblue4", size = 2) +
  geom_ribbon(aes(ymin=effective_ICER_kaz+0.02, ymax=effective_ICER_kaz-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="Kazahstan", ylim=c(0.75, 1), x = "Willingness-to-pay (WTP) thresholds in $", y ="Percentage of simulations mSTR regimen is cost-effective (ICER<WTP)")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 10080), breaks=c(0, 1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)) +
  scale_y_continuous(limits=c(.75, 0.92), breaks=c(0.75,0.8, 0.85, 0.9, 0.95), labels = c("75","80", "85", "90", "95")) +
  geom_vline(xintercept = 10073, linetype="dotted", color = "blue", size=1.5)
dev.off()



#(2):Belarus ##########
n_c<- 801
psa_results_mSTR_bel <- t(sapply(
  X = split(params_psa_mSTR_bel, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
psa_results_CON_bel <- t(sapply(
  X = split(params_psa_CON_bel, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
#Calculations

theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="belarus_psa.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(psa_results_mSTR_bel[,2], psa_results_mSTR_bel[,1], type='p', cex.lab = 1,  cex.axis = 0.8, main ="mSTR treatment, Belarus", 
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(6, 14), ylim=c(5000, 30000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_mSTR_bel[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_mSTR_bel[,1])),col="dark green", lty=8)
points(mean(psa_results_mSTR_bel[,2]), mean(psa_results_mSTR_bel[,1]), col="dark green", cex=2)

plot(psa_results_CON_bel[,2], psa_results_CON_bel[,1], type='p', cex.lab = 1,  cex.axis = 0.8,  main ="Conventional treatment, Belarus",
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(3, 14), ylim=c(10000, 70000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_CON_bel[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_CON_bel[,1])),col="dark green", lty=8)
points(mean(psa_results_CON_bel[,2]), mean(psa_results_CON_bel[,1]), col="dark green", cex=2)
dev.off()

#GRAPH ICERS
ICER_psa_Bel <- (psa_results_mSTR_bel[,1] - psa_results_CON_bel[,1])/( psa_results_mSTR_bel[,2] -psa_results_CON_bel[,2])
ICER_psa_bel1 <- psa_results_mSTR_bel[,1] - psa_results_CON_bel[,1]
ICER_psa_bel2 <- psa_results_mSTR_bel[,2] -psa_results_CON_bel[,2]
#view(ICER_psa_Bel)

tiff(file="belarus_psaICER.tiff",
     width=10, height=7, units="in", res=400)
par(mfrow=c(1,1))
plot(ICER_psa_bel1, ICER_psa_bel2, pch = 21, main= "Incremental outcomes/costs between mSTR and conevntional treatment for MDR/RR TB, Belarus",
     bg = "bisque",   # Fill color
     col = "darkgoldenrod1", # Border color
     cex = 2,      # Symbol size
     lwd = 3, 
     type='p', cex.lab = 1,  cex.axis = 0.8,
     ylab = "Incremental QALYs gained",
     xlab = "Incremental economic costs in 2022 USD", las = 1, xlim=c(-60000, 20000), ylim=c(-5, 8))      
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_bel1),mean(ICER_psa_bel2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-36000, -40, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
dev.off()


#Graph2;effectiveness as percentage depending on scenario(4 scenarios)
name <- list("+QALYs, - Costs", "-QALYs, -Costs", "+QALYs, +Costs", "-QALYs, +Costs")
value <- list(sum(ICER_psa_bel1<0.01 & ICER_psa_bel2>0.01), sum(ICER_psa_bel1<0.01 & ICER_psa_bel2<0.01), sum(ICER_psa_bel1>0.01 & ICER_psa_bel2>0.01), sum(ICER_psa_bel1>0.01 & ICER_psa_bel2<0.01))
databelnodf <- cbind( name, value)
datacebel <- data.frame(databelnodf)
#Use values from var values below!
datacebel[1,3]= 636/1000*100
datacebel[2,3]= 298/1000*100
datacebel[3,3]= 41/1000*100
datacebel[4,3]= 25/1000*100
# Uniform color GRAPH!
tiff(file="belarus_psaICER_varA.tiff",
     width=10, height=9, units="in", res=400)
par(mfrow=c(1,1))
barplot(height=datacebel$V3, names=datacebel$name, 
        col= "bisque1", cex.axis=1.2, cex.names=1.4,
        horiz=T, las=1, yaxt="n", xlab = "Percentage", ylab = "", main= "belarus", xlim=c(0, 100))
text(60,1:4,name,cex=1.5)
mtext("Variation in Costs and QALYs values between mSTR and conventional treatment", side=2, line=2.2, cex=1.3)
dev.off()
#Graph3; effectiveness using WTP-thresholds
effective_ICER_bel <- matrix(0, nrow= 1000, ncol=1)
wtp_bel <- matrix(0, nrow= 1000, ncol=1)
wtp_bel[1,1] <-0
effective_ICER_bel[1,1] <- (sum(ICER_psa_Bel < 0, na.rm=TRUE))/1000
for (i in 2:1000) {
  effective_ICER_bel[i,1] <- (sum(ICER_psa_Bel < (i-1)*10, na.rm=TRUE))/1000
  wtp_bel[i,1] <- (i-1)*10
}
# Create DataFrame for Plotting
DF <- data.frame(wtp_bel, effective_ICER_bel)
tiff(file="belarus_psaICER_varB.tiff",
     width=10, height=8, units="in", res=400)
par(mfrow=c(1,1))
ggplot(DF, aes(wtp_bel, effective_ICER_bel)) +                                     
  geom_line(color = "dodgerblue4", linewidth = 2) +
  geom_ribbon(aes(ymin=effective_ICER_bel+0.02, ymax=effective_ICER_bel-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="Belarus", ylim=c(0.80, 1), x = "Willingness-to-pay (WTP) thresholds in $", y ="Percentage of simulations mSTR regimen is cost-effective (ICER<WTP)")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 10000), breaks=c(0, 1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)) +
  scale_y_continuous(limits=c(0.45, 0.82), breaks=c(0.5,0.6, 0.65, 0.70, 0.75, 0.8), labels = c("50","60", "65", "70", "75", "80")) +
  geom_vline(xintercept = 6102, linetype="dotted", color = "blue", size=1.5)
dev.off()


#(3):Moldova ##########
n_c<- 593
psa_results_mSTR_mol <- t(sapply(
  X = split(params_psa_mSTR_mol, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
psa_results_CON_mol <- t(sapply(
  X = split(params_psa_CON_mol, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
#Calculations

theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="Moldova_psa.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(psa_results_mSTR_mol[,2], psa_results_mSTR_mol[,1], type='p', cex.lab = 1,  cex.axis = 0.8, main ="mSTR treatment, Moldova", 
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(7, 13), ylim=c(3000, 30000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_mSTR_mol[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_mSTR_mol[,1])),col="dark green", lty=8)
points(mean(psa_results_mSTR_mol[,2]), mean(psa_results_mSTR_mol[,1]), col="dark green", cex=2)

mean(psa_results_CON_mol[,2])
mean(psa_results_CON_mol[,1])
plot(psa_results_CON_mol[,2], psa_results_CON_mol[,1], type='p', cex.lab = 1,  cex.axis = 0.8,  main ="Conventional treatment, Moldova",
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(3, 12), ylim=c(5000, 35000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_CON_mol[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_CON_mol[,1])),col="dark green", lty=8)
points(mean(psa_results_CON_mol[,2]), mean(psa_results_CON_mol[,1]), col="dark green", cex=2)
dev.off()

#GRAPH ICERS
ICER_psa_mol <- (psa_results_mSTR_mol[,1] - psa_results_CON_mol[,1])/( psa_results_mSTR_mol[,2] -psa_results_CON_mol[,2])
ICER_psa_mol1 <- psa_results_mSTR_mol[,1] - psa_results_CON_mol[,1]
ICER_psa_mol2 <- psa_results_mSTR_mol[,2] -psa_results_CON_mol[,2]
#view(ICER_psa_Kaz)

tiff(file="Moldova_psaICER.tiff",
     width=10, height=7, units="in", res=400)
par(mfrow=c(1,1))
plot(ICER_psa_mol1, ICER_psa_mol2, pch = 21, main= "Incremental outcomes/costs between mSTR and conevntional treatment for MDR/RR TB, Moldova",
     bg = "bisque",   # Fill color
     col = "darkgoldenrod1", # Border color
     cex = 2,      # Symbol size
     lwd = 3, 
     type='p', cex.lab = 1,  cex.axis = 0.8,
     ylab = "Incremental QALYs gained",
     xlab = "Incremental economic costs in 2022 USDs", las = 1, xlim=c(-40000, 15000), ylim=c(-3, 8))      
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_mol1),mean(ICER_psa_mol2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-36000, -40, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
dev.off()

#mean(ICER_psa_mol1)
#mean(ICER_psa_mol2)

#Graph2;effectiveness as percentage depending on scenario(4 scenarios)
name <- list("+QALYs, - Costs", "-QALYs, -Costs", "+QALYs, +Costs", "-QALYs, +Costs")
value <- list(sum(ICER_psa_mol1<0.01 & ICER_psa_mol2>0.01), sum(ICER_psa_mol1<0.01 & ICER_psa_mol2<0.01), sum(ICER_psa_mol1>0.01 & ICER_psa_mol2>0.01), sum(ICER_psa_mol1>0.01 & ICER_psa_mol2<0.01))
datamolnodf <- cbind( name, value)
datacemol <- data.frame(datamolnodf)
#CHANGE THE NUMbERS bELOW DEPending ON COuNTry!
datacemol[1,3]= 854/1000*100
datacemol[2,3]= 60/1000*100
datacemol[3,3]= 84/1000*100
datacemol[4,3]= 2/1000*100
# Uniform color GRAPH!
tiff(file="Moldova_psaICER_varA.tiff",
     width=10, height=9, units="in", res=400)
par(mfrow=c(1,1))
barplot(height=datacemol$V3, names=datacemol$name, 
        col= "bisque1", cex.axis=1.2, cex.names=1.4,
        horiz=T, las=1, yaxt="n", xlab = "Percentage", ylab = "", main= "Moldova", xlim=c(0, 100))
text(60,1:4,name,cex=1.5)
mtext("Variation in Costs and QALYs values between mSTR and conventional treatment", side=2, line=2.2, cex=1.3)
dev.off()
#Graph3; effectiveness using WTP-thresholds
effective_ICER_mol <- matrix(0, nrow= 1000, ncol=1)
wtp_mol <- matrix(0, nrow= 1000, ncol=1)
wtp_mol[1,1] <-0
effective_ICER_mol[1,1] <- (sum(ICER_psa_mol < 0, na.rm=TRUE))/1000
for (i in 2:1000) {
  effective_ICER_mol[i,1] <- (sum(ICER_psa_mol < (i-1)*10, na.rm=TRUE))/1000
  wtp_mol[i,1] <- (i-1)*10
}
# Create DataFrame for Plotting
DF <- data.frame(wtp_mol, effective_ICER_mol)
tiff(file="Moldova_psaICER_varB.tiff",
     width=10, height=8, units="in", res=400)
par(mfrow=c(1,1))
ggplot(DF, aes(wtp_mol, effective_ICER_mol)) +                                     
  geom_line(color = "dodgerblue4", size = 2) +
  geom_ribbon(aes(ymin=effective_ICER_mol+0.02, ymax=effective_ICER_mol-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="Moldova", ylim=c(0.80, 1), x = "Willingness-to-pay (WTP) thresholds in $", y ="Percentage of simulations mSTR regimen is cost-effective (ICER<WTP)")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 3000), breaks=c(0, 250, 500,750,1000, 1250, 1500,1750,2000,2250,2500, 2750,3000)) +
  scale_y_continuous(limits=c(.68, 0.92), breaks=c(0.7,0.75,0.80, 0.85, 0.9), labels = c("70","75","80","85", "90")) +
  geom_vline(xintercept = 1446, linetype="dotted", color = "blue", size=1.5)
dev.off()




#(4):Georgia ##########
n_c<- 187
psa_results_mSTR_geo <- t(sapply(
  X = split(params_psa_mSTR_geo, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
psa_results_CON_geo <- t(sapply(
  X = split(params_psa_CON_geo, 1:N_psa),
  FUN = modelpsa,
  simplify = TRUE
))
#Calculations
theme(legend.spacing.x = unit(0.3, 'cm'))
tiff(file="Georgia_psa.tiff",
     width=6.5, height=9, units="in", res=400)
par(mfrow=c(2,1))
plot(psa_results_mSTR_geo[,2], psa_results_mSTR_geo[,1], type='p', cex.lab = 1,  cex.axis = 0.8, main ="mSTR treatment, Georgia", 
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(6, 14), ylim=c(2500, 20000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_mSTR_geo[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_mSTR_geo[,1])),col="dark green", lty=8)
points(mean(psa_results_mSTR_geo[,2]), mean(psa_results_mSTR_geo[,1]), col="dark green", cex=2)

plot(psa_results_CON_geo[,2], psa_results_CON_geo[,1], type='p', cex.lab = 1,  cex.axis = 0.8,  main ="Conventional treatment, Georgia",
     xlab = "QALYs gained",
     ylab = "Economic costs in 2022 USDs", las = 1, xlim=c(4, 14), ylim=c(5000, 50000),
     col = transparent("salmon", .5),
     #bg = "deepskyblue",   # Fill color
     #col = "blue", # Border color
     cex = 1,      # Symbol size
     lwd = 1, pch=19)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
abline(v=c(mean(psa_results_CON_geo[,2])),col="dark green", lty=8)
abline(h=c(mean(psa_results_CON_geo[,1])),col="dark green", lty=8)
points(mean(psa_results_CON_geo[,2]), mean(psa_results_CON_geo[,1]), col="dark green", cex=2)
dev.off()

#GRAPH ICERS
ICER_psa_Geo <- (psa_results_mSTR_geo[,1] - psa_results_CON_geo[,1])/( psa_results_mSTR_geo[,2] -psa_results_CON_geo[,2])
ICER_psa_geo1 <- psa_results_mSTR_geo[,1] - psa_results_CON_geo[,1]
ICER_psa_geo2 <- psa_results_mSTR_geo[,2] -psa_results_CON_geo[,2]
#view(ICER_psa_Geo)

tiff(file="Georgia_psaICER.tiff",
     width=10, height=7, units="in", res=400)
par(mfrow=c(1,1))
a31 <- plot(ICER_psa_geo1, ICER_psa_geo2, pch = 21, main= "Incremental outcomes/costs between mSTR and conevntional treatment for MDR/RR TB, Georgia",
            bg = "bisque",   # Fill color
            col = "darkgoldenrod1", # Border color
            cex = 2,      # Symbol size
            lwd = 3, 
            type='p', cex.lab = 1,  cex.axis = 0.8,
            ylab = "Incremental QALYs gained",
            xlab = "Incremental economic costs in 2022 USD", las = 1, xlim=c(-50000, 20000), ylim=c(-5, 10))      
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_geo1),mean(ICER_psa_geo2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-36000, -40, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
dev.off()


#Graph2;effectiveness as percentage depending on scenario(4 scenarios)
name <- list("+QALYs, - Costs", "-QALYs, -Costs", "+QALYs, +Costs", "-QALYs, +Costs")
value <- list(sum(ICER_psa_geo1<0.01 & ICER_psa_geo2>0.01), sum(ICER_psa_geo1<0.01 & ICER_psa_geo2<0.01), sum(ICER_psa_geo1>0.01 & ICER_psa_geo2>0.01), sum(ICER_psa_geo1>0.01 & ICER_psa_geo2<0.01))
datageonodf <- cbind( name, value)
datacegeo <- data.frame(datageonodf)
datacegeo[1,3]= 555/1000*100
datacegeo[2,3]= 280/1000*100
datacegeo[3,3]= 104/1000*100
datacegeo[4,3]= 61/1000*100
# Uniform color GRAPH!
tiff(file="Georgia_psaICER_varA.tiff",
     width=10, height=9, units="in", res=400)
par(mfrow=c(1,1))
barplot(height=datacegeo$V3, names=datacegeo$name, 
        col= "bisque1", cex.axis=1.2, cex.names=1.4,
        horiz=T, las=1, yaxt="n", xlab = "Percentage", ylab = "", main= "Georgia", xlim=c(0, 100))
text(60,1:4,name,cex=1.5)
mtext("Variation in Costs and QALYs values between mSTR and conventional treatment", side=2, line=2.2, cex=1.3)
dev.off()
#Graph3; effectiveness using WTP-thresholds
effective_ICER_geo <- matrix(0, nrow= 1000, ncol=1)
wtp_geo <- matrix(0, nrow= 1000, ncol=1)
wtp_geo[1,1] <-0
effective_ICER_geo[1,1] <- (sum(ICER_psa_Geo < 0, na.rm=TRUE))/1000
for (i in 2:1000) {
  effective_ICER_geo[i,1] <- (sum(ICER_psa_Geo < (i-1)*10, na.rm=TRUE))/1000
  wtp_geo[i,1] <- (i-1)*10
}
# Create DataFrame for Plotting
DF <- data.frame(wtp_geo, effective_ICER_geo)
tiff(file="Georgia_psaICER_varB.tiff",
     width=10, height=8, units="in", res=400)
par(mfrow=c(1,1))
ggplot(DF, aes(wtp_geo, effective_ICER_geo)) +                                     
  geom_line(color = "dodgerblue4", size = 2) +
  geom_ribbon(aes(ymin=effective_ICER_geo+0.02, ymax=effective_ICER_geo-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="Georgia", ylim=c(0.75, 1), x = "Willingness-to-pay (WTP) thresholds in $", y ="Percentage of simulations mSTR regimen is cost-effective (ICER<WTP)")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 6000), breaks=c(0, 1000,2000,3000,4000,5000,6000)) +
  scale_y_continuous(limits=c(.60, 0.95), breaks=c(0.6, 0.65, 0.7, 0.75,0.8, 0.85, 0.9, 0.95), labels = c("60", "65", "70","75","80", "85", "90", "95")) +
  geom_vline(xintercept = 2324, linetype="dotted", color = "blue", size=1.5)
dev.off()


#SavingICERScountries#####

tiff(file="bel_mol_kaz_geo_univSA_ICER.tiff",
     width=12, height=10, units="in", res=600)
par(mfrow = c(2, 2))
plot(ICER_psa_kaz1, ICER_psa_kaz2, pch = 21, main= "(A) Kazakhstan      ",cex.main = 1.8,
            bg = "bisque",   # Fill color
            col = "darkgoldenrod1", # Border color
            cex = 2,      # Symbol size
            lwd = 3, 
            type='p', cex.lab = 1,  cex.axis = 1.3,
            ylab = "Incremental QALYs gained", cex.lab = 1.3,
            xlab = "Incremental economic costs in 2022 USD", las = 1, xlim=c(-50000, 20000), ylim=c(-5, 10))
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_kaz1),mean(ICER_psa_kaz2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-27000, -5, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
plot(ICER_psa_bel1, ICER_psa_bel2, pch = 21, main= "(B) Belarus      ", cex.main = 1.8,
     bg = "bisque",   # Fill color
     col = "darkgoldenrod1", # Border color
     cex = 2,      # Symbol size
     lwd = 3, 
     type='p', cex.lab = 1,  cex.axis = 1.3,
     ylab = "Incremental QALYs gained", cex.lab = 1.3,
     xlab = "Incremental economic costs in 2022 USD", las = 1, xlim=c(-60000, 20000), ylim=c(-5, 8))      
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_bel1),mean(ICER_psa_bel2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-36000, -40, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
plot(ICER_psa_mol1, ICER_psa_mol2, pch = 21, main= "(C) Moldova      ",cex.main = 1.8,
     bg = "bisque",   # Fill color
     col = "darkgoldenrod1", # Border color
     cex = 2,      # Symbol size
     lwd = 3, 
     type='p', cex.lab = 1,  cex.axis = 1.3,
     ylab = "Incremental QALYs gained", cex.lab = 1.3,
     xlab = "Incremental economic costs in 2022 USDs", las = 1, xlim=c(-40000, 15000), ylim=c(-3, 8))      
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_mol1),mean(ICER_psa_mol2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)
text(-36000, -200, "▲ = Average ICER for mSTR compared to conventional treatment", cex = 0.9)
plot(ICER_psa_geo1, ICER_psa_geo2, pch = 21, main= "(D) Georgia      ",cex.main = 1.8,
     bg = "bisque",   # Fill color
     col = "darkgoldenrod1", # Border color
     cex = 2,      # Symbol size
     lwd = 3, 
     type='p', cex.lab = 1,  cex.axis = 1.3,
     ylab = "Incremental QALYs gained", cex.lab = 1.3,
     xlab = "Incremental economic costs in 2022 USD", las = 1, xlim=c(-50000, 20000), ylim=c(-5, 10))
abline(v=c(0),col="black", lty=12, lwd=2)
abline(h=c(0),col="black", lty=12, lwd=2)
points(mean(ICER_psa_geo1),mean(ICER_psa_geo2), col="black", cex=3,  bg = "black", lwd = 2, type = "p", pch = 17)

dev.off()
#WTP per COUNTRY####

library(patchwork)
p1<- ggplot(DF, aes(wtp_kaz, effective_ICER_kaz)) +                                     
  geom_line(color = "dodgerblue4", size = 2) +
  geom_ribbon(aes(ymin=effective_ICER_kaz+0.02, ymax=effective_ICER_kaz-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="(A) Kazahstan     ", ylim=c(0.70, 1), x = "WTP thresholds in $", y ="% of simulations indicating ICER<WTP")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 10080), breaks=c(0, 2000,4000,6000,8000,10000)) +
  scale_y_continuous(limits=c(.70, 0.95), breaks=c(0.7,0.75,0.8, 0.85, 0.9, 0.95), labels = c("70","75","80", "85", "90", "95")) +
  geom_vline(xintercept = 10073, linetype="dotted", color = "blue", size=1.5)
p2<- ggplot(DF, aes(wtp_bel, effective_ICER_bel)) +                                     
  geom_line(color = "dodgerblue4", linewidth = 2) +
  geom_ribbon(aes(ymin=effective_ICER_bel+0.02, ymax=effective_ICER_bel-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="(B) Belarus     ", ylim=c(0.50, 0.85), x = "WTP thresholds in $", y ="% of simulations indicating ICER<WTP")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 10000), breaks=c(0, 2000,4000,6000,8000,10000)) +
  scale_y_continuous(limits=c(0.50, 0.85), breaks=c(0.5,0.55, 0.6, 0.65, 0.70, 0.75, 0.8, 0.85), labels = c("50","55","60", "65", "70", "75", "80", "85")) +
  geom_vline(xintercept = 6102, linetype="dotted", color = "blue", size=1.5)
p3<-ggplot(DF, aes(wtp_mol, effective_ICER_mol)) +                                     
  geom_line(color = "dodgerblue4", size = 2) +
  geom_ribbon(aes(ymin=effective_ICER_mol+0.02, ymax=effective_ICER_mol-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="(C) Moldova     ", ylim=c(0.75, 0.95), x = "WTP thresholds in $", y ="% of simulations indicating ICER<WTP")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 3000), breaks=c(0, 500,1000, 1500,2000,2500,3000)) +
  scale_y_continuous(limits=c(.75, 0.95), breaks=c(0.75,0.80, 0.85, 0.9, 0.95), labels = c("75","80","85", "90", "95")) +
  geom_vline(xintercept = 1446, linetype="dotted", color = "blue", size=1.5)
p4<- ggplot(DF, aes(wtp_geo, effective_ICER_geo)) +                                     
  geom_line(color = "dodgerblue4", size = 2) +
  geom_ribbon(aes(ymin=effective_ICER_geo+0.02, ymax=effective_ICER_geo-0.02), alpha=0.1, fill = "deepskyblue", 
              color = "deepskyblue", linetype = "solid") +
  labs(title="(D) Georgia     ", ylim=c(0.6, 0.85), x = "WTP thresholds in $", y ="% of simulations indicating ICER<WTP")+
  theme_linedraw(base_size = 16 ) +
  theme(title = element_text(face="bold"))+
  scale_x_continuous(limits=c(0, 6000), breaks=c(0, 1000,2000,3000,4000,5000,6000)) +
  scale_y_continuous(limits=c(.60, 0.85), breaks=c(0.6, 0.65, 0.7, 0.75,0.8, 0.85), labels = c("60", "65", "70","75","80", "85")) +
  geom_vline(xintercept = 2324, linetype="dotted", color = "blue", size=1.5)

tiff(file="WTP_combo.tiff",
     width=12, height=10, units="in", res=600)
p1+p2+p3+p4 + plot_layout(ncol = 2, nrow = 2) 
dev.off()

################################################################################
 #Graph of the costs per country.
library(readxl)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
library(cowplot)
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/Bloomsbury_PolicyLab/0_Article/Figures")
# Read an Excel file into R
data <- read_excel("Costs_Book1.xlsx")
df <- as.data.frame(data)


p1<-ggplot(df, aes(x = CountryT, y=Value, fill = Cost_t)) +
  geom_bar(stat = "identity", position = "stack", color='black') +
  coord_flip()+
  scale_fill_brewer(palette="Set3") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 16),     # x-axis labels
        axis.text.y = element_text(size = 18),     # y-axis labels
        axis.title.x = element_text(size = 23),    # x-axis title
        axis.title.y = element_text(size = 23),    # y-axis title
        legend.title = element_text(size = 19),    # Legend title
        legend.text = element_text(size = 18)) +  # Legend text)
  labs(y = "Average total treatment cost per patient (USD)", fill = "Cost type", x = "Country, treatment arm") +
  scale_y_continuous(breaks = c(0, 2500, 5000, 7500, 10000, 12500, 15000)) +   # Custom x-axis breaks
  scale_fill_brewer(palette="Set3", labels = c("Drugs", "Inpatient", "Monitoring", "Observation", "Other", "Outpatient"))

ggsave(filename = "Country_arm_cost.tiff", plot = p1, device = "tiff", 
       width = 12, height = 9, dpi = 600) 


p2<-ggplot(df, aes(x = CountryT, y=Value, fill = Cost_t)) +
  geom_bar(stat = "identity", position = "fill", color='black') +
  coord_flip()+
  scale_fill_brewer(palette="Set3") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 16),     # x-axis labels
        axis.text.y = element_text(size = 18),     # y-axis labels
        axis.title.x = element_text(size = 23),    # x-axis title
        axis.title.y = element_text(size = 23),    # y-axis title
        legend.title = element_text(size = 19),    # Legend title
        legend.text = element_text(size = 18)) +  # Legend text)
  labs(y = "Percentage (%)", fill = "Cost type", x = "Country, treatment arm") +
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1), labels = c("0", "20", "40", "60", "80", "100")) +   # Custom x-axis breaks
  scale_fill_brewer(palette="Set3", labels = c("Drugs", "Inpatient", "Monitoring", "Observation", "Other", "Outpatient"))

  
ggsave(filename = "Country_arm_cost_p.tiff", plot = p2, device = "tiff", 
       width = 12, height = 9, dpi = 600) 

combined_plot <- plot_grid(p1, p2, ncol=1, labels= c("(A)", "(B)"), label_size=24)

ggsave(filename = "Country_arm_cost_FULL.tiff", plot = combined_plot, device = "tiff", 
       width = 13, height = 13, dpi = 600) 

################################################################################
# MAPS                                                                          #
##############################################################################
#Armenia ARM
#Republic of Moldova MDA
#Kazakhstan KAZ
#Azerbaijan AZE
#Belarus BLR
#Ukraine UKR
#Uzbekistan UZB
#Tajikistan TJK
#Turkmenistan TKM
#Georgia GEO
#Kyrgyzstan KGZ
library(dplyr)
require(maps)
require(viridis)
theme_set(
  theme_void()
)

world_map <- map_data("world")
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")

some.countries <- c(
  "Armenia", "Moldova", "Kazakhstan", "Azerbaijan", "Belarus",
  "Ukraine", "Uzbekistan", "Tajikistan", "Turkmenistan",
  "kyrgyzstan", "Georgia")

somemaps <- map_data("world", region = some.countries)

region.lab.data <- somemaps %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

ggplot(somemaps, aes(x = long, y = lat)) +
  geom_polygon(aes( group = group, fill = region))+
  geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none") 
  

#https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html
#https://www.datanovia.com/en/blog/how-to-create-a-map-using-ggplot2/