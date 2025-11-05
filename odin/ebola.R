#   #########
#  ###########
# #############
# # VHF MODEL #
# #############
#  ###########
#   #########



#######################################################
# Core equations for transitions between compartments #
#######################################################

## Note the use of [i] syntax for indexing patches
##
## Compartments include:
## - S: susceptible
## - E: exposed (not infectious yet)
## - I: infectious
## - F: dead causing funeral exposures
## - C: convalescent, displaying lesser infectivity (e.g. sexual transmission);
##   only a fraction of infected individuals who survive will get to this state,
##   the others getting directly 'recovered'
## - V: vaccinated, who are no longer infectious; no waning immunity, so this
##   compartment is a dead end
## - R: recovered, no longer infectious
## - D: dead, no longer infectious
##
## Note that in addition, the intervention status `interv` will be added to
## outputs.

update(S[]) <- S[i] - n_SE[i] - n_SV[i]
update(V[]) <- V[i] + n_SV[i]
update(E[]) <- E[i] + n_SE[i] - n_EI[i]
update(I[]) <- I[i] + n_EI[i] - n_IR[i] - n_IC[i] - n_ID[i] - n_IF[i]
update(F[]) <- F[i] + n_IF[i] - n_FD[i]
update(D[]) <- D[i] + n_FD[i] + n_ID[i] 
update(C[]) <- C[i] + n_IC[i] - n_CR[i]
update(R[]) <- R[i] + n_IR[i] + n_CR[i]

## N counts all living individuals, so all compartments except F and D.
N[] <- S[i] + E[i] + I[i] + C[i] + V[i] + R[i]




###########################
# Monitoring intervention #
###########################

## We monitor intervention for each patch at the current step. There are 3
## different states, recorded at the current time step by by `interv_state`:
##
## - interv_state:
##      + 1: 'naive', i.e. no intervention
##      + 2: 'alert', i.e. some intervention due to neighbouring patches having
##      cases
##      + 3: 'response', i.e. full intervention due to cases in the patch
##
## The following user-provided inputs are used to handle these states:
##
## - alert_threshold: an integer of length 2; alert_threshold[1] is the
##   minimum (inclusive) average number of infected coming from neighboring
##   patches for intervention to kick in; alert_threshold[2] is the number of
##   infected below which intervention will stop, after a number of days of this
##   condition being true; note that number of infected needs to be strictly 
##   less than alert_threshold[2] for stoppage
##
## - alert_release: the number of days below the alert threshold after
##   which the intervention stops
##
## - response_threshold: an integer of length 2; response_threshold[1] is the
##   minimum (inclusive) number of infected inside the patch for intervention to
##   kick in; response_threshold[2] is the number of infected below which
##   intervention will stop, after a number of days of this condition being
##   true; note that number of infected needs to be strictly less than
##   response_threshold[2] for stoppage
##
## - response_release: the number of days below the response threshold after
##   which the intervention stops
##
## These variables are used locally (not user-provided) as timers for the 
## release of alert and response states
##
## - days_below_alert_threshold: counts the number of days I[] has been below 
##   the intervention threshold
##
## - days_below_response_threshold: counts the number of days I[] has been below
##   the intervention threshold
##
##
## Examples of response trigger and release:
##
## 1) response_threshold = c(1, 0) and response_release = 10: interventions will
## start as soon as there is one infected in a patch, and will stop after 10
## days without cases
##
## 2) response_threshold = c(50, 10) and response_release = 21: interventions
## will start as soon as there are at least 50 infected in a patch, and will
## stop if there are strictly less than 10 cases every day for 21 days

### calculation of average incoming I for alert triggers
mat_I_neighbours[, ] <- delta[i, j] * I[i]
sum_I_neighbours[] <- sum(mat_I_neighbours[, i]) - delta[i, i] * I[i]

### initialize interv_state
interv_state[] <- if (as.integer(step) == 1) 1 else interv_state[i]

### trigger alert/response; response takes precedence over alert
interv_state[] <- if (sum_I_neighbours[i] >= alert_threshold[1]) 2 else interv_state[i]
interv_state[] <- if (I[i] >= response_threshold[1]) 3 else interv_state[i]

### timers of alert/response release
days_below_alert_threshold[] <- if 
  (sum_I_neighbours[i] < alert_threshold[2]) days_below_alert_threshold[i] + 1 else 0
days_below_response_threshold[] <- if 
  (I[i] < response_threshold[2]) days_below_response_threshold[i] + 1 else 0


### release alert/response where timers expired
###
### - 'alert' patches that hit the release timer go from 2 to 1
### - 'response' patches that hit the release timer go from 3 to 2 or 1
### - 'response' needs to be demoted first, so that situation like 3->2 and 2->1 
###  could happen (whenever both release timers expire at the same time)

### demotion 3 -> 2
interv_state[] <- if (interv_state[i] == 3 &&  # was in response
                      (days_below_response_threshold[i] >= response_release) && # hit the timer
                      (sum_I_neighbours[i] >= alert_threshold[1]) # alert is triggered
                     ) 2 else interv_state[i]

### demotion 3 -> 1
interv_state[] <- if (interv_state[i] == 3 && # was in response
                      (days_below_response_threshold[i] >= response_release) && # hit the timer
                      (sum_I_neighbours[i] < alert_threshold[1]) # no alert
                     ) 1 else interv_state[i]

### demotion 2 -> 1
interv_state[] <- if (interv_state[i] == 2 && # was in alert
                      (days_below_alert_threshold[i] >= alert_release) # hit the timer
                     ) 1 else interv_state[i]

## We also output the intervention state, as it will be needed for
## re-calculating FOIs when post-processing simulation results
output(interv[]) <- interv_state[i]



############################
# Transition S->E and S->V #
############################

# Note on individual probabilities of transition: see file sir.odin for detailed
# explanations on the calculations of the individual force of infection and
# disease dispersal across patches.

## Individuals can leave S for two reasons:
## 1. they become exposed (S->E); this is determined by a FOI with different
##   components
## 2. they get vaccinated (S->V); this happens with a given probability, if
## intervention has been triggered


############################
## Rate calculation: S->E ##
############################

## Notes on `foi_mode`
##
## An integer indicating the mode of calculation of the FOI; currently can be:
## - 1: (default) frequency dependent: FOI = beta * S * I / N
## - 2: density dependent: FOI = beta * S * I
##
## This applies to both FOIs, from infected individuals as well as funeral
## exposures.

## Notes on the terms for S->E
##
## Calculations a similar to the transition S->I in the SIR model (file
## `sir.odin`) but with an extra term for funeral exposure.

## The FOI has four sources:
## 1. infected individuals with a rate beta_I; can be lessened during
##    intervention by a factor `response_efficacy`
## 2. funeral exposure with a rate beta_F
## 3. convalescent individuals with a rate beta_C
## 4. zoonotic introduction with a rate eta (changing in time and space)

## Note on scaling FOIs
##
## Most transitions between compartments are handled via the following process:
##
## 1. calculate patch-wise FOI
## 2. convert into *individual* probas of transition
## 3. binomial draws using these probas
##
## We need to be very careful to rescale patch-wise FOI to per-capita rates
## (i.e. per susceptible individual) before converting to probabilities; this is
## done by ommitting the '* S' term in all the FOI calculations below.

### FOI from I - can be frequency or density dependent
lambda_prod_infec[, ] <- if
  (as.integer(foi_mode_I) == 1 && N[i] > 0) (delta[i, j] * I[i] / N[i]) / area[i] else (delta[i, j] * I[i]) / area[i]
lambda_infec[] <- sum(lambda_prod_infec[, i]) * beta_I # sum by column

lambda_infec[] <- if 
  (interv_state[i] == 3) lambda_infec[i] * (1 - response_efficacy) else if # response mode
  (interv_state[i] == 2) lambda_infec[i] * (1 - alert_efficacy) else # alert mode
   lambda_infec[i] # naive mode

### FOI from F - can be frequency or density dependent
lambda_prod_funer[, ] <- if (as.integer(foi_mode_F) == 1 && (F[i] + N[i]) > 0)
                           (delta[i, j] * F[i] / (F[i] + N[i])) / area[i] else (delta[i, j] * F[i]) / area[i]
lambda_funer[] <- sum(lambda_prod_funer[, i]) * beta_F # sum by column

### FOI from C - always frequency-dependent
lambda_prod_conval[, ] <- if (N[i] > 0) (delta[i, j] * C[i] / N[i]) / area[i] else 0
lambda_conval[] <- sum(lambda_prod_conval[, i]) * beta_C # sum by column

### Total FOI (adding FOI 4)
rate_SE[] <- lambda_infec[i] + lambda_funer[i] + lambda_conval[i] + eta[as.integer(step), i]


############################
## Rate calculation: S->V ##
############################
##
## This part is a bit convoluted as we need to express this transition as a rate
## first, so that we can determine the number of people leaving S, and then
## decide if they go to E or V. However, user inputs for vaccination (S->V) are
## probabilities (i.e. vaccination coverage), which we need to convert this back
## to a rate.

p_SV[] <- if
  (interv_state[i] == 3) vacc_coverage_response * vacc_efficacy else if # response mode
  (interv_state[i] == 2) vacc_coverage_alert * vacc_efficacy else # alert mode
  0 # naive mode

rate_SV[] <- -log(1 - p_SV[i])


##################################
## Binomial draws S->E and S->V ## 
##################################

## There are competing events to leave S: either go to E or V. We use a binomial
## distribution to decide where individuals go, based on the ratio of the
## respective rates. However we need to be careful when no individual leaves to
## avoid 0/0, as rbinom(0, 0/0) outputs NaN. See issue 24 on gitlab:
## https://gitlab.geomatys.com/produits/geomatys/librairies/r_/epirs/-/issues/24

rate_S_out[] <- rate_SE[i] + rate_SV[i]
p_S_out[] <- 1 - exp(-rate_S_out[i])
n_S_out[] <- rbinom(S[i], p_S_out[i])
n_SE[] <- if (n_S_out[i] > 0) rbinom(n_S_out[i], rate_SE[i] / rate_S_out[i]) else 0
n_SV[] <- n_S_out[i] - n_SE[i]


###################
# Transition E->I #
###################

p_EI <- 1 - exp(-alpha) # E to I
n_EI[] <- rbinom(E[i], p_EI)


#####################################
# Transition I->F, I->C, I->R, I->D #
#####################################

## General process:
##
## These two transitions are linked and must be handled together. The strategy
## is as follows:
##
## 1. Determine the rate at which individuals leave I, including both possible
##    outcomes, and the corresponding individual probability of leaving I
## 2. Draw the number of individuals leaving I from a Binomial distribution
## 3. Determine how many individuals survive or die.
## 4. Determine how many deaths result in funeral exposure.
## 5. Determine how many survivors move to C or R from the prop of convalescent.

### those leaving I
p_I_out <- 1 - exp(-gamma)
n_I_out[] <- rbinom(I[i], p_I_out)

### those who die (F and C)
n_I_dead[] <- rbinom(n_I_out[i], cfr)
p_IF[] <- if 
  (interv_state[i] == 3) (1 - p_sdb_response) else if # response mode
  (interv_state[i] == 2) (1 - p_sdb_alert) else # alert mode
  (1 - p_sdb) # naive mode

n_IF[] <- rbinom(n_I_dead[i], p_IF[i]) # funeral exposures
n_ID[] <- n_I_dead[i] - n_IF[i]

### those who survive (C and R)
n_survivors[] <- n_I_out[i] - n_I_dead[i]
n_IC[] <- rbinom(n_survivors[i], pi) # convalescents
n_IR[] <- n_survivors[i] - n_IC[i] # recovered


###################
# Transition F->D #
###################

## This one is simply governed by the rate at which infectious bodies are
## burried (`tau`), which is constant over time and space.

p_FD <- 1 - exp(-tau)
n_FD[] <- rbinom(F[i], p_FD)


###################
# Transition C->R #
###################

## Individuals become fully recovered after a period of reduced infectiousness
## (typically sexual transmission) at a constant rate.

p_CR <- 1 - exp(-rho)
n_CR[] <- rbinom(C[i], p_CR)


##################
# Initial states #
##################

initial(S[]) <- ini_S[i]
initial(V[]) <- ini_V[i]
initial(E[]) <- ini_E[i]
initial(I[]) <- ini_I[i]
initial(F[]) <- ini_F[i]
initial(D[]) <- ini_D[i]
initial(C[]) <- ini_C[i]
initial(R[]) <- ini_R[i]


# User defined parameters - default in parentheses
## Initialisations

n_patches <- user(1)

ini_S[] <- user()
ini_V[] <- user()
ini_E[] <- user()
ini_I[] <- user()
ini_F[] <- user()
ini_D[] <- user()
ini_C[] <- user()
ini_R[] <- user()

area[] <- user()

alpha <- user()
beta_I <- user()
beta_F <- user()
beta_C <- user()
gamma <- user()
eta[,] <- user()
delta[,] <- user()
tau <- user()
pi <- user()
rho <- user()
cfr <- user()
p_sdb <- user()

alert_threshold[] <- user()
alert_release <- user()
alert_efficacy <- user()
p_sdb_alert <- user()
vacc_coverage_alert <- user()

response_threshold[] <- user()
response_release <- user()
response_efficacy <- user()
p_sdb_response <- user()
vacc_coverage_response <- user()

vacc_efficacy <- user()
foi_mode_I <- user(1)
foi_mode_F <- user(1)


## Dimensions of the variables

dim(ini_S) <- n_patches
dim(ini_V) <- n_patches
dim(ini_E) <- n_patches
dim(ini_I) <- n_patches
dim(ini_F) <- n_patches
dim(ini_D) <- n_patches
dim(ini_C) <- n_patches
dim(ini_R) <- n_patches

dim(S) <- n_patches
dim(E) <- n_patches
dim(I) <- n_patches
dim(F) <- n_patches
dim(D) <- n_patches
dim(C) <- n_patches
dim(V) <- n_patches
dim(R) <- n_patches
dim(N) <- n_patches

dim(area) <- n_patches

dim(interv_state) <- n_patches
dim(interv) <- n_patches
dim(days_below_alert_threshold) <- n_patches
dim(days_below_response_threshold) <- n_patches

dim(n_SE) <- n_patches
dim(n_SV) <- n_patches
dim(n_S_out) <- n_patches
dim(n_EI) <- n_patches
dim(n_I_out) <- n_patches
dim(n_survivors) <- n_patches
dim(n_I_dead) <- n_patches
dim(n_IF) <- n_patches
dim(n_ID) <- n_patches
dim(n_IC) <- n_patches
dim(n_IR) <- n_patches
dim(n_CR) <- n_patches
dim(n_FD) <- n_patches

dim(delta) <- c(n_patches, n_patches)
dim(mat_I_neighbours) <- c(n_patches, n_patches)
dim(sum_I_neighbours) <- n_patches
dim(lambda_prod_infec) <- c(n_patches, n_patches)
dim(lambda_infec) <- n_patches
dim(lambda_prod_funer) <- c(n_patches, n_patches)
dim(lambda_funer) <- n_patches
dim(lambda_prod_conval) <- c(n_patches, n_patches)
dim(lambda_conval) <- n_patches

dim(alert_threshold) <- 2 # threshold to [1] activate [2] deactivate alert
dim(response_threshold) <- 2 # threshold to [1] activate [2] deactivate response

dim(eta) <- user()
dim(rate_SE) <- n_patches
dim(p_SV) <- n_patches
dim(rate_SV) <- n_patches
dim(rate_S_out) <- n_patches
dim(p_S_out) <- n_patches
dim(p_IF) <- n_patches

