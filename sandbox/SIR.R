# https://staff.math.su.se/hoehle/blog/2020/03/16/flatteningthecurve.html


######################################################################
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list with one component being a vector of length two containing
#  dS(t)/dt and dI(t)/dt
######################################################################

sir <- function(t, y, parms) {
  beta <- parms[1]
  gamma <- parms[2]
  S <- y[1]
  I <- y[2]
  return(list(c(S = -beta * S * I, I = beta * S * I - gamma * I)))
}

# Population size
N <- 1e6
# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
gamma <- 1/5
# Infectious contact rate - beta = R0/N*gamma and when R0 \approx 2.25 then  2.25/N*gamma
beta <- 4.5e-07
# R0 for the beta and gamma values
R0 <- beta*N/gamma

# Load package to numerically solve ODEs
suppressPackageStartupMessages(library(deSolve))

# Grid where to evaluate
max_time <- 150
times <- seq(0, max_time, by=0.01)

# Solve ODE system using Runge-Kutta numerical method.
ode_solution <- rk4(y = c(N - 10, 10), times = times,
                    func = sir, parms = c(beta, gamma)) %>%
  as.data.frame() %>%
  setNames(c("t", "S", "I")) %>%
  mutate(beta = beta, gama = gamma, R0 = N * beta / gamma,
         s = S / N, i = I / N, type = "without_intervention")

