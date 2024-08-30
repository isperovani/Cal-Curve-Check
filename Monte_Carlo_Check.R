# This R script aims to perform Monte Carlo simulations for evaluating the probability of getting an out of specification result when checking the preparation of the calibration curve standard solutions

# In the pharmaceutical industry, this process is performed when the calibration curve is prepared by series dilutions of a standard stock solution (solution 1) to prepare the work standard solutions usually in a range from 80 - 120% of the API concentration in the drug product

# In order to check the preparation of the first standard work solution another standard stock solution (solution 2) and another work solution is prepared at a concentration close to the API concentration in the drug product

# Defining Variables 

mass = c() # Vector containing the values drawn from the Monte Carlo simulation (mg)

vpip = c() # Vector containing all the Monte Carlo simulation values of vpip (µL)

vinj = c() # Vector containing all the Monte Carlo simulation values of vinj (µL)

vflask = c() # Vector containing all the Monte Carlo simulation values of vflask (mL)

check = c() # Vector containing all the values calculated using the Monte Carlo simulation values

m_1 = c() # Vector containing values for the mass weighted for making the calibration curve (mg)

m_2 = c() # Vector containing values for the mass weighted to check the first weighting (mg)

vinj_1 = c() # Vector containing values for the injection volumes of the first injection (µL)

vinj_2 = c() # Vector containing values for the injection volumes of the check injection (µL)

vflask_s1 = c() # Vector containing values for the volumetric flask used to make the stock solution 1 (mL)

vflask_d1 = c() # Vector containing values for the volumetric flask used to make the series dilution of solution 1 (mL)

vflask_s2 = c() # Vector containing values for the volumetric flask used to make the stock solution 2 (mL)

vflask_d2 = c() # Vector containing values for the volumetric flask used to make the series dilution of solution 2 (mL)

vpip_1 = c() # Vector containing values pipetted for the calibration curve (µL)

vpip_2 = c() # Vector containing values pipetted for the check standard (µL)

mass_avg = 10 # Target mass to be weighted for calibration and check standards (mg)

sd_mass = 0.01 # Standard deviation of the weighting process (mg)

vinj_avg = 20 # Target injection volume of the method (µL)

sd_vinj = 0.1 # Standard deviation of the injection process (µL)

vpip_avg = 800 # Target pipeting volume for the second standard and check standard (µL)

sd_vpip = 0.96 # Standard deviation of the pipetting process (µL)

vflask_avg = 10 # Target size of the volumetric flasks for the calibration and check standards (mL)

sd_vflask = 0.02 # Standard deviation of the volumetric flask used CLASS A (mL)

k = 100 # Proportionality constant used to calculate integrated area of the peak (1/pg)

n_sim = 20000 # Number of Monte Carlo Simulations

drawings = 10000

# Monte Carlo Simulation

mass = rnorm(n_sim,
             mean = mass_avg,
             sd = sd_mass)

vpip = rnorm(n_sim,
             mean = vpip_avg,
             sd = sd_vpip)

vflask = rnorm(n_sim,
             mean = vflask_avg,
             sd = sd_vflask)

vinj = rnorm(n_sim,
               mean = vinj_avg,
               sd = sd_vinj)

# Drawing values for m_1, m_2, vpip_1, vpip_2, vinj_1, vinj_2, vflask_1 and vflask_2

m_1 = sample(mass,
             size = drawings,
             replace = F)

m_2 = sample(mass,
             size = drawings,
             replace = F)

vinj_1 = sample(vinj,
             size = drawings,
             replace = F)

vinj_2 = sample(vinj,
                size = drawings,
                replace = F)

vflask_s1 = sample(vflask,
                size = drawings,
                replace = F)

vflask_w1 = sample(vflask,
                   size = drawings,
                   replace = F)

vflask_s2 = sample(vflask,
                size = drawings,
                replace = F)

vflask_w2 = sample(vflask,
                   size = drawings,
                   replace = F)

vpip_1 = sample(vpip,
                size = drawings,
                replace = F)

vpip_2 = sample(vpip,
                size = drawings,
                replace = F)

# Calculating the check variable

for (i in seq_along(vflask_1)){
  
  check[i] = (k*m_1[i]*vpip_1[i]*vinj_1[i]*m_2[i]*vflask_w1[i]*vflask_s1[i])/(vflask_s2[i]*vflask_w2[i]*k*m_2[i]*vpip_2[i]*vinj_2[i]*m_1[i])
}

p_out = 100*(sum(check < 0.98) + sum(check >1.02))/length(check)

# Creates the histogram

hist(check,                                                 
     col = "blue",
     border = "black",
     main = "Figure B - Distribution of Check Standard Values",
     xlab = "Check Standard Values",
     ylab = "Frequency")

abline(v = 0.98,     # v is for vertical line, lty = 2 makes it dashed, lwt = 3 makes the line thicker
       col = "red",
       lty = 2,
       lwd = 3)  
abline(v = 1.02, 
       col = "red",
       lty = 2,
       lwd = 3)
