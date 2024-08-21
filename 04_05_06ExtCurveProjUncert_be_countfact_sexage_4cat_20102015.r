################################################################################
# This code is based on the article:
#   "A hands-on tutorial on a modelling framework for projections of climate 
#    change impacts on health"
#   Ana M. Vicedo-Cabrera, Francesco Sera, Antonio Gasparrini.

# The supplement data can be found at (https://github.com/huychanchan98/....)
# ADAPTED BY: HUY TRAN DUC - IUPWARE - KULEUVEN - MASTER THESIS PROJECT
################################################################################

################################################################################
# 04 EXTRAPOLATION OF THE EXPOSURE-RESPONSE CURVE
# 05 PROJECTION & QUANTIFICATION OF THE IMPACT
# 06 ENSEMBLE ESTIMATES & QUANTIFICATION OF THE UNCERTAINTY
################################################################################

# The three last steps of the analysis (extrapolation of the curve, 
#   impact projections and quantification of the uncertainty) can be performed 
#   sequentially using the following code.

# In brief, once we extrapolate the curve, we estimate the daily number of attributable 
#   deaths (AN) in each scenario, GCM and temperature range. 
# Then, we compute the sum of ANs per each period for the temperature range. 
# By dividing between the total mortality, we estimate the corresponding 
#   attributable fractions (AFs).
# Uncertainty of the estimated impacts is expressed in terms of empirical 
#   confidence intervals, defined as the 2.5th and 97.5th percentiles of the 
#   empirical distribution of the impacts across coefficients samples and GCMs. 
#   The distribution is obtained through Monte Carlo simulation of the coefficients.

library(zoo)
# EXTRACT COEFFICIENTS AND VCOV FROM THE MODEL IN STAGE 1
# With the crossreduce function we can reduce the fit of the bidimensional DLNM
#   (of the stage 1 using the observed temperature-mortality series)
#   to summaries defined in the exposure-response dimension. We can then 
#   extract the coefficients and the covariance matrix defining the overall 
#   exposure-response association.
setwd("D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/code/tutorial/be2020")

# DEFINE PROJECTED MORTALITY SERIES
# It is computed as the average mortality for each day of the year 
#   from daily observed deaths, then repeated along the same projection period
#   of the modelled temperature series.

# Define a function to compute the lagged values
Lag <- function(x, k) {
  sapply(k, function(ki) {
    if (ki > 0) {
      c(rep(NA, ki), x[1:(length(x) - ki)])
    } else if (ki < 0) {
      c(x[(-ki + 1):length(x)], rep(NA, -ki))
    } else {
      x
    }
  })
}

deathdoy <- tapply(obs$all,as.numeric(format(obs$date,"%j")),
                  mean,na.rm=T)[seq(365)]
while(any(isna <- is.na(deathdoy)))
 deathdoy[isna] <- rowMeans(Lag(deathdoy,c(-1,1)),na.rm=T)[isna]
deathproj <- rep(deathdoy,length=nrow(obs)) #Check length=nrow() obs or rcp8p5 !

# STORE THE MORTALITY BY PERIOD
deathperiod <- sum(deathdoy)*10
##########################

# Age groups list
m_category<-list(m0_64 ,m65_74,m75_84,m85plus)

# DIMENSION - RANGE OF TEMPERATURES
#temprange <- c("tot","cold","heat")

# DIMENSION - ABSOLUTE AN/DIFFERENCE IN AN
#absrel <- c("abs","rel")
# DIMENSION - GENERAL CIRCULATION MODELS
#gcm_be <- c("RMI"="tmean_rmi")

################################################################################
####### ADAPTATION WITH THE COUNTERFACTUALS CODE ###############################
# DIMENSION - NUMBER OF ITERATION IN THE MONTE-CARLO SIMULATION 
nsim <- 1000

# DEFINE THE ARRAY TO STORE RESULTS
ansim_attr <- array(NA,dim=list(length(m_category),5,nsim+1),
                    dimnames=list(c("<65y","65-74y","75-84y",">85y"),
                                  c("obs","mod1","mod2","mod3","mod4"),
                                  c("est",paste0("sim",seq(nsim)))))

# RUN ACROSS AGE GROUPS
for (i in seq(m_category)){
  
  m_cat=m_category[[i]]
  red <- crossreduce(cb,m_cat,cen=cen)
  coef <- coef(red)
  vcov <- vcov(red)

  # FACTORS TO DEFINE SCENARIOS (FACTUAL + 4 COUNTERFACTUALS)
  # NEED TO CHECK ADJUST FOR BELGIUM
  scen_fact <- c(0,2.06,2.65,1.09,2.17) #2.65 (HADCrut5), (-0.1 compare to Swiss)
  
    for (g in seq(length(scen_fact))){
    
      # DEFINE THE TMEAN (SUBSTRACT FACTOR SCENARIO)
      tmeanmod <- obs[,"tmean"] - scen_fact[g]
        
      # DERIVE THE CENTERED BASIS
      bvar <- do.call(onebasis,c(list(x=tmeanmod),argvar))
      cenvec <- do.call(onebasis,c(list(x=cen),argvar))
      bvarcen <- scale(bvar,center=cenvec,scale=F)
    
      # INDICATOR FOR HEAT DAYS
      indheat <- tmeanmod>cen 
    
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- (1-exp(-bvarcen%*%coef))*deathproj
    
      # SUM 
      ansim_attr[i,g,"est"] <- sum(an[indheat], na.rm=T) #Check [] number of subscript

      # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
      set.seed(13041975+g)
      coefsim <- mvrnorm(nsim,coef,vcov)
    
    # LOOP ACROSS ITERATIONS
      for(s in seq(nsim)) {
      
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
        an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathproj
      
      # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
      # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
        ansim_attr[i,g,s+1] <- sum(an[indheat], na.rm=T)
    }
  }}


################################################################################
# SUMMARIZE THE RESULTS 
# COMPUTE AN/AF (95%CI) IN THE ENSEMBLE, BY RANGE & PERIOD & RCP

# CREATE NEW OBJECTS TO STORE RESULTS
# We now create 4 new arrays (2 arrays to store the abs and rel AN, and another 2
#   to store the estimated rel and abs AF) with 4 dimensions each to store the 
#   ensemble estimates (average impacts across GCMs) with the empirical 
#   95% confidence intervals. 
# In this case, the 4 dimensions correspond to the 10-year period, the point estimate 
#   and the CI, the temperature range and the scenario. 

estci <- c("est","ci.l","ci.u")

anabs <- afabs <- anabs_mod <- afabs_mod <- array(NA,dim=c(length(m_category),length(estci)),
                                          dimnames=list(c("<65y","65-74y","75-84y",">85y"),estci))

# ATTRIBUTABLE NUMBERS 
# ABSOLUTE
## Factual (Obs data in ansim_attr) 
anabs[,"est"] <- apply(ansim_attr[,1,],1,mean)
anabs[,"ci.l"] <- apply(ansim_attr[,1,],1,quantile,0.025)
anabs[,"ci.u"] <- apply(ansim_attr[,1,],1,quantile,0.975)

## Counterfactual (mean of mod 1 2 3 4 in ansim_attr)
anabs_mod[,"est"] <- apply(ansim_attr[,2:5,], 1, mean)
anabs_mod[,"ci.l"] <- apply(ansim_attr[,2:5,], 1, quantile, 0.025)
anabs_mod[,"ci.u"] <- apply(ansim_attr[,2:5,], 1, quantile, 0.975)

# ATTRIBUTABLE FRACTION
afabs[,] <- anabs[,]/deathperiod*100
afabs_mod[,] <- anabs_mod[,]/deathperiod*100

# REMOVE ansim
#rm(ansim)

#----------------------------- BARPLOT -----------------------------------------
# Create afabs_data data frame
afabs_data <- data.frame(
  age_group = factor(c("<65y", "65-74y", "75-84y", ">85y"), levels = c("<65y", "65-74y", "75-84y", ">85y")),
  est = afabs[, "est"],
  ci.l = afabs[, "ci.l"],
  ci.u = afabs[, "ci.u"]
)

# Create obs_abs_f_data data frame
afabs_mod_data <- data.frame(
  age_group = factor(c("<65y", "65-74y", "75-84y", ">85y"), levels = c("<65y", "65-74y", "75-84y", ">85y")),
  est = c(afabs_mod[,"est"]),
  ci.l = c(afabs_mod[,"ci.l"]),
  ci.u = c(afabs_mod[,"ci.u"])
)
# Add a column to identify the dataset
afabs_data$dataset <- "afabs"
afabs_mod_data$dataset <- "afabs_mod"

# Combine both datasets
combined_data <- rbind(afabs_data, afabs_mod_data)

# Create the barplot
barplot <- ggplot(combined_data, aes(x=age_group, y=est, fill=dataset)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  geom_errorbar(aes(ymin=ci.l, ymax=ci.u), width=0.2, position=position_dodge(width=0.7), size=0.7) +
  scale_fill_manual(values = c("red","blue"),
                    labels=c("Factual","Counterfactual")) +
  labs(x = "Age Group", y = "Heat-related Mortality (%)",
       title = "Heat-related Mortality by Age Group") +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.spacing.x= unit(0.3, 'cm'),
        legend.text = element_text(size=11),
        plot.margin = unit(c(2, 1, 1, 0), units = "lines"),
        axis.text.y = element_text(size=10), 
        axis.title.y = element_text(size=12, vjust = 0.5),
        panel.grid.major.y = element_line(linewidth=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

# Display the plot
print(barplot)
# Save the plot
#ggsave("Heatmort_bel_male_age_20102015.pdf", plot = barplot, width=10, height=6)
#ggsave("Heatmort_bel_female_age_20102015.pdf", plot = barplot, width=10, height=6)



# #####################################################
# # FIGURE DENSITY BETWEEN OBS AND COUNTERFACTUAL SCENARIOS
# #pdf("Density_bel_summer2020.pdf",height=5,width=9)
# 
# # Initialize necessary variables
# scen_fact <- c(0, 2.06, 2.65, 1.09, 2.17) # Factor scenarios
# 
# # Initialize a list to store data frames for each scenario
# all_data_long <- list()
# 
# for (i in seq(m_category)) {
#   
#   m_cat = m_category[[i]]
#   red <- crossreduce(cb, m_cat, cen = cen)
#   coef <- coef(red)
#   vcov <- vcov(red)
#   
#   for (g in seq(length(scen_fact))) {
#     
#     # Define the modified temperature (subtract scenario factor)
#     tmeanmod <- obs$tmean - scen_fact[g]
#     
#     # Prepare the data frame
#     data <- data.frame(
#       obs_tmean = obs$tmean,
#       tmeanmod = tmeanmod,
#       scenario = paste("Scenario", scen_fact[g])
#     )
#     
#     # Reshape the data using pivot_longer and update Type names
#     data_long <- data %>%
#       pivot_longer(cols = c(obs_tmean, tmeanmod), names_to = "Type", values_to = "Temperature") %>%
#       mutate(Type = ifelse(Type == "obs_tmean", "Observed", scenario))
#     
#     # Append to the list
#     all_data_long[[g]] <- data_long
#   }
# }
# 
# # Combine all the data frames into one
# combined_data_long <- bind_rows(all_data_long)
# 
# # Ensure 'Observed' is plotted last
# combined_data_long$Type <- factor(combined_data_long$Type, levels = c("Scenario 0", "Scenario 2.06", "Scenario 2.65", "Scenario 1.09", "Scenario 2.17", "Observed"))
# 
# # Define colors for the plots with increasing opacity and lighter shades of yellow
# colors <- c("Observed" = "red", 
#             "Scenario 0" = "yellow", 
#             "Scenario 2.06" = "lightgoldenrodyellow", 
#             "Scenario 2.65" = "khaki", 
#             "Scenario 1.09" = "gold", 
#             "Scenario 2.17" = "orange")
# 
# # Create the density plot
# ggplot(combined_data_long, aes(x = Temperature, color = Type, fill = Type)) +
#   geom_density(alpha = 0.5, size = 1, color = "white") +
#   scale_color_manual(values = colors) +
#   scale_fill_manual(values = colors) +
#   labs(title = "Observed vs. Counterfactual Temperature",
#        x = "Mean Temperature (°C)",
#        y = "Density") +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#   plot.title = element_text(hjust = 0.5) # Center the title
# ) +
#   scale_x_continuous(limits = c(5, 40), breaks = seq(5, 40, by = 5)) # Set x-axis limits and breaks
# 
# #dev.off()
# # NB: uncomment pdf() and dev.off() for saving the plot in pdf format.

