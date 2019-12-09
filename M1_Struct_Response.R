# M1 implementation using lateral performance data and limit
# state functions for structural response (maximum transient
# interstory drift ratio at DBE and MCE earthquake Sa)
# RSB Overview Illustration
# Madeleine Flint, 2019-07 to 2019-12

library(mvtnorm)

# define names, IDA slopes, and other parameters of interest -----
# parameters from Multi-Hazard Performance Data Inventory, for
# mid-rise structures in moderate-to-high seismic zones
df.lat <- data.frame(Lateral = c("RCMF_Goul",
                                 "RCMF_Pit",
                                 "SCBF_Ak",
                                 "SCBF_Uriz",
                                 "SBRB_Uriz",
                                 "SBRT_Wong",
                                 "SMF_Lig",
                                 "SMF_MinW_Tehr",
                                 "SMF_Int_Tehr",
                                 "SMF_UnifE_Tehr",
                                 "SPSW_Purba"),
                     b =  c(38.90, # g/in/in, Sa = b*IDR
                            43.49,
                            62.89,
                            45.79,
                            63.90,
                            33.33,
                            23.86,
                            25.50,
                            30.54,
                            43.15,
                            47.37),
                     T = c(1.32, # period, seconds
                           0.66,
                           0.46,
                           0.5,
                           0.5,
                           0.94,
                           1.32,
                           0.69,
                           0.69,
                           0.69,
                           0.36),
                     Sa.C = c(2.80,
                              1.767,
                              2.53,
                              2.36,
                              3.07,
                              1.50,
                              1.67,
                              2.01,
                              2.31,
                              3.03,
                              3.6),
                     beta = c(0.34, # lognormal std. dev
                              0.3,  # or dispersion of Sa|IDR
                              0.37, # also equals beta_IDR|Sa
                              0.15, # as obtained from orthogonal
                              0.15, # regression
                              0.29,
                              0.39,
                              0.40,
                              0.18,
                              0.30,
                              0.40),
                     stringsAsFactors = FALSE)
save(df.lat, file = file.path("Data","M1.lateral.systems.RData"))
write.table(df.lat, file = file.path("Data","M1.lateral.systems.txt"), 
            sep = "\t", row.names = FALSE)


# Load hazard curve and prepare for numerical convolution -------
load(list.files("Data",pattern = "interp.RData", full.names = TRUE))
Sa <- Sa.interp
Sa.Exp <- c(0.0024, Sa, 3.1)
dSa <- Sa.Exp[2:length(Sa.Exp)] - Sa.Exp[1:(length(Sa.Exp)-1)]

# set up limit state functions for allowable drift at DBE and MCE events ----
MAFE.interp <- c(1/475,1/2475) # DBE, MCE
g1 <- function(x){ 0.01 -x}
g2 <- function(x){ 0.07 -x}

# find probability of exceeding allowable drift and system failure probability ---
# will use 0 and 1 as correlation bounds between limit state surfaces
df.lat[,c("Sa_DBE","Sa_MCE","med_IDR_DBE","med_IDR_MCE","P_C_MCE","MAF_C")] <- NA_real_
for(s in 1:nrow(df.lat)){
  Period <- paste0("Sa",df.lat$T[s])
  MAFE.Sa <- haz.all.melt$MAFE[haz.all.melt$Period==Period]
  df.lat[s,"Sa_DBE"] <- exp(approx(x=log(MAFE.Sa),y=log(Sa),xout=log(1/475))$y)
  df.lat[s,"Sa_MCE"] <- exp(approx(x=log(MAFE.Sa),y=log(Sa),xout=log(1/2475))$y)
  df.lat[s,"med_IDR_DBE"] <- 1/df.lat$b[s]*df.lat[s,"Sa_DBE"]
  df.lat[s,"med_IDR_MCE"] <- 1/df.lat$b[s]*df.lat[s,"Sa_MCE"]
  df.lat[s,"P_C_MCE"] <- pnorm((log(df.lat[s,"Sa_MCE"])-log(df.lat$Sa.C[s]))/df.lat$beta[s])
  df.lat[s,"beta1"] <- -(log(0.005)-log(df.lat[s,"med_IDR_DBE"]))/df.lat$beta[s]
  df.lat[s,"beta2"] <- -(log(0.03)-log(df.lat[s,"med_IDR_MCE"]))/df.lat$beta[s]
  df.lat[s,"pf12_cor0"] <- pmvnorm(lower=rep(-Inf,2), upper = c(df.lat[s,"beta1"],df.lat[s,"beta2"]),
                                 corr = matrix(c(1,0,0,1),nrow=2,ncol=2))
  df.lat[s,"pf12_cor1"] <- pmvnorm(lower=rep(-Inf,2), upper = c(df.lat[s,"beta1"],df.lat[s,"beta2"]),
                                      corr = matrix(1,nrow=2,ncol=2))
}

# view and save results ------
df.lat$Lateral[order(df.lat$pf12_cor0)]
df.lat$Lateral[order(df.lat$pf12_cor1)]
save(df.lat,file=file.path("Data","M1.lateral.systems.ranked.RData"))
write.table(df.lat, file = file.path("Data","M1.lateral.systems.ranked.txt"), 
            sep = "\t", row.names = FALSE)
