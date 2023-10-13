statusarray = function(df_frame, ages, nstates) {
  ar = array(
    data = 0,
    dim = c(length(ages), nrow(df_frame), length(nstates)+1),
    dimnames = list(ages, as.numeric(as.character(
      df_frame$identifier
    )), c(nstates,"out"))
  )
  
  df_frame$cad.prs = scale(df_frame$cad.prs)
  
  
  # from Health (1) to Health directly alone
  
  for (i in 1:length(ages)) {
    age = ages[i]
    agename = as.character(age)
    nx = age + 1
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age &
                        age < Death_Censor_Age,]
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Health"] =
      1
    
    rm(atrisk)
    
    
    
    # From HT
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age &
                      age < Death_Censor_Age,]
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Ht"] = 1
    
    rm(atrisk)
    ##  from Hyperlip
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < Dm_0_censor_age
                        & age < Death_Censor_Age,]
    
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "HyperLip"] =
      1
    
    rm(atrisk)
    
    
    
    
    ## from Dm2
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Death_Censor_Age,]
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Dm"] = 1
    
    rm(atrisk)
    
    atrisk = df_frame[age > Cad_0_censor_age &
                        Cad_0_Any == 2 & Death_Censor_Age > age,]
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Cad"] = 1
    
    rm(atrisk)
    
    
    ### 2 risk states
    ## from Ht and DM
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age & Dm_0_Any == 2 &
                        age > Ht_0_censor_age & Ht_0_Any == 2
                      & age < HyperLip_0_censor_age&age < Death_Censor_Age,]
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Ht&Dm"] = 1
    
    rm(atrisk)
    
    ## from HyperLip and DM to CAD marginal
    ## from HyperLip and Dm to Death directly alone (but could die with CAD) ###### BAD
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age &age < Death_Censor_Age,]
    
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "HyperLip&Dm"] =
      1
    
    rm(atrisk)
    
    ## from  HyperLip and Ht to HyperLip and Ht
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2
                      & age < Dm_0_censor_age&age < Death_Censor_Age, ]
    
                      
   ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Ht&HyperLip"] =
      1
    
    rm(atrisk)
    
    ## three risk states
    ## from HyperLip and Ht $ Dmto CAD marginal
    
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2
                      & age > Dm_0_censor_age & Dm_0_Any == 2&age < Death_Censor_Age,]
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "Ht&HyperLip&Dm"] =
      1
    
    
    rm(atrisk)
    
    ## from all 3 to death (but could die with CAD)
    
    
    atrisk = df_frame[Death_Censor_Age < age &
                        Death_censor_Any == 2,]
    
    if(nrow(atrisk)>0){ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "death"] =1}
    rm(atrisk)
    
  
    ### if you
    atrisk = df_frame[age > Death_Censor_Age & Death_censor_Any ==1,]
    
    ar[agename, rownames(ar[i, , ]) %in% atrisk$identifier, "out"] =1
  }
  return(ar)
}
