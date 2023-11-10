### here we provide the code for creating the model fit and a sample script for using.
## access to training data from UKB available by request, surbut@mgh.harvard.edu
require(data.table)


fitfunc = function(df_frame, ages, nstates, mode,covariates) {
  nar = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  event = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )

  mean = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  endlist = startlist = vector("list", length = length(nstates))
  names(endlist) = names(startlist) = nstates

  agelist = vector("list", length = length(ages))
  timeatrisk = agelist

  names(agelist) = names(timeatrisk) = ages
  for (i in 1:length(agelist)) {
    agelist[[i]] = startlist
    for (j in 1:length(startlist)) {
      agelist[[i]][[j]] = endlist
    }
  }

  for (i in 1:length(timeatrisk)) {
    timeatrisk[[i]] = startlist
  }






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
                        age < Death_Censor_Age, ]
    atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                 atrisk$statin_age <= nx, 1, 0)
    atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                  atrisk$htn_age <= nx, 1, 0)

    NAR = dim(atrisk)[1]

    #atrisk=atrisk%>%mutate(yearsinstate=age) all same years in state

    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Health"]] = NA

      agelist[[agename]][["Health"]][["Health"]] = NA

      agelist[[agename]][["death"]][["Health"]] = NA

      agelist[[agename]][["Ht"]][["Health"]] = NA


      agelist[[agename]][["HyperLip"]][["Health"]] = NA

      agelist[[agename]][["Dm"]][["Health"]] = NA
    } else{

      my_list <- covariates

      # Remove the word "banana"
      new_list <- gsub("\\+yearsinstate", "", my_list)

      ## from Health (1) to Health
      censored = dim(atrisk[which(
        Cad_0_censor_age > nx &
          nx < Ht_0_censor_age &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]

      nar[i, "Health", 1] = NAR
      event[i, "Health", 1] = censored
      mean[i, "Health", 1] = censored / NAR
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age > nx &
            nx < Ht_0_censor_age &
            nx < HyperLip_0_censor_age &
            nx < Dm_0_censor_age & nx < Death_Censor_Age,
          1,
          0
        ) ~", new_list)),
        data = atrisk
      )
      agelist[[agename]][["Health"]][["Health"]] =summary(fit2)$coefficients

      rm(censored)

      ## from Health (1) to CAD (or CAD death) directly alone

      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 #&
                                    #nx < Death_Censor_Age
                                    ), ])[1]
      nar[i, "Cad", 1] = NAR
      event[i, "Cad", 1] = censored
      mean[i, "Cad", 1] = censored / NAR

      rm(censored)
      # fit2 = glm(
      #   family = mode,as.formula(paste0("
      #   ifelse(
      #     Cad_0_censor_age <= nx &
      #       Cad_0_Any == 2 & nx < Death_Censor_Age,
      #     1,
      #     0
      #   ) ~", new_list)),
      #   ## here age represents time in state
      #   data = atrisk
      # )

      fit2 = glm(
        family = mode,as.formula(paste0("
        ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2,
          1,
          0
        ) ~", new_list)),
        ## here age represents time in state
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Health"]] =summary(fit2)$coefficients

      ## from Health (1) to Death directly alone (but could die with CAD, make CAD and death marginal,
      ## i.e., don't make conditions that other conditions don't develop that year)

      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]

      nar[i, "death", 1] = NAR
      event[i, "death", 1] = censored
      mean[i, "death", 1] = censored / NAR
      rm(censored)
      fit2 = glm(as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0)~",new_list)),
                 family=mode,data=atrisk)
      agelist[[agename]][["Death"]][["Health"]]=summary(fit2)$coefficients
      #rm(atrisk)



    ### from health to HT

      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx & Ht_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      nar[i, "Ht", 1] = NAR
      event[i, "Ht", 1] = censored
      mean[i, "Ht", 1] = censored / NAR
      rm(censored)

      ### from health to HyperLip


      censored = dim(atrisk[which(
        HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Ht_0_censor_age & nx < Dm_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]

      nar[i, "HyperLip", 1] = NAR
      event[i, "HyperLip", 1] = censored
      mean[i, "HyperLip", 1] = censored / NAR
      rm(censored)


      ## from health to Dm

      censored = dim(atrisk[which(
        Dm_0_censor_age <= nx & Dm_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Ht_0_censor_age & nx < HyperLip_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]

      nar[i, "Dm", 1] = NAR
      event[i, "Dm", 1] = censored
      mean[i, "Dm", 1] = censored / NAR
      rm(censored)
      rm(atrisk)
    }

    ## From one risk states to 2,3,4,CAD or death

    # From HT
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]
    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Ht"]] = NA

      agelist[[agename]][["Ht"]][["Ht"]] = NA

      agelist[[agename]][["death"]][["Ht"]] = NA

      agelist[[agename]][["Ht&HyperLip"]][["Ht"]] = NA

      agelist[[agename]][["Ht&Dm"]][["Ht"]] = NA

      agelist[[agename]][["Ht&HyperLip&Dm"]][["Ht"]] = NA



      #se_high[i,"death",2]=NA
    } else{

      atrisk = atrisk %>% mutate(yearsinstate = (age - Ht_0_censor_age))
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)

      timeatrisk[[agename]][["Ht"]] = atrisk$yearsinstate

      ## censored for CAD

      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 2] = NAR

      event[i, "Cad", 2] = censored
      mean[i, "Cad", 2] = censored / NAR
      rm(censored)


      fit2 = glm(
        family = mode,
        ## in order to parse the formulat
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 & nx < Death_Censor_Age,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Ht"]] = summary(fit2)$coefficients

      ## from HT (2) to Ht directly alone
      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx & Ht_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      nar[i, "Ht", 2] = NAR
      event[i, "Ht", 2] = censored
      mean[i, "Ht", 2] = censored / NAR
      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Ht_0_censor_age <= nx & Ht_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < HyperLip_0_censor_age &
            nx < Dm_0_censor_age & nx < Death_Censor_Age,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )

      agelist[[agename]][["Ht"]][["Ht"]] =summary(fit2)$coefficients
      #

      ## from HT (2) to Death directly alone (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any  == 2 &
                                    Ht_0_censor_age < nx &
                                    Ht_0_Any == 2), ])[1]
      nar[i, "death", 2] = NAR
      event[i, "death", 2] = censored
      mean[i, "death", 2] = censored / NAR
      rm(censored)


      fit2 = glm(as.formula(paste0("ifelse(Death_Censor_Age <= nx &Death_censor_Any == 2, 1, 0) ~",covariates)),
        family=mode,data=atrisk)
      agelist[[agename]][["death"]][["Ht"]] =summary(fit2)$coefficients

      ## from HT (2) to Ht and Hyperlip directly
      censored = dim(atrisk[which(HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                   Cad_0_censor_age > nx &
                                    nx < Dm_0_censor_age), ])[1]
      nar[i, "Ht&HyperLip", 2] = NAR
      event[i, "Ht&HyperLip", 2] = censored
      mean[i, "Ht&HyperLip", 2] = censored / NAR

      ## from HT (2) to Ht and Dm directly
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < HyperLip_0_censor_age), ])[1]
      nar[i, "Ht&Dm", 2] = NAR
      event[i, "Ht&Dm", 2] = censored
      mean[i, "Ht&Dm", 2] = censored / NAR

      ## from HT (2) to Ht and Dm and Hyperip directly
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i,"Ht&HyperLip&Dm", 2] = NAR ## not recording
      event[i,"Ht&HyperLip&Dm", 2] = censored
      mean[i,"Ht&HyperLip&Dm", 2] = censored / NAR


      rm(atrisk)

    }



    ##  from Hyperlip
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]





    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["HyperLip"]] = NA

      agelist[[agename]][["HyperLip"]][["HyperLip"]] = NA

      agelist[[agename]][["death"]][["HyperLip"]] = NA

      agelist[[agename]][["HyperLip&Dm"]][["HyperLip"]] = NA

      agelist[[agename]][["Ht&HyperLip"]][["HyperLip"]] = NA

      agelist[[agename]][["Ht&HyperLip&Dm"]][["HyperLip"]] = NA

    } else{
      ## from Hyperlip (3) to CAD marginal (i.e., could develop other states in the interim year)
      atrisk = atrisk %>% mutate(yearsinstate = age - HyperLip_0_censor_age)
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)

      timeatrisk[[agename]][["HyperLip"]] = atrisk$yearsinstate

      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2
                                  & nx < Death_Censor_Age), ])[1]


      nar[i, "Cad", 3] = NAR
      event[i, "Cad", 3] = censored
      mean[i, "Cad", 3] = censored / NAR

      rm(censored)


      fit2 = glm(
        as.formula(paste0("ifelse(Cad_0_censor_age <= nx &Cad_0_Any == 2, 1, 0) ~",
          covariates)),
        data = atrisk,family=mode)

      agelist[[agename]][["Cad"]][["HyperLip"]] =summary(fit2)$coefficients


      ## from Hyperlip (2) to Hyperlip directly alone


      censored = dim(atrisk[which(
        HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Ht_0_censor_age & nx < Dm_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]

      nar[i, "HyperLip", 3] = NAR
      event[i, "HyperLip", 3] = censored
      mean[i, "HyperLip", 3] = censored / NAR
      rm(censored)


      fit2 = glm(
        family = mode,as.formula(paste0("
        ifelse(
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          & Cad_0_censor_age > nx &
            nx < Ht_0_censor_age &
            nx < Dm_0_censor_age &
            nx < Death_Censor_Age,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )

      agelist[[agename]][["HyperLip"]][["HyperLip"]] =summary(fit2)$coefficients

      ## from HyperLip (2) to Death directly alone (but could die with CAD)

      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      nar[i, "death", 3] = NAR
      event[i, "death", 3] = censored
      mean[i, "death", 3] = censored / NAR
      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0)~",covariates)),
        data = atrisk
      )

      agelist[[agename]][["death"]][["HyperLip"]] =summary(fit2)$coefficients



      ## now do the uninteresting ones
      ## from HyperLip (2) to Ht and Hyperlip directly
      censored = dim(atrisk[which(Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Dm_0_censor_age), ])[1]

      nar[i, "Ht&HyperLip", "HyperLip"] = NAR
      event[i, "Ht&HyperLip", "HyperLip"] = censored
      mean[i, "Ht&HyperLip", "HyperLip"] = censored / NAR

      ## from Hyperlip (2) to Hyperlip and Dm directly
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Ht_0_censor_age), ])[1]

      nar[i, "HyperLip&Dm", "HyperLip"] = NAR
      event[i, "HyperLip&Dm", "HyperLip"] = censored
      mean[i, "HyperLip&Dm", "HyperLip"] = censored / NAR

      ## from Hyperlip to Ht and Dm and Hyperip directly
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i, "Ht&HyperLip&Dm", "HyperLip"] = NAR
      event[i, "Ht&HyperLip&Dm", "HyperLip"] = censored
      mean[i, "Ht&HyperLip&Dm", "HyperLip"] = censored / NAR




      rm(atrisk)

    }



    ## from Dm2

    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age, ]
    NAR = dim(atrisk)[1]


    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Dm"]] = NA

      agelist[[agename]][["Dm"]][["Dm"]] = NA

      agelist[[agename]][["death"]][["Dm"]] = NA

      agelist[[agename]][["Ht&Dm"]][["Dm"]] = NA

      agelist[[agename]][["HyperLip&Dm"]][["Dm"]] = NA

      agelist[[agename]][["Ht&HyperLip&Dm"]][["Dm"]] = NA


    } else{
      # from Dm (4) to CAD directly alone

      atrisk = atrisk %>% mutate(yearsinstate = (age - Dm_0_censor_age))
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)

      timeatrisk[[agename]][["Dm"]] = atrisk$yearsinstate

      censored = dim(atrisk[which(
        Cad_0_censor_age <= nx & Cad_0_Any == 2 &
          Dm_0_censor_age <= nx & Dm_0_Any == 2
        #&nx<HyperLip_0_censor_age&nx<Ht_0_censor_age
        & nx < Death_Censor_Age
      ), ])[1]

      nar[i, "Cad", 4] = NAR
      event[i, "Cad", 4] = censored
      mean[i, "Cad", 4] = censored / NAR

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 &
            Dm_0_censor_age <= nx & Dm_0_Any == 2
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )

      agelist[[agename]][["Cad"]][["Dm"]] =summary(fit2)$coefficients



      rm(censored)

      ## from DM (4) to Dm directly alone
      censored = dim(atrisk[which(
        Dm_0_censor_age <= nx &
          Dm_0_Any == 2 &
          Cad_0_censor_age > nx &
          nx < HyperLip_0_censor_age &
          nx < Ht_0_censor_age & nx < Death_Censor_Age
      ), ])[1]

      nar[i, "Dm", 4] = NAR
      event[i, "Dm", 4] = censored
      mean[i, "Dm", 4] = censored / NAR
      agelist[[agename]][["Dm"]][["Dm"]] =summary(fit2)$coefficients

      #rm(fit)
      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Dm_0_censor_age <= nx &
            Dm_0_Any == 2 &
            Cad_0_censor_age > nx &
            nx < HyperLip_0_censor_age & nx < Ht_0_censor_age &
            nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )

      agelist[[agename]][["Dm"]][["Dm"]] =summary(fit2)$coefficients

      ## from Dm to Death directly alone (but could die with CAD)

      censored = dim(atrisk[which(
        Death_Censor_Age <= nx &
          Death_censor_Any == 2), ])[1]

      nar[i, "death", 4] = NAR
      event[i, "death", 4] = censored
      mean[i, "death", 4] = censored / NAR
      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Dm"]] =summary(fit2)$coefficients


     ## from Dm (4) to Ht and Dm directly
      censored = dim(atrisk[which(Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < HyperLip_0_censor_age), ])[1]
      nar[i, "Ht&Dm", "Dm"] = NAR
      event[i, "Ht&Dm", "Dm"] = censored
      mean[i, "Ht&Dm", "Dm"] = censored / NAR

      ## from Dm (4) to Hyperlip and Dm directly
      censored = dim(atrisk[which(HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Ht_0_censor_age), ])[1]
      nar[i, "HyperLip&Dm",4] = NAR
      event[i, "HyperLip&Dm", 4] = censored
      mean[i, "HyperLip&Dm", 4] = censored / NAR

      ## from Dm (2) to Ht and Dm and Hyperip directly
      censored = dim(atrisk[which(HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                    Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i, "Ht&HyperLip&Dm", 4] = NAR
      event[i, "Ht&HyperLip&Dm", 4] = censored
      mean[i, "Ht&HyperLip&Dm", 4] = censored / NAR




      rm(atrisk)
    }

    ## from CAD to death (to do, should we add RF plus CAD as a starting RF? (no, because that would add the complexity of having to add RF + CAD as an ending state and right now we consider CAD as an 'absorbing' state i.e., if you are diagnosed with an additional RF in the interin  year between diagnosis, not considered, ending with CAD is the same as ending with CAD + ...))
    ## from CAD marginal to Death directly

    atrisk = df_frame[age > Cad_0_censor_age &
                        Cad_0_Any == 2 & Death_Censor_Age > age, ]
    NAR = dim(atrisk)[1]

    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Cad"]] = NA



      agelist[[agename]][["death"]][["Cad"]] = NA
    } else{
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]

      atrisk = atrisk %>% mutate(yearsinstate = age - Cad_0_censor_age)
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)

      timeatrisk[[agename]][["Cad"]] = atrisk$yearsinstate

      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      # ,1,0)~cad.prs.lev,data = atrisk)

      nar[i, "death", 5] = NAR
      event[i, "death", 5] = censored
      mean[i, "death", 5] = censored / NAR
      rm(censored)
      #rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )


      censored = dim(atrisk[which(Death_Censor_Age > nx &
                                    Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2), ])[1]

      nar[i, "Cad", 5] = NAR
      event[i, "Cad", 5] = censored
      mean[i, "Cad", 5] = censored / NAR
      agelist[[agename]][["Cad"]][["Cad"]] =summary(fit2)$coefficients


      rm(censored)
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Death_Censor_Age > nx &
            Cad_0_censor_age <= nx & Cad_0_Any == 2,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Cad"]] =summary(fit2)$coefficients


      rm(atrisk)
    }

    ### 2 risk states
    ## from Ht and DM to CAD marginal

    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age & Dm_0_Any == 2 &
                        age > Ht_0_censor_age & Ht_0_Any == 2
                      & age < HyperLip_0_censor_age, ]
    NAR = dim(atrisk)[1]

    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["Ht&Dm"]] = NA

      agelist[[agename]][["Ht&Dm"]][["Ht&Dm"]] = NA
      agelist[[agename]][["death"]][["Ht&Dm"]] = NA


    } else{
      atrisk$yearsinstate = apply(atrisk[, c("Ht_0_censor_age", "Dm_0_censor_age")], 1, function(x) {
        age - max(x[1], x[2])
      })
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)


      timeatrisk[[agename]][["Ht&Dm"]] = atrisk$yearsinstate

      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 9] = NAR
      event[i, "Cad", 9] = censored
      mean[i, "Cad", 9] = censored / NAR
      rm(censored)
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Death_Censor_Age > nx & Cad_0_censor_age <= nx & Cad_0_Any == 2
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Ht&Dm"]] =summary(fit2)$coefficients
      ## from Ht and Dm to Death directly alone (but could die with CAD)

      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      #                ,1,0)~cad.prs.lev,data = atrisk)
      nar[i, "death", 9] = NAR
      event[i, "death", 9] = censored
      mean[i, "death", 9] = censored / NAR


      rm(censored)
      #rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Ht&Dm"]] =summary(fit2)$coefficients


      ## from Ht and Dm to Ht and Dm directly

      censored = dim(atrisk[which(
        nx < Cad_0_censor_age &
          nx >= Dm_0_censor_age & Dm_0_Any == 2 &
          nx >= Ht_0_censor_age & Ht_0_Any == 2
        & nx < HyperLip_0_censor_age
      ), ])[1]

      nar[i, "Ht&Dm", 9] = NAR
      event[i, "Ht&Dm", 9] = censored
      mean[i, "Ht&Dm", 9] = censored / NAR
      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          nx < Cad_0_censor_age &
            nx >= Dm_0_censor_age & Dm_0_Any == 2 &
            nx >= Ht_0_censor_age & Ht_0_Any == 2
          & nx < HyperLip_0_censor_age
          ,
          1,
          0
        ) ~", covariates)),
        data = atrisk
      )

      agelist[[agename]][["Ht&Dm"]][["Ht&Dm"]] =summary(fit2)$coefficients
      rm(atrisk)

    }



    ## from HyperLip and DM to CAD marginal
    ## from HyperLip and Dm to Death directly alone (but could die with CAD) ###### BAD

    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age, ]
    NAR = dim(atrisk)[1]

    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["HyperLip&Dm"]] = NA

      agelist[[agename]][["HyperLip&Dm"]][["HyperLip&Dm"]] = NA
      agelist[[agename]][["death"]][["HyperLip&Dm"]] = NA


    } else{
      atrisk$yearsinstate = apply(atrisk[, c("HyperLip_0_censor_age", "Dm_0_censor_age")], 1, function(x) {
        age - max(x[1], x[2])
      })
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)

      timeatrisk[[agename]][["HyperLip&Dm"]] = atrisk$yearsinstate

      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]

      nar[i, "death", 8] = NAR
      event[i, "death", 8] = censored
      mean[i, "death", 8] = censored / NAR

      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["HyperLip&Dm"]] =summary(fit2)$coefficients


      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    Death_Censor_Age > nx), ])[1]
      nar[i, "Cad", 8] = NAR
      event[i, "Cad", 8] = censored
      mean[i, "Cad", 8] = censored / NAR

      rm(censored)
      #rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 & Death_Censor_Age > nx
          ,
          1,
          0
        ) ~",covariates)),data = atrisk
      )
      agelist[[agename]][["Cad"]][["HyperLip&Dm"]] =summary(fit2)$coefficients


      censored = dim(atrisk[which(
        nx < Cad_0_censor_age &
          nx > Dm_0_censor_age &
          Dm_0_Any == 2 &
          nx > HyperLip_0_censor_age &
          HyperLip_0_Any == 2 & nx < Ht_0_censor_age
      ), ])[1]
      nar[i, "HyperLip&Dm", 8] = NAR
      event[i, "HyperLip&Dm", 8] = censored
      mean[i, "HyperLip&Dm", 8] = censored / NAR

      #rm(fit)

      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          nx < Cad_0_censor_age &
            nx > Dm_0_censor_age &
            Dm_0_Any == 2 &
            nx > HyperLip_0_censor_age &
            HyperLip_0_Any == 2 & nx < Ht_0_censor_age
          ,
          1,
          0
        ) ~", covariates)),
        data = atrisk
      )
      agelist[[agename]][["HyperLip&Dm"]][["HyperLip&Dm"]] =summary(fit2)$coefficients
      rm(atrisk)

    }



    ## from HyperLip and Ht to CAD marginal

    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2
                      & age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]


    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["Ht&HyperLip"]] = NA

      agelist[[agename]][["Ht&HyperLip"]][["Ht&HyperLip"]] = NA
      agelist[[agename]][["death"]][["Ht&HyperLip"]] = NA

    } else{
      atrisk$yearsinstate = apply(atrisk[, c("HyperLip_0_censor_age", "Ht_0_censor_age")], 1, function(x) {
        age - max(x[1], x[2])
      })
      timeatrisk[[agename]][["Ht&HyperLip"]] = atrisk$yearsinstate
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)


      censored = dim(atrisk[which(
        Cad_0_censor_age <= nx &
          Cad_0_Any == 2 &
          Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx &
          HyperLip_0_Any == 2 & nx < Death_Censor_Age
      ), ])[1]




      nar[i, "Cad", 7] = NAR
      event[i, "Cad", 7] = censored
      mean[i, "Cad", 7] = censored / NAR


      rm(censored)
      ##rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 &
            Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2 &
            nx < Death_Censor_Age,1,0) ~",covariates)),
        data = atrisk
      )


      agelist[[agename]][["Cad"]][["Ht&HyperLip"]] =summary(fit2)$coefficients


      ## from  HyperLip and Ht to HyperLip and Ht

      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]



      nar[i, "Ht&HyperLip", 7] = NAR
      event[i, "Ht&HyperLip", 7] = censored
      mean[i, "Ht&HyperLip", 7] = censored / NAR

      rm(censored)
      #rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < Dm_0_censor_age & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )



      agelist[[agename]][["Ht&HyperLip"]][["Ht&HyperLip"]] =summary(fit2)$coefficients
      ## from HyperLip and HT to Death directly alone (but could die with CAD)

      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]


      nar[i, "death", 7] = NAR
      event[i, "death", 7] = censored
      mean[i, "death", 7] = censored / NAR


      rm(censored)
      #rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Ht&HyperLip"]] =summary(fit2)$coefficients
      rm(atrisk)

    }


    ## three risk states
    ## from HyperLip and Ht $ Dmto CAD marginal


    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2
                      & age > Dm_0_censor_age & Dm_0_Any == 2, ]
    NAR = dim(atrisk)[1]


    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["Ht&HyperLip&Dm"]] = NA


      agelist[[agename]][["Ht&HyperLip&Dm"]][["Ht&HyperLip&Dm"]] = NA

      agelist[[agename]][["death"]][["death"]] = NA



    } else{
      atrisk$yearsinstate = apply(atrisk[, c("HyperLip_0_censor_age",
                                             "Ht_0_censor_age",
                                             "Dm_0_censor_age")], 1, function(x) {
                                               age - max(x[3], max(x[1], x[2]))
                                             })
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)

      timeatrisk[[agename]][["Ht&HyperLip&Dm"]] = atrisk$yearsinstate

      censored = nrow(atrisk[which(Cad_0_censor_age <= nx &
                                     Cad_0_Any == 2 &
                                     nx < Death_Censor_Age), ])


      nar[i, "Cad", 10] = NAR
      event[i, "Cad", 10] = censored
      mean[i, "Cad", 10] = censored / NAR
      #rm(fit)
      rm(censored)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )

      agelist[[agename]][["Cad"]][["Ht&HyperLip&Dm"]] =summary(fit2)$coefficients


      # STAY
      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          Dm_0_censor_age <= nx &
          Dm_0_Any == 2 & nx < Death_Censor_Age
      ), ])[1]


      nar[i, "Ht&HyperLip&Dm", 10] = NAR
      event[i, "Ht&HyperLip&Dm", 10] = censored
      mean[i, "Ht&HyperLip&Dm", 10] = censored / NAR


      rm(censored)
      #rm(fit)

      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          &
            Cad_0_censor_age > nx &
            Dm_0_censor_age <= nx &
            Dm_0_Any == 2 & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Ht&HyperLip&Dm"]][["Ht&HyperLip&Dm"]] =summary(fit2)$coefficients


      ## from all 3 to death (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]

      nar[i, "death", 10] = NAR
      event[i, "death", 10] = censored
      mean[i, "death", 10] = censored / NAR


      rm(censored)
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Ht&HyperLip&Dm"]] =summary(fit2)$coefficients


    }

    ##Dead

    atrisk = df_frame[Death_Censor_Age < age &
                        Death_censor_Any == 2,]

    nar[i,"death","death"]=nrow(atrisk)


   # print(i)

  }


  mylist = list(
    "events" = event,
    "rates" = mean,
    "AR" = nar,
    "model_list" = agelist,
    "yearsinstate" = timeatrisk

  )
  return(mylist)
}

fitfunc2=function(df_frame,ages,nstates,mode,covariates){
  suppressWarnings(fitfunc(df_frame,ages,nstates,mode,covariates))
}


# ### to use
#
#
#
# source("~/multistate2//code/utils.R")
# source("~/multistate2//code/smoothtest.R")
# source("~/multistate2//code/newsmooth.R")
# source("~/multistate2/code/fitarray.R")
# library("reshape2")
# source("~/multistate2/code/arrayindicate.R")
#
# ## training data: ask for access
#
# load("~/Library/CloudStorage/Dropbox-Personal///pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
# train=dfh[1:(nrow(dfh)*0.80),]
# ages=c(40:80)
# nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")
#
# modelfit = fitfunc2(
# data.table(train),
# ages = ages,
# nstates = nstates,
# mode = "binomial",
# covariates = "cad.prs+f.31.0.0+smoke+antihtn_now+statin_now")
#
# ### now compute projection from healthy state for example, for a sample 40 year old
#
# a = coefplotsmooth2(
#   ages = ages,
#   start = nstates[1],
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
# healthy_coefs = a$custom_smooth
#
# ## create sample matrix
#
# intercept=1
# cad.prs=c(-2,-1,0,1,2)
# sex=c(0,1)
# smoke=c(0,1)
# antihtn_now=c(0,1)
# statin_now=c(0,1)
# atrisk=expand.grid(intercept,cad.prs,sex,smoke,antihtn_now,statin_now)
#
#
# mso = compute_prediction_product_matrix(
# coefmat = healthy_coefs,
# atrisk = atrisk,
# agepredinterval = c(40:80)
# )
#
# # remaining lifetime risk matrix starting from each 40
#
# mso$PredictedIntervalrisk
#
# # remaining lifetime risk matrix under treatment starting from 40
#
# mso$Hazard_treated
#
#
#
# ### can be repeated over all states by changing start state and extracting coefficients,
### can be computed for any age by changing age predinterval
#
#

# b = coefplotsmooth2(
#   ages = ages,
#   start = "Ht",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 30,
#   span = 0.75,
#   degree = 2
# )
# ht_coefs = b$custom_smooth
#
# c = coefplotsmooth2(
#   ages = ages,
#   start = "HyperLip",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
# hyperlip_coefs = c$custom_smooth
#
# d = coefplotsmooth2(
#   ages = ages,
#   start = "Dm",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
# dm_coefs = d$custom_smooth
#
#
# e = coefplotsmooth2(
#   ages = ages,
#   start = "Ht&HyperLip",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
# HH_coefs = e$custom_smooth
#
# f = coefplotsmooth2(
#   ages = ages,
#   start = "HyperLip&Dm",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
# HD_coefs = f$custom_smooth
#
# g = coefplotsmooth2(
#   ages = ages,
#   start = "Ht&Dm",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
# TD_coefs = g$custom_smooth
#
# h = coefplotsmooth2(
#   ages = ages,
#   start = "Ht&HyperLip&Dm",
#   stop = "Cad",
#   modelfit = modelfit,
#   window_width = 20,
#   span = 0.75,
#   degree = 2
# )
#
# HHD_coefs = h$custom_smooth
#
# coeflist = list(
#   healthy_coefs,
#   ht_coefs,
#   hyperlip_coefs,
#   dm_coefs,
#   HH_coefs,
#   HD_coefs,
#   TD_coefs,
#   HHD_coefs
# )
#


