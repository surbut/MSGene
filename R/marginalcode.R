
modelfitfun=function(ages,transit_mat,covariates,df,statusarray,mode="binomial"){

  states=colnames(transit_mat)
  results_by_age=list()

  nar = array(
    data = NA,
    dim = c(length(ages), length(states), length(states)),
    dimnames = list(ages, states, states)
  )
  event = array(
    data = NA,
    dim = c(length(ages), length(states), length(states)),
    dimnames = list(ages, states, states)
  )

  mean = array(
    data = NA,
    dim = c(length(ages), length(states), length(states)),
    dimnames = list(ages, states, states)
  )

  # Loop through each age



  for (age in ages) {

    agename=as.character(age)
    nx=age+1
    nxname=as.character(nx)
    s=statusarray
    # Create a new list for this age
    age_list <- list()



    # Loop through each possible ending state, we need to do it this way because the previous modelfit worked this way
    for (end_state in states) {

      end_state_list <- list()
      # Loop through each possible starting state, excluding the two absorbing states
      for (start_state in states[-c(5,6)]) {

        # Create a list for this starting state
        start_state_list <- list()
        ar=which(s[agename,,start_state]==1)
        at_risk = s[agename,ar,]
        if(transit_mat[start_state,end_state] == 1){
          # print(age)
          # print(agename)
          # print(start_state)
          # print(end_state)
          censored=which(s[nxname,ar,end_state]==1)

          ## now for model fit

          fitmat=df[df$identifier%in%names(ar),]
          fitmat$y=0
          fitmat[fitmat$identifier%in%names(censored),"y"]=1
          fitmat$statin_now = ifelse(fitmat$statin == 1 &
                                       fitmat$statin_age <= nx, 1, 0)
          fitmat$antihtn_now = ifelse(fitmat$antihtn == 1 &
                                        fitmat$htn_age <= nx, 1, 0)

          modfit=glm(
            family = mode,
            as.formula(paste0("y~", covariates)),
            data = fitmat
          )
          nar[agename,start_state,end_state]=dim(at_risk)[1]
          mean[agename,start_state,end_state]=length(censored)/dim(at_risk)[1]
          event[agename,start_state,end_state]=length(censored)



          end_state_list[[start_state]] <- summary(modfit)$coefficients
        }
      }

      # Store the 'start_state' list in the 'age_list' using the state's name as the key
      age_list[[end_state]] <- end_state_list

    }
    results_by_age[[as.character(age)]] <- age_list

    # Store each 'age_list' in the 'results_by_age' list using the age as the key

  }

  modelfit=list()
  modelfit$model_list=results_by_age
  modelfit$AR=nar
  modelfit$events=event
  modelfit$rates=mean

  return(modelfit)}



## input transition matrix and
coef_array_func=function(ages,transit_mat,modelfit){
  ncovariates=dim(modelfit$model_list[[1]][[1]][[1]])[1]
  covariates=rownames(modelfit$model_list[[1]][[1]][[1]])
  transit_coef=array(NA,dim=c(length(ages),dim(transit_mat)[1],dim(transit_mat)[2],ncovariates))
  dimnames(transit_coef)=list(ages,rownames(transit_mat),colnames(transit_mat),covariates)

  states=rownames(transit_mat)
  for (start_state in states[-c(5,6)]) {

    # Create a list for this starting state
    start_state_list <- list()

    # Loop through each possible ending state
    for (end_state in states) {

      if(transit_mat[start_state,end_state] == 1){
        # start state age < current age & next_state age < current age + 1


        a = coefplotsmooth2(
          ages = ages,
          start = start_state,
          stop = end_state,
          modelfit = modelfit,
          window_width = 20,
          span = 0.75,
          degree = 2
        )
        transit_coef[,start_state,end_state,]=a$custom_smooth
      }}}
  return(transit_coef)}


generate_single_transition_matrix <- function(age, coeff_array, covariates, absorbing_states = NULL) {
  # Extract the dimensions
  n_states <- dim(coeff_array)[2]  # Assuming state x state matrix

  # Create an empty matrix
  transition_matrix <- matrix(0, nrow = n_states, ncol = n_states)

  # Populate the transition matrix
  for (start_state in 1:n_states) {
    for (end_state in 1:n_states) {
      current_coeffs <- coeff_array[as.character(age), start_state, end_state, ]

      # Check if the length of covariates matches the number of coefficients
      if (length(current_coeffs) != length(covariates)) {
        stop("The length of covariates must match the number of coefficients.")
      }

      if (!all(is.na(current_coeffs))) {  # Check if coefficients are not all NA
        # Calculate the log odds ratio and convert to probability
        logOR <- sum(current_coeffs * covariates)
        probability <- exp(logOR) / (1 + exp(logOR))

        transition_matrix[start_state, end_state] <- probability
      }
    }
  }

  # Set the absorbing states
  if (!is.null(absorbing_states)) {
    for (state in absorbing_states) {
      transition_matrix[state, ] <- 0  # No transitions from an absorbing state
      transition_matrix[state, state] <- 1  # Set self-transition to 1
    }
  }

  # Adjust diagonal to make row sums to 1
  for (row in 1:n_states) {
    if (!row %in% absorbing_states) { # Skip absorbing states
      transition_probs <- transition_matrix[row, ]
      non_absorbing <- which(!is.na(transition_probs) & 1:n_states != row) # Exclude self-transition
      remaining_prob <- 1 - sum(transition_probs[non_absorbing], na.rm = TRUE)
      transition_matrix[row, row] <- max(0, remaining_prob) # Ensure it's not negative due to floating-point precision
    }
  }

  colnames(transition_matrix)=rownames(transition_matrix)=dimnames(coeff_array)[[2]]


  return(transition_matrix)
}


mat.test = function(coefficients, covariates, age_start, age_end, absorbing_states) {
  # Pre-allocate the list with the required size
  matlist = vector("list", length = age_end - age_start + 1)
  names(matlist) = age_start:age_end

  for(age in age_start:age_end) {
    current_matrix <- generate_single_transition_matrix(age = age,
                                                        coeff_array = coefficients,
                                                        covariates = covariates,
                                                        absorbing_states = absorbing_states)
    matlist[[as.character(age)]] = current_matrix
  }

  return(matlist)
}




# # Usage of the function:
# age_specific_matrix <- generate_single_transition_matrix(age = 50, coeff_array, covariates, absorbing_states = c(5, 6))
# print(age_specific_matrix)

calculate_matrix_product_over_interval <- function(coefficients, covariates, age_start, age_end,absorbing_states) {
  # 'age_start' and 'age_end' define the interval of ages.

  # Generate the first transition matrix.
  product_matrix <- generate_single_transition_matrix(age = age_start,coeff_array = coefficients,covariates = covariates,absorbing_states = absorbing_states)

  # Multiply the transition matrices for each subsequent year in the interval.
  for (age in (age_start + 1):(age_end-1)) {
    current_matrix <- generate_single_transition_matrix(age = as.character(age),coeff_array = coefficients,covariates = covariates,absorbing_states = absorbing_states)
    product_matrix <- product_matrix %*% current_matrix  # Matrix multiplication.
  }

  return(product_matrix)
}

