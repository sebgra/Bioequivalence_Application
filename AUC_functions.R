library(ggplot2)

total_AUC <-function(times,concentrations){
  AUC = 0
  for (i in 2:length(times)) {
    AUC = AUC + 0.5*(times[i]-times[i-1])*(concentrations[i]+concentrations[i-1])
    #print(AUC)
  }
  
  df = as.data.frame(cbind(times,concentrations))
  
  #p <- ggplot(data = df, aes(x = times,y = concentrations)) + geom_area(fill = "pink") + geom_point() + ggtitle("Total AUC")
  #print(p)
  return(AUC)
}

# iAUC <- function(times,concentrations){
#   inc_AUC = 0
#   for (i in 2:length(times)) {
#     inc_AUC = inc_AUC + 0.5*((concentrations[i]-concentrations[1])+(concentrations[i-1]-concentrations[1]))*(times[i]-times[i-1])
#   }
#   
#   df = as.data.frame(cbind(times,concentrations))
#   
#     
#     return(inc_AUC)
# }

iAUC <-function(times,concentrations){
  areas = vector()
  
  for (i in 1:(length(concentrations)-1)) {
    #print(i)
    A = 0
    
    if ((concentrations[i]>= concentrations[1])&&(concentrations[i+1]>= concentrations[1])) {
      A = (times[i+1]-times[i])*(concentrations[i+1]+concentrations[i]-2*concentrations[1])*0.5
      areas[i] = A
      
    }
    
    if ((concentrations[i] >= concentrations[1])&&(concentrations[i+1]) < concentrations[1]) {
      A = (((concentrations[i]-concentrations[1])^2) * (times[i+1]-times[i])) / (2*(concentrations[i]-concentrations[i+1]))
      areas[i] = A
      
    }
    
    if ((concentrations[i] < concentrations[1]) && (concentrations[i+1] >= concentrations[1])) {
      
      A = (((concentrations[i+1]-concentrations[1])^2) * (times[i+1]-times[i])) / (2*(concentrations[i+1]-concentrations[i]))
      areas[i] = A
    }
    
    if((concentrations[i] < concentrations[1]) && (concentrations[i+1]<concentrations[1]))
      
      A = 0
    areas[i] = A
    
  }
  return(sum(areas))
  
}

logAUC <-function(times,concentrations){
  logarithmic_AUC = 0
  
  for (i in 2:length(concentrations)) {
    logarithmic_AUC = logarithmic_AUC + ((concentrations[i-1]-concentrations[i])/(log2(concentrations[i-1])-log2(concentrations[i]))*(times[i]-times[i-1]))
  }
  return(logarithmic_AUC)
}

linearlogAUC <-function(times, concentrations){
  linlogAUC = 0
  
  for (i in 2:length(concentrations)) {
    
    if (((concentrations[i]-concentrations[i-1])/(times[i]-times[i-1]))>0) {
      
      # print("increasing")
      linlogAUC = linlogAUC + 0.5*(concentrations[i]+concentrations[i-1])*(times[i]-times[i-1])
      
    }
    else if (((concentrations[i]-concentrations[i-1])/(times[i]-times[i-1]))<0) {
      # print("decreasing")
      linlogAUC = linlogAUC + (((concentrations[i-1]-concentrations[i])/(log2(concentrations[i-1])-log2(concentrations[i])))*(times[i]-times[i-1]))
      
    }
  }
  return(linlogAUC)
}


iAUCcut<-function(times,concentrations){
  cut_iAUC = 0
  i = 1
  
  first_variation = (concentrations[2]-concentrations[1])/(times[2]-times[1])
  
  if (first_variation < 0) {
    
    # print("Variation negative")
    #return("Variation negative")
    return(0)
    
  }
  
  else if(first_variation > 0){
    
    while ((concentrations[i+1]-concentrations[i])/(times[i+1]-times[i])>0) {
      i = i+1
    }
    
    times_subset = times[1:i]
    concentrations_subset = concentrations[1:i]
    cut_iAUC = iAUC(times_subset,concentrations_subset)
    # print(i)
    return(cut_iAUC)
    
  }
}

cut_tail_AUC<-function(times,concentrations){
  i = 1
  ct_AUC = 0
  
  first_variation = (concentrations[2]-concentrations[1])/(times[2]-times[1])
  
  if (first_variation < 0) {
    
    return("Variation negative")
    
  }
  
  else if(first_variation > 0){
    
    while ((concentrations[i+1]-concentrations[i])/(times[i+1]-times[i])>0) {
      i = i+1
    }
    
    times_subset = times[1:i]
    concentrations_subset = concentrations[1:i]
    ct_AUC = total_AUC(times_subset,concentrations_subset)
    # print(i)
    return(ct_AUC)
    
  }
  
  return(ct_AUC)
}

iAUCmin<-function(times,concentrations){
  
  mini_iAUC = 0
  baseline_position = which.min(concentrations)
  
  mini_iAUC = total_AUC(times,concentrations) - (times[length(times)]*concentrations[baseline_position])
  
  return(mini_iAUC)
}

Simpsons <-function(times,concentrations){ # integrate even number forcing ?
  
  if (length(times) %% 2 == 1 ) {
    return("Number of times and concentrations must be even ")
  }
  
  Area = 0
  concentrations_copy = concentrations
  first = concentrations_copy[1]
  last = concentrations_copy[length(concentrations_copy)]
  concentrations_copy = concentrations_copy[-1]
  concentrations_copy = concentrations_copy[-length(concentrations_copy)]
  even <- vector()
  odd <- vector()
  
  for (i in concentrations_copy) {
    
    if (i %% 2 == 0){
      
      even <- append(even,i)
    }
    
    else {
      odd <- append(odd,i)
    }
  }
  
  Area = ((times[2]-times[1])/3)*(first + last + 4*sum(odd)+2*sum(even))
  return(Area)
}

LPS_AUC<-function(times,concentrations){ ## need to be generalised to n measurement points
  
  estimated_AUC = 0
  degree = length(concentrations) - 1
  fit = lm(concentrations~poly(times,degree,raw=TRUE))
  coeff = fit$coefficients
  
  
  integrand <-function(x){coeff[1]+ coeff[2]*x + coeff[3]*(x^2)+coeff[4]*(x^3)+coeff[5]*(x^4)+coeff[6]*(x^5)+coeff[7]*(x^6)}
  
  LPS_AUC = integrate(integrand, lower = 0, upper = times[length(times)])
  return(LPS_AUC$value)
  
  
}

logScaleAUC <- function(times,concentrations){
  LS_AUC = 0
  LS_AUC = total_AUC(log10(times + 1), concentrations)
  return(LS_AUC)
}

standard_scale_AUC<-function(concentrations){
  sd_scale_AUC = 0
  sd_scale_AUC = total_AUC(c(1:length(concentrations)), concentrations)
  return(sd_scale_AUC)
  
}


first_cut_iAUC <-function(times,concentrations){
  
  if (length(which(concentrations<concentrations[1]))==0) {
    return(iAUC(times,concentrations))
    
  }
  baseline = concentrations[1]
  flag = TRUE
  
  for (i in c(1:length(concentrations)-1)) {
    
    if (concentrations[i+1] < baseline && flag == TRUE) {
      
      a = ( (concentrations[i]-concentrations[i+1])/ (times[i]-times[i+1]) )
      b = concentrations[i] - a*(times[i])
      
      position_to_cut = i
      
      time_to_cross = ( (concentrations[1]-b) / a)
      flag = FALSE
    }
    
  }
  
  new_times = c(times[c(1:position_to_cut)],time_to_cross)
  
  new_concentrations = c(concentrations[c(1:position_to_cut)],concentrations[1])
  # print(new_times)
  # print(new_concentrations)
  
  Area = iAUC(new_times,new_concentrations)
  
  return(Area)
}




cut_different_phases <-function(times,concentrations){
  
  
  
  baseline = concentrations[1]
  flag = TRUE
  flag_triphasic = FALSE
  
  for (i in c(1:length(concentrations)-1)) {
    
    if (concentrations[i+1] < baseline && flag == TRUE) {
      
      a = (concentrations[i] - concentrations[i+1]) / (times[i]-times[i+1])
      
      b = concentrations[i] - a*times[i]
      
      position_to_cut = i
      
      times_to_cross = (baseline - b) / a
      
      flag = FALSE
      
    }
    
    
  }
  
  times_phase_1 = c(times[c(1:position_to_cut)],times_to_cross)
  concentrations_phase_1 = c(concentrations[c(1:position_to_cut)],baseline)
  
  # print(times_phase_1)
  # print(concentrations_phase_1)
  
  times_phase_2 = c(times_to_cross,times[-c(1:position_to_cut)])
  concentrations_phase_2 = c(baseline,concentrations[-c(1:position_to_cut)])
  
  # times_phase_2 = c(times[-c(1:position_to_cut)])
  # concentrations_phase_2 = c(concentrations[-c(1:position_to_cut)])
  
  # print(times_phase_2)
  # print(concentrations_phase_2)
  
  
  flag_phase_2 = TRUE
  
  for (j in c(2:length(concentrations_phase_2)-1)) {
    
    if (concentrations_phase_2[j+1] > baseline && flag_phase_2 ==TRUE) {
      
      c = (concentrations_phase_2[j]-concentrations_phase_2[j+1]) / (times_phase_2[j]-times_phase_2[j+1])
      d = concentrations_phase_2[j] - c*times_phase_2[j]
      
      position_to_cut_2 = j
      
      times_to_cross_2 = (concentrations_phase_2[1] - d) / c
      # print("time to cross 2")
      # print(times_to_cross_2)
      
      flag_phase_2 = FALSE
      
      flag_triphasic = TRUE
      
    }
    
  }
  
  if (flag_triphasic == TRUE) {
    # print(times_to_cross_2)
    
    times_phase_inter = c(times_phase_2[c(1:position_to_cut_2)],times_to_cross_2)
    concentrations_phase_inter = c(concentrations_phase_2[c(1:position_to_cut_2)],baseline)
    
    # print("temps phase inter")
    # print(times_phase_inter)
    # print("concentrations phase inter")
    # print(concentrations_phase_inter)
    
    
    times_phase_finale = c(times_to_cross_2,times_phase_2[-c(1:position_to_cut_2)])
    concentrations_phase_finale = c(baseline,concentrations_phase_2[-c(1:position_to_cut_2)])
    
    # print(times_phase_finale)
    # print(concentrations_phase_finale)
    
    phase_beg = data.frame(times_phase_1,concentrations_phase_1)
    phase_inter = data.frame(times_phase_inter,concentrations_phase_inter)
    phase_finale = data.frame(times_phase_finale,concentrations_phase_finale)
    
    colnames(phase_beg) = c("Time","Concentrations")
    colnames(phase_inter) = c("Time","Concentrations")
    colnames(phase_finale) = c("Time","Concentrations")
    
    list_df = list(phase_beg,phase_inter,phase_finale)
    
    # print("test")
    # print(length(list_df))
    
    return(list_df)
    
  }
  
  if (flag_triphasic == FALSE) {
    
    phase_beg = data.frame(times_phase_1,concentrations_phase_1)
    phase_finale = data.frame(times_phase_2,concentrations_phase_2)
    
    colnames(phase_beg) = c("Time","Concentrations")
    colnames(phase_finale) = c("Time","Concentrations") 
    list_df = list(phase_beg,phase_finale)
    return(list_df)
    
  }
  
}



reverse_area <- function(times,concentrations){
  
  baseline = concentrations[1]
  reverse_concentration = vector()
  
  for (i in c(1:length(concentrations))) {
    
    reverse_concentration[i] = baseline + abs(baseline-concentrations[i])
  }
  # print(reverse_concentration)
  return(reverse_concentration)
}


net_iAUC <-function(times,concentrations){
  
  if (length(which(concentrations<concentrations[1]))==0) {
    return(iAUC(times,concentrations))
    
  }
  
  global = cut_different_phases(times,concentrations)
  
  
  if (length(global) == 3) {
    
    P_1_times = global[[1]]$Time
    P_2_times = global[[2]]$Time
    P_3_times = global[[3]]$Time
    
    P_1_concentrations = global[[1]]$Concentrations
    P_2_concentrations = global[[2]]$Concentrations
    P_3_concentrations = global[[3]]$Concentrations
    
    
    
    P_2_concentrations_reversed = reverse_area(P_2_times,P_2_concentrations)
    # print(P_2_concentrations_reversed)
    
    Area = iAUC(P_1_times,P_1_concentrations) + iAUC(P_3_times,P_3_concentrations) - iAUC(P_2_times,P_2_concentrations_reversed)
    
    return(Area)
    
  }
  
  if (length(global) ==2) {
    
    P_1_times = global[[1]]$Time
    P_2_times = global[[2]]$Time
    
    P_1_concentrations = global[[1]]$Concentrations
    P_2_concentrations = global[[2]]$Concentrations
    
    
    
    P_2_concentrations_reversed = reverse_area(P_2_times,P_2_concentrations)
    # print(P_2_concentrations_reversed)
    
    Area = iAUC(P_1_times,P_1_concentrations) - iAUC(P_2_times,P_2_concentrations_reversed)
    
    return(Area)
    
  }
  
  
}
