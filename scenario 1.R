library(fBasics)
library(expm)
library(IMIFA)
library(gtools)
library(ggplot2)
library(plotrix)
library(Cairo)
library(vegan)
library(faux)


source("resistant_codes_LMS_LTS.txt")


set.seed(9)
rho_err = 0.4
error_frac = 0.1
population_size = 1000
sd_P = 1
Population_P = rnorm_multi(population_size, 2, 0, 1, rho_err,empirical = TRUE)
Population_P = as.matrix(Population_P)
Population_Q = Population_P
error_matrix = matrix(rnorm(population_size*2,mean=0,sd=error_frac*sd_P),nrow=population_size,ncol=2)
Population_P = Population_P + error_matrix
total_set = cbind(Population_P,Population_Q)#generate correlated population

n_landmark = 15
total_sim = 100
scale_fix = 2
translation_fix = 1
n_outliers = round(0.13*n_landmark)
#n_outliers = round(0.4*n_landmark)
outlier_factor = -2
reflection_factor = diag(c(-1,1)) #matrix(c(-1, 0, 0, 1), nrow = 2, ncol = 2) # y axis
######################

MM2_save = matrix(NA,nrow=total_sim,ncol=5)           #column 1: 0--no reflection; 1--with reflection
Huber_save = matrix(NA,nrow=total_sim,ncol=5)         #column 2: angle
biweight_save = matrix(NA,nrow=total_sim,ncol=5)      #column 3: scale
procrustes_save = matrix(NA,nrow=total_sim,ncol=5)    #column 4: x translation 
S_est_save = matrix(NA,nrow=total_sim,ncol=5)         #column 5: y translation
LMS_est_save = matrix(NA,nrow=total_sim,ncol=5)
Prof_est_save = matrix(NA,nrow=total_sim,ncol=5)
LMS_LTS_est_save = matrix(NA,nrow=total_sim,ncol=5)


# Res_MM2_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim) 
# Res_Huber_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim)
# Res_Biweight_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim)
# Res_procrustes_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim)
# Res_S_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim)
# Res_LMS_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim)
# Res_Prof_save = matrix(NA,nrow=n_landmark,ncol=2*total_sim)

MM2_ind_error = 0
Huber_ind_error = 0
Biweight_ind_error = 0
Procrustes_ind_error = 0
S_ind_error = 0
LMS_ind_error = 0
Prof_ind_error = 0
LMS_LTS_ind_error = 0

Res_MM2_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_Huber_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_Biweight_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_procrustes_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_S_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_LMS_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_Prof_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)
Res_LMS_LTS_dis_save = matrix(NA,nrow=n_landmark,ncol=total_sim)


ind_save = matrix(NA,nrow=total_sim,ncol=n_outliers)

# Pearson_matrix_1 = matrix(NA,nrow=total_sim,ncol=21)
# colnames(Pearson_matrix_1) = c("MM2/HB","MM2/BW","MM2/Procrustes","MM2/S","MM2/LMS","MM2/RS",
#                                "HB/BW","HB/Procrustes","HB/S","HB/LMS","HB/RS",
#                                "BW/Procrustes","BW/S","BW/LMS","BW/RS",
#                                "Procrustes/S","Procrustes/LMS","Procrustes/RS",
#                                "S/LMS","S/RS",
#                                "LMS/RS")
# Pearson_matrix_2 = matrix(NA,nrow=total_sim,ncol=21)
# colnames(Pearson_matrix_2) = c("MM2/HB","MM2/BW","MM2/Procrustes","MM2/S","MM2/LMS","MM2/RS",
#                                "HB/BW","HB/Procrustes","HB/S","HB/LMS","HB/RS",
#                                "BW/Procrustes","BW/S","BW/LMS","BW/RS",
#                                "Procrustes/S","Procrustes/LMS","Procrustes/RS",
#                                "S/LMS","S/RS",
#                                "LMS/RS")
# 
# Spearman_matrix_1 = matrix(NA,nrow=total_sim,ncol=21)
# colnames(Spearman_matrix_1) = c("MM2/HB","MM2/BW","MM2/Procrustes","MM2/S","MM2/LMS","MM2/RS",
#                                 "HB/BW","HB/Procrustes","HB/S","HB/LMS","HB/RS",
#                                 "BW/Procrustes","BW/S","BW/LMS","BW/RS",
#                                 "Procrustes/S","Procrustes/LMS","Procrustes/RS",
#                                 "S/LMS","S/RS",
#                                 "LMS/RS")
# 
# Spearman_matrix_2 = matrix(NA,nrow=total_sim,ncol=21)
# colnames(Spearman_matrix_2) = c("MM2/HB","MM2/BW","MM2/Procrustes","MM2/S","MM2/LMS","MM2/RS",
#                                 "HB/BW","HB/Procrustes","HB/S","HB/LMS","HB/RS",
#                                 "BW/Procrustes","BW/S","BW/LMS","BW/RS",
#                                 "Procrustes/S","Procrustes/LMS","Procrustes/RS",
#                                 "S/LMS","S/RS",
#                                 "LMS/RS")


k_4 = rep(1,n_landmark)

for (i in 1:total_sim) {
  sub_population = total_set[sample.int(1000,size=n_landmark),] #default replace=FALSE
  P_mat = sub_population[,1:2]
  Q_mat = sub_population[,3:4]
  rot_mat = diag(rep(cos(pi/6),2))
  rot_mat[2,1] = sin(pi/6)
  rot_mat[1,2] = -rot_mat[2,1]
  rot_mat = rot_mat%*%reflection_factor
  Q_mat = Q_mat%*%t(rot_mat)
  ind = sort(sample.int(n_landmark,n_outliers)) 
  P_mat[ind,] = outlier_factor*P_mat[ind,]
  Q_mat = scale_fix*Q_mat
  Q_mat = Q_mat + translation_fix
  
  #######################################################
  ##MM2##
  result_1= MM2_estimator_with_reflection(Q_mat,P_mat) # when the data has reflection, use MM2_estimator_with_reflection
  rot_now = result_1$rot                               # when the data has no reflection, use MM2_estimator_no_reflection
  MM2_save[i,1] = rot_check(rot_now)
  if(MM2_save[i,1]==1){
    rot_now = reflection_factor%*%rot_now
    MM2_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
    MM2_save[i,3] = result_1$scal
    MM2_save[i,4] = result_1$trans_1
    MM2_save[i,5] = result_1$trans_2
  }  
  Res_MM2 = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_1,result_1$trans_2))
  Res_MM2_dis = sqrt(rowSums(Res_MM2*Res_MM2))
  MM2_ind = tail(order(Res_MM2_dis), n = n_outliers)
  if(any(sort(ind)!=sort(MM2_ind))){
    MM2_ind_error = MM2_ind_error+1
  }
  
  ##Huber##
  result_1= majorization_tuned(P_mat,Q_mat,1)
  rot_now = result_1$rot
  Huber_save[i,1] = rot_check(rot_now)
  if(result_1$converge==0){
    if(Huber_save[i,1]==1){
      rot_now = reflection_factor%*%rot_now
      Huber_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
      Huber_save[i,3] = result_1$scal
      Huber_save[i,4] = result_1$trans_P[1,1]
      Huber_save[i,5] = result_1$trans_P[2,1]
    } 
  }
  Res_Huber = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_Huber_dis = sqrt(rowSums(Res_Huber*Res_Huber))
  Huber_ind = tail(order(Res_Huber_dis), n = n_outliers)
  if(any(sort(ind)!=sort(Huber_ind))){
    Huber_ind_error = Huber_ind_error+1
  }
  
  
  
  ##biweight##
  result_1= majorization_tuned(P_mat,Q_mat,2)
  rot_now = result_1$rot
  biweight_save[i,1] = rot_check(rot_now)
  if(result_1$converge==0){
    if(biweight_save[i,1]==1){
      rot_now = reflection_factor%*%rot_now
      biweight_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
      biweight_save[i,3] = result_1$scal
      biweight_save[i,4] = result_1$trans_P[1,1]
      biweight_save[i,5] = result_1$trans_P[2,1]
    } 
  }
  Res_Biweight = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_Biweight_dis = sqrt(rowSums(Res_Biweight*Res_Biweight))
  Biweight_ind = tail(order(Res_Biweight_dis), n = n_outliers)
  if(any(sort(ind)!=sort(Biweight_ind))){
    Biweight_ind_error = Biweight_ind_error+1
  }
  
  
  ##procrustes##
  result_1= Procrustes_fun(P_mat,Q_mat)
  rot_now = result_1$rot
  procrustes_save[i,1] = rot_check(rot_now)
  if(procrustes_save[i,1]==1){
    rot_now = reflection_factor%*%rot_now
    procrustes_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
    procrustes_save[i,3] = result_1$scal
    procrustes_save[i,4] = result_1$trans_P[1,1]
    procrustes_save[i,5] = result_1$trans_P[2,1]
  } 
  Res_procrustes = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- 
    k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_procrustes_dis = sqrt(rowSums(Res_procrustes*Res_procrustes))
  Procrustes_ind = tail(order(Res_procrustes_dis), n = n_outliers)
  if(any(sort(ind)!=sort(Procrustes_ind))){
    Procrustes_ind_error = Procrustes_ind_error+1
  }
  
  
  
  ##S-estimator##
  result_1= S_estimator(P_mat,Q_mat,c_value,A_value,0,s_l,s_u)
  rot_now = result_1$rot
  S_est_save[i,1] = rot_check(rot_now)
  if(S_est_save[i,1]==1){
    rot_now = reflection_factor%*%rot_now
    S_est_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
    S_est_save[i,3] = result_1$scal
    S_est_save[i,4] = result_1$trans_P[1,1]
    S_est_save[i,5] = result_1$trans_P[2,1]
  }  
  Res_S = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_S_dis = sqrt(rowSums(Res_S*Res_S))
  S_ind = tail(order(Res_S_dis), n = n_outliers)
  if(any(sort(ind)!=sort(S_ind))){
    S_ind_error = S_ind_error+1
  }
  
  
  ##LMS-estimator##
  result_1= LMS_estimator(P_mat,Q_mat,0,s_l_LMS,s_u_LMS)
  rot_now = result_1$rot
  LMS_est_save[i,1] = rot_check(rot_now)
  if(LMS_est_save[i,1]==1){
    rot_now = reflection_factor%*%rot_now
    LMS_est_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
    LMS_est_save[i,3] = result_1$scal
    LMS_est_save[i,4] = result_1$trans_P[1,1]
    LMS_est_save[i,5] = result_1$trans_P[2,1]
  }  
  Res_LMS = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_LMS_dis = sqrt(rowSums(Res_LMS*Res_LMS))
  LMS_ind = tail(order(Res_LMS_dis), n = n_outliers)
  if(any(sort(ind)!=sort(LMS_ind))){
    LMS_ind_error = LMS_ind_error+1
  }
  
  
  ##LTS-estimator##
  result_1= LTS_estimator(P_mat,Q_mat,0)
  rot_now = result_1$rot
  Prof_est_save[i,1] = rot_check(rot_now)
  if(Prof_est_save[i,1]==1){
    rot_now = reflection_factor%*%rot_now
    Prof_est_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
    Prof_est_save[i,3] = result_1$scal
    Prof_est_save[i,4] = result_1$trans_P[1,1]
    Prof_est_save[i,5] = result_1$trans_P[2,1]
  }  
  Res_Prof = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_Prof_dis = sqrt(rowSums(Res_Prof*Res_Prof))
  Prof_ind = tail(order(Res_Prof_dis), n = n_outliers)
  if(any(sort(ind)!=sort(Prof_ind))){
    Prof_ind_error = Prof_ind_error+1
  }
  
  
  ##LMS_LTS-estimator##
  result_1= LMS_LTS_estimator(P_mat,Q_mat,0,s_l_LMS,s_u_LMS)
  rot_now = result_1$rot
  LMS_LTS_est_save[i,1] = rot_check(rot_now)
  if(LMS_LTS_est_save[i,1]==1){
    rot_now = reflection_factor%*%rot_now
    LMS_LTS_est_save[i,2] = atan2(rot_now[1,2],rot_now[1,1])*180/pi
    LMS_LTS_est_save[i,3] = result_1$scal
    LMS_LTS_est_save[i,4] = result_1$trans_P[1,1]
    LMS_LTS_est_save[i,5] = result_1$trans_P[2,1]
  }  
  Res_LMS_LTS = Q_mat - result_1$scal*P_mat%*%(result_1$rot)- k_4%*%t(c(result_1$trans_P[1,1],result_1$trans_P[2,1]))
  Res_LMS_LTS_dis = sqrt(rowSums(Res_LMS_LTS*Res_LMS_LTS))
  LMS_LTS_ind = tail(order(Res_LMS_LTS_dis), n = n_outliers)
  if(any(sort(ind)!=sort(LMS_LTS_ind))){
    LMS_LTS_ind_error = LMS_LTS_ind_error+1
  }
  
  
  ##################################################################
  ####save residual
  # Res_MM2_save[,c(2*i-1,2*i)] = Res_MM2
  # Res_Huber_save[,c(2*i-1,2*i)] = Res_Huber
  # Res_Biweight_save[,c(2*i-1,2*i)] = Res_Biweight
  # Res_procrustes_save[,c(2*i-1,2*i)] = Res_procrustes
  # Res_S_save[,c(2*i-1,2*i)] = Res_S 
  # Res_LMS_save[,c(2*i-1,2*i)] = Res_LMS
  # Res_Prof_save[,c(2*i-1,2*i)] = Res_Prof
  #######################################################################
  #####save residual distance
  Res_MM2_dis_save[,i] = Res_MM2_dis
  Res_Huber_dis_save[,i] = Res_Huber_dis
  Res_Biweight_dis_save[,i] = Res_Biweight_dis
  Res_procrustes_dis_save [,i] = Res_procrustes_dis
  Res_S_dis_save[,i] = Res_S_dis
  Res_LMS_dis_save[,i] = Res_LMS_dis
  Res_Prof_dis_save[,i] = Res_Prof_dis
  Res_LMS_LTS_dis_save[,i] = Res_LMS_LTS_dis
  #######################################################################
  ##save ind
  ind_save[i,] = ind
  
  #######################################################################
  # data_res_1 = data.frame(Res_MM2[,1],Res_Huber[,1],Res_Biweight[,1],Res_procrustes[,1],Res_S[,1],Res_LMS[,1],Res_Prof[,1])
  # cor_p1 = round(cor(data_res_1[,1:7]),2)
  # Pearson_matrix_1[i,1:6] = cor_p1[1,2:7] 
  # Pearson_matrix_1[i,7:11] = cor_p1[2,3:7] 
  # Pearson_matrix_1[i,12:15] = cor_p1[3,4:7] 
  # Pearson_matrix_1[i,16:18] = cor_p1[4,5:7] 
  # Pearson_matrix_1[i,19:20] = cor_p1[5,6:7] 
  # Pearson_matrix_1[i,21] = cor_p1[6,7] 
  # data_res_2 = data.frame(Res_MM2[,2],Res_Huber[,2],Res_Biweight[,2],Res_procrustes[,2],Res_S[,2],Res_LMS[,2],Res_Prof[,2])
  # cor_p2 = round(cor(data_res_2[,1:7]),2)
  # Pearson_matrix_2[i,1:6] = cor_p2[1,2:7] 
  # Pearson_matrix_2[i,7:11] = cor_p2[2,3:7] 
  # Pearson_matrix_2[i,12:15] = cor_p2[3,4:7] 
  # Pearson_matrix_2[i,16:18] = cor_p2[4,5:7] 
  # Pearson_matrix_2[i,19:20] = cor_p2[5,6:7] 
  # Pearson_matrix_2[i,21] = cor_p2[6,7]
  # 
  # 
  # cor_s1 = round(cor(data_res_1[,1:7], method = "spearman"),2)
  # Spearman_matrix_1[i,1:6] = cor_s1[1,2:7] 
  # Spearman_matrix_1[i,7:11] = cor_s1[2,3:7] 
  # Spearman_matrix_1[i,12:15] = cor_s1[3,4:7] 
  # Spearman_matrix_1[i,16:18] = cor_s1[4,5:7] 
  # Spearman_matrix_1[i,19:20] = cor_s1[5,6:7] 
  # Spearman_matrix_1[i,21] = cor_s1[6,7]
  # 
  # cor_s2 = round(cor(data_res_2[,1:7], method = "spearman"),2)
  # Spearman_matrix_2[i,1:6] = cor_s2[1,2:7] 
  # Spearman_matrix_2[i,7:11] = cor_s2[2,3:7] 
  # Spearman_matrix_2[i,12:15] = cor_s2[3,4:7] 
  # Spearman_matrix_2[i,16:18] = cor_s2[4,5:7] 
  # Spearman_matrix_2[i,19:20] = cor_s2[5,6:7] 
  # Spearman_matrix_2[i,21] = cor_s2[6,7]
  # 
  
  
}

write.csv(ind_save,file = "simulated ind.csv")

# write.csv(Res_MM2_save,file = "MM2 residual.csv")
# write.csv(Res_Huber_save,file = "Huber residual.csv")
# write.csv(Res_Biweight_save,file = "Biweight residual.csv")
# write.csv(Res_procrustes_save,file = "Procrustes residual.csv")
# write.csv(Res_S_save,file = "S residual.csv")
# write.csv(Res_LMS_save,file = "LMS residual.csv")
# write.csv(Res_Prof_save,file = "LTS residual.csv")

write.csv(Res_MM2_dis_save,file = "MM2 residual dis.csv")
write.csv(Res_Huber_dis_save,file = "Huber residual dis.csv")
write.csv(Res_Biweight_dis_save,file = "Biweight residual dis.csv")
write.csv(Res_procrustes_dis_save,file = "Procrustes residual dis.csv")
write.csv(Res_S_dis_save,file = "S residual dis.csv")
write.csv(Res_LMS_dis_save,file = "LMS residual dis.csv")
write.csv(Res_Prof_dis_save,file = "LTS residual dis.csv")
write.csv(Res_LMS_LTS_dis_save,file = "LMS LTS residual dis.csv")


save(Res_MM2_dis_save,file="Res_MM2_dis_save.RData")
save(Res_Huber_dis_save,file="Res_Huber_dis_save.RData")
save(Res_Biweight_dis_save,file="Res_Biweight_dis_save.RData")
save(Res_procrustes_dis_save,file="Res_procrustes_dis_save.RData")
save(Res_S_dis_save,file="Res_S_dis_save.RData")
save(Res_LMS_dis_save,file="Res_LMS_dis_save.RData")
save(Res_Prof_dis_save,file="Res_Prof_dis_save.RData")
save(Res_LMS_LTS_dis_save,file="Res_LMS_LTS_dis_save.RData")



# write.csv(Pearson_matrix_1,file = "Pearson_V1_with_outliers.csv")
# write.csv(Pearson_matrix_2,file = "Pearson_V2_with_outliers.csv")
# write.csv(Spearman_matrix_1,file = "Spearman_V1_with_outliers.csv")
# write.csv(Spearman_matrix_2,file = "Spearman_V2_with_outliers.csv")

#identify reflection
MM2_error= length(which(MM2_save[,1]==0))
Huber_error= length(which(Huber_save[,1]==0))
biweight_error= length(which(biweight_save[,1]==0))
procrustes_error= length(which(procrustes_save[,1]==0))
S_est_error = length(which(S_est_save[,1]==0))
LMS_est_error = length(which(LMS_est_save[,1]==0))
Prof_est_error = length(which(Prof_est_save[,1]==0))
LMS_LTS_est_error = length(which(LMS_LTS_est_save[,1]==0))
reflection_identify = matrix(NA,nrow=1,ncol=8)
colnames(reflection_identify) = c("MM2","Huber","Biweight","Procrustes","S","LMS","LTS","ILMS")
reflection_identify[,1:8] = c(100-MM2_error,100-Huber_error,100-biweight_error,100-procrustes_error,
                              100-S_est_error,100-LMS_est_error,100-Prof_est_error,100-LMS_LTS_est_error)
write.csv(reflection_identify,file = "reflection identification accuracy.csv")


#identify outliers
outlier_identify = matrix(NA,nrow=1,ncol=8)
colnames(outlier_identify) = c("MM2","Huber","Biweight","Procrustes","S","LMS","LTS","ILMS")
outlier_identify[,1:8] = c(100-MM2_ind_error,100-Huber_ind_error,100-Biweight_ind_error,100-Procrustes_ind_error,
                           100-S_ind_error,100-LMS_ind_error,100-Prof_ind_error,100-LMS_LTS_ind_error)
write.csv(outlier_identify,file = "outlier identification accuracy.csv")






