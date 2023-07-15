#Download our dataset 'logD_final.xlsx'
#Set the location where the database was downloaded as working directory 
#if your terminal does not have the required packages, remove the '#' on the next lines:
#install.packages('rJava')
#install.packages('rcdk')
#install.packages('readxl')
#install.packages('caret')
#install.packages('dplyr')
#install.packages('ggplot2')
#install.packages('Metrics')
#install.packages('bestglm')
#install.packages('pROC')
#install.packages('glmnet')
#install.packages('randomForest')
#install.packages('e1071')

#Run from row 20 to 304 to train the ML models


library("rJava")
library(rcdk)
library(readxl)
library(caret)
library(dplyr)
library(ggplot2)

data <- read_excel("logD_final.xlsx")
data$cond <- as.factor(ifelse(data$d3 > 0.2, 0,1))
acids <- data %>% filter(type=='Acid',is.na(type2))
bases <- data %>% filter(type=='Base',is.na(type2))
zw <- data %>% filter(type2=='Z')

#generate smiles objects
smilesA <- acids$SMILES
smilesB <- bases$SMILES
molsA <- parse.smiles(smilesA)
molsB <- parse.smiles(smilesB)

#generate descriptor matrix
descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
descsA <- eval.desc(molsA, descNames)
descsB <- eval.desc(molsB, descNames)

#agregar descriptores experimentales
descsA$pKa <- acids$pKa
descsA$delta <- acids$delta
descsA$logPN <- acids$logPN
descsA$CH_strength <- acids$CH_strength
descsA$XH_strength <- acids$XH_strength
descsA$HBA_strength <- acids$HBA_strength
descsA$Hyd_Apolar <- acids$Hyd_Apolar
descsA$Hyd_Polar <- acids$Hyd_Polar
descsA$Hyd <- acids$Hyd


descsB$pKa <- bases$pKa
descsB$delta <- bases$delta
descsB$logPN <- bases$logPN
descsB$CH_strength <- bases$CH_strength
descsB$XH_strength <- bases$XH_strength
descsB$HBA_strength <- bases$HBA_strength
descsB$Hyd_Apolar <- bases$Hyd_Apolar
descsB$Hyd_Polar <- bases$Hyd_Polar
descsB$Hyd <- bases$Hyd

#reduction step#
descsA <- descsA[, !apply(descsA, 2, function(x) any(is.na(x)) )] #remove NAs
descsA <- descsA[, !apply( descsA, 2, function(x) length(unique(x)) == 1 )] #remove constant columns

descsB <- descsB[, !apply(descsB, 2, function(x) any(is.na(x)) )] #remove NAs
descsB <- descsB[, !apply( descsB, 2, function(x) length(unique(x)) == 1 )] #remove constant columns



#remove correlated descriptors
r2 <- which(cor(descsA)^2 > .6, arr.ind=TRUE)
r2 <- r2[ r2[,1] > r2[,2] , ]
descsA <- descsA[, -unique(r2[,2])]


r2 <- which(cor(descsB)^2 > .6, arr.ind=TRUE)
r2 <- r2[ r2[,1] > r2[,2] , ]
descsB <- descsB[, -unique(r2[,2])]





#agregar $cond a descs
descsA$cond <- acids$cond
descsB$cond <- bases$cond



#Dividir los dfs con los descriptores entre los que tienen cond de 0 o 1
descsA0 <- descsA %>% filter(cond==0) %>% select(-cond)
descsA1 <- descsA %>% filter(cond==1) %>% select(-cond)

descsB0 <- descsB %>% filter(cond==0) %>% select(-cond)
descsB1 <- descsB %>% filter(cond==1) %>% select(-cond)



#Creando un DF con los promedios de cada descriptor y la diferencia relativa de los promedios
meanA <- data.frame(
  'descs' = c(colnames(descsA0))
)
for (i in 1:length(meanA[,1])){
  meanA$mean0[i] = mean(descsA0[,i])
  meanA$mean1[i] = mean(descsA1[,i])
  meanA$diff_m[i] = abs(meanA$mean0[i]-meanA$mean1[i])
  meanA$rsd0[i] = sd(descsA0[,i])/sqrt(length(descsA0[,i]))
  meanA$rsd1[i] = sd(descsA1[,i])/sqrt(length(descsA1[,i]))
  meanA$diff_rsd[i] = sqrt(meanA$rsd0[i]^2+meanA$rsd1[i]^2)
}


meanB <- data.frame(
  'descs' = c(colnames(descsB0))
)
for (i in 1:length(meanB[,1])){
  meanB$mean0[i] = mean(descsB0[,i])
  meanB$mean1[i] = mean(descsB1[,i])
  meanB$diff_m[i] = abs(meanB$mean0[i]-meanB$mean1[i])
  meanB$rsd0[i] = sd(descsB0[,i])/sqrt(length(descsB0[,i]))
  meanB$rsd1[i] = sd(descsB1[,i])/sqrt(length(descsB1[,i]))
  meanB$diff_rsd[i] = sqrt(meanB$rsd0[i]^2+meanB$rsd1[i]^2)
}





meanA <- meanA %>% filter(diff_rsd<diff_m)
meanB <- meanB %>% filter(diff_rsd<diff_m)







acids_d <- data.frame(matrix(NA,
                             nrow = length(smilesA),
                             ncol = length(descsA)))

for (i in 1:length(descsA)){
  for (n in 1:nrow(meanA)){
    if (colnames(descsA)[i] == meanA$descs[n]) {
      acids_d[,i] <- descsA[,i]
    }
  }
}
acids_d <- acids_d %>% select_if(~ !any(is.na(.)))
colnames(acids_d) <- meanA$descs
acids_d$cond <- acids$cond





bases_d <- data.frame(matrix(NA,
                             nrow = length(smilesB),
                             ncol = length(descsB)))

for (i in 1:length(descsB)){
  for (n in 1:nrow(meanB)){
    if (colnames(descsB)[i] == meanB$descs[n]) {
      bases_d[,i] <- descsB[,i]
    }
  }
}
bases_d <- bases_d %>% select_if(~ !any(is.na(.)))
colnames(bases_d) <- meanB$descs
bases_d$cond <- bases$cond


#_______________________________WELCH'S T-TEST_________________________________
#acids
pvalue_a <- c()
for (i in 1:ncol(acids_d)){
  pvalue_a <- c(pvalue_a,t.test(ifelse(acids_d[,ncol(acids_d)]==1,acids_d[,i],NA),ifelse(acids_d[,ncol(acids_d)]==0,acids_d[,i],NA))[[3]])
}
meanA$welchs_p <- pvalue_a
meanA <- meanA %>% filter(welchs_p<0.05)
meanA <- meanA[-8,]

acids_d <- data.frame(matrix(NA,
                             nrow = length(smilesA),
                             ncol = length(descsA)))

for (i in 1:length(descsA)){
  for (n in 1:nrow(meanA)){
    if (colnames(descsA)[i] == meanA$descs[n]) {
      acids_d[,i] <- descsA[,i]
    }
  }
}
acids_d <- acids_d %>% select_if(~ !any(is.na(.)))
colnames(acids_d) <- meanA$descs
acids_d$cond <- acids$cond




#bases
pvalue_b <- c()
for (i in 1:ncol(bases_d)){
  pvalue_b <- c(pvalue_b,t.test(ifelse(bases_d[,ncol(bases_d)]==1,bases_d[,i],NA),ifelse(bases_d[,ncol(bases_d)]==0,bases_d[,i],NA))[[3]])
}
meanB$welchs_p <- pvalue_b
meanB <- meanB %>% filter(welchs_p<0.05)

bases_d <- data.frame(matrix(NA,
                             nrow = length(smilesB),
                             ncol = length(descsB)))

for (i in 1:length(descsB)){
  for (n in 1:nrow(meanB)){
    if (colnames(descsB)[i] == meanB$descs[n]) {
      bases_d[,i] <- descsB[,i]
    }
  }
}
bases_d <- bases_d %>% select_if(~ !any(is.na(.)))
colnames(bases_d) <- meanB$descs
bases_d$cond <- bases$cond

#_____________END WELCHS T-TEST______________________________





#------------------------------------------------------------


##CREAR TRARINING Y TEST SETS
library(Metrics)


#Crear el training set y test set de los ACIDOS y bases
set.seed(1234)
sample.index <- sample(1:nrow(acids_d),nrow(acids_d)*0.8,replace=FALSE)
training.set_A <- acids_d[sample.index,]
test.set_A <- acids_d[-sample.index,]

sample.index <- sample(1:nrow(bases_d),nrow(bases_d)*0.8,replace=FALSE)
training.set_B <- bases_d[sample.index,]
test.set_B <- bases_d[-sample.index,]



#--------------LOGISTIC REGRESSION-----------------------------------
library(bestglm)
library(pROC)
library(glmnet)
best.fit_A <- bestglm(training.set_A, IC = "AIC", family = binomial, method = "exhaustive")
best.fit_A$BestModel$coefficients
acids_dLR <- acids_d[,c('ALogp2','MDEC.11','khs.sCH3','HybRatio','C1SP3','C2SP3','delta','cond')]
training.set_ALR <- training.set_A[,c('ALogp2','MDEC.11','khs.sCH3','HybRatio','C1SP3','C2SP3','delta','cond')]
test.set_ALR <- test.set_A[,c('ALogp2','MDEC.11','khs.sCH3','HybRatio','C1SP3','C2SP3','delta','cond')]

fit_A_LR <- train(cond~.,data=training.set_A,method='glm',family='binomial')


best.fit_B <- bestglm(training.set_B, IC = "AIC", family = binomial, method = "exhaustive")
best.fit_B$BestModel$coefficients
bases_dLR <- bases_d[,c('ALogp2','delta','cond')]
training.set_BLR <- training.set_B[,c('ALogp2','delta','cond')]
test.set_BLR <- test.set_B[,c('ALogp2','delta','cond')]

fit_B_LR <- train(cond~.,data=training.set_BLR,method='glm',family='binomial')

#--------RANDOM FOREST----------------------
library(randomForest)
fit_A_RF <- randomForest(cond~.,data=training.set_A)


mtry <- tuneRF(training.set_A[-ncol(training.set_A)],training.set_A$cond,
               stepFactor=1.2,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
fit_A_RF <- randomForest(cond~.,data=training.set_A,mtry=best.m,importance=TRUE)


fit_B_RF <- randomForest(cond~.,data=training.set_B)

mtry <- tuneRF(training.set_B[-ncol(training.set_B)],training.set_B$cond,
               stepFactor=1.2,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
fit_B_RF <- randomForest(cond~.,data=training.set_B,mtry=best.m,importance=TRUE)


#--------SVM----------
library(e1071)

#-------linear kernel----------------------
###MODELOS PARA LOS ACIDOS

fit_A_SVML <- svm(cond~.,data=training.set_A, type='C-classification',kernel="linear")
fit_A_SVMP <- svm(cond~.,data=training.set_A, type='C-classification',kernel="polynomial")

fit_B_SVML <- svm(cond~.,data=training.set_B, type='C-classification',kernel="linear")
fit_B_SVMP <- svm(cond~.,data=training.set_B, type='C-classification',kernel="polynomial")



#-------------------PROGRAM-------------------------
####INSTRUCTIONS
#Create a .xlsx file that has the following columns:
#- SMILES: Use the canonical SMILES Code for each molecule 
#- type (acid or base) 
#- pKa (it can be calculated with ChemAxon) 
#- pH: The desired pH at which you want to classificate your molecule




#Use the script 'jazzy_descriptors_v3.sh' to calculate the bond strengths and hydration energies with Jazzy. 
#Follow ths instructions from the script and include the results into in the .xlsx file.

predictor <- read_excel('test.xlsx')
predictor$delta <- ifelse(predictor$type=='acid',predictor$pH-predictor$pKa,predictor$pKa-predictor$pH)


#calcular los descriptores diferenciados entre los elegidos para ácidos y bases
smilesP <- predictor$SMILES
molsP <- parse.smiles(smilesP)
descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
descsP<- eval.desc(molsP, descNames)
descsP$delta <- predictor$delta



if(predictor$type[1]=='acid'){
  descsP$CH_strength <- predictor$CH_strength
   acids_P <- data.frame(matrix(NA,
                                nrow = length(smilesP),
                                ncol = length(descsP)))
  
  for (i in 1:ncol(descsP)){
    for (n in 1:nrow(meanA)){
      if (colnames(descsP)[i] == colnames(acids_d)[n]) {
        acids_P[,i] <- descsP[,i]
      }
    }}
  acids_P <- acids_P %>% select_if(~ !any(is.na(.)))
  colnames(acids_P) <- meanA$descs
  } else {
    bases_P <- data.frame(matrix(NA,
                                 nrow = length(smilesP),
                                 ncol = length(descsP)))
    
    for (i in 1:ncol(descsP)){
      for (n in 1:nrow(meanB)){
        if (colnames(descsP)[i] == colnames(bases_d)[n]) {
          bases_P[,i] <- descsP[,i]
        }
      }}
      bases_P <- bases_P %>% select_if(~ !any(is.na(.)))
      
      
      colnames(bases_P) <- meanB$descs
  }


#evaluar modelos de ML
if(predictor$type=='acid'){
  predictor$LR <- ifelse(predict(fit_A_LR,acids_P)==0,'Eq.2','Eq.1')
  predictor$RFC <- ifelse(predict(fit_A_RF,acids_P)==0,'Eq.2','Eq.1')
  predictor$SVML <- ifelse(predict(fit_A_SVML,acids_P)==0,'Eq.2','Eq.1')
  predictor$SVMP <- ifelse(predict(fit_A_SVMP,acids_P)==0,'Eq.2','Eq.1')
} else {
  predictor$LR <- ifelse(predict(fit_B_LR,bases_P)==0,'Eq.2','Eq.1')
  predictor$RFC <- ifelse(predict(fit_B_RF,bases_P)==0,'Eq.2','Eq.1')
  predictor$SVML <- ifelse(predict(fit_B_SVML,bases_P)==0,'Eq.2','Eq.1')
  predictor$SVMP <- ifelse(predict(fit_B_SVMP,bases_P)==0,'Eq.2','Eq.1')
}



library(writexl)
write_xlsx(predictor,path ='results_ML.xlsx')


