library(tidyr)
library(BBmisc)
library(arm)
#wczytanie zbioru do ramki data
data <- read.csv("pulsar_stars.csv",header = TRUE, sep=",")
attach(data)
#wyœwietlenie pierwszych 6 wierszy danych
head(data)
#wyœwietlenie zmiennych i typu danych
str(data)
#podstawowa statystyka dla zbioru data
summary(data)
#sprawdzenie czy s¹ brakuj¹ce dane
sapply(data, function(x) sum(is.na(x)))


#Analiza czêstotliwoœci wyników
hist(Mean.of.the.integrated.profile, ylim = c(0,4000))
hist(Standard.deviation.of.the.integrated.profile, ylim = c(0,6000))
hist(Excess.kurtosis.of.the.integrated.profile)
hist(Skewness.of.the.integrated.profile)
hist(Mean.of.the.DM.SNR.curve)
hist(Standard.deviation.of.the.DM.SNR.curve)
hist(Excess.kurtosis.of.the.DM.SNR.curve, ylim = c(0,5000))
hist(Skewness.of.the.DM.SNR.curve, ylim = c(0,10000))
hist(target_class,breaks = 4, ylim = c(0,20000))

detach(data)

data_norm <- normalize(data, method = "range", range = c(0, 1), margin = 2)


#podzia³ znormalizowanego zbioru na testowy i treningowy
set.seed(666)
indeks <- floor(0.7 * nrow(data_norm))
sample <- sample(seq_len(nrow(data_norm)), size = indeks)
train <- data_norm[sample, ]
test <- data_norm[-sample, ]


#utworzenie modelu na podstawie zbioru treningowego
model.1 <- glm(target_class~., data=train, family="binomial")
summary(model.1)

#predykcja przy wykorzystaniu model.1 na podstawie zbioru testowego
przewidywanie1 <- predict(model.1,test,interval="confidence",type = "response", level=0.95)
przewidywane.zmienne1 <- ifelse(przewidywanie1 > 0.5, "1", "0")
head(przewidywane.zmienne1)

#obliczenie procentu dobrze przewidzianych wyników
CM1<-table(przewidywane.zmienne1, test$target_class)
accuracy1<-(CM1[1,1]+CM1[2,2])/nrow(test)
accuracy1


#spierwiastkowanie zmiennej Excess.kurtosis.of.the.DM.SNR.curve i sprawdzenie czy jej znaczenie sie nie zwiekszy
model.1sqrt <- glm(target_class~.-Excess.kurtosis.of.the.DM.SNR.curve+sqrt(Excess.kurtosis.of.the.DM.SNR.curve),data=train,family="binomial")
summary(model.1sqrt)
przewidywanie1sq <- predict(model.1sqrt,test,interval="confidence",type = "response", level=0.95)
przewidywane.zmienne1sq <- ifelse(przewidywanie1sq > 0.5, "1", "0")
head(przewidywane.zmienne1sq)

CM1sq<-table(przewidywane.zmienne1sq, test$target_class)
accuracy1sq<-(CM1sq[1,1]+CM1sq[2,2])/nrow(test)
accuracy1sq


#model pomniejszony o zmienn¹ o najmniejszym znaczeniu Excess.kurtosis.of.the.DM.SNR.curve
model.2 <- glm(target_class~.-Excess.kurtosis.of.the.DM.SNR.curve,data=train,family="binomial")
summary(model.2)
przewidywanie2 <- predict(model.2,test,interval="confidence",type = "response", level=0.95)
przewidywane.zmienne2 <- ifelse(przewidywanie2 > 0.5, "1", "0")
head(przewidywane.zmienne2)


#Diagnostyka modelu
CM2<-table(przewidywane.zmienne2, test$target_class)
accuracy2<-(CM2[1,1]+CM2[2,2])/nrow(test)
accuracy2
plot(residuals.glm(model.2,type="response"), xlab = "Observation index", ylab = "Response residuals", main="model.2")
plot(residuals.glm(model.2, type="deviance"), xlab = "Observation index", ylab = "Deviance residuals", main="model.2")
plot(residuals.glm(model.2, type="pearson"), xlab = "Observation index", ylab = "Pearson residuals", main="model.2")
binnedplot(fitted(model.2), 
           residuals(model.2, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot model 2", 
           cex.pts = 0.5, 
           col.pts = 1, 
           col.int = "gray")
plot(przewidywanie2, xlab = "Observation index")
plot(model.2, which = 4, id.n = 3)


cooksd<-cooks.distance(model.2)
influential <- as.numeric(names(cooksd)[(cooksd > (4/nrow(train)))])
top_x_outlier <- 3
influential <- as.numeric(names(sort(cooksd, decreasing = TRUE)[1:top_x_outlier]))
#Usuniecie trzech obserwacji wp³ywowych ze zbioru train
trainv2 <- train[!row.names(train) %in% influential,]


#Utworzenie model.3 wykorzystujac zbior testv2  usunieta zmienna
model.3 <- glm(target_class~.-Excess.kurtosis.of.the.DM.SNR.curve, data=trainv2, family="binomial")
summary(model.3)
przewidywanie3 <- predict(model.3,test,interval="confidence",type = "response", level=0.95)
przewidywane.zmienne3 <- ifelse(przewidywanie3 > 0.5, "1", "0")
head(przewidywane.zmienne3)


#Diagnostyka modelu
CM3<-table(przewidywane.zmienne3, test$target_class)
accuracy3<-(CM3[1,1]+CM3[2,2])/nrow(test)
accuracy3
plot(residuals.glm(model.3,type="response"), xlab = "Observation index", ylab = "Response residuals", main="model.3")
plot(residuals.glm(model.3, type="deviance"), xlab = "Observation index", ylab = "Deviance residuals", main="model.3")
plot(residuals.glm(model.3, type="pearson"), xlab = "Observation index", ylab = "Pearson residuals", main="model.3")
binnedplot(fitted(model.3), 
           residuals(model.3, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot model 3", 
           cex.pts = 0.5, 
           col.pts = 1, 
           col.int = "gray")
plot(przewidywanie3, xlab = "Observation index")

#McFadden R^2 aby oszacowaæ dopasowanie modelu funkcja dla modeli logistycznych
library(pscl)
pR2(model.1)
pR2(model.1sqrt)
pR2(model.2)
pR2(model.3)






