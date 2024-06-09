####################################################
####################################################
#######                                      #######
####### SIMULAZIONE LOBLOLLY PINE POPULATION #######
#######                                      #######
####################################################
####################################################

################################
# GENERAZIONE DI DISTRIBUZIONI #
################################

# Imposto i parametri della distribuzione che andro' a generare successivamente, in modo da averne una simile a quella dell'area basale dei pini Loblolly, analizzata nell'articolo.
N <- 14443 
mu <- 0.8 
sigma <- 0.72
# Media e Deviazione standard li ho scelti a priori per poter generare una distribuzione abbastanza simile a quella del dataset.

# Calcolo di k (forma) e theta (scala) basato sui parametri forniti
k <- (mu^2) / (sigma^2)
theta <- (sigma^2) / mu

## DISTRIBUZIONE GAMMA
set.seed(524)
gamma_distr <- rgamma(N, shape = k, scale = theta)
summary(gamma_distr)

## DISTRIBUZIONE NORMALE
# Utilizzo gli stessi parametri (mu e sigma)
set.seed(524)
norm_distr <- rnorm(N, mean = mu, sd = sigma)
summary(norm_distr)

##########################
# ISTOGRAMMI A CONFRONTO #
##########################

par(mfrow=c(1,2))
## ISTOGRAMMA DELLA DISTRIBUZIONE GAMMA
hist(gamma_distr,breaks = seq(min(gamma_distr),max(gamma_distr), length.out=50), main = "Distribuzione asimmetrica\nsimulata con una\nDistribuzione Gamma", xlab = "valori",col="#d692df",border = "#c071d0")
abline(v=mu, lwd=2)
text(1.1,1000, cex=1.5, labels=expression(paste(mu)))

## ISTOGRAMMA DELLA DISTRIBUZIONE NORMALE
hist(norm_distr,breaks = seq(min(norm_distr),max(norm_distr), length.out=50), main = "Approssimazione Normale", xlab = "valori",col="#d692df",border = "#c071d0")
abline(v=mu, lwd=2)
text(1.1,500, cex=1.5, labels=expression(paste(mu)))

######################################
# FORMULA USUALE PER IL CALCOLO DI n #
######################################

## CALCOLO LA NUMEROSITA' CAMPIONARIA CON LA FORMULA 2

# Definisco i parametri mancanti
delta <- 0.3 # La precisione desiderata (margine di errore)
alpha <- 0.1 # Il livello di significatività per un livello di confidenza del 90%

# Calcolo il valore critico t 
t_value <- qt(1 - alpha/2, df = N-1)

# Calcolo la numerosità campionaria
n_pilota <- round((t_value * sigma / delta)^2)
n_pilota

# Itero una seconda volta per portare i gradi di libertà ad uno in meno
t_value <- qt(1 - alpha/2, df = n_pilota-1)
n_pilota <- round((t_value * sigma / delta)^2)
n_pilota

# Ora estrarrò dalle popolazioni Gamma e Gaussiana rispettivamente due campioni con Campionamento Casuale Semplice con un pacchetto specifico
# install.packages("sampling")
# library(sampling)
# indices <- srswor(n_pilota, N) # 'indices' contiene gli indici degli elementi selezionati dalla popolazione
# gamma_camp <- gamma_distr[indices]
# norm_camp <- norm_distr[indices]


# Poiche' e' impossibile nella realta' che la deviazione standard sia data a priori, la stimiamo attraverso il calcolo della deviazione standard dal campione pilota.

# Estraggo i campioni con Campionamento Casuale Semplice fissando il seme
set.seed(524)
gamma_camp_pilota <- sample(gamma_distr, n_pilota, replace = FALSE)
norm_camp_pilota <- sample(norm_distr, n_pilota,replace=FALSE)
sigma_hat <-  sd(gamma_camp_pilota)
sigma_hat_norm <- sd(norm_camp_pilota)

# Calcolo nuovamente la Formula (2) iterando due volte
t_value <- qt(1 - alpha/2, df = n_pilota-1)
n <- ceiling((t_value * sigma_hat / delta)^2)
n
t_value <- qt(1 - alpha/2, df = n-1)
n <- ceiling((t_value * sigma_hat / delta)^2)
n # Dimensione del campione desiderata

# Estraggo i campioni con la nuova numerosità campionaria
gamma_camp <- sample(gamma_distr,n,replace = FALSE)
norm_camp <- sample(norm_distr,n,replace = FALSE)

cat("La numerosita' campionaria calcolata con la formula usuale e' di", n, "data dalla seconda iterazione, per portare i gradi di libertà ad uno in meno ripetto alla dimensione del campione stimata. I due campioni estratti con Campionamento Casuale Semplice Senza Ripetizione sono: \nGAMMA\n", gamma_camp, "\nNORMALE\n", norm_camp)

###############
# SIMULAZIONI #
###############

# Ora che e' chiaro il procedimento eseguo 50.000 simulazioni per ottenere una rappresentazione realistica delle distribuzioni di campionamento
M <- 50000

## ITERAZIONE PER LA DISTRIBUZIONE GAMMA

# Inizializzo un vettore per l'allocazione dei risultati
estimated_n <- numeric(M)
set.seed(524)
for (i in 1:M) {
  # Simulo il prelievo di un campione pilota
  campione_ssr <- sample(gamma_distr, size = n_pilota, replace = FALSE)
  
  # Calcolo la deviazione standard del campione pilota
  deviazione_standard_campione <- sd(campione_ssr)
  
  # Calcolo la dimensione del campione desiderata n tramite la formula usuale
  n <- ceiling((qt(1 - alpha / 2, length(campione_ssr) - 1) * deviazione_standard_campione / delta)^2)
  
  # Memorizzo il risultato nel vettore preallocato
  estimated_n[i] <- max(ceiling(n),2)
}
# Visualizzo i primi 10 risultati
head(estimated_n)
cat("Mean estimate n from (2): ", ceiling(mean(estimated_n)), "\nMaximum estimate n: ",max(estimated_n))

# ITERAZIONE PER LA DISTIBUZIONE NORMALE
est_n_norm <- numeric(M)
set.seed(524)
for (i in 1:M) {
  campione_ssr<-sample(norm_distr,size = n_pilota,replace = FALSE)
  deviazione_standard_campione <-sd(campione_ssr)
  n <- ceiling((qt(1-alpha/2,length(campione_ssr)-
             1)*deviazione_standard_campione/delta)^2)
  est_n_norm[i]<-max(ceiling(n),2)
}
head(est_n_norm)
cat("Mean estimate n from (2): ", ceiling(mean(est_n_norm)), "\nMaximum estimate n: ",max(est_n_norm))

#######################################
# MEAN ESTIMATED MARGIN ERROR USING n #
#######################################

# RICAVO DELTA DALLA DISTRIBUZIONE GAMMA
delta_gamma <- round(qnorm(1 - alpha/2) * sigma_hat / sqrt(estimated_n),2) # Utilizzo sigma_hat
summary(delta_gamma)
cat("Mean estimated margin error using n from Gamma ditribution: ",round(mean(delta_gamma),2))

# RICAVO DELTA DALLA DISTRIBUZIONE NORMALE
delta_norm <- round(qnorm(1 - alpha/2) * sigma_hat_norm / sqrt(est_n_norm),2)
cat("Mean estimate margin error using n from gaussian ditribution: ",round(mean(delta_norm),2))

par(mfrow= c(2,2))
# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n PER LA GAMMA
hist(estimated_n, breaks=seq(min(estimated_n),max(estimated_n),length.out=70), main = "Istogramma della distribuzione\ndegli n con la (2) per la distribuzione Gamma", ylim = c(0,0.07),col = "#7dbedf",border = "#7aa9ce", xlab = "n", freq = F)
abline(v=ceiling(mean(estimated_n)), lwd=3)
text(20, 0.05,labels=expression(paste(bar(n))), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta PER LA GAMMA
hist(delta_gamma, breaks = seq(min(delta_gamma),max(delta_gamma),length.out=55), main = "Istogramma della distribuzione\ndei delta corrispondente", col = "#f8c583",border = "#fa8282", xlab = expression(paste(delta)), freq = F)
abline(v=mean(delta_gamma),lwd=3)
text(0.24,6,labels=expression(paste(bar(delta))), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n PER LA NORMALE
hist(est_n_norm, breaks=seq(min(est_n_norm),max(est_n_norm),length.out=50), main = "Istogramma della distribuzione\ndegli n con la (2) per la Normale",  col = "#7dbedf",border = "#7aa9ce", xlab = "n",freq = F)
abline(v=ceiling(mean(est_n_norm)), lwd=3)
text(19, 0.04,labels=expression(paste(bar(n))), cex=1.5)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta PER LA NORMALE
hist(delta_norm, breaks = seq(min(delta_norm),max(delta_norm),length.out=60), col = "#f8c583",border = "#fa8282", main = "Istogramma della distribuzione\ndei delta corrispondente", xlab = expression(paste(delta)),freq = F)
abline(v=mean(delta_norm), lwd=3)
text(0.30,10,labels=expression(paste(bar(delta))), cex=1.5, lwd=2)

############################
# INTERVALLI DI CONFIDENZA #
############################

## PER CAMPIONE SINGOLO

# LIVELLO DI CONFIDENZA AL 90% PER LA GAMMA

# Calcolo del margine di errore
margine_errore <- qnorm(1-alpha/2) * (sigma_hat / sqrt(ceiling(mean(estimated_n)))) 
mu_camp_G <- mean(gamma_camp)
# Calcolo l'intervallo di confidenza
inf <- mu_camp_G - margine_errore
sup <- mu_camp_G + margine_errore
cat("Intervallo di confidenza al 90% per la media: [", inf, ",", sup, "]\n")

# LIVELLO DI CONFIDENZA AL 90% PER LA NORMALE

# Eseguo lo stesso procedimento utilizzato sopra, cambiando la deviazione standard e la stima di n
margin_error <- qnorm(1-alpha/2) * (sigma_hat_norm / sqrt(ceiling(mean(est_n_norm)))) 
mu_camp_N <- mean(norm_camp)
inf_norm <- mu_camp_N - margin_error
sup_norm <- mu_camp_N + margin_error
cat("Intervallo di confidenza della Normale al 90% per la media: [", inf_norm, ",", sup_norm, "]\n")

####################################################
# INTERVALLI DI  CONFIDENZA PER 50.000 SIMULAZIONI #
####################################################

## LIVELLI DI CONFIDENZA AL 90% PER LA GAMMA

# Anche qui, ora che il procedimento è chiaro eseguo M simulazioni
# Inizializzo un vettore per memorizzare i risultati degli intervalli di confidenza
IC <- matrix(NA, nrow = M, ncol = 2)

# Imposto il livello di confidenza al 90%
livello_confidenza <- 1 - alpha

# Itero il processo M volte
for (i in 1:M) {

  # Seleziono un campione SSR di dimensione n (dato dalla media delle stime di n)
  campione_ssr <- sample(gamma_distr, size = ceiling(mean(estimated_n)), replace = FALSE)
  
  # Calcolo l'intervallo di confidenza per la media
  intervallo_confidenza <- t.test(campione_ssr, conf.level = livello_confidenza)$conf.int
  
  # Memorizzo l'intervallo di confidenza
  IC[i, ] <- intervallo_confidenza
}

#Visualizzo i primi 10 intervalli di confidenza
head(IC)

## LIVELLI DI CONFIDENZA AL 90% PER LA NORMALE

# Utilizzo lo stesso procedimento della Gamma
I_C <- matrix(NA, nrow = M, ncol = 2)

livello_confidenza <- 1 - alpha

for (i in 1:M) {
  
  campione_ssr <- sample(norm_distr, size = ceiling(mean(est_n_norm)), replace = FALSE)
  
  
  intervallo_confidenza <- t.test(campione_ssr, conf.level = livello_confidenza)$conf.int
  
  I_C[i, ] <- intervallo_confidenza
}

head(I_C)

#######################################
# PROBABILITA' DI TOLLERANZA EMPIRICA #
#######################################

## CALCOLO LA PROBABILITA' DI TOLLERANZA EMPIRICA PER LA DISTRIBUZIONE GAMMA

# Inizializzo il contatore per intervalli di confidenza che coprono la media della popolazione
copertura <- 0

# Controllo ogni intervallo di confidenza per vedere se copre la media della popolazione
for (i in 1:M) {
  if (IC[i, 1] <= mu && IC[i, 2] >= mu) {
    copertura <- copertura + 1
  }
}
# Calcola la probabilità di tolleranza empirica
ETP <- round(copertura / M, 2)

# ESEGUO LO STESSO PROCEDIMENTO PER LA NORMALE
copertura_norm <- 0

# Controlla ogni intervallo di confidenza per vedere se copre la media vera
for (i in 1:M) {
  if (I_C[i, 1] <= mu && I_C[i, 2] >= mu) {
    copertura_norm <- copertura_norm + 1
  }
}
# Calcolo la probabilità di tolleranza empirica
norm_ETP <- round(copertura_norm / M, 2)

cat("La probabilià di tolleranza empirica per la distribuzione Gamma (1-gamma) e' di:",ETP,"\nLa probabilià di tolleranza empirica per la distribuzione Normale (1-gamma) e' di:",norm_ETP)
# Al contrario di quanto descritto nell'articolo il tasso di copertura è molto alto per la stima della dimensione campionaria con la formula usuale, ma soprattutto è più alto per la normale.

#################################
# MARGINE DI ERRORE TOLLERABILE #
#################################

## GAMMA

# Inizializzo un vettore per memorizzare le medie campionarie per la Gamma
mu_camp_G <- numeric(M)

# Estraggo i campioni e calcolo le medie
for (i in 1:M) {

    campione_ssr <- sample(gamma_distr, ceiling(mean(estimated_n)), replace = TRUE)
  
    mu_camp_G[i] <- mean(campione_ssr)
}

# Calcolo quante medie campionarie sono entro delta dalla media della popolazione
medie_entro_delta_G <- sum(abs(mu_camp_G - mu) <= delta)

percentuale_entro_delta_G <- ceiling((medie_entro_delta_G / M) * 100)

cat("Percentage of sample means within delta of the Gamma population mean:", percentuale_entro_delta_G, "%\n")

# UTILIZZO LO STESSO PROCEDIMENTO PER LA NORMALE
# Inizializzo un vettore per memorizzare le medie campionarie per la Gamma
mu_camp_N <- numeric(M)

# Estraggo i campioni e calcolo le medie
for (i in 1:M) {
  n <- ceiling((qt(1 - alpha / 2, length(campione_ssr) - 1) * deviazione_standard_campione / delta)^2)
  
  campione_ssr <- sample(norm_distr, ceiling(mean(est_n_norm)), replace = TRUE)
  
  mu_camp_N[i] <- mean(campione_ssr)
}

# Calcolo quante medie campionarie sono entro delta dalla media della popolazione
medie_entro_delta_N <- sum(abs(mu_camp_N - mu) <= delta)

percentuale_entro_delta_N <- ceiling((medie_entro_delta_N / M) * 100)

cat("Percentage of sample means within delta of the Normal population mean:", percentuale_entro_delta_N, "%\n")

##############################################################################

#####################################
# FORMULA PER n* DI KUPPER E HAFNER #
#####################################

# Ricordo i parametri che mi servono
alpha <- 0.1 # Livello di confidenza al 90%
gamma <- 0.1 # Livello di tolleranza al 95%

# Funzione per calcolare n*
calculate_n_star <- function(n, alpha, gamma) {
  n_star <- n # Inizializzo n_star con n calcolato con la formula usuale
  while(TRUE) {
    chi_sq_n_star_minus_1 <- qchisq(1 - gamma, df = n_star - 1)
    f_1_n_star_minus_1 <- qf(1 - alpha, df1 = 1, df2 = n_star - 1)
    chi_sq_1 <- qchisq(1 - alpha, df = 1)
    right_side <- n * chi_sq_n_star_minus_1 * f_1_n_star_minus_1 / chi_sq_1
    left_side <- n_star * (n_star - 1)
    
    if(left_side >= right_side) {
      break
    } else {
      n_star <- n_star + 1
    }
  }
  return(n_star)
}

## FORMULA (5) PER LA DISTRIBUZIONE GAMMA

# Inizializzo un vettore per memorizzare i risultati di n*
n_star_risult <- numeric(length(estimated_n))

# Calcolo n* per ogni simulazione
for (i in 1:length(estimated_n)) {
  n_star_risult[i] <- calculate_n_star(estimated_n[i], alpha, gamma)
}

# A questo punto, n_star_risult contiene i valori di n* calcolati per ogni n in estimated_n
cat("Mean estimate n* from (5): ", ceiling(mean(n_star_risult)), "\nMaximum estimate n*: ", max(n_star_risult))

## FORMULA (5) PER LA DISTRIBUZIONE NORMALE

# Utilizzando la funzione  gia' definita, ri-eseguo il ciclo for per est_n_norm
# Inizializzo un vettore per memorizzare i risultati di n*
n_star_risult_norm <- numeric(length(est_n_norm))

for (i in 1:length(est_n_norm)) {
  n_star_risult_norm[i] <- calculate_n_star(est_n_norm[i], alpha, gamma)
}

cat("Mean estimate n* from (5): ", ceiling(mean(n_star_risult_norm)), "\nMaximum estimate n*: ", max(n_star_risult_norm))


########################################
# MEAN ESTIMATED MARGIN ERROR USING n* #
########################################

## RICAVO DELTA STAR DALLA DISTRIBUZIONE GAMMA
delta_star_gamma <- round(qnorm(1 - alpha/2) * sigma_hat / sqrt(n_star_risult),2) # Utilizzo sigma o sigma_hat?????????
summary(delta_star_gamma)
cat("Mean estimated margin error using n* from Gamma ditribution: ",round(mean(delta_star_gamma),2))

## RICAVO DELTA STAR DALLA DISTRIBUZIONE NORMALE
delta_star_norm <- round(qnorm(1 - alpha/2) * sigma_hat_norm / sqrt(n_star_risult_norm),2)
summary(delta_star_norm)
cat("Mean estimate margin error using n* from gaussian ditribution: ",round(mean(delta_star_norm),2))

par(mfrow=c(2,2))
# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n* PER LA GAMMA
hist(n_star_risult, breaks=seq(min(n_star_risult),max(n_star_risult),length.out=40), main = "Istogramma della distribuzione\ndegli n* con la (5) per la Gamma", col = "#d0d38f",border = "#a6a872", xlab = "n", freq = F)
abline(v=ceiling(mean(n_star_risult)), lwd=3)
text(30, 0.04,labels=expression(paste(bar(n)),"   *"), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta star PER LA GAMMA
hist(delta_star_gamma, breaks = seq(min(delta_star_gamma),max(delta_star_gamma),length.out=40), main = "Istogramma della distribuzione\ndei delta star corrispondente", col = "#fe91be",border = "#f567a4", xlab = expression(paste(delta)), freq = F)
abline(v=mean(delta_star_gamma), lwd=3)
text(0.18,15,labels=expression(paste(bar(delta)),"   *"), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n* PER LA NORMALE
hist(n_star_risult_norm, breaks=seq(min(n_star_risult_norm),max(n_star_risult_norm),length.out=40), main = "Istogramma della distribuzione\ndegli n* con la (5) per la Normale",  col = "#d0d38f",border = "#a6a872", xlab = "n*",freq = F)
abline(v=ceiling(mean(n_star_risult_norm)), lwd=3)
text(29, 0.06,labels=expression(paste(bar(n)),"   *"), cex=1.5)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta star PER LA NORMALE
hist(delta_star_norm, breaks = seq(min(delta_star_norm),max(delta_star_norm),length.out=35), col = "#fe91be",border = "#f567a4", main = "Istogramma della distribuzione\ndei delta star corrispondente", xlab = expression(paste(delta)),freq = F)
abline(v=mean(delta_star_norm), lwd=3)
text(0.26,14,labels=expression(paste(bar(delta)),"   *"), cex=1.5, lwd=2)

####################################################
# INTERVALLI DI  CONFIDENZA PER 50.000 SIMULAZIONI #
####################################################

## LIVELLI DI CONFIDENZA AL 90% PER LA GAMMA

# Secondo vettore inizializzato per allocare gli intervalli di confidenza per le stime di n*
II_IC <- matrix(NA, nrow = M, ncol = 2)

livello_confidenza <- 1 - alpha

for (i in 1:M) {
  
  # Seleziono un campione SSR di dimensione n (dato dalla media delle stime di n*)
  campione_ssr <- sample(gamma_distr, size = ceiling(mean(n_star_risult)), replace = FALSE)
  
  intervallo_confidenza <- t.test(campione_ssr, conf.level = livello_confidenza)$conf.int
  
  II_IC[i, ] <- intervallo_confidenza
}

#Visualizzo i primi 10 intervalli di confidenza
head(II_IC)

## LIVELLI DI CONFIDENZA AL 90% PER LA NORMALE

# Utilizzo lo stesso procedimento sopra
II_I_C <- matrix(NA, nrow = M, ncol = 2)

livello_confidenza <- 1 - alpha

for (i in 1:M) {
  
  campione_ssr <- sample(norm_distr, size = ceiling(mean(n_star_risult_norm)), replace = FALSE)
  
  
  intervallo_confidenza <- t.test(campione_ssr, conf.level = livello_confidenza)$conf.int
  
  II_I_C[i, ] <- intervallo_confidenza
}

head(II_I_C)

#######################################
# PROBABILITA' DI TOLLERANZA EMPIRICA #
#######################################

##  LA DISTRIBUZIONE GAMMA

# Inizializzo il contatore per intervalli di confidenza che coprono la media della popolazione
copert_nst <- 0

# Controllo ogni intervallo di confidenza per vedere se copre la media della popolazione
for (i in 1:M) {
  if (II_IC[i, 1] <= mu && II_IC[i, 2] >= mu) {
    copert_nst <- copert_nst + 1
  }
}
# Calcola la probabilità di tolleranza empirica
ETP_nst <- round(copert_nst / M, 2)

# ESEGUO LO STESSO PROCEDIMENTO PER LA NORMALE
copertura_norm_nst <- 0

# Controlla ogni intervallo di confidenza per vedere se copre la media vera
for (i in 1:M) {
  if (II_I_C[i, 1] <= mu && II_I_C[i, 2] >= mu) {
    copertura_norm_nst <- copertura_norm_nst + 1
  }
}
# Calcolo la probabilità di tolleranza empirica
norm_ETP_nst <- round(copertura_norm_nst / M, 2)

cat("In riferimento al secondo calcolo della numerosità campionaria con la formula di Kupper e Hafner:\n","La probabilià di tolleranza empirica per la distribuzione Gamma (1-gamma) e' di:",ETP_nst,"\nLa probabilià di tolleranza empirica per la distribuzione Normale (1-gamma) e' di:",norm_ETP_nst)
# Al contrario di quanto descritto nell'articolo il tasso di copertura è molto alto anche per la stima della dimensione campionaria con la formula di Kupper e Hafner, ma soprattutto è più alto per la normale.

#################################
# MARGINE DI ERRORE TOLLERABILE #
#################################

## GAMMA

# Inizializzo un vettore per memorizzare le medie campionarie per la Gamma
mu_camp_G <- numeric(M)

# Estraggo i campioni e calcolo le medie
for (i in 1:M) {
  
  campione_ssr <- sample(gamma_distr, ceiling(mean(n_star_risult)), replace = TRUE)
  
  mu_camp_G[i] <- mean(campione_ssr)
}

# Calcolo quante medie campionarie sono entro delta dalla media della popolazione
mu_w_delta_G <- sum(abs(mu_camp_G - mu) <= delta)

perc_entro_delta_G <- ceiling((mu_w_delta_G / M) * 100)

cat("Percentage of sample means within delta of the Gamma population mean:", perc_entro_delta_G, "%\n")

## UTILIZZO LO STESSO PROCEDIMENTO PER LA NORMALE
# Inizializzo un vettore per memorizzare le medie campionarie per la Gamma
mu_camp_N <- numeric(M)

# Estraggo i campioni e calcolo le medie
for (i in 1:M) {

  campione_ssr <- sample(norm_distr, ceiling(mean(n_star_risult_norm)), replace = TRUE)
  
  mu_camp_N[i] <- mean(campione_ssr)
}

# Calcolo quante medie campionarie sono entro delta dalla media della popolazione
mu_w_delta_N <- sum(abs(mu_camp_N - mu) <= delta)

perc_w_delta_N <- ceiling((mu_w_delta_N / M) * 100)

cat("Percentage of sample means within delta of the Normal population mean:", perc_w_delta_N, "%\n")

########################################
# COEFFICIENTE DI ASIMMETRIA DI FISHER #
########################################

# Installo il pacchetto e1071 e lo carico
install.packages("e1071")
library(e1071)


# La funzione skewness calcola il coefficiente di asimmetria di Fisher (G1), ovvero il il momento terzo della media del campione, diviso la deviazione standard al cubo
G1_gamma <- skewness(gamma_distr,type = 3)
G1_norm <- skewness(norm_distr, type = 3)

# Elevo al quadrato (per considerare solamente la magnitudine dell'asimmetria, ignorando la direzione) e moltiplico per 25 (per adeguare l'effetto dell'asimmetria sulla dimensione del campione necessaria per raggiungere un determinato livello di precisione)
G1_Gamma_ajusted <- round(25 * G1_gamma^2)
G1_Gamma_ajusted
G1_Norm_adjusted <- round(25*G1_norm^2)
G1_Norm_adjusted

tabella1 <- data.frame(
  Population=c("Mean estimate n from (2)", 
               "Maximum estimate n",
               "Mean estimate margin error using n",
               "Empirical Tollerance Probability",
               "Percentage of sample means within delta of the population mean:",
               "Mean estimate n* from (2)", 
               "Maximum estimate n*",
               "Mean estimate margin error using n*",
               "Empirical Tollerance Probability",
               "Percentage of sample means within delta of the population mean:",
               "25G1^2"), # da calcolare con n*
  Normale= round(c(ceiling(mean(est_n_norm)),
             max(est_n_norm),
             mean(delta_norm),
             norm_ETP,
             percentuale_entro_delta_N,
             ceiling(mean(n_star_risult_norm)), 
             max(n_star_risult_norm),
             mean(delta_star_norm),
             norm_ETP_nst,
             perc_w_delta_N,
             G1_Norm_adjusted),2),
  Gamma= round(c(ceiling(mean(estimated_n)), 
           max(estimated_n),
           mean(delta_gamma),
           ETP,
           percentuale_entro_delta_G,
           ceiling(mean(n_star_risult)), 
           max(n_star_risult),
           mean(delta_star_gamma),
           ETP_nst,
           perc_entro_delta_G,
           G1_Gamma_ajusted),2))


View(tabella1)


