---
title: "Simulazione con Gamma del dataset 'Alabama Loblolly Pine'"
author: "Federica Prosperini"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
N <- 14443 
mu <- 0.8 
sigma <- 0.72
k <- (mu^2) / (sigma^2)
theta <- (sigma^2) / mu
#DISTRIBUZIONE GAMMA
set.seed(524)
gamma_distr <- rgamma(N, shape = k, scale = theta)
summary(gamma_distr)
# DISTRIBUZIONE NORMALE
set.seed(524)
norm_distr <- rnorm(N, mean = mu, sd = sigma)
summary(norm_distr)
delta <- 0.3 # La precisione desiderata (margine di errore)
alpha <- 0.1 # Il livello di significatività per un livello di confidenza del 90%
gamma <- 0.05 # Livello di tolleranza al 95%

```


## Plotto Le due distribuzioni

Ho impostato i parametri  (media e deviazione standard) che andro' a generare con una Distribuzione Gamma, in modo da averne una simile a quella dell'area basale dei pini Loblolly, analizzata nell'articolo.
Successivamente, allo stesso modo, andrò a generare una Distribuzione Gaussiana per avere un'approssimazione normale.

```{r pressure, echo=FALSE}
par(mfrow=c(1,2))
## ISTOGRAMMA DELLA DISTRIBUZIONE GAMMA
hist(gamma_distr,breaks = seq(min(gamma_distr),max(gamma_distr), length.out=50), main = "Distribuzione asimmetrica\nsimulata con una\nDistribuzione Gamma", xlab = "valori",col="#ffc7f4",border = "#f683c8")
abline(v=mu, col="#b71336", lwd=3)
text(2,1000, col="black", "Media")

## ISTOGRAMMA DELLA DISTRIBUZIONE NORMALE
hist(norm_distr,breaks = seq(min(norm_distr),max(norm_distr), length.out=50), main = "Approssimazione Normale", xlab = "valori",col="#ffc7f4",border = "#f683c8")
abline(v=mu, col="#b71336", lwd=3)
text(2,600,  "Media")
```

## CALCOLO LA NUMEROSITA' CAMPIONARIA CON LA FORMULA 2
Stimo la dimensione campionaria con la formula usuale:$$
n = \left( \frac{t_{1-\frac{\alpha}{2}, n-1} \hat{\sigma}}{\delta} \right)^2
$$in modo da trovare il campione pilota di dimensione n.

```{r formula usuale, echo=TRUE}
# Calcolo il valore critico t 
t_value <- qt(1 - alpha/2, df = N-1)

# Calcolo la numerosità campionaria
n_pilota <- round((t_value * sigma / delta)^2)
n_pilota

# Itero una seconda volta per portare i gradi di libertà ad uno in meno
t_value <- qt(1 - alpha/2, df = n_pilota-1)
n_pilota <- round((t_value * sigma / delta)^2)
n_pilota

```
Poiche' e' impossibile nella realta' che la deviazione standard sia data a priori, la stimiamo attraverso il calcolo della deviazione standard dal campione pilota, estratto con Campionamento Casuale Semplice.
```{r campionamento}
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
```
```{r risultati, echo=FALSE}
cat("La numerosita' campionaria calcolata con la formula usuale e' di", n, "data dalla seconda iterazione, per portare i gradi di libertà ad uno in meno ripetto alla dimensione del campione stimata. I due campioni estratti con Campionamento Casuale Semplice Senza Ripetizione sono: \nGAMMA\n", gamma_camp, "\nNORMALE\n", norm_camp)
```
Ora che e' chiaro il procedimento eseguo 50.000 simulazioni per ottenere una rappresentazione realistica delle distribuzioni di campionamento
```{r simulazioni per la Gamma}
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
```

```{r gamma risultati, echo=FALSE}
cat("Mean estimate n from (2): ", ceiling(mean(estimated_n)), "\nMaximum estimate n: ",max(estimated_n))
```

Allo stesso modo svolgo le simulazioni per la distribuzione normale ed ottengo:
```{r simulazione per la normale, include=FALSE}
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
```

```{r risultati simulazione normale, echo=FALSE}
cat("Mean estimate n from (2): ", ceiling(mean(est_n_norm)), "\nMaximum estimate n: ",max(est_n_norm))
```

## MEAN ESTIMATED MARGIN ERROR USING n
$$
\delta = \frac{Z_{1-\frac{\alpha}{2}} \sigma}{\sqrt{n}}
$$
Indica la metà del margine di un intervallo di confidenza basato sul disegno al 100 × (1 - α)%.
```{r delta}
delta_gamma <- round(qnorm(1 - alpha/2) * sigma / sqrt(estimated_n),2)
delta_norm <- round(qnorm(1 - alpha/2) * sigma / sqrt(est_n_norm),2)
```
```{r delta risultati, echo=FALSE}
cat("Mean estimated margin error using n from Gamma ditribution: ",round(mean(delta_gamma),2))

cat("Mean estimate margin error using n from gaussian ditribution: ",round(mean(delta_norm),2))
```

## CONFRONTI GRAFICI
Con questi istogrammi voglio rappresentare le distribuzioni empiriche delle 50000 simulazioni svolte, sia per il calcolo della dimensione campionaria, sia per il calcolo della metà del margine di errore. Sia per la distribuzione Gamma, sia per la Normale.
```{r confronto plot, echo=FALSE}
par(mfrow=c(2,2))
# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n PER LA GAMMA
hist(estimated_n, breaks=seq(min(estimated_n),max(estimated_n),length.out=600), main = "Istogramma della distribuzione\ndegli n con la (2) per la\ndistribuzione GAMMA", col = "#FF90BC",border = "#EC4C72", xlab = "n", freq = F)
abline(v=ceiling(mean(estimated_n)), col="#9f33b1", lwd=3)
text(20, 0.04,labels=expression(paste(bar(n))), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta PER LA GAMMA
hist(delta_gamma, breaks = seq(min(delta_gamma),max(delta_gamma),length.out=100), main = "Istogramma della distribuzione\ndei delta corrispondente", col = "#f8c583",border = "#fa8282", xlab = expression(paste(delta)), freq = F)
abline(v=mean(delta_gamma), col="darkorange", lwd=3)
text(0.36,4,labels=expression(paste(bar(delta))), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n PER LA NORMALE
hist(est_n_norm, breaks=seq(min(est_n_norm),max(est_n_norm),length.out=100), main = "Istogramma della distribuzione\ndegli n con la (2) per la normale",  col = "#FF90BC",border = "#EC4C72", xlab = "n",freq = F)
abline(v=ceiling(mean(est_n_norm)), col="#9f33b1", lwd=3)
text(20, 0.04,labels=expression(paste(bar(n))), cex=1.5)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta PER LA NORMALE
hist(delta_norm, breaks = seq(min(delta_norm),max(delta_norm),length.out=100), col = "#f8c583",border = "#fa8282", main = "Istogramma della distribuzione\ndei delta corrispondente", xlab = expression(paste(delta)),freq = F)
abline(v=mean(delta_norm), col="darkorange", lwd=3)
text(0.32,8,labels=expression(paste(bar(delta))), cex=1.5, lwd=2)
```

E' evidante che le distribuzioni di campionamento sono fortemente asimmetriche, soprattutto quella della dimensione campionaria calcolata sulla gamma, in quanto all'interno della formula hanno una Chi quadro.
Ma notiamo che anche le distribuzioni dei margini di errore sono leggermente asimmetriche.

## INTERVALLI DI  CONFIDENZA PER 50.000 SIMULAZIONI

Nel seguito procederò a calcolare gli intervalli di confidenza per le due distribuzioni in modo tale da andarne a valutare l'affidabilità attraverso la probabilità di tolleranza empirica

```{r IC}
## LIVELLI DI CONFIDENZA AL 90% PER LA GAMMA

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
```

## PROBABILITA' DI TOLLERANZA
La probabilità di tolleranza empirica è la capacità di coprire una certa proporzione della popolazione. 

```{r prob toll emp}
# PROCEDIMENTO PER LA GAMMA

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

for (i in 1:M) {
  if (I_C[i, 1] <= mu && I_C[i, 2] >= mu) {
    copertura_norm <- copertura_norm + 1
  }
}
norm_ETP <- round(copertura_norm / M, 2)
```

```{r risultati etp, echo=FALSE}
cat("La probabilià di tolleranza empirica per la distribuzione Gamma (1-gamma) e' di:",ETP,"\nLa probabilià di tolleranza empirica per la distribuzione Normale (1-gamma) e' di:",norm_ETP)
```
Al contrario di quanto descritto nell'articolo il tasso di copertura è molto alto per la stima della dimensione campionaria con la formula usuale, ma soprattutto è più alto per la normale.

## MARGINE DI ERRORE TOLLERABILE 

```{r margin error}
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



# UTILIZZO LO STESSO PROCEDIMENTO PER LA NORMALE
mu_camp_N <- numeric(M)

for (i in 1:M) {
  n <- ceiling((qt(1 - alpha / 2, length(campione_ssr) - 1) * deviazione_standard_campione / delta)^2)
  
  campione_ssr <- sample(norm_distr, ceiling(mean(est_n_norm)), replace = TRUE)
  
  mu_camp_N[i] <- mean(campione_ssr)
}

medie_entro_delta_N <- sum(abs(mu_camp_N - mu) <= delta)

percentuale_entro_delta_N <- ceiling((medie_entro_delta_N / M) * 100)


```

```{r % mean within delta, echo=FALSE}
cat("Percentage of sample means within delta of the Gamma population mean:", percentuale_entro_delta_G, "%\n")
cat("Percentage of sample means within delta of the Normal population mean:", percentuale_entro_delta_N, "%\n")

```
Le medie stimate basate su campioni di dimensione n riescono ad essere entro δ dalla media della popolazione leggermente più del (1 - α) × 100% = 90% del tempo.


## FORMULA PER n* DI KUPPER E HAFNER 
*"Sulla base del loro studio di simulazione,Kupper and Hafner (1989, p. 102) hanno affermato che la "formula della dimensione
del campione" usuale "porta sempre a una sottostima grave della dimensione del campione richiesta", anzi hanno consigliato
che il suo utilizzo dovrebbe essere evitato"*

$$
n^*(n^* - 1) \geq n\chi^2_{n^*-1,1-\gamma} F_{1,n^*-1,1-\alpha}/\chi^2_{1,1-\alpha}
$$
La formula sopra viene utilizzata per un calcolo meno approssimativo della precedente. Il valore che ci serve n* verrà ricavato risolvendo l'equazione quadratica implicita.

Nel seguito saranno resi noti solamenti i risultati più interessanti in quanto i procedimenti saranno gli stessi, con valori diversi.
```{r n*}
# Definisco la funzione per calcolare n*
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

## FORMULA (5) PER LA DISTRIBUZIONE NORMALE

# Utilizzando la funzione  gia' definita, ri-eseguo il ciclo for per est_n_norm
n_star_risult_norm <- numeric(length(est_n_norm))

for (i in 1:length(est_n_norm)) {
  n_star_risult_norm[i] <- calculate_n_star(est_n_norm[i], alpha, gamma)
}

```

```{r risultati n*, echo=FALSE}
cat("Mean estimate n* from (5): ", ceiling(mean(n_star_risult)), "\nMaximum estimate n*: ", max(n_star_risult))
cat("\nMean estimate n* from (5): ", ceiling(mean(n_star_risult_norm)), "\nMaximum estimate n*: ", max(n_star_risult_norm))
```


```{r iterazione per n*, include=FALSE}
## RICAVO DELTA STAR DALLA DISTRIBUZIONE GAMMA
delta_star_gamma <- round(qnorm(1 - alpha/2) * sigma / sqrt(n_star_risult),2) # Utilizzo sigma o sigma_hat?????????

## RICAVO DELTA STAR DALLA DISTRIBUZIONE NORMALE
delta_star_norm <- round(qnorm(1 - alpha/2) * sigma / sqrt(n_star_risult_norm),2)
```

```{r principali risultati, echo=FALSE}
cat("Mean estimated margin error using n* from Gamma ditribution: ",round(mean(delta_star_gamma),2))
cat("Mean estimate margin error using n* from gaussian ditribution: ",round(mean(delta_star_norm),2))

```

## RAPPRESENTAZIONE GRAFICA
# Del campionamento con la nuova formula e del calcolo di delta con i nuovi valori di n*

```{r plot, echo=FALSE}
par(mfrow=c(2,2))
# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n* PER LA GAMMA
hist(n_star_risult, breaks=seq(min(n_star_risult),max(n_star_risult),length.out=50), main = "Istogramma della distribuzione\ndegli n* con la (5) per la GAMMA", col = "#c0ffee",border = "#76ffe2", xlab = "n", freq = F)
abline(v=ceiling(mean(n_star_risult)), col="#209fba", lwd=3)
text(33, 0.04,labels=expression(paste(bar(n)),"   *"), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta star PER LA GAMMA
hist(delta_star_gamma, breaks = seq(min(delta_star_gamma),max(delta_star_gamma),length.out=50), main = "Istogramma della distribuzione\ndei delta star corrispondente", col = "#996699",border = "#a64d8c", xlab = expression(paste(delta)), freq = F)
abline(v=mean(delta_star_gamma), col="#bf1a73", lwd=3)
text(0.25,10,labels=expression(paste(bar(delta)),"   *"), cex=1.5, lwd=2)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEGLI n* PER LA NORMALE
hist(n_star_risult_norm, breaks=seq(min(n_star_risult_norm),max(n_star_risult_norm),length.out=50), main = "Istogramma della distribuzione\ndegli n* con la (5) per la NORMALE",  col = "#c0ffee",border = "#76ffe2", xlab = "n*",freq = F)
abline(v=ceiling(mean(n_star_risult_norm)), col="#209fba", lwd=3)
text(31, 0.06,labels=expression(paste(bar(n)),"   *"), cex=1.5)

# ISTOGRAMMA DELLA DISTRIBUZIONE DEI delta star PER LA NORMALE
hist(delta_star_norm, breaks = seq(min(delta_star_norm),max(delta_star_norm),length.out=50), col = "#996699",border = "#a64d8c", main = "Istogramma della distribuzione\ndei delta star corrispondente", xlab = expression(paste(delta)),freq = F)
abline(v=mean(delta_star_norm), col="#bf1a73", lwd=3)
text(0.25,20,labels=expression(paste(bar(delta)),"   *"), cex=1.5, lwd=2)

```

```{r risult n*, include=FALSE}
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


# UTILIZZO LO STESSO PROCEDIMENTO PER LA NORMALE
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

```


```{r princip risult, echo=FALSE}

cat("In riferimento al secondo calcolo della numerosità campionaria con la formula di Kupper e Hafner:\n","La probabilià di tolleranza empirica per la distribuzione Gamma (1-gamma) e' di:",ETP_nst,"\nLa probabilià di tolleranza empirica per la distribuzione Normale (1-gamma) e' di:",norm_ETP_nst)


cat("Percentage of sample means within delta of the Gamma population mean:", perc_entro_delta_G, "%\n")

cat("Percentage of sample means within delta of the Normal population mean:", perc_w_delta_N, "%\n")
```


## TABELLA
# Visualizzo i principali risultati in una tabella che distingue il campionamento usuale da quello di Kupper e Hafner e distingue anche i risultati per la distribuzione Gamma e per la Gaussiana
```{r tab, echo=FALSE, cache=FALSE}
library(knitr)

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
             0),2),
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
           85),2))

kable(tabella1, caption = "Tabella riassuntiva")


```
