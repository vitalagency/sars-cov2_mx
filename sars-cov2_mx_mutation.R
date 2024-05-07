
# Mexico Secuencias 2020 - 2024

## ------------------------------------------------------------------------------------------------

## Librerias

library(seqinr)
library(dplyr)
library(ggplot2)

## ------------------------------------------------------------------------------------------------

## Traductor

trad =    c(UUU="F", UUC="F", UUA="L", UUG="L",
            UCU="S", UCC="S", UCA="S", UCG="S",
            UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
            UGU="C", UGC="C", UGA="STOP", UGG="W",
            CUU="L", CUC="L", CUA="L", CUG="L",
            CCU="P", CCC="P", CCA="P", CCG="P",
            CAU="H", CAC="H", CAA="Q", CAG="Q",
            CGU="R", CGC="R", CGA="R", CGG="R",
            AUU="I", AUC="I", AUA="I", AUG="M",
            ACU="T", ACC="T", ACA="T", ACG="T",
            AAU="N", AAC="N", AAA="K", AAG="K",
            AGU="S", AGC="S", AGA="R", AGG="R",
            GUU="V", GUC="V", GUA="V", GUG="V",
            GCU="A", GCC="A", GCA="A", GCG="A",
            GAU="D", GAC="D", GAA="E", GAG="E",
            GGU="G", GGC="G", GGA="G", GGG="G")

## ------------------------------------------------------------------------------------------------

## To ARN Función

ToARN = function(input) {
  newARN = as.vector(input)
  newARN[which(newARN=="t")] = "u"
  newARN = toupper(newARN)
  return (newARN)
}

## ------------------------------------------------------------------------------------------------

## Data Frames

df2020_2021 = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

df2021_2022 = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

df2022_2023 = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

df2023_2024 = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------------------------------------------

## Lectura de Secuencias

ref_2020 = read.fasta("ref_2020.txt")
length(ref_2020)

ref_2021 = read.fasta("ref_2021.txt")
length(ref_2021)

ref_2022 = read.fasta("ref_2022.txt")
length(ref_2022)

ref_2023 = read.fasta("ref_2023.txt")
length(ref_2023)

ref_2024 = read.fasta("ref_2024.txt")
length(ref_2024)


## ------------------------------------------------------------------------------------------------
##
## 2020 --> 2021
##

cat("Procesando ", as.integer(length(ref_2020)/12), " genomas de Primera Secuencia 2020 y Primera Secuencia 2021 \n")

nObs2020_2021 = 1
for (i in seq(1,length(ref_2020),1)){
  if (i==2) next
  anotaciones = attr(ref_2020[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2020[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2021), 12)){
    
    genref_2021 = ToARN( ref_2021[[k]] )  
    cat(i, k, length(genRef), length(genref_2021), "\n")
    
    if (length(genRef) == length(genref_2021)){
      dif = which(genRef != genref_2021)
      cat("length",length(dif))
      
      if (length(dif) > 0){
        for (x in dif){
          muta = paste(genRef[x],"to",genref_2021[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2021 = paste(genref_2021[inicioCodon], genref_2021[inicioCodon+1], genref_2021[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonref_2021, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonref_2021], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2021]) && trad[codonOri]!=trad[codonref_2021]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df2020_2021[nObs2020_2021,] = obs 
            nObs2020_2021 = nObs2020_2021+1
          }
        }
      }
    }else{
    }
  }
}

##
## 2020 --> 2021 END
##
## ------------------------------------------------------------------------------------------------
##
##
## 2021 --> 2022 START
##
cat("Procesando ", as.integer(length(ref_2021)/12), " genomas de Primera Secuencia 2021 y Primera Secuencia 2022 \n")

nObs2021_2022 = 1
for (i in seq(1,length(ref_2021),1)){
  if (i==2) next
  anotaciones = attr(ref_2021[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2021[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2022), 12)){
    
    genref_2022 = ToARN( ref_2022[[k]] )  
    cat(i, k, length(genRef), length(genref_2022), "\n")
    
    if (length(genRef) == length(genref_2022)){
      dif = which(genRef != genref_2022)
      cat("length",length(dif))
      
      if (length(dif) > 0){
        for (x in dif){
          muta = paste(genRef[x],"to",genref_2022[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2022 = paste(genref_2022[inicioCodon], genref_2022[inicioCodon+1], genref_2022[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonref_2022, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonref_2022], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2022]) && trad[codonOri]!=trad[codonref_2022]){
            obs2 = list(muta,codonChange,aminoChange,geneName)
            df2021_2022[nObs2021_2022,] = obs2
            nObs2021_2022 = nObs2021_2022+1
          }
        }
      }
    }else{
    }
  }
}
##
## 2021 - 2022 END
##
## ------------------------------------------------------------------------------------------------
##
## 2022 - 2023 Start
##
nObs2022_2023 = 1
for (i in seq(1,length(ref_2022),1)){
  if (i==2) next
  anotaciones = attr(ref_2022[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2022[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2023), 12)){
    
    genref_2023 = ToARN( ref_2023[[k]] )  
    cat(i, k, length(genRef), length(genref_2023), "\n")
    
    if (length(genRef) == length(genref_2023)){
      dif = which(genRef != genref_2023)
      cat("length",length(dif))
      
      if (length(dif) > 0){
        for (x in dif){
          muta = paste(genRef[x],"to",genref_2023[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2023 = paste(genref_2023[inicioCodon], genref_2023[inicioCodon+1], genref_2023[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonref_2023, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonref_2023], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2023]) && trad[codonOri]!=trad[codonref_2023]){
            obs3 = list(muta,codonChange,aminoChange,geneName)
            df2022_2023[nObs2022_2023,] = obs3
            nObs2022_2023 = nObs2022_2023+1
          }
        }
      }
    }else{
    }
  }
}
##
## 2022 - 2023 END
##
## ------------------------------------------------------------------------------------------------
##
## 2023 - 2024 Start
##
nObs2023_2024 = 1
for (i in seq(1,length(ref_2023),1)){
  if (i==2) next
  anotaciones = attr(ref_2023[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2023[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2024), 12)){
    
    genref_2024 = ToARN( ref_2024[[k]] )  
    cat(i, k, length(genRef), length(genref_2024), "\n")
    
    if (length(genRef) == length(genref_2024)){
      dif = which(genRef != genref_2024)
      cat("length",length(dif))
      
      if (length(dif) > 0){
        for (x in dif){
          muta = paste(genRef[x],"to",genref_2024[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2024 = paste(genref_2024[inicioCodon], genref_2024[inicioCodon+1], genref_2024[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonref_2024, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonref_2023], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2024]) && !is.na(trad[codonOri]) && trad[codonOri]!=trad[codonref_2024]){
            obs4 = list(muta,codonChange,aminoChange,geneName)
            df2023_2024[nObs2023_2024,] = obs4
            nObs2023_2024 = nObs2023_2024+1
          }
        }
      }
    }else{
    }
  }
}
##
## 2023 - 2024 END
##
## ------------------------------------------------------------------------------------------------
##
## 2020 --> 2021 Plot
##
p = ggplot(df2020_2021)
p = p + aes(x=Mutation, fill=Mutation, label=after_stat(count))
p = p + ggtitle("Mutaciones de sustitución México 2020 - 2021")
p = p + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
p = p + geom_bar(stat = "count")
p = p + geom_text(stat = "count", vjust=1.5)
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

df2020_2021gr = filter(
  summarise(
    select(
      group_by(df2020_2021, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

df2020_2021gr = df2020_2021 %>%
  filter(!is.na(Amino)) %>%
  group_by(Amino) %>%
  summarise(
    mutation2 = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ) %>%
  filter(Cuenta > 1) %>%
  arrange(desc(Cuenta)) %>%
  slice(1:20)

head(df2020_2021gr)
nrow(df2020_2021gr)

p2 = ggplot(df2020_2021gr)
p2 = p2 + aes(x = Amino, y = Cuenta, fill = Amino, label = Cuenta)
p2 = p2 + ggtitle("Cambio de Aminoácidos México 2020 - 2021")
p2 = p2 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p2 = p2 + geom_bar(stat = "identity")
p2 = p2 + geom_text(stat = "identity", vjust = 1.5)
p2 = p2 + facet_wrap(~ Gene, scales = "free", space = "free_x")
p2 = p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

##
## 2020 --> 2021 Plot END
##
## ------------------------------------------------------------------------------------------------
##
## 2021 --> 2022 Plot Start
##
## 
p2021_2022 = ggplot(df2021_2022)
pp2021_2022 = pp2021_2022 + aes(x=Mutation, fill=Mutation, label=after_stat(count))
pp2021_2022 = pp2021_2022 + ggtitle("Mutaciones de sustitución México 2021 - 2022")
pp2021_2022 = pp2021_2022 + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
pp2021_2022 = pp2021_2022 + geom_bar(stat = "count")
p2021_2022 = p2021_2022 + geom_text(stat = "count", vjust=1.5)
p2021_2022 = p2021_2022 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2021_2022

df2021_2022gr = filter(
  summarise(
    select(
      group_by(df2021_2022, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

df2021_2022gr = df2021_2022 %>%
  filter(!is.na(Amino)) %>%
  group_by(Amino) %>%
  summarise(
    mutation2 = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ) %>%
  filter(Cuenta > 1) %>%
  arrange(desc(Cuenta)) %>%
  slice(1:20)

head(df2021_2022gr)
nrow(df2021_2022gr)

p3 = ggplot(df2021_2022gr)
p3 = p3 + aes(x = Amino, y = Cuenta, fill = Amino, label = Cuenta)
p3 = p3 + ggtitle("Cambio de Aminoácidos México 2021 - 2022")
p3 = p3 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p3 = p3 + geom_bar(stat = "identity")
p3 = p3 + geom_text(stat = "identity", vjust = 1.5)
p3 = p3 + facet_wrap(~ Gene, scales = "free", space = "free_x")
p3 = p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3
##
## 2021 --> 2022 Plot END
##
## ------------------------------------------------------------------------------------------------
##
## 2022 --> 2023 Plot Start
##
## 
p2022_2023 = ggplot(df2022_2023)
pp2022_2023 = pp2021_2022 + aes(x=Mutation, fill=Mutation, label=after_stat(count))
pp2022_2023 = pp2021_2022 + ggtitle("Mutaciones de sustitución México 2022 - 2023")
pp2022_2023 = pp2021_2022 + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
pp2022_2023 = pp2021_2022 + geom_bar(stat = "count")
p2022_2023 = p2021_2022 + geom_text(stat = "count", vjust=1.5)
p2022_2023 = p2021_2022 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2022_2023

df2022_2023gr = filter(
  summarise(
    select(
      group_by(df2022_2023, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

df2022_2023gr = df2022_2023 %>%
  filter(!is.na(Amino)) %>%
  group_by(Amino) %>%
  summarise(
    mutation2 = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ) %>%
  filter(Cuenta > 1) %>%
  arrange(desc(Cuenta)) %>%
  slice(1:20)

head(df2022_2023gr)
nrow(df2022_2023gr)

p4 = ggplot(df2022_2023gr)
p4 = p4 + aes(x = Amino, y = Cuenta, fill = Amino, label = Cuenta)
p4 = p4 + ggtitle("Cambio de Aminoácidos México 2022 - 2023")
p4 = p4 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p4 = p4 + geom_bar(stat = "identity")
p4 = p4 + geom_text(stat = "identity", vjust = 1.5)
p4 = p4 + facet_wrap(~ Gene, scales = "free", space = "free_x")
p4 = p4 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4
##
## 2022 --> 2023 Plot END
##
## ------------------------------------------------------------------------------------------------
##
## 2023 --> 2024 Plot Start
##
## 
p2023_2024 = ggplot(df2023_2024)
pp2023_2024 = pp2023_2024 + aes(x=Mutation, fill=Mutation, label=after_stat(count))
pp2023_2024 = pp2023_2024 + ggtitle("Mutaciones de sustitución México 2023 - 2024")
pp2023_2024 = pp2023_2024 + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
pp2023_2024 = pp2023_2024 + geom_bar(stat = "count")
p2023_2024 = p2023_2024 + geom_text(stat = "count", vjust=1.5)
p2023_2024 = p2023_2024 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2023_2024

df2023_2024gr = filter(
  summarise(
    select(
      group_by(df2023_2024, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

df2023_2024gr = df2023_2024 %>%
  filter(!is.na(Amino)) %>%
  group_by(Amino) %>%
  summarise(
    mutation2 = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ) %>%
  filter(Cuenta > 1) %>%
  arrange(desc(Cuenta)) %>%
  slice(1:20)

head(df2023_2024gr)
nrow(df2023_2024gr)

p5 = ggplot(df2023_2024gr)
p5 = p5 + aes(x = Amino, y = Cuenta, fill = Amino, label = Cuenta)
p5 = p5 + ggtitle("Cambio de Aminoácidos México 2023 - 2024")
p5 = p5 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p5 = p5 + geom_bar(stat = "identity")
p5 = p5 + geom_text(stat = "identity", vjust = 1.5)
p5 = p5 + facet_wrap(~ Gene, scales = "free", space = "free_x")
p5 = p5 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p5
## ------------------------------------------------------------------------------------------------
##
## 2023 --> 2024 Plot END
##
## 