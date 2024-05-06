library(seqinr)

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

ToARN = function(input) {
  newARN = as.vector(input)
  newARN[which(newARN=="t")] = "u"
  newARN = toupper(newARN)
  return (newARN)
}

df = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

df_mut1 = data.frame(
  Mutation1 = character(),
  Codon1 = character(),
  Amino1 = character(),
  Gene1 = character(),
  stringsAsFactors = FALSE
)


df_mut2 = data.frame(
  Mutation2 = character(),
  Codon2 = character(),
  Amino2 = character(),
  Gene2 = character(),
  stringsAsFactors = FALSE
)

df_mut3 = data.frame(
  Mutation3 = character(),
  Codon3 = character(),
  Amino3 = character(),
  Gene3 = character(),
  stringsAsFactors = FALSE
)

df_mut4 = data.frame(
  Mutation4 = character(),
  Codon4 = character(),
  Amino4 = character(),
  Gene4 = character(),
  stringsAsFactors = FALSE
)

df_mut5 = data.frame(
  Mutation5 = character(),
  Codon5 = character(),
  Amino5 = character(),
  Gene5 = character(),
  stringsAsFactors = FALSE
)


# Primera secuencia México
ref_2020 = read.fasta("mx1_sequence")
length(ref_2020)

# Mutaciones México 2 : 5
ref_2021 = read.fasta("mx2_sequence")
length(ref_2021)

ref_2022 = read.fasta("mx3_sequence")
length(ref_2022)

ref_2023 = read.fasta("mx4_sequence")
length(ref_2024)

ref_2024 = read.fasta("mx5_sequence")
length(fMx5)

resultados = lista_ref()

for (year in c(2020, 2021, 2022, 2023, 2024)){
  ref_actual = get(paste("ref_", year, sep=""))
}

  for (compyear in c(2020, 2021, 2022, 2023, 2024)){
    if (compyear >= year) next
    
    ref_comp = get(paste("ref_", compyear, sep=""))
    
    df_mut1 = data.frame(
      Mutation1 = character(),
      Codon1 = character(),
      Amino1 = character(),
      Gene1 = character(),
      stringsAsFactors = FALSE
    )
    
    df_mut2 = data.frame(
      Mutation2 = character(),
      Codon2 = character(),
      Amino2 = character(),
      Gene2 = character(),
      stringsAsFactors = FALSE
    )
    
    df_mut3 = data.frame(
      Mutation3 = character(),
      Codon3 = character(),
      Amino3 = character(),
      Gene3 = character(),
      stringsAsFactors = FALSE
    )
    
    df_mut4 = data.frame(
      Mutation4 = character(),
      Codon4 = character(),
      Amino4 = character(),
      Gene4 = character(),
      stringsAsFactors = FALSE
    )
    
    df_mut5 = data.frame(
      Mutation5 = character(),
      Codon5 = character(),
      Amino5 = character(),
      Gene5 = character(),
      stringsAsFactors = FALSE
    )
    
    
    
  }


cat("Procesando ", as.integer(length(ref_2020)/12), " genomas \n")
cat("Procesando ", as.integer(length(ref_2021)/12), " genomas \n")
cat("Procesando ", as.integer(length(ref_2022)/12), " genomas \n")
cat("Procesando ", as.integer(length(ref_2024)/12), " genomas \n")
cat("Procesando ", as.integer(length(fMx5)/12), " genomas \n")

nObs1 = 1
nObs2 = 1
nObs3 = 1
nObs4 = 1
nObs5 = 1



# mutacion 1 
for (i in seq(1,length(ref_2020),1)){
  if (i==2) next
  anotaciones = attr(ref_2020[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN(ref_2020[[i]])  
  cat("#",geneName)
  for (k in seq(i, length(ref_2020), 12)){
    
    genfRMx = ToARN( ref_2020[[k]] )  
    cat(i, k, length(genRef), length(genfRMx), "\n")
    
    if (length(genRef) == length(genfRMx)){
      dif = which(genRef != genfRMx)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genfRMx[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO1 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2020 = paste(genfRMx[inicioCodon], genfRMx[inicioCodon+1], genfRMx[inicioCodon+2],sep="")
          codonChange = paste(codonO1,"to",codonfRMx, sep="")
          aminoChange = paste(trad[codonO1],numCodon,trad[codonfRMx], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfRMx]) && trad[codonO1]!=trad[codonfRMx]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df_mut1[nObs1,] = obs 
            nObs1 = nObs1+1
          }
        }
      }
    }else{
    }
  }
}


#mutacion 2 
for (i in seq(1,length(ref_2021),1)){
  if (i==2) next
  anotaciones = attr(ref_2021[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2021[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2021), 12)){
    
    genref_2021 = ToARN( ref_2021[[k]] )  
    cat(i, k, length(genRef), length(genref_2021), "\n")
    
    if (length(genRef) == length(genref_2021)){
      dif = which(genRef != genref_2021)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genref_2021[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO2 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2021 = paste(genref_2021[inicioCodon], genref_2021[inicioCodon+1], genref_2021[inicioCodon+2],sep="")
          codonChange = paste(codonO2,"to",codonref_2021, sep="")
          aminoChange = paste(trad[codonO2],numCodon,trad[codonref_2021], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2021]) && trad[codonO2]!=trad[codonref_2021]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df_mut2[nObs2,] = obs 
            nObs2 = nObs2+1
          }
        }
      }
    }else{
    }
  }
}


#mutacion 3 
for (i in seq(1,length(ref_2022),1)){
  if (i==2) next
  anotaciones = attr(ref_2022[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2022[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2022), 12)){
    
    genref_2022 = ToARN( ref_2022[[k]] )  
    cat(i, k, length(genRef), length(genref_2022), "\n")
    
    if (length(genRef) == length(genref_2022)){
      dif = which(genRef != genref_2022)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genref_2022[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO3 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2022 = paste(genref_2022[inicioCodon], genref_2022[inicioCodon+1], genref_2022[inicioCodon+2],sep="")
          codonChange = paste(codonO3,"to",codonref_2022, sep="")
          aminoChange = paste(trad[codonO3],numCodon,trad[codonref_2022], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2022]) && trad[codonO3]!=trad[codonref_2022]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df_mut3[nObs3,] = obs 
            nObs3 = nObs3+1
          }
        }
      }
    }else{
    }
  }
}


#mutacion 4
for (i in seq(1,length(ref_2024),1)){
  if (i==2) next
  anotaciones = attr(ref_2024[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( ref_2020[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(ref_2024), 12)){
    
    genref_2024 = ToARN( ref_2024[[k]] )  
    cat(i, k, length(genRef), length(genref_2024), "\n")
    
    if (length(genRef) == length(genref_2024)){
      dif = which(genRef != genref_2024)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genref_2024[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO4 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonref_2024 = paste(genref_2024[inicioCodon], genref_2024[inicioCodon+1], genref_2024[inicioCodon+2],sep="")
          codonChange = paste(codonO4,"to",codonref_2024, sep="")
          aminoChange = paste(trad[codonO4],numCodon,trad[codonref_2024], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonref_2024]) && trad[codonO4]!=trad[codonref_2024]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df[nObs4,] = obs 
            nObs4 = nObs4+1
          }
        }
      }
    }else{
    }
  }
}


head(df)
nrow(df)

library(ggplot2)
#


# Graficar Mutaciones 1..2..3..4


p = ggplot(df)
p = p + aes(x=Mutation, fill=Mutation, label=after_stat(count))
p = p + ggtitle("Mutaciones de sustitución")
p = p + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
p = p + geom_bar(stat = "count")
p = p + geom_text(stat = "count", vjust=1.5)
#p = p + facet_grid(~Gene)
p



# Graficas:

library(dplyr)
dfgraph = filter(
  summarise(
    select(
      group_by(df, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

dfgraph = dfgraph[order(-dfgraph$Cuenta), ]
dfgraph = dfgraph[1:20, ]
head(dfgraph)
nrow(dfgraph)



p2 = ggplot(dfgraph)
p2 = p2 + aes(x=Amino, y=Cuenta, fill=Amino, label=Cuenta)
p2 = p2 + ggtitle("Cambio de Aminoácidos")
p2 = p2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
p2 = p2 + geom_bar(stat = "identity")
p2 = p2 + geom_text(stat = "identity", vjust=1.5)
p2 = p2 + facet_grid(~Gene, scales="free", space="free_x")
p2