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
fMx1 = read.fasta("mx1_sequence")
length(fMx1)

# Mutaciones México 2 : 5
fMx2 = read.fasta("mx2_sequence")
length(fMx2)

fMx3 = read.fasta("mx3_sequence")
length(fMx3)

fMx4 = read.fasta("mx4_sequence")
length(fMx4)

fMx5 = read.fasta("mx5_sequence")
length(fMx5)


cat("Procesando ", as.integer(length(fMx1)/12), " genomas \n")
cat("Procesando ", as.integer(length(fMx2)/12), " genomas \n")
cat("Procesando ", as.integer(length(fMx3)/12), " genomas \n")
cat("Procesando ", as.integer(length(fMx4)/12), " genomas \n")
cat("Procesando ", as.integer(length(fMx5)/12), " genomas \n")

nObs1 = 1
nObs2 = 1
nObs3 = 1
nObs4 = 1
nObs5 = 1



# mutacion 1 
for (i in seq(1,length(fMx1),1)){
  if (i==2) next
  anotaciones = attr(fMx1[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN(fMx1[[i]])  
  cat("#",geneName)
  for (k in seq(i, length(fMx1), 12)){
    
    genfRMx = ToARN( fMx1[[k]] )  
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
          codonfMx1 = paste(genfRMx[inicioCodon], genfRMx[inicioCodon+1], genfRMx[inicioCodon+2],sep="")
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
for (i in seq(1,length(fMx2),1)){
  if (i==2) next
  anotaciones = attr(fMx2[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fMx2[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(fMx2), 12)){
    
    genfMx2 = ToARN( fMx2[[k]] )  
    cat(i, k, length(genRef), length(genfMx2), "\n")
    
    if (length(genRef) == length(genfMx2)){
      dif = which(genRef != genfMx2)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genfMx2[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO2 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfMx2 = paste(genfMx2[inicioCodon], genfMx2[inicioCodon+1], genfMx2[inicioCodon+2],sep="")
          codonChange = paste(codonO2,"to",codonfMx2, sep="")
          aminoChange = paste(trad[codonO2],numCodon,trad[codonfMx2], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfMx2]) && trad[codonO2]!=trad[codonfMx2]){
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
for (i in seq(1,length(fMx3),1)){
  if (i==2) next
  anotaciones = attr(fMx3[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fMx3[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(fMx3), 12)){
    
    genfMx3 = ToARN( fMx3[[k]] )  
    cat(i, k, length(genRef), length(genfMx3), "\n")
    
    if (length(genRef) == length(genfMx3)){
      dif = which(genRef != genfMx3)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genfMx3[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO3 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfMx3 = paste(genfMx3[inicioCodon], genfMx3[inicioCodon+1], genfMx3[inicioCodon+2],sep="")
          codonChange = paste(codonO3,"to",codonfMx3, sep="")
          aminoChange = paste(trad[codonO3],numCodon,trad[codonfMx3], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfMx3]) && trad[codonO3]!=trad[codonfMx3]){
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
for (i in seq(1,length(fMx4),1)){
  if (i==2) next
  anotaciones = attr(fMx4[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fMx4[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(fMx4), 12)){
    
    genfMx4 = ToARN( fMx4[[k]] )  
    cat(i, k, length(genRef), length(genfMx4), "\n")
    
    if (length(genRef) == length(genfMx4)){
      dif = which(genRef != genfMx4)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genfMx4[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonO4 = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfMx4 = paste(genfMx4[inicioCodon], genfMx4[inicioCodon+1], genfMx4[inicioCodon+2],sep="")
          codonChange = paste(codonO4,"to",codonfMx4, sep="")
          aminoChange = paste(trad[codonO4],numCodon,trad[codonfMx4], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfMx4]) && trad[codonO4]!=trad[codonfMx4]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df_mut4[nObs4,] = obs 
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