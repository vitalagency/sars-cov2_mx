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

# Sequence Wuhan
fRef = read.fasta("mx1_sequence")
length(fRef)

# Mutaciones México 
fMx1 = read.fasta("mx2_sequence")
length(fMx1)

fMx2 = read.fasta("mx3_sequence")
length(fMx2)

fMx3 = read.fasta("mx4_sequence")
length(fMx3)

fMx4 = read.fasta("mx5_sequence")
length(fMx4)


cat("Procesando ", as.integer(length(fRef)/12), " genomas \n")

nObs = 1


# mutacion 1 
for (i in seq(1,length(fRef),1)){
  if (i==2) next
  anotaciones = attr(fRef[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fRef[[i]] )  
  genMut = ToARN(fMut1[[i]]) # fMut2[[i]], fMut3[[i]], fMut4[[i]]) Para las demás mutaciones
  cat("#",geneName)
  for (k in seq(i, length(fB117), 12)){
    
    genfB117 = ToARN( fRef[[k]] )  
    cat(i, k, length(genRef), length(genfB117), "\n")
    
    if (length(genRef) == length(genfB117)){
      dif = which(genRef != genfB117)
      cat("length",length(dif))
      print(dif)
      
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genfB117[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfB117 = paste(genfB117[inicioCodon], genfB117[inicioCodon+1], genfB117[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonfB117, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonfB117], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfB117]) && trad[codonOri]!=trad[codonfB117]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df[nObs,] = obs 
            nObs = nObs+1
          }
        }
      }
    }else{
    }
  }
}


#mutacion 2 
for (i in seq(1,length(fRef),1)){
  if (i==2) next
  anotaciones = attr(fRef[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fRef[[i]] )  
  genMut = ToARN(fMut2[[i]]) # fMut2[[i]], fMut3[[i]], fMut4[[i]]) Para las demás mutaciones
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
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfB117 = paste(genfMx2[inicioCodon], genfMx2[inicioCodon+1], genfMx2[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonfMx2, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonfMx2], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfMx2]) && trad[codonOri]!=trad[codonfMx2]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df[nObs,] = obs 
            nObs = nObs+1
          }
        }
      }
    }else{
    }
  }
}


#mutacion 3 
for (i in seq(1,length(fRef),1)){
  if (i==2) next
  anotaciones = attr(fRef[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fRef[[i]] )  
  genMut = ToARN(fMut3[[i]]) # fMut2[[i]], fMut3[[i]], fMut4[[i]]) Para las demás mutaciones
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
          muta = paste(genRef[x],"to",genfBMx3[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfB117 = paste(genfMx3[inicioCodon], genfMx3[inicioCodon+1], genfMx3[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonfMx3, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonfMx3], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfMx3]) && trad[codonOri]!=trad[codonfMx3]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df[nObs,] = obs 
            nObs = nObs+1
          }
        }
      }
    }else{
    }
  }
}


#mutacion 4
for (i in seq(1,length(fRef),1)){
  if (i==2) next
  anotaciones = attr(fRef[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fRef[[i]] )  
  genMut = ToARN(fMut4[[i]]) # fMut2[[i]], fMut3[[i]], fMut4[[i]]) Para las demás mutaciones
  cat("#",geneName)
  for (k in seq(i, length(fB117), 12)){
    
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
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfB117 = paste(genfMx4[inicioCodon], genfMx4[inicioCodon+1], genfMx4[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonfMx4, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonfMx4], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfMx4]) && trad[codonOri]!=trad[codonfMx4]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df[nObs,] = obs 
            nObs = nObs+1
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