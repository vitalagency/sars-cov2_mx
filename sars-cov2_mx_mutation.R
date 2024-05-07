
# Mexico Secuencias 2020 - 2024

library(seqinr)
library(dplyr)
library(ggplot2)

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

df2023_2024 = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)


ref_2020 = read.fasta("ref_2020.txt")
length(ref_2020)

ref_2021 = read.fasta("ref_2022.txt")
length(ref_2021)

ref_2022 = read.fasta("ref_2022")
length(ref_2022)

ref_2023 = read.fasta("ref_2023")
length(ref_2023)

ref_2024 = read.fasta("ref_2024")
length(ref_2024)
##


cat("Procesando ", as.integer(length(ref_2020)/12), " genomas de Primera Secuencia 2020 y Primera Secuencia 2021 \n")

nObs = 1
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
            df2020_2021[nObs,] = obs 
            nObs = nObs+1
          }
        }
      }
    }else{
    }
  }
}

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


