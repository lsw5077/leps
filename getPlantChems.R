library(pdftools)
library(stringr)
library(qdapRegex)
library(dplyr)
library(data.table)

### TODOs
# 1. dblchk that everything in caps is a genera or chem_class
# 2. many chem_classes are hierarchical, ie. ALKALOIDS and *specific sub-class* ALKALOID;
#           We could just lump everything into overall chem_class
# 3. ALso, need to handle singular (ALKALOID) vs plural (ALKALOIDS)

# read in data from pdf
# each page is one element in the list "text"
text <- pdf_text("data/dictionary_plant_metabolites.pdf")  
text <- text[8:1290] # drop preface and index of book

# All genera names and chemical classes are CAPITALIZED in the pdf. We can use this pattern to
# extract them

# pulling out caps like this includes a lot more than just chem class genera
# using ex_caps also splits up phrases like "AMINO ACID" into "AMINO", "ACID"
caps <- ex_caps(text[[4]])
caps 

# Use ex_caps_phrase to preserve phrases like "AMINO ACID"
caps <- ex_caps_phrase(text[[4]])
caps

### Clean up cap phrase data
# ------
# get all cap phrases
caps <- unlist(sapply(1:length(text), function(x) ex_caps_phrase(text[[x]])))

# however, there are still lots of elements that are not genera or chem_class
caps_length <- nchar(caps)         
caps[head(order(caps_length),500)] # most short words are nonsense

# there are also NAs
caps[which(is.na(caps))] 

# check in book to see where NAs are introduced
caps[which(is.na(caps))-1]
caps[which(is.na(caps))+1] # NAs are introduced whenever there is a blank page

# so we can safely drop NAs from caps_vec
caps <- caps[-c(which(is.na(caps)))]
caps_length <- nchar(caps)

# look at short cap words
which(caps_length < 4) # there are lot. Most are not genera or chem_class

unique(caps[which(caps_length == 2)]) # no valid words are nchar(caps)==2

# drop nchar(caps) < 3
caps <- caps[which(caps_length > 2 )]
caps_length <- nchar(caps)
unique(caps[which(caps_length == 2)]) # all removed

unique(caps[which(caps_length == 3)]) # only ASA, IVA, ZEA is valid (someone else dbl check)

word3 <- unique(caps[which(caps_length == 3)])
word3 <- word3[which(!(word3 %in% c("ASA","IVA","ZEA")))]

caps <- caps[-which((caps %in% word3))] 
caps_length <- nchar(caps)

caps[which(caps_length==3)] # confirm nchar(caps)==3 are correct

unique(caps[which(caps_length == 4)]) # 
word4 <- rm_non_words(unique(caps[which(caps_length == 4)])) # removes hyphen
word4 <- unlist(rm_nchar_words(word4,n=4))   # get all words less than 4 char
word4 <- append(word4[which(word4 != "")],c("USSR","SSSR","IRCS","ATCC","XVII", "XLII", "XLIV", "XLVI", 
           "LIII","VIII","XIII","IIIB")) # get all unwanted nchar(caps)==4
word4
caps <- caps[-which((caps %in% word4))] 
caps_length <- nchar(caps)
caps[which(caps_length==4)] # still have non words

word4 <- caps[which(caps_length==4)][which(grepl("-",caps[which(caps_length==4)]))]
caps <- caps[-which((caps %in% word4))] 
caps_length <- nchar(caps)
unique(caps[which(caps_length==4)]) # all nchar(caps)==4 are fixed

unique(caps[which(caps_length == 5)]) # all unwanted words have hyphens

word5 <- caps[which(caps_length==5)][which(grepl("-",caps[which(caps_length==5)]))] 
word5 # everything with hyphens

caps <- caps[-which((caps %in% word5))] 
caps_length <- nchar(caps)
unique(caps[which(caps_length == 5)]) # should dblcheck all are genus or chem

unique(caps[which(caps_length == 6)]) # still something with hyphens
word6 <- caps[which(caps_length==6)][which(grepl("-",caps[which(caps_length==6)]))]
word6

caps <- caps[-which((caps %in% word6))] 
caps_length <- nchar(caps)
unique(caps[which(caps_length == 6)]) # should dblcheck

# havent seen unwanted words after nchar(caps) >= 7 (I've checked by eye for words up to 10)
unique(caps[which(caps_length == 7)]) ### someone could do more checking here

# (tentatively) caps ONLY contains genera and chem_class info

# get caps, chemical classes, and genera from every page in entire book
caps <- setDT(list(caps))
chem_class <- setDT(list(unlist(lapply(1:length(text), 
                                       function(x) ex_caps_phrase(ex_between(text[[x]],
                                                                             "\r\n",":"))))))
chem_class <- chem_class[-which(is.na(chem_class))] # rm NAs caused by empty pages in book
genera <- caps[!caps$V1 %in% chem_class$V1] 

# get genera using plant genera reference 
read.csv("data/plantGenera.csv", header = TRUE, stringsAsFactors = FALSE) %>%
  transmute(Genus = toupper(genus)) -> plantGenera

genera2 <- (caps[(caps$V1 %in% plantGenera$Genus)]) # some genera missing in plantGenera

# 935 genera not in plantGenera; all appear to be valid genera names
genDiff <- setdiff(genera[[1]],genera2[[1]]) 

# Some differences are typos, e.g., ABELOM*O*SCHUS instead of ABELMOSCHUS
# make df that has all Genera from dictionary + "difference" check if also in plantGenera
genera %>% mutate(different = ifelse(genera$V1 %in% plantGenera$Genus, "no","yes")) -> generaDiff
names(generaDiff)[1] <- "genus"

write.csv(generaDiff, "data/generaDiff.csv", row.names = FALSE)


### after cleaning up caps, chem_class, and genera 
# create dataframe with chemical classes as columns and genera as rows

defense <- data.frame(matrix(0, nrow = nrow(genera), ncol = nrow(unique(chem_class))+1))
names(defense) <- c("Genus", unique(chem_class$V1)) 

defense[,1] <- genera$V1
defense <- data.table(defense)

#This is clunky, but not too slow. It a loop that categories each capitalized word in book as a genus or a chemical family and fills out the matrix accordingly
i = 0
for(n in 1:length(caps$V1)){
  if(is.element(caps$V1[n], genera$V1)){
    i = i +1}
  if(is.element(caps$V1[n], genera$V1) == FALSE){
    defense[i, caps$V1[n]:=1]}
}

defense<-as.data.frame(defense)
defense[,2:ncol(defense)]<-sapply(defense[,2:ncol(defense)], as.numeric)
write.csv(defense, "data/PlantChems.csv", row.names = FALSE)
