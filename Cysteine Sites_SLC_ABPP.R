#Cys data analysis for Miljan K's
# Nat Biotech paper

library(tidyverse)
library(httr)
library(jsonlite)
library(xml2)
setwd("C:\\Users\\Mintaka\\Dropbox (Personal)\\UW\\Writing\\Miljan_NatBiotech")
cysDf = read.table("natb_mk_rev_proteins.txt",header=T)

full.list = list()
baseUrl = "https://www.ebi.ac.uk/proteins/api/proteins/"
# baseUrl.features = "https://www.ebi.ac.uk/proteins/api/features/"

for(i in 1:nrow(cysDf)){
  accession = cysDf[i,1]
  print(accession)
  requestURL = paste0(baseUrl,accession)
  proteinInfo = GET(requestURL, accept("application/json"))
  
  warn_for_status(proteinInfo)

  json = toJSON(content(proteinInfo))
  full.list[[accession]] = fromJSON(json)
}


cys.sites = read.csv("Cysteine Sites_SLC_ABPP.csv",header = T)

cys.sites %>%
  filter(Site.Position == 205,Accession == "Q9Y5T5")

checkFeatureRes = function(input, residue = 205, 
                           feature.types = c("ACT_SITE","SITE","MOD_RES","DISULFID","CROSSLNK"),
                           return.types = TRUE,
                           return.desc = TRUE,...){
  if(is.null(input)){
    return("")
  }
  if(return.desc){
    return.types =FALSE
  }
  
  if(length(feature.types) <=1){
    temp.input = input %>%
      filter(begin == residue,end == residue)
  } else {
    temp.input = input %>%
      filter(begin == residue,end == residue, type %in% feature.types)
  }

  if(return.desc){
    return(paste(temp.input$description,collapse = ","))
  }else if(return.types){
    return(paste(temp.input$type,collapse = ","))
  }else{
    return(temp.input)
  }
}

cys.sites.ann = cys.sites
cys.sites.ann$site.desc = apply(cys.sites,1,FUN = function(x){
  acc = x["Accession"]
  res = x["Site.Position"]
  full.list[[acc]]$features %>% 
    checkFeatureRes(res)
  })

cys.sites.ann$site.type = apply(cys.sites,1,FUN = function(x){
  acc = x["Accession"]
  res = x["Site.Position"]
  full.list[[acc]]$features %>% 
    checkFeatureRes(res,return.desc = F)
})

cys.sites.ann$site.desc.any = apply(cys.sites,1,FUN = function(x){
  acc = x["Accession"]
  res = x["Site.Position"]
  full.list[[acc]]$features %>% 
    checkFeatureRes(res,feature.types=vector())
})

cys.sites.ann$site.type.any = apply(cys.sites,1,FUN = function(x){
  acc = x["Accession"]
  res = x["Site.Position"]
  full.list[[acc]]$features %>% 
    checkFeatureRes(res,return.desc = F,feature.types=vector())
})
 
cys.table = as.data.frame(table(cys.sites.ann$site.desc))
sum(cys.table[-1,"Freq"]) / sum(cys.table[,"Freq"])

cys.type.table = as.data.frame(table(cys.sites.ann$site.type))
sum(cys.type.table[-1,"Freq"]) / sum(cys.type.table[,"Freq"])

cys.table.any = as.data.frame(table(cys.sites.ann$site.desc.any))
sum(cys.table.any[-1,"Freq"]) / sum(cys.table.any[,"Freq"])

cys.type.table.any = as.data.frame(table(cys.sites.ann$site.type.any))
sum(cys.type.table.any[-1,"Freq"]) / sum(cys.type.table.any[,"Freq"])

write.csv(cys.sites.ann,"output.csv")

# sub cellular localization
full.list$P00533$comments[full.list$P00533$comments$type == "SUBCELLULAR_LOCATION","locations"]

checkComments = function(input, query = "membrane", comment.type = "SUBCELLULAR_LOCATION",...){
  if(is.null(input)){
    return("")
  }
  x=unlist(input$comments[input$comments$type == "SUBCELLULAR_LOCATION","locations"][[1]]$location$value)
  return(any(str_detect(tolower(x),query)))
}

mem.df = apply(cys.sites,1,FUN = function(x){
  acc = x["Accession"]
  full.list[[acc]] %>% 
    checkComments()
})

cys.sites.wMem = cys.sites.ann
cys.sites.wMem$is.mem.prot = mem.df


ggplot(cys.sites.wMem,aes(is.mem.prot,fill = is.mem.prot)) +
  geom_bar(stat = "count") + xlab("Site on Membrane Protein") +
  geom_text(stat = "count",aes(label=..count..),vjust=-1) +
  scale_fill_manual(values = c("darkgrey","darkred")) +
  theme_bw()

write.csv(cys.sites.wMem,"output_withMem.csv")
