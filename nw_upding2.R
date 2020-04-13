# 20191212
library(tidyverse)

# producing urls to dowload equations
make_urls <- readr::read_lines("~/Documents/5gatom_draft/gatom_network_updating/data/current_KEGG_reaction_20191212.txt")
head(make_urls)
new_lines <- c()
for(line in make_urls){
  get_id <- stringr::str_extract(line, "rn:R[:digit:]{5}")
  id <- gsub("rn:", "", get_id)
  # new_lines <- c(new_lines, id)
  new_lines <- c(new_lines, paste0("http://rest.kegg.jp/get/", id))
}
head(new_lines)
readr::write_lines(new_lines, path = "~/Documents/5gatom_draft/gatom_network_updating/data/urls")


# cheking which one were not dowloaded
make_urls <- readr::read_lines("~/Documents/5gatom_draft/gatom_network_updating/data/current_KEGG_reaction_20191212.txt")
head(make_urls)
new_lines <- c()
for(line in make_urls){
  get_id <- stringr::str_extract(line, "rn:R[:digit:]{5}")
  id <- gsub("rn:", "", get_id)
  new_lines <- c(new_lines, id)
}
head(new_lines)
downloaded <- list.files("~/Documents/5gatom_draft/gatom_network_updating/data/reactions/")
head(downloaded)
diff <- setdiff(new_lines, downloaded); diff
cat(paste0("http://rest.kegg.jp/get/", diff, "\n"))


# making eautions.RDa
lines <- readr::read_lines("~/Documents/5gatom_draft/gatom_network_updating/data/equations.txt")
head(lines)
equations <- lines %>%
  str_split_fixed(":", 2) %>% {
    setNames(.[,1], .[,2])
  }
glimpse(equations)
# cheking which one were not dowloaded
diff <- setdiff(new_lines, unname(equations)); diff
length(names(equations))
length(new_lines)
which(nchar(names(equations)) < 6)
readr::write_rds(equations, "~/Documents/5gatom_draft/gatom_network_updating/data/equations.rds")


# downloading compounds url
make_urls <- readr::read_lines("~/Documents/5gatom_draft/gatom_network_updating/data/current_KEGG_compound_20191212.txt")
head(make_urls)
new_lines <- c()
for(line in make_urls){
  get_id <- stringr::str_extract(line, "cpd:C[:digit:]{5}")
  id <- gsub("cpd:", "", get_id)
  new_lines <- c(new_lines, id)
  # new_lines <- c(new_lines, paste0("http://rest.kegg.jp/get/", id, "/mol"))
}
head(new_lines)
write_lines(new_lines, path = "~/Documents/5gatom_draft/gatom_network_updating/data/compound_list_20191212.txt")


# checking whether all compounds from compound_list_20191212.txt   are presentes as files in KEGG_MOL_files 
downloaded <- list.files("~/Documents/5gatom_draft/gatom_network_updating/data/KEGG_MOL_files/")
head(downloaded)
diff <- setdiff(new_lines, downloaded); diff
new_lines[17888]
downloaded[17888]


################## produce .rxn #########################
eqs <- readRDS("~/Documents/5gatom_draft/gatom_network_updating/data/equations.rds")
head(eqs)
length(eqs) # 11334
pattern <- "([:digit:]*)\\s*([A-za-z][:digit:]{5})"
pattern_mult <- "([:digit:]+)\\s*(C[:digit:]{5})"
reaction_error <- c()

# rr <- 6831
for(rr in seq_along(eqs)){
  
  reaction <- names(eqs[rr]); reaction
  reaction_id <- unname(eqs[rr]); reaction_id
  
  
  # ------------ replacing "(n+1)" with "2" -----------------
  if(stringr::str_detect(reaction, "\\(n\\+1\\)")){
    reaction <- gsub("\\(n\\+1\\)", 2, reaction)
    reaction <- gsub(" n", "", reaction)
  }
  print(paste(reaction_id, "________", reaction))
  
  # ------------ dealing with multiple compounds -----------------
  m <- stringr::str_extract_all(reaction, pattern) #; m
  # checking <- stringr::str_detect(m, "([:digit:]+)\\s*(C[:digit:]{5})")
  # checking <- stringr::str_detect(unlist(m), "([:digit:]+)\\s*(C[:digit:]{5})")
  need_to_unmult <- stringr::str_extract(unlist(m), pattern_mult) #; need_to_mult
  
  if(stringr::str_detect(m, pattern)){
    
    if(any(!is.na(need_to_unmult))){
      for(not_na in need_to_unmult[!is.na(need_to_unmult)]){
        # print(not_na)
        get_elemens <- stringr::str_match(not_na, pattern_mult); get_elemens
        reaction <- gsub(not_na,
                         paste(rep(get_elemens[3], get_elemens[2]), collapse = " + "),
                         reaction) #; reaction
      }
    }
  }
  # reaction
  
  # ------ extract compound IDs -------
  comp_ids <- unlist(stringr::str_extract_all(reaction, "[A-za-z][:digit:]{5}")) #; comp_ids
  
  if(all(file.exists(
    paste0("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/KEGG_MOL_files/",
           comp_ids)))){
    
    # ------ detect how many reactants and products -------
    arrow <- unlist(strsplit(reaction, "<=>"))
    n <- stringr::str_count(arrow[1], pattern = "C[:digit:]{5}")
    m <- stringr::str_count(arrow[2], pattern = "C[:digit:]{5}")
    
    # ------ traslate reaction in RXN format and save -----
    content_of_output <- c("$RXN", "", "  WHATEVER     blabla", "",
                           paste0("  ", n, "  ", m))
    
    for(comp in comp_ids){
      
      # try to extract all mol files, if not --> next reaction
      comp_mol <- readr::read_lines(
        paste0("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/KEGG_MOL_files/", comp))
      
      where_stop <- which(grepl("> <ENTRY>", comp_mol)) - 1
      id <- gsub("cpd:", "", comp_mol[grepl("cpd:C", comp_mol)])
      
      content_of_output <- c(content_of_output,
                             "$MOL", 
                             id, "  WHATEVER  0000000000", 
                             comp_mol[-c(1, 2, where_stop:length(comp_mol))])
    }
    # ------ save file ------ 
    readr::write_lines(content_of_output, 
                       path = paste0("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/KEGG_IN_RXN_20191212/", 
                                     reaction_id, ".rxn"))  }
  else {
    reaction_error <<- c(reaction_error, reaction_id)
  }
}
length(reaction_error) # 1 626 from 
head(reaction_error)
readr::write_lines(reaction_error, 
                   path = 
                     paste0("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/reactionsNotIn_KEGG_IN_RXN_20191212.txt"))

'''
Warning messages:
  1: In stri_detect_regex(string, pattern, negate = negate, opts_regex = opts(pattern)) :
  argument is not an atomic vector; coercing
'''
################## produce .rxn #########################


input_list <- list.files(path = "/home/octopus/Documents/5gatom_draft/my/KEGG_IN_RXN_20191212", pattern = "*.rxn")
input_list <- gsub(".rxn", "", input_list)
length(input_list) # 9708
output_list <- list.dirs("/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/", full.names = F)
output_list <- output_list[-1]
length(output_list) # 9708
diff <- setdiff(input_list, output_list); diff

# they all are without files inside, no foldersw with x=1,2 files
# 59 filed
failed_RDT1 <- read.table("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/failed_RDT_20191212.txt")
failed_RDT <- read.table("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/failed_RDT_2_20191212.txt")
failed_RDT1 <- stringr::str_sort(gsub("/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/", "", failed_RDT1$V1))
failed_RDT1
length(failed_RDT) # 59
readr::write_lines(failed_RDT, path = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/failed_RDT_20191212.txt")



### analyse RDT output
# work with those where not less than 3 files
length(list.dirs("/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/")) # 9709

RDT_outputs <- list.dirs("/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/")
head(RDT_outputs)
RDT_outputs <- RDT_outputs[-1]
head(RDT_outputs)


############# #######################
# address <- "/home/octopus/Documents/5gatom_draft/RDT_outputs/R00001/ECBLAST_R00001_AAM.rxn"
address <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs/R05314/ECBLAST_R05314_AAM.rxn" # single hydrogen
address <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs/R00001/ECBLAST_R00001_AAM.rxn"
address <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs/R01492/ECBLAST_R01492_AAM.rxn" # 99 

address <- "/home/octopus/Documents/5gatom_draft/my/ECBLAST_R01492_AAM.txt"
address <- "/home/octopus/Documents/5gatom_draft/my/ECBLAST_R05314_AAM.txt"

# R00001 with double thing
# R00006 R00008 R00014 
# C00022

rdt_lines <- readr::read_lines(address)

rdt <- ChemmineR::read.SDFindex(address, 
                                index = data.frame(
                                  "A" = which(stringr::str_detect(rdt_lines, "C[:digit:]{5}")), 
                                  "B" = which(stringr::str_detect(rdt_lines, "END"))
                                ))

ttt <- ChemmineR::read.SDFindex(address,
                                index=data.frame(A=1, 
                                                 B=length(readr::read_lines(address))))
ChemmineR::atomblock(ttt)


ChemmineR::bondblock(rdt)
ChemmineR::header(rdt)
ChemmineR::atomblock(rdt)
ChemmineR::sdfid(rdt)
View(ChemmineR::atomblock(rdt)$CMP1)
View(ChemmineR::atomblock(rdt)[[1]])
############# #######################

library(ChemmineR)
# 


RDT_outputs <- list.dirs("/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212")
RDT_outputs <- RDT_outputs[-1]
RDT_foutputs <- readr::read_lines(
  "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/failed_RDT_2_20191212.txt")
# RDT_outputs[1]
# length(RDT_foutputs)
# length(RDT_outputs)
# head(RDT_outputs)
RDT_outputs <- RDT_outputs[-which(RDT_outputs %in% RDT_foutputs)]
length(RDT_outputs) # 9651


dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R00029"
# dir <- RDT_outputs[11] # 16: R00023
# dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R02480"
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R02480"
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R05314" # одиночный водород
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R00011" # одиночное всё
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R09339" 

dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R05314" # как получается NA в неодиночном водороде?
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R10000" # звёздочки
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R05110" # больше 100
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R00075" 
dir <- "/home/octopus/Documents/5gatom_draft/my/RDT_outputs_20191212/R00078" # одиночным бывает не только Н


final_mapping <- dplyr::tibble()

ncolnot15 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(ncolnot15) <- c("react_id", "number_of_not")

withHydros <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(withHydros) <- c("react_id", "number_of_pos")

hydro_true <- FALSE


for(dir in RDT_outputs){
  # print(dir)
  # >>1<< not three files
  rxn <- list.files(dir, pattern = ".rxn$", full.names = T)
  rdt_lines <- readr::read_lines(rxn)
  rdt <- ChemmineR::read.SDFindex(rxn, 
                                  index = data.frame(
                                    "A" = which(stringr::str_detect(rdt_lines, "[CMG][:digit:]{5}")), 
                                    "B" = which(stringr::str_detect(rdt_lines, "END"))
                                    # type = "file"
                                  ))
  # atomblock(rdt)
  
  react_id <- stringr::str_extract(rdt_lines[3], "R[:digit:]{5}")
  print(react_id)
  
  if(any(sapply(atomblock(rdt), dim)[2,] == 2)){
    hydro_pos <- which(sapply(atomblock(rdt), dim)[2,] == 2)
    hydro_true <- TRUE
    withHydros <- rbind(withHydros, 
                        data.frame("react_id" = react_id, "number_of_pos" = length(hydro_pos)))
    rdt <- ChemmineR::read.SDFindex(rxn,
                                    index = data.frame(
                                    "A" = which(stringr::str_detect(rdt_lines, "[CMG][:digit:]{5}"))[-hydro_pos],                                         
                                    "B" = which(stringr::str_detect(rdt_lines, "END"))[-hydro_pos]
                                    # type = "file"
                                    ))
  }
    
  # if(any(!sapply(atomblock(rdt), dim)[2,] %in% c(2, 15))){  # :(
  if(any(sapply(atomblock(rdt), dim)[2,] != 15)){ # instead of "<"
    ncolnot15 <- rbind(ncolnot15, 
                       data.frame("react_id" = react_id, 
                                  "number_of_not" = length(which(sapply(atomblock(rdt), dim)[2,] != 15))))
    hydro_true <- FALSE
  } else {
    nn <- stringr::str_match_all(rdt_lines[5], "[:digit:]+"); nn 
    n_reactants <- nn[[1]][1]
    # n_products <- nn[[1]][2]
    list_of_reactants <- data.frame()
    list_of_products <- data.frame()
    
    if(hydro_true){
      hydro_true <- FALSE
      pos_reactants <- seq(as.numeric(n_reactants))
      pos_to_remove <- c()
      for(pos in hydro_pos){
        # print(pos)
        if(pos <= n_reactants){
          pos_to_remove <- c(pos_to_remove, pos)
        }
      }
      if(length(pos_to_remove) > 0){
        pos_reactants <- pos_reactants[-pos_to_remove]
        n_reactants <- length(pos_reactants)
      }
    }
    
    if(n_reactants == 0 | n_reactants == length(rdt)){
      withHydros <- rbind(withHydros, 
                          data.frame("react_id" = react_id, "number_of_pos" = "turned zero"))
      next
    }
    
    for(el in seq_along(rdt)){
      # el <- 10
      curr <- atomblock(rdt)[[el]]
      doing_cid <- sdfid(rdt)[[el]]
      
      doing_aid <- paste(doing_cid,
                         gsub("_.*", "", rownames(curr)), curr[, 1], curr[, 2], 
                         sep = "_")
      
      if(el <= n_reactants){
        list_of_reactants <- rbind(list_of_reactants,
                                   data.frame("reaction" = react_id,
                                              "reactant" = doing_aid,
                                              "num" = unname(curr[, 13])))
      } else{
        list_of_products <- rbind(list_of_products,
                                  data.frame("reaction" = react_id,
                                             "product" = doing_aid,
                                             "num" = unname(curr[, 13])))
      }
    }
    # list_of_reactants
    # list_of_products
    final_mapping <- dplyr::bind_rows(final_mapping,
                                      dplyr::full_join(list_of_reactants, 
                                                       list_of_products, 
                                                       by = c("reaction", "num")))
    # View(final_mapping)
  }
}



View(final_mapping)

rUtils::write.tsv(final_mapping, 
                  file = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/final_mapping_20200311.tsv")
save(final_mapping, 
     file = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/final_mapping_20200311.Rda")
load("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/final_mapping_20200311.Rda")


# rUtils::write.tsv(not_successful_reactions, 
#                   file = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/not_successful_reactions_20191212.tsv")
# save(not_successful_reactions, 
#      file = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/not_successful_reactions_20191212.Rda")
# load("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/not_successful_reactions_20191212.Rda")
# >>>>> delete NAs <<<<<, decide with !=  15


############ looking at my lists ###########
head(ncolnot15) # R00017
head(withHydros) # check for NA presence when hydros mapping is excluded from mapping 

dim(ncolnot15)  # 307
dim(withHydros) # 2179
View(ncolnot15)  
View(withHydros)

# H + 15 = 2311 (the same which is when "<")

length(intersect(withHydros$react_id, ncolnot15$react_id)) # 172
############ looking at my lists ###########


# final_mapping
dim(final_mapping[complete.cases(final_mapping), ]) # 469 780      4
dim(final_mapping[!complete.cases(final_mapping), ]) # 14 337      4

dim(final_mapping) # 484 117      4
View(final_mapping)
View(head(final_mapping, 40))
View(final_mapping[complete.cases(final_mapping), ])

length(unique(final_mapping$reaction)) # 9 341 - take those with are with enzyme, just rigth here (to have all mapping)


files <- list.files("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/reactions", full.names = T)
# keggReaction <- readLines("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/reactions/R00062")
reacEnzymeMap <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(reacEnzymeMap) <- c("reaction", "enzyme")

for(keggEntry in files){
  keggReaction <- readLines(keggEntry)
  # print(keggReaction)
  reaction <- stringr::str_extract(keggReaction[grep("ENTRY", keggReaction)], "R[0-9]{5}")
  pre_enzyme <- stringr::str_extract(keggReaction[grep("ENZYME", keggReaction)], "[0-9\\.]+")
  enzyme <- ifelse(identical(pre_enzyme, character(0)), NA, pre_enzyme)
  reacEnzymeMap <- rbind(reacEnzymeMap, 
                         data.frame("reaction" = reaction, 
                                    "enzyme" = enzyme))
}

rUtils::write.tsv(reacEnzymeMap, 
                  file = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/reacEnzymeMap_20200311.tsv")
save(reacEnzymeMap, 
     file = "/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/reacEnzymeMap_20200311.Rda")


dim(reacEnzymeMap) # 11332     2
dim(reacEnzymeMap[complete.cases(reacEnzymeMap), ]) # 10000     2
View(head(reacEnzymeMap, 4000)) # R00645, R00887 (много таких)
library(dplyr)
final_mapping <- tibble::as.tibble(final_mapping)
reacEnzymeMap <- tibble::as.tibble(reacEnzymeMap)

length(intersect(unique(final_mapping$reaction),
                 unique(reacEnzymeMap$reaction))) # 9 339 of 9 341 in final_mapping
length(unique(reacEnzymeMap$reaction)) # 11 332 (1 993 left)
