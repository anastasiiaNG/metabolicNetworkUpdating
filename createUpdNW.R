load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
dplyr::glimpse(network)
length(unique(network$reactions$reaction)) # 10 489
head(network$reaction2rpair)
length(unique(network$reaction2rpair$reaction)) # 9 726

head(unique(network$reactions$reaction))
dim(network$rpair2align) # 119 278      3

rabs <- setdiff(unique(network$reaction2rpair$reaction), unique(final_mapping$reaction))
length(rabs) # in network
sample(rabs, 40)
length(setdiff(unique(final_mapping$reaction), unique(network$reactions$reaction))) # in final_mapping
length(intersect(unique(final_mapping$reaction), unique(network$reactions$reaction)))

length(setdiff(unique(network$reactions$reaction), new_lines)) # in network
length(setdiff(new_lines, unique(network$reactions$reaction))) # in new_lines # 896
length(intersect(new_lines, unique(network$reactions$reaction))) # 10 438

make_network.R # from gatom - slide 23 - table "was / now"


load("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/final_mapping_20191212.Rda")
# NA and *
dim(final_mapping) # 331 976      4
length(unique(final_mapping$reaction)) # 331 976      4
check_na <- final_mapping[union(which(is.na(final_mapping[,2])), which(is.na(final_mapping[,4]))),]
dim(check_na)
View(check_na)
length(unique(check_na$reaction))
check_h <- check_na[union(grep(".*H.*", check_na$reactant), grep(".*H.*", check_na$product)),]
dim(check_h)
View(check_h)
length(unique(check_h$reaction))
check_rest <- check_na[-union(grep(".*H.*", check_na$reactant), grep(".*H.*", check_na$product)),]
View(check_rest)




dplyr::glimpse(final_mapping)
length(unique(final_mapping$reaction)) # 7 340
head(unique(final_mapping$reaction))

network_upd <- network[names(network) %in% c("reactions",
                                             "enzyme2reaction",
                                             "rpairs", # ? (reaction and metabolites)
                                             "atoms",
                                             "metabolite2atom")]
dplyr::glimpse(network_upd)
network_upd$reaction2align <- final_mapping[, -3] # instead of "rpair2align" and "reaction2rpair"

load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
dplyr::glimpse(org.Mm.eg.gatom.anno)




# getting reaction attributes:
equations <- readRDS("/home/octopus/Documents/5gatom_draft/gatom_network_updating/data/equations.rds")
dplyr::glimpse(equations)
equations[grep(".*C00194.*", names(equations))]
head(equations)

