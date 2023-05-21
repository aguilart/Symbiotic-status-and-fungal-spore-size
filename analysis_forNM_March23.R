## Review comments PNAS

library(rcartocolor)


To_Analysis$simpleFunct <- 
gsub("Plant Pathogen Mildew", "Plant pathogenic biotroph", To_Analysis$simpleFunct,
     ignore.case = T)

To_Analysis$simpleFunct <- 
  gsub("Plant Pathogen Rust", "Plant pathogenic biotroph", To_Analysis$simpleFunct,
       ignore.case = T)

To_Analysis$simpleFunct <- 
  gsub("Plant Pathogen Smut", "Plant pathogenic biotroph", To_Analysis$simpleFunct,
       ignore.case = T)

To_Analysis$simpleFunct <- 
  gsub("Plant Pathogen undefined", "Plant Pathogen Necrotroph", To_Analysis$simpleFunct,
       ignore.case = T)

Li_tree_names$simpleFunct <- 
  gsub("Plant Pathogen Mildew", "Plant pathogenic biotroph", Li_tree_names$simpleFunct,
       ignore.case = T)

Li_tree_names$simpleFunct <- 
  gsub("Plant Pathogen Rust", "Plant pathogenic biotroph", Li_tree_names$simpleFunct,
       ignore.case = T)

Li_tree_names$simpleFunct <- 
  gsub("Plant Pathogen Smut", "Plant pathogenic biotroph", Li_tree_names$simpleFunct,
       ignore.case = T)

Li_tree_names$simpleFunct <- 
  gsub("Plant Pathogen undefined", "Plant Pathogen Necrotroph", Li_tree_names$simpleFunct,
       ignore.case = T)

# Functions to use later

limpia <- function(x){
  #names(x)[1] <- "group" 
  x$group[grep("tercept" ,x$group)] <- "01_Asymbiotic saprotrophs"
  x$group[grep("human" ,x$group, ignore.case = T)] <- "02_Human pathogens"
  x$group[grep("endoph" ,x$group, ignore.case = T)] <- "03_Plant endophytes"
  x$group[grep("necro" ,x$group, ignore.case = T)] <- "04_Plant pathogenic necrotrophs"
  #x$group[grep("undefined" ,x$group, ignore.case = T)] <- "05_Undefined plant pathogens"
  x$group[grep("Lichen" ,x$group, ignore.case = T)] <- "05_Lichen fungi"
  x$group[grep("insect" ,x$group, ignore.case = T)] <- "06_Insect pathogens"
  x$group[grep("Ecto" ,x$group, ignore.case = T)] <- "07_Ectomycorrhizal"
  x$group[grep("biotroph" ,x$group, ignore.case = T)] <- "08_Plant pathogenic biotrophs"
  # x$group[grep("rust" ,x$group, ignore.case = T)] <- "10_Rust fungi(plant pathogens)"
  # x$group[grep("smut" ,x$group, ignore.case = T)] <- "11_Smut fungi (plant pathogens)"
  
  
  # x$group[grep("Life_style1" ,x$group, ignore.case = T)] <- "Asymbionts vs Symbionts"
  # x$group[grep("Life_style2" ,x$group, ignore.case = T)] <- "Facultative symbionts vs Obligate symbionts"
  
  x
}

making_tables <- function(modelo){
  
  muestra <- modelo$n
  modelo <- summary(modelo)
  
  step1 <- 
    data.frame(
      model = paste(modelo$call$data, modelo$call$phy, sep = "_"),
      Pagels_lambda = modelo$optpar,
      r2 = modelo$r.squared,
      n = muestra)
  
  step1 <- rapply(step1, classes = "numeric",
                  how = "replace", function(x){round(x, 2)})
  
  step2 <- as.data.frame(modelo$coefficients, optional = T)
  step2$group <- rownames(step2); step2 <- step2[,c(5, 1:4) ]
  step2$StdErr <- NULL
  rownames(step2) <- NULL
  step2 <- rapply(step2, classes = "numeric",
                  how = "replace", function(x){round(x, 2)})
  
  step2 <- limpia(step2)
  step2 <- step2[order(step2[,1]),]
  step2[,1] <- gsub("\\d\\d_", "",step2[,1])
  
  tabla <-
    rbind(matrix(names(step1), nrow = 1, ncol = 4), as.matrix(step1),
          matrix(names(step2), nrow = 1, ncol = 4), as.matrix(step2))
  
  tabla[which(tabla[ ,4]== "0.00"|tabla[ ,4]== "0.0"|tabla[ ,4]== "0"), 4] <- "<0.01"
  tabla
  
}

cleaning_functions <- function(x){
  
  x$simpleFunct[grep("Saprotroph|Intercept" ,x$simpleFunct)] <- "01_Asymbiotic\nsaprotrophs"
  x$simpleFunct[grep("human" ,x$simpleFunct, ignore.case = T)] <- "02_Human pathogens"
  x$simpleFunct[grep("Animal" ,x$simpleFunct, ignore.case = T)] <- "02.5_Animal endoparasite"
  x$simpleFunct[grep("Fungi" ,x$simpleFunct, ignore.case = T)] <- "03_Fungal Parasite"
  x$simpleFunct[grep("Lichen" ,x$simpleFunct, ignore.case = T)] <- "04_Lichen fungi"
  x$simpleFunct[grep("insect" ,x$simpleFunct, ignore.case = T)] <- "05_Insect pathogens"
  x$simpleFunct[grep("endoph" ,x$simpleFunct, ignore.case = T)] <- "06_Plant endophyte"
  #x$simpleFunct[grep("undefined" ,x$simpleFunct, ignore.case = T)] <- "07_Undefined disease\n(plant pathogens)"
  x$simpleFunct[grep("necro" ,x$simpleFunct, ignore.case = T)] <- "07_Necrotroph\n(plant pathogens)"
  x$simpleFunct[grep("Ecto" ,x$simpleFunct, ignore.case = T)] <- "08_Ectomycorrhizal"
  x$simpleFunct[grep("biotroph" ,x$simpleFunct, ignore.case = T)] <- "09_Biotroph\n(plant pathogens)"
  #x$simpleFunct[grep("rust" ,x$simpleFunct, ignore.case = T)] <- "11_Rust fungi\n(plant pathogens)"
  #x$simpleFunct[grep("smut" ,x$simpleFunct, ignore.case = T)] <- "12_Smut fungi\n(plant pathogens)"
  
  x
}

get_predicted <- function(modelo, phy, mod) {
  
  modelo <- summary(modelo)
  step2 <- as.data.frame(modelo$coefficients)
  step2$simpleFunct <- rownames(step2); rownames(step2) <- NULL
  
  step2$Estimate[2:length(step2$Estimate)] <- 
    step2$Estimate[2:length(step2$Estimate)] + step2$Estimate[1]
  
  step2$upper <- step2$Estimate + step2$StdErr
  step2$lower <- step2$Estimate - step2$StdErr
  
  step2 <- rapply(step2, classes = "numeric",
                  how = "replace", function(x){round(x, 2)})
  step2 <- cleaning_functions(step2)
  step2$phylum <- phy
  step2$model <- mod
  
  step2
}

# They requested to split the analysis into phyla. I am using the genaral groupings from the
# Mycota


######  Ascospores  #####

# Li tree dataset

as <- which(Li_tree_names$SporeType=="Meiospores"&
              Li_tree_names$phylum=="Ascomycota")

ascospores <- Li_tree_names[as,]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(ascospores$simpleFunct)<3))
ascospores <- ascospores[-which(ascospores$simpleFunct %in% few), ]
ascospores <- ascospores[-which(ascospores$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
ascospores <- ascospores[which(ascospores$Life_style %in% styles), ]

rownames(ascospores)<-ascospores$orginal_tree_names

mod_Li_ascos <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = ascospores,
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_ascos)

# tax tree dataset

ascospores2 <- To_Analysis[which(grepl("eiospores",To_Analysis$SporeType)&
                                 To_Analysis$phylum == "Ascomycota"),]

ascospores2 <- 
  ascospores2[which(ascospores2$names_to_use%in%tax_tree2$tip.label), ]

ascospores2 <- ascospores2[-which(ascospores2$simpleFunct == "Unknown"), ]
ascospores2 <- ascospores2[which(ascospores2$Life_style %in% styles), ]
few <- names(which(table(ascospores2$simpleFunct)<5))
ascospores2 <- ascospores2[-which(ascospores2$simpleFunct %in% few), ]
rownames(ascospores2) <- ascospores2$names_to_use

mod_tax_ascos <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = ascospores2,
          
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_tax_ascos)

##### Basidiospores ######

bas <- which(Li_tree_names$SporeType=="Meiospores"&
               Li_tree_names$phylum=="Basidiomycota")

basidiospores <- Li_tree_names[bas,]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(basidiospores$simpleFunct)<4))
basidiospores <- basidiospores[-which(basidiospores$simpleFunct %in% few), ]
#basidiospores <- basidiospores[-which(basidiospores$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
basidiospores <- basidiospores[which(basidiospores$Life_style %in% styles), ]
rownames(basidiospores)<-basidiospores$orginal_tree_names

mod_Li_basidios <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = basidiospores,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_basidios)

#### basidiospores again but using the tax tree

basidiospores2<-To_Analysis[which(grepl("eiospores",To_Analysis$SporeType)&
                                    To_Analysis$phylum == "Basidiomycota"),]
basidiospores2 <- 
  basidiospores2[which(basidiospores2$names_to_use%in%tax_tree2$tip.label), ]
basidiospores2 <- basidiospores2[which(basidiospores2$Life_style %in% styles), ]
few <- names(which(table(basidiospores2$simpleFunct)<5))
basidiospores2 <- basidiospores2[-which(basidiospores2$simpleFunct %in% few), ]
#basidiospores2 <- basidiospores2[-which(basidiospores2$simpleFunct == "Unknown"), ]

rownames(basidiospores2) <- basidiospores2$names_to_use

mod_tax_basidios <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = basidiospores2,
          
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_tax_basidios)

###################  MITOSPORES  ##############################

mas <- which(Li_tree_names$SporeType=="Mitospores"&
               Li_tree_names$phylum=="Ascomycota")
#mitospores_as
mitospores_as <- Li_tree_names[mas,]
mitospores_as <- mitospores_as[-which(mitospores_as$simpleFunct == "Unknown"), ]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(mitospores_as$simpleFunct)<3))
mitospores_as <- mitospores_as[-which(mitospores_as$simpleFunct %in% few), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
mitospores_as <- mitospores_as[which(mitospores_as$Life_style %in% styles), ]
rownames(mitospores_as)<-mitospores_as$orginal_tree_names


mod_Li_mitos_as <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_as,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_mitos_as)

#### mitospores_as again but using the tax tree

mitospores_as2<-To_Analysis[which(grepl("itospores",To_Analysis$SporeType)&
                                    To_Analysis$phylum == "Ascomycota"),]
mitospores_as2 <- mitospores_as2[-which(mitospores_as2$simpleFunct == "Unknown"), ]
mitospores_as2 <- mitospores_as2[which(mitospores_as2$Life_style %in% styles), ]
mitospores_as2 <- 
  mitospores_as2[which(mitospores_as2$names_to_use%in%tax_tree2$tip.label), ]

few <- names(which(table(mitospores_as2$simpleFunct)<5))
mitospores_as2 <- mitospores_as2[-which(mitospores_as2$simpleFunct %in% few), ]
rownames(mitospores_as2) <- mitospores_as2$names_to_use

mod_tax_mitos_as <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_as2,
          
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_tax_mitos_as)



# Now forw mitospores_bas

bas <- which(Li_tree_names$SporeType=="Mitospores"&
               Li_tree_names$phylum=="Basidiomycota")

mitospores_bas <- Li_tree_names[bas,]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(mitospores_bas$simpleFunct)<3))
mitospores_bas <- mitospores_bas[-which(mitospores_bas$simpleFunct %in% few), ]
#mitospores_bas <- mitospores_bas[-which(mitospores_bas$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
mitospores_bas <- mitospores_bas[which(mitospores_bas$Life_style %in% styles), ]
rownames(mitospores_bas)<-mitospores_bas$orginal_tree_names


mod_Li_mitos_bas <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_bas,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_mitos_bas)

#### mitospores_bas again but using the tax tree


mitospores_bas2<-To_Analysis[which(grepl("itospores",To_Analysis$SporeType)&
                                     To_Analysis$phylum == "Basidiomycota"),]
mitospores_bas2 <- 
  mitospores_bas2[which(mitospores_bas2$names_to_use%in%tax_tree2$tip.label), ]

mitospores_bas2 <- mitospores_bas2[-which(mitospores_bas2$simpleFunct == "Unknown"), ]
mitospores_bas2 <- mitospores_bas2[which(mitospores_bas2$Life_style %in% styles), ]
few <- names(which(table(mitospores_bas2$simpleFunct)<5))
mitospores_bas2 <- mitospores_bas2[-which(mitospores_bas2$simpleFunct %in% few), ]
rownames(mitospores_bas2) <- mitospores_bas2$names_to_use

#To check: 4 taxa not in the tree
mod_tax_mitos_bas <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_bas2,
          
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_tax_mitos_bas)

### For Zygomycetous

# Li tree

zy <- which(Li_tree_names$SporeType=="Mitospores"&
              grepl("Mucoromycota|Zoopagomycota|Zygomycetous", Li_tree_names$phylum))

mitospores_zy <- Li_tree_names[zy,]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(mitospores_zy$simpleFunct)<3))
mitospores_zy <- mitospores_zy[-which(mitospores_zy$simpleFunct %in% few), ]
#mitospores_zy <- mitospores_zy[-which(mitospores_zy$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
mitospores_zy <- mitospores_zy[which(mitospores_zy$Life_style %in% styles), ]
rownames(mitospores_zy)<-mitospores_zy$orginal_tree_names

mod_Li_mitos_zy <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_zy,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_mitos_zy)

# tax tree
# To check: there are 10 taxa that are not in the tree... this is weird....
# Actually this makes sense, the tree is composed of species with complete taxonomy

mitospores_zy2 <- To_Analysis[which(To_Analysis$SporeType=="Mitospores"&
                                      grepl("Mucoromycota|Zoopagomycota|Zygomycetous",
                                            To_Analysis$phylum)),]
mitospores_zy2 <- 
  mitospores_zy2[which(mitospores_zy2$names_to_use%in%tax_tree2$tip.label), ]

mitospores_zy2 <- mitospores_zy2[which(mitospores_zy2$Life_style %in% styles), ]
few <- names(which(table(mitospores_zy2$simpleFunct)<5))
mitospores_zy2 <- mitospores_zy2[-which(mitospores_zy2$simpleFunct %in% few), ]
#mitospores_zy2 <- mitospores_zy2[-which(mitospores_zy2$simpleFunct == "Unknown"), ]

rownames(mitospores_zy2) <- mitospores_zy2$names_to_use

mod_tax_mitos_zy2 <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_zy2,
          
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_tax_mitos_zy2)

# The "lower" fungi: For zoosporic fungi and microporidan fungi it makes no 
# sense to use the taxonomy tree as almost all the entries are in Li tree so
# I can just use the genome level tree
##### Zoosporic fungi

# Li tree
zoo <- which(Li_tree_names$SporeType=="Mitospores"&
               grepl("Blastocladiomycota|Neocallimastigomycota|Chytridiomycota",
                     Li_tree_names$phylum))

mitospores_zoo <- Li_tree_names[zoo,]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(mitospores_zoo$simpleFunct)<3))
mitospores_zoo <- mitospores_zoo[-which(mitospores_zoo$simpleFunct %in% few), ]
#mitospores_zoo <- mitospores_zoo[-which(mitospores_zoo$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
mitospores_zoo <- mitospores_zoo[which(mitospores_zoo$Life_style %in% styles), ]
rownames(mitospores_zoo)<-mitospores_zoo$orginal_tree_names

mod_Li_mitos_zoo <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_zoo,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_mitos_zoo)

#checking things:
# zoo_in_li <- gsub("_", " ", mitospores_zoo$orginal_tree_names)
# zoo_in_tax <- mitospores_zoo2$names_to_use
# 
# zoo_in_li[which(!zoo_in_li%in%zoo_in_tax)]# three species
# 
# zoo_in_tax[which(!zoo_in_tax%in%zoo_in_li)] # three species too and
# # the functions of those three are precisely the functions not present in the li dataset
# mitospores_zoo2$simpleFunct[which(mitospores_zoo2$names_to_use%in%zoo_in_tax[which(!zoo_in_tax%in%zoo_in_li)])]
# 
# # the functions reported in the li tree are only saprotrophs (12) and Animal (3)
# # the functions reported in the tax tree are saprotrophs (11), Algae, Animal, Fungi and plant, each with 1
# 
# mitospores_zoo2$simpleFunct[which(mitospores_zoo2$names_to_use%in%zoo_in_tax[which(!zoo_in_tax%in%zoo_in_li)])]

# Tax tree: Update this does not work because there is not enough species with functional
# groups. There are 12 saprotrophs and there are around 3-4 animal pathogens, however, out
# these ones only one has a fully resolved taxonomy (no missing entries or "incerti sedis")
# leaving nothing to compare really.

# mitospores_zoo2 <- To_Analysis[which(To_Analysis$SporeType=="Mitospores"&
#                                        grepl("Blastocladiomycota|Neocallimastigomycota|Chytridiomycota",
#                                              To_Analysis$phylum)),]
# mitospores_zoo2 <- 
#   mitospores_zoo2[which(mitospores_zoo2$names_to_use%in%tax_tree2$tip.label), ]

# mitospores_zoo2 <- mitospores_zoo2[which(mitospores_zoo2$Life_style %in% styles), ]
# few <- names(which(table(mitospores_zoo2$simpleFunct)<5))
# mitospores_zoo2 <- mitospores_zoo2[-which(mitospores_zoo2$simpleFunct %in% few), ]
# #mitospores_zoo2 <- mitospores_zoo2[-which(mitospores_zoo2$simpleFunct == "Unknown"), ]
# rownames(mitospores_zoo2) <- mitospores_zoo2$names_to_use
# 
# mod_tax_mitos_zoo2 <-
#   phylolm(log10(SporeVolume) ~ simpleFunct,
#           data = mitospores_zoo2,
#           
#           phy = tax_tree2,
#           #phy = Li_tree2,
#           model = "lambda")
# 
# summary(mod_tax_mitos_zoo2)

##### Microsporidia #########

mic <- which(Li_tree_names$SporeType=="Mitospores"&
               grepl("Microsporidia",
                     Li_tree_names$phylum))

mitospores_mic <- Li_tree_names[mic,]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(mitospores_mic$simpleFunct)<3))
#mitospores_mic <- mitospores_mic[-which(mitospores_mic$simpleFunct %in% few), ]
#mitospores_mic <- mitospores_mic[-which(mitospores_mic$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
mitospores_mic <- mitospores_mic[which(mitospores_mic$Life_style %in% styles), ]
rownames(mitospores_mic)<-mitospores_mic$orginal_tree_names

mod_Li_mitos_mic <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = mitospores_mic,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_mitos_mic)

# # Tax tree
# 
# mitospores_mic2 <- To_Analysis[which(To_Analysis$SporeType=="Mitospores"&
#                                        grepl("Microsporidia",
#                                              To_Analysis$phylum)),]
# 
# mitospores_mic2 <- 
#   mitospores_mic2[which(mitospores_mic2$names_to_use%in%tax_tree2$tip.label), ]
# 
# mitospores_mic2 <- mitospores_mic2[which(mitospores_mic2$Life_style %in% styles), ]
# few <- names(which(table(mitospores_mic2$simpleFunct)<3))
# #mitospores_mic2 <- mitospores_mic2[-which(mitospores_mic2$simpleFunct %in% few), ]
# #mitospores_mic2 <- mitospores_mic2[-which(mitospores_mic2$simpleFunct == "Unknown"), ]
# rownames(mitospores_mic2) <- mitospores_mic2$names_to_use
# 
# mod_tax_mitos_mic2 <-
#   phylolm(log10(SporeVolume) ~ simpleFunct,
#           data = mitospores_mic2,
#           
#           phy = tax_tree2,
#           #phy = Li_tree2,
#           model = "lambda")
# 
# summary(mod_tax_mitos_mic2)

# It is interesnting that when using the two trees for analyzign the same data, the 
# tax tree actually is less conservative than the phylo tree... something to discuss
# with others

#### Multinucleate spores of the zygomycetous fungi

mu <- which(Li_tree_names$SporeType=="Multinucleate sexual spores"&
              !grepl("Basidiomycota",
                     Li_tree_names$phylum))

multi_sex_mu <- Li_tree_names[mu, ]
#Removing the groups for which we have less than 3 species per function and "unknown" functions
few <- names(which(table(multi_sex_mu$simpleFunct)<3))
multi_sex_mu <- multi_sex_mu[-which(multi_sex_mu$simpleFunct %in% few), ]
#multi_sex_mu <- multi_sex_mu[-which(multi_sex_mu$simpleFunct == "Unknown"), ]
styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association") 
multi_sex_mu <- multi_sex_mu[which(multi_sex_mu$Life_style %in% styles), ]
rownames(multi_sex_mu)<-multi_sex_mu$orginal_tree_names

mod_Li_multi_sex <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = multi_sex_mu,
          
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_Li_multi_sex)

# Tax tree

multi_sex_mu2 <- To_Analysis[which(To_Analysis$SporeType=="Multinucleate sexual spores"&
                                     !grepl("Basidiomycota",
                                            To_Analysis$phylum)),]

multi_sex_mu2 <- 
  multi_sex_mu2[which(multi_sex_mu2$names_to_use%in%tax_tree2$tip.label), ]

multi_sex_mu2 <- multi_sex_mu2[which(multi_sex_mu2$Life_style %in% styles), ]
few <- names(which(table(multi_sex_mu2$simpleFunct)<5))
multi_sex_mu2 <- multi_sex_mu2[-which(multi_sex_mu2$simpleFunct %in% few), ]
#multi_sex_mu2 <- multi_sex_mu2[-which(multi_sex_mu2$simpleFunct == "Unknown"), ]
rownames(multi_sex_mu2) <- multi_sex_mu2$names_to_use

mod_tax_multi_sex <-
  phylolm(log10(SporeVolume) ~ simpleFunct,
          data = multi_sex_mu2,
          
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

summary(mod_tax_multi_sex)

##### saving the results



# Saving the data for meiospores groups

ascos_basidios <- 
  rbind(
    making_tables(mod_Li_ascos),
    making_tables(mod_tax_ascos),
    making_tables(mod_Li_basidios),
    making_tables(mod_tax_basidios)
  )

ascos_basidios[, 1] <- gsub("Li_tree", "genome_tree", ascos_basidios[, 1])
ascos_basidios[, 1] <- gsub("tax_tree2", "taxonomy_tree", ascos_basidios[, 1])

write.csv(ascos_basidios,
          "model_outputs/simple_funct_models_ascos_basidios.csv",
          row.names = F)

# Saving data for mitospores groups

mitos_as_bas_zy <- 
  rbind(
    making_tables(mod_Li_mitos_as),
    making_tables(mod_tax_mitos_as),
    making_tables(mod_Li_mitos_bas),
    making_tables(mod_tax_mitos_bas),
    making_tables(mod_Li_mitos_zy),
    making_tables(mod_tax_mitos_zy2),
    making_tables(mod_Li_mitos_zoo),
    #making_tables(mod_tax_mitos_zoo2),
    making_tables(mod_Li_mitos_mic)#,
    #making_tables(mod_tax_mitos_mic2)
    )

mitos_as_bas_zy[, 1] <- gsub("Li_tree", "genome_tree", mitos_as_bas_zy[, 1])
mitos_as_bas_zy[, 1] <- gsub("tax_tree2", "taxonomy_tree", mitos_as_bas_zy[, 1])

write.csv(mitos_as_bas_zy,
          "model_outputs/simple_funct_models_mitos_all.csv",
          row.names = F)

# write.csv(mitos_as_bas_zy,
#           "model_outputs/simple_funct_models_mitos_all_no_chlam.csv",
#           row.names = F)

# ###############################################################################
# 

# 
# 
To_Analysis_pa <- To_Analysis[grep("Sapro|Plant", To_Analysis$simpleFunct), ]
Li_tree_names_pa <- Li_tree_names[grep("Sapro|Plant", Li_tree_names$simpleFunct), ]


# meiospores
Li_tree_names_pa_me <-
  Li_tree_names_pa[which(Li_tree_names_pa$SporeType == "Meiospores"), ]
Li_tree_names_pa_me <-
  Li_tree_names_pa_me[-which(Li_tree_names_pa_me$phylum == "Blastocladiomycota"), ]#removing the only species from this phylum


#specifying that only life styles in the hypothesis are including

styles <- c("AFree living", "B_Opportunistic association",
            "Facultative association", "Obligate association")

Li_tree_names_pa_me <-
  Li_tree_names_pa_me[which(Li_tree_names_pa_me$Life_style %in% styles), ]
rownames(Li_tree_names_pa_me) <- Li_tree_names_pa_me$orginal_tree_names

# mitospores
Li_tree_names_pa_mi <-
  Li_tree_names_pa[which(Li_tree_names_pa$SporeType == "Mitospores"), ]

#Removing the groups for which we have less than 3 species per function and "unknown" functions
Li_tree_names_pa_mi <-
  Li_tree_names_pa_mi[which(Li_tree_names_pa_mi$Life_style %in% styles), ]
rownames(Li_tree_names_pa_mi)<-Li_tree_names_pa_mi$orginal_tree_names

# For the taxonomy tree data subset

# meiospores
To_Analysis_pa_me <- To_Analysis_pa[grep("eiospores",To_Analysis_pa$SporeType),]
To_Analysis_pa_me <- To_Analysis_pa_me[which(To_Analysis_pa_me$Life_style %in% styles), ]
To_Analysis_pa_me <- 
  To_Analysis_pa_me[-which(To_Analysis_pa_me$phylum=="Blastocladiomycota"), ]
To_Analysis_pa_me <- 
  To_Analysis_pa_me[which(To_Analysis_pa_me$names_to_use%in%tax_tree2$tip.label), ]

rownames(To_Analysis_pa_me) <- To_Analysis_pa_me$names_to_use

# mitospores
To_Analysis_pa_mi <-
  To_Analysis_pa[grep("itospores", To_Analysis_pa$SporeType), ]
To_Analysis_pa_mi <-
  To_Analysis_pa_mi[which(To_Analysis_pa_mi$Life_style %in% styles), ]

# Updating the entry of Fusarium anthophilum: On the Spore functions database, this fungus is listed as Opportunistic
# because it has been reported indeed as a human pathogen. However, it also listed as endophyte and plant pathogen.
# as thhis analysis is restricted only to plant pathogens (human pathogens are out) it is necessary to update the
# entry for Life style
To_Analysis_pa_mi$Life_style[which(To_Analysis_pa_mi$names_to_use == "Fusarium anthophilum")] <- "Facultative association"
To_Analysis_pa_mi <- 
  To_Analysis_pa_mi[which(To_Analysis_pa_mi$names_to_use%in%tax_tree2$tip.label), ]
rownames(To_Analysis_pa_mi)<-To_Analysis_pa_mi$names_to_use


#Meiospores: Ascomycota

ascospores_pa <- 
 ascospores[grep("Plant", ascospores$simpleFunct), ]

mod_Li_pa_me_symb_as <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = Li_tree_names_pa_me[which(Li_tree_names_pa_me$phylum == "Ascomycota"), ],
          data = ascospores_pa,
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")


ascospores2_pa <- ascospores2[grep("Plant", ascospores2$simpleFunct), ]

mod_tax_pa_me_symb_as <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = To_Analysis_pa_me[which(To_Analysis_pa_me$phylum == "Ascomycota"), ],
          data = ascospores2_pa,
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

#Meiospores: Basidiomycota

basidiospores_pa <- basidiospores[grep("Plant", basidiospores$simpleFunct),]

mod_Li_pa_me_symb_bas <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = Li_tree_names_pa_me[which(Li_tree_names_pa_me$phylum == "Basidiomycota"), ],
          data = basidiospores_pa,
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")

basidiospores2_pa <- basidiospores2[grep("Plant", basidiospores2$simpleFunct), ]

mod_tax_pa_me_symb_bas <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = To_Analysis_pa_me[which(To_Analysis_pa_me$phylum == "Basidiomycota"), ],
          data = basidiospores2_pa,
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")


#Mitospores: Ascomycota

mitospores_as_pa <- 
 mitospores_as[grep("Plant", mitospores_as$simpleFunct), ]

mod_Li_pa_mi_symb_as <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = Li_tree_names_pa_mi[which(Li_tree_names_pa_mi$phylum == "Ascomycota"), ],
          data = mitospores_as_pa,
          phy = Li_tree,
          #phy = Li_tree2,
          model = "lambda")



mitospores_as2$Life_style[mitospores_as2$names_to_use == "Fusarium anthophilum"] <-
  "Facultative association"
mitospores_as2_pa <- 
mitospores_as2[grep("Plant", mitospores_as2$simpleFunct), ]


mod_tax_pa_mi_symb_as <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = To_Analysis_pa_mi[which(To_Analysis_pa_mi$phylum == "Ascomycota"), ],
          data = mitospores_as2_pa,
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

#Mitospores: Basidiomycota does not work for Li because there are not enough species
# mod_Li_pa_mi_symb_bas <- # Low number of species
#   phylolm(log10(SporeVolume) ~ Life_style,
#           #data = Li_tree_names_pa_mi[which(Li_tree_names_pa_mi$phylum == "Basidiomycota"), ],
#           data = mitospores_bas[grep("Plant", mitospores_bas$simpleFunct), ],
#           phy = Li_tree,
#           #phy = Li_tree2,
#           model = "lambda")

# I think I will also justify not using these ones either because the only
# obligate biotrophs I have here are the smut fungi
mod_tax_pa_mi_symb_bas <-
  phylolm(log10(SporeVolume) ~ Life_style,
          #data = To_Analysis_pa_mi[which(To_Analysis_pa_mi$phylum == "Basidiomycota"), ],
          data = mitospores_bas2[grep("Plant", mitospores_bas2$simpleFunct), ],
          phy = tax_tree2,
          #phy = Li_tree2,
          model = "lambda")

#Mitospores: non-dykaria
# I cannot make the comparison because there is only 1 Obligate association in this group
# for both the genome tree and the taxonomy tree
table(Li_tree_names_pa_mi$Life_style[-grep("Asco|Basidio",
                                           Li_tree_names_pa_mi$phylum)])

table(To_Analysis_pa_mi$Life_style[-grep("Asco|Basidio",
                                         To_Analysis_pa_mi$phylum)])

table(basidiospores$simpleFunct[grep("Plant",
                                         basidiospores$simpleFunct)])


table(mitospores_as2$simpleFunct[grep("Plant",
                                      mitospores_as2$simpleFunct)])

# # Saving the data
# 

making_tables2 <- function(modelo){

tabla <-
cbind(
as.data.frame(summary(modelo)$coefficients)[-1,],
data.frame(
  model = paste(modelo$call$data, modelo$call$phy, sep = "_"),
  Pagels_lambda = modelo$optpar,
  r2 = modelo$r.squared,
  n = modelo$n)
)
rownames(tabla) <- NULL
#tabla$Contrast <- "Oblibate vs Facultative symbiosis"
tabla <- tabla[,c(5,1:4,6:8)]
names(tabla)[2] <- "Diff. Facultative vs Obligate symbiosis"
tabla <- rapply(tabla, classes = "numeric",
                how = "replace", function(x){round(x, 2)})
tabla$p.value[which(tabla$p.value < 0.01)] <- "<0.01"
tabla
}

all_contrasts <- bind_rows(
making_tables2(mod_Li_pa_me_symb_as), #matrix(c("","","",""), ncol = 9),
making_tables2(mod_tax_pa_me_symb_as),# matrix(c("","","",""), ncol = 9),
making_tables2(mod_Li_pa_me_symb_bas), #matrix(c("","","",""), ncol = 9),
making_tables2(mod_tax_pa_me_symb_bas),# matrix(c("","","",""), ncol = 9),
making_tables2(mod_Li_pa_mi_symb_as), #matrix(c("","","",""), ncol = ),
making_tables2(mod_tax_pa_mi_symb_as)#, matrix(c("","","",""), ncol = 4)#,
#making_tables(mod_Li_pa_mi_symb_bas), matrix(c("","","",""), ncol = 4),
#making_tables(mod_tax_pa_mi_symb_bas), matrix(c("","","",""), ncol = 4)
)

all_contrasts$model <- 
  gsub("Li_tree", "(genome tree model)", all_contrasts$model)

all_contrasts$model <- 
  gsub("tax_tree2", "(taxonomy tree model)", all_contrasts$model)

all_contrasts$model <- 
  gsub("_pa_", " ", all_contrasts$model)

all_contrasts$model <- 
  gsub("\\d", "", all_contrasts$model)

write.csv(all_contrasts,
          "model_outputs/contrasts_plant_associated.csv",
          row.names = F)

# write.csv(all_contrasts,
#           "model_outputs/contrasts_plant_associated_nochlam.csv",
#           row.names = F)
# 
