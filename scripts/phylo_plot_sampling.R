# phylo_plot_sampling.R

require(ape);packageVersion("ape")#5.5
require(dplyr);packageVersion("dplyr")#1.1.4
require(ggplot2);packageVersion("ggplot2")#3.4.4
require(ggtree);packageVersion("ggtree")#3.0.4
require(cowplot);packageVersion("cowplot")#1.1.1
require(caper); packageVersion("caper")#1.0.1

# host phylogeny
tree <- read.nexus("../raw_data/upham_tree_666.nex")
tax <- read.csv("../raw_data/taxonomy_mamPhy_5911species.csv")

hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)

#subset to only host reported and upham
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()


# order-level trees
bats <- tax[tax$ord=="CHIROPTERA",]
setdiff(bats$Species_Name, tree$tip.label)# all there

rodents <- tax[tax$ord=="RODENTIA",]
setdiff(rodents$Species_Name, tree$tip.label)# all there

bat_tree <- drop.tip(tree, setdiff(tree$tip.label, bats$Species_Name))
# 1287 species in phylogeny

rodent_tree <- drop.tip(tree, setdiff(tree$tip.label, rodents$Species_Name))
# 2392 species in phylogeny


# virus taxonomy
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
vtax <- as.data.frame(unclass(vtax), stringsAsFactors=TRUE)

viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- dplyr::left_join(viruses, vtax)

# make tree
frm <- ~superkingdom/realm/kingdom/phylum/class/order/family/genus/Virus_ICTV
vtree <- as.phylo(frm, data = vtax, collapse=FALSE)
vtree$edge.length <- rep(1, nrow(vtree$edge))


# experimental data
primary <- read.csv("../raw_data/primary_screen.csv")
# new column to keep count of primary papers / hosts
primary$primary_paper <- 1

speciesonly <- read.csv("../raw_data/species_only_data.csv")
indiv <- read.csv("../raw_data/individual_data.csv")

# add names
primary <- dplyr::left_join(primary, viruses)
primary <- dplyr::left_join(primary, hosts)
speciesonly <- dplyr::left_join(speciesonly, viruses)
speciesonly <- dplyr::left_join(speciesonly, hosts)
indiv <- dplyr::left_join(indiv, viruses)
indiv <- dplyr::left_join(indiv, hosts)

# merge speciesonly with individual data
dat <- full_join(indiv, speciesonly)
dat <- full_join(dat, primary)

# removing hosts and viruses with NA names
dat <- dat[!is.na(dat$Virus_ICTV),]
dat <- dat[!is.na(dat$Host_Upham),]


# susceptible hosts and viruses with pathology reported
path_hosts <- unique(dat$Host_Upham[dat$Susceptible_primary_YN=="Y" & dat$PathologyReported_YN=="Y"])
path_viruses <- unique(dat$Virus_ICTV[dat$Susceptible_primary_YN=="Y" & dat$PathologyReported_YN=="Y"])
length(unique(path_hosts[!is.na(path_hosts)]))# 93
length(unique(path_viruses[!is.na(path_viruses)]))# 56

# susceptible hosts and viruses with pathology reported AND individual level data
indiv_hosts <- unique(dat$Host_Upham[dat$Susceptible_primary_YN=="Y" & dat$PathologyReported_YN=="Y" & dat$IndividualData_YN=="Y"])
indiv_viruses <- unique(dat$Virus_ICTV[dat$Susceptible_primary_YN=="Y" & dat$PathologyReported_YN=="Y" & dat$IndividualData_YN=="Y"])
length(unique(indiv_hosts[!is.na(indiv_hosts)]))# 84
length(unique(indiv_viruses[!is.na(indiv_viruses)]))# 47

# plot only species with susceptibility and pathology reported
df1 <- rbind(bats, rodents)
df1$sampled <- df1$Species_Name%in%path_hosts*1
rownames(df1) <- df1$Species_Name



# Identifying family level nodes for labelling
# Bats
unique(df1$fam[df1$sampled==1 & df1$ord=="CHIROPTERA"])
unique(df1$fam[df1$Species_Name %in% path_hosts & df1$ord=="CHIROPTERA"])
n_bats <- length(unique(df1$Species_Name[df1$sampled==1 & df1$ord=="CHIROPTERA"]))

bat_families <- tax$fam[tax$ord=="CHIROPTERA"] %>% unique()
bat_families

# which are most speciose?
bat_sp_by_fam <- tax[tax$ord=="CHIROPTERA",] %>% select(Species_Name, fam) %>% 
				unique() %>% group_by(fam) %>% summarise(n_spec=n_distinct(Species_Name))

arrange(bat_sp_by_fam, desc(n_spec))
# let's plot sampled families, plus unsampled bat families with >10 species

# represented
vesps <- tax$Species_Name[tax$fam=="VESPERTILIONIDAE"]
phyllos <- tax$Species_Name[tax$fam=="PHYLLOSTOMIDAE"]
pteros <- tax$Species_Name[tax$fam=="PTEROPODIDAE"]
molossids <- tax$Species_Name[tax$fam=="MOLOSSIDAE"]

#unrepresented
hippos <- tax$Species_Name[tax$fam=="HIPPOSIDERIDAE"]
rhinolophids <- tax$Species_Name[tax$fam=="RHINOLOPHIDAE"]
# one rhinolphid was not reported to species level and other did not have pathology data reported
emballos <- tax$Species_Name[tax$fam=="EMBALLONURIDAE"]
nycterids <- tax$Species_Name[tax$fam=="NYCTERIDAE"]
natalids <- tax$Species_Name[tax$fam=="NATALIDAE"]

# rhinos_MRCA <- MRCA(.data=bat_tree, .node1=rhinolophids)
vesp_MRCA <- MRCA(.data=bat_tree, .node1=vesps)
phyllos_MRCA <- MRCA(.data=bat_tree, .node1=phyllos)
pteros_MRCA <- MRCA(.data=bat_tree, .node1=pteros)
molossids_MRCA <- MRCA(.data=bat_tree, .node1=molossids)
hippos_MRCA <- MRCA(.data=bat_tree, .node1=hippos)
rhinolophids_MRCA <- MRCA(.data=bat_tree, .node1=rhinolophids)
emballos_MRCA <- MRCA(.data=bat_tree, .node1=emballos)
nycterids_MRCA <- MRCA(.data=bat_tree, .node1=nycterids)
natalids_MRCA <- MRCA(.data=bat_tree, .node1=natalids)


b <- ggtree(bat_tree, layout='circular', open.angle=35, color="darkgrey")

sum(df1$sampled)#93 species
df1 %>% filter(ord=="CHIROPTERA") %>% summarise(sum(sampled))# 34
df1 %>% filter(ord=="RODENTIA") %>% summarise(sum(sampled))# 59

df1$sampled <- as.character(df1$sampled)
df2 <- select(df1, c(sampled))


colours_BR <- c("#1B85BF", "#AB1808")
colours_YB <-c("#FFC20A","#0C7BDC")

heatmap.colours <- c("white",colours_BR[1])
# heatmap.colours <- c("white","darkred")
names(heatmap.colours) <- 0:1


p_bats <- gheatmap(b, df2, offset = -3, color=NULL, width=0.12, 
				colnames=FALSE) +
  		 		scale_fill_manual(values=heatmap.colours, breaks=0:1, name="Sampled") + 
  		 		theme(legend.position="none") + 
  		 		geom_cladelabel(node=vesp_MRCA, label="Vespertilionidae", color=colours_BR[1], offset=10, offset.text=2) +
  		 		geom_hilight(node=vesp_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=phyllos_MRCA, label="Phyllostomidae", color=colours_BR[1], offset=10, offset.text=17, hjust=0.45) +
  		 		geom_hilight(node=phyllos_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=molossids_MRCA, label="Molossidae", color=colours_BR[1], offset=10, offset.text=16, hjust=0.35) +
  		 		geom_hilight(node=molossids_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=pteros_MRCA, label="Pteropodidae", color=colours_BR[1], offset=10, offset.text=2) + 
  		 		geom_hilight(node=pteros_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
  		 		# unsampled
  		 		geom_cladelabel(node=rhinolophids_MRCA, label="Rhinolophidae", color="gray", offset=10, offset.text=1.2, vjust=-1) + 
  		 		geom_cladelabel(node=hippos_MRCA, label="Hipposideridae", color="gray", offset=10, offset.text=1.5) +
  		 		geom_cladelabel(node=nycterids_MRCA, label="Nycteridae", color="gray", offset=10, offset.text=20.5) +
  		 		geom_cladelabel(node=natalids_MRCA, label="Natalidae", color="gray", offset=10, offset.text=20.5) +
  		 		geom_cladelabel(node=emballos_MRCA, label="Emballonuridae", color="gray", offset=10, offset.text=28.5) 

# p_bats 


# Rodents
unique(df1$fam[df1$sampled==1 & df1$ord=="RODENTIA"])
unique(df1$fam[df1$Species_Name %in% path_hosts & df1$ord=="RODENTIA"])
n_rodents <- length(unique(df1$Species_Name[df1$sampled==1 & df1$ord=="RODENTIA"]))

rodent_families <- tax$fam[tax$ord=="RODENTIA"] %>% unique()
rodent_families

# which are most speciose?
rodent_sp_by_fam <- tax[tax$ord=="RODENTIA",] %>% select(Species_Name, fam) %>% 
				unique() %>% group_by(fam) %>% summarise(n_spec=n_distinct(Species_Name))

print(arrange(rodent_sp_by_fam, desc(n_spec)),n=20)
# let's plot represented families, plus unrepresented rodent families with >20 species

# represented
murids <- tax$Species_Name[tax$fam=="MURIDAE"]
sciurids <- tax$Species_Name[tax$fam=="SCIURIDAE"]
cricetids <- tax$Species_Name[tax$fam=="CRICETIDAE"]
castors <- tax$Species_Name[tax$fam=="CASTORIDAE"]
nesomyds <- tax$Species_Name[tax$fam=="NESOMYIDAE"]
dasys <- tax$Species_Name[tax$fam=="DASYPROCTIDAE"]
heteromys <- tax$Species_Name[tax$fam=="HETEROMYIDAE"]
erethiz <- tax$Species_Name[tax$fam=="ERETHIZONTIDAE"]
glirids <- tax$Species_Name[tax$fam=="GLIRIDAE"]
hystrics <- tax$Species_Name[tax$fam=="HYSTRICIDAE"]
echimys <- tax$Species_Name[tax$fam=="ECHIMYIDAE"]

# not represented with >60 species
ctenomys <- tax$Species_Name[tax$fam=="CTENOMYIDAE"]
dipods <- tax$Species_Name[tax$fam=="DIPODIDAE"]
geomys <- tax$Species_Name[tax$fam=="GEOMYIDAE"]
# geomys looks to be within heteromys? - don't plot?
spalacs <- tax$Species_Name[tax$fam=="SPALACIDAE"]
cavys <- tax$Species_Name[tax$fam=="CAVIIDAE"]


murids_MRCA <- MRCA(.data=rodent_tree, .node1=murids)
sciurids_MRCA <- MRCA(.data=rodent_tree, .node1=sciurids)
cricetids_MRCA <- MRCA(.data=rodent_tree, .node1=cricetids)
castors_MRCA <- MRCA(.data=rodent_tree, .node1=castors)
nesomyds_MRCA <- MRCA(.data=rodent_tree, .node1=nesomyds)
dasys_MRCA <- MRCA(.data=rodent_tree, .node1=dasys)
heteromys_MRCA <- MRCA(.data=rodent_tree, .node1=heteromys)
erethiz_MRCA <- MRCA(.data=rodent_tree, .node1=erethiz)
glirids_MRCA <- MRCA(.data=rodent_tree, .node1=glirids)
hystrics_MRCA <- MRCA(.data=rodent_tree, .node1=hystrics)
echimys_MRCA <- MRCA(.data=rodent_tree, .node1=echimys)

ctenomys_MRCA <- MRCA(.data=rodent_tree, .node1=ctenomys)
dipods_MRCA <- MRCA(.data=rodent_tree, .node1=dipods)
geomys_MRCA <- MRCA(.data=rodent_tree, .node1=geomys)
spalacs_MRCA <- MRCA(.data=rodent_tree, .node1=spalacs)
cavys_MRCA <- MRCA(.data=rodent_tree, .node1=cavys)


r <- ggtree(rodent_tree, layout='circular', open.angle=35, color="darkgrey")

heatmap.colours <- c("white",colours_BR[2])
names(heatmap.colours) <- 0:1

p_rodents <- gheatmap(r, df2, offset = -3, color=NULL, width=0.12,
				colnames=FALSE) +
  		 		scale_fill_manual(values=heatmap.colours, breaks=0:1, name="Sampled") + 
  		 		theme(legend.position="none") + 
  		 		geom_cladelabel(node=murids_MRCA, label="Muridae", color=colours_BR[2], offset=10, offset.text=3) +
  		 		geom_hilight(node=murids_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=sciurids_MRCA, label="Sciuridae", color=colours_BR[2], offset=10, offset.text=3) +
  		 		geom_hilight(node=sciurids_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=cricetids_MRCA, label="Cricetidae", color=colours_BR[2], offset=10, offset.text=20) +
  		 		geom_hilight(node=cricetids_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=castors_MRCA, label="Castoridae", color=colours_BR[2], offset=10, offset.text=4, hjust=0.2) +
  		 		geom_hilight(node=castors_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=nesomyds_MRCA, label="Nesomyidae", color=colours_BR[2], offset=10, offset.text=5, hjust=0.9) +
  		 		geom_hilight(node=nesomyds_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=dasys_MRCA, label="Dasyproctidae", color=colours_BR[2], offset=10, offset.text=2, vjust=0.8) +
  		 		geom_hilight(node=dasys_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=heteromys_MRCA, label="Heteromyidae", color=colours_BR[2], offset=10, offset.text=6, hjust=0.9) +
  		 		geom_hilight(node=heteromys_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=glirids_MRCA, label="Gliridae", color=colours_BR[2], offset=10, offset.text=2) +
  		 		geom_hilight(node=glirids_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=hystrics_MRCA, label="Hystricidae", color=colours_BR[2], offset=10, offset.text=2) +
  		 		geom_hilight(node=hystrics_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=erethiz_MRCA, label="Erethizontidae", color=colours_BR[2], offset=10, offset.text=2) +
  		 		geom_hilight(node=erethiz_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) +
  		 		geom_cladelabel(node=echimys_MRCA, label="Echimyidae", color=colours_BR[2], offset=10, offset.text=2) +
  		 		geom_hilight(node=echimys_MRCA, fill=colours_BR[2], color="white", alpha=0.6, extend=0.9) + 

  		 		#not represented  		 		
				geom_cladelabel(node=ctenomys_MRCA, label="Ctenomyidae", color="gray", offset=10, offset.text=2) +
  		 		# geom_hilight(node=ctenomys_MRCA, fill="gray", color="white", alpha=0.6, extend=0.9) +
				# geom_cladelabel(node=geomys_MRCA, label="Geomyidae", color="gray", offset=10, offset.text=5, hjust=0.8) +
  		 		# geom_hilight(node=geomys_MRCA, fill="gray", color="white", alpha=0.6, extend=0.9) +
				geom_cladelabel(node=spalacs_MRCA, label="Spalacidae", color="gray", offset=10, offset.text=4, hjust=0.8) +
  		 		# geom_hilight(node=spalacs_MRCA, fill="gray", color="white", alpha=0.6, extend=0.9) +
				geom_cladelabel(node=cavys_MRCA, label="Cavviidae", color="gray", offset=10, offset.text=2, vjust=0.6) +
  		 		# geom_hilight(node=cavys_MRCA, fill="gray", color="white", alpha=0.6, extend=0.9) +
				geom_cladelabel(node=dipods_MRCA, label="Dipodidae", color="gray", offset=10, offset.text=5, hjust=0.8) #+
  		 		# geom_hilight(node=dipods_MRCA, fill="gray", color="white", alpha=0.6, extend=0.9) 
# p_rodents 


# Joint plots

trans_bg <- theme(panel.background = element_rect(fill = "transparent", colour = NA_character_), # necessary to avoid drawing panel outline
                  plot.background = element_rect(fill = "transparent", colour = NA_character_))      

# vertical orientation

heatmap.colours <- c("white",colours_BR[1])
names(heatmap.colours) <- 0:1

p_bats <- gheatmap(b, df2, offset = -2.5, color=NULL, width=0.11, 
        colnames=FALSE) +
          scale_fill_manual(values=heatmap.colours, breaks=0:1, name="Sampled") + 
          theme(legend.position="none") + 
          geom_cladelabel(node=vesp_MRCA, label="Vespertilionidae", color=colours_BR[1], offset=10, offset.text=2) +
          geom_hilight(node=vesp_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
          geom_cladelabel(node=phyllos_MRCA, label="Phyllostomidae", color=colours_BR[1], offset=10, offset.text=14, hjust=0.45) +
          geom_hilight(node=phyllos_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
          geom_cladelabel(node=molossids_MRCA, label="Molossidae", color=colours_BR[1], offset=10, offset.text=13, hjust=0.35) +
          geom_hilight(node=molossids_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
          geom_cladelabel(node=pteros_MRCA, label="Pteropodidae", color=colours_BR[1], offset=10, offset.text=2) + 
          geom_hilight(node=pteros_MRCA, fill=colours_BR[1], color="white", alpha=0.6, extend=0.9) +
          # unrepresented
          geom_cladelabel(node=rhinolophids_MRCA, label="Rhinolophidae", color="gray", offset=10, offset.text=1.2, vjust=-1) + 
          geom_cladelabel(node=hippos_MRCA, label="Hipposideridae", color="gray", offset=10, offset.text=1.5) +
          geom_cladelabel(node=nycterids_MRCA, label="Nycteridae", color="gray", offset=10, offset.text=19) +
          geom_cladelabel(node=natalids_MRCA, label="Natalidae", color="gray", offset=10, offset.text=19) +
          geom_cladelabel(node=emballos_MRCA, label="Emballonuridae", color="gray", offset=10, offset.text=20.5) 


joint_plot <- plot_grid(
            p_rodents + theme(plot.margin=grid::unit(c(-60,-10,-80,-10), "mm")) + trans_bg, 
            NULL,
            p_bats + theme(plot.margin=grid::unit(c(-70,-18,-80,-18), "mm")) + trans_bg, 
            nrow=3, rel_heights = c(0.99, -0.000001, 0.99), 
        # labels = c('A) Bats', "", 'B) Rodents'), align="l", greedy=FALSE)
        labels = c('', "", ''), align="l", greedy=FALSE)
# joint_plot

joint_plot_imgs <- ggdraw() + 
              draw_plot(joint_plot) + 
              draw_image(rat, width=0.09, hjust = -1.3, vjust = -0.45) +
              draw_image(bat, width=0.09, hjust = -1.1, vjust = 0.04) 

joint_plot_imgs

ggsave("../plots_tables/phylos_sampled_species_vert.pdf", joint_plot_imgs, height=15, width=10, bg=NULL)
ggsave("../plots_tables/phylos_sampled_species_vert.png", joint_plot_imgs, height=15, width=10, bg = 'white', dpi=900)




### PHYLOGENETIC SIGNAL IN SAMPLING ACROSS PHYLOGENY

df1_bat <- filter(df1, ord=="CHIROPTERA")
df1_rodent <- filter(df1, ord=="RODENTIA")

d_bat <- phylo.d(df1_bat, bat_tree, Species_Name, sampled, permut = 1000, rnd.bias=NULL)
d_rodent <- phylo.d(df1_rodent, rodent_tree, Species_Name, sampled, permut = 1000, rnd.bias=NULL)

# 0 = phylogenetically conserved as expected under a Brownian threshold model 
# 1 = random

summary(d_bat)
# plot(d_bat, bw=0.02)
# estimated D=0.99 (prob of resulting from random structure 0.464, p brownian=0)

summary(d_rodent)
# plot(d_rodent, bw=0.02)
# estimated D=0.805 (prob of resulting from random structure 0.001, p brownian=0)

# Note this is phylogenetic signal based on binary presence/absence, 
# not on the counts of sampling / number of papers per species



# Number of hosts per virus 
# (susceptible / susceptible via experiment / with individual pathology data)

hr <- read.csv("../clean_data/host_range_data.csv")
dat <- left_join(dat, hr)

# subset to relevant columns
df1 <- dat %>% dplyr::select(c(PaperID, Virus_ICTV, 
				primary_paper, batHosts, rodentHosts, 
				Susceptible_primary_YN, IndividualData_YN, 
				Host_Upham)) %>% unique()


df2 <- df1 %>%  group_by(Virus_ICTV) %>% 
						summarise( n_bats_rodents = batHosts+rodentHosts,
								   n_hosts_primary = n_distinct(Host_Upham[primary_paper==1 & Susceptible_primary_YN=="Y"]),
								   n_species_indiv = n_distinct(Host_Upham[Susceptible_primary_YN=="Y" & IndividualData_YN=="Y"]), 
								   .groups = 'drop') %>% unique()

dim(df2)#74...

# remove viruses with no susceptibe bat or rodent hosts
df2 <- df2[df2$n_bats_rodents>0,]

dim(df2)#70...

# remove viruses with no susceptibe bat or rodent hosts
vtree <- drop.tip(vtree, setdiff(vtree$tip.label, df2$Virus_ICTV))

df2 <- as.data.frame(df2)
rownames(df2) <- df2$Virus_ICTV 
dim(df2)# 70

df2 <- df2 %>% dplyr::select(-Virus_ICTV)


# shorten virus names except for ebolaviruses
rownames(df2) <- gsub(" virus$", "", rownames(df2))
rownames(df2)[grep("Avian ortho",rownames(df2))] <- "Newcastle Disease"
rownames(df2)[grep("Middle East",rownames(df2))] <- "MERS coronavirus"
rownames(df2)[grep("Severe acute",rownames(df2))] <- "SARS coronavirus"
rownames(df2)[grep("ebolavirus", rownames(df2), invert=T)] <- gsub(" [a-z]*virus$", "", rownames(df2)[grep("ebolavirus", rownames(df2), invert=T)])
rownames(df2) <- gsub("Venezuelan equine encephalitis", "VEE", rownames(df2))
rownames(df2) <- gsub("Western equine encephalitis", "WEE", rownames(df2))
rownames(df2) <- gsub("Eastern equine encephalitis", "EEE", rownames(df2))
rownames(df2) <- gsub(" disease$", "", rownames(df2))

vtree$tip.label <- gsub(" virus$", "", vtree$tip.label)
vtree$tip.label[grep("Avian ortho",vtree$tip.label)] <- "Newcastle Disease"
vtree$tip.label[grep("Middle East",vtree$tip.label)] <- "MERS coronavirus"
vtree$tip.label[grep("Severe acute",vtree$tip.label)] <- "SARS coronavirus"
vtree$tip.label[grep("ebolavirus", vtree$tip.label, invert=T)] <- gsub(" [a-z]*virus$", "", vtree$tip.label[grep("ebolavirus", vtree$tip.label, invert=T)])
vtree$tip.label <- gsub("Venezuelan equine encephalitis", "VEE", vtree$tip.label)
vtree$tip.label <- gsub("Western equine encephalitis", "WEE", vtree$tip.label)
vtree$tip.label <- gsub("Eastern equine encephalitis", "EEE", vtree$tip.label)
vtree$tip.label <- gsub(" disease$", "", vtree$tip.label)

v <- ggtree(ape::rotateConstr(vtree, rev(vtree$tip.label)), ladderize=FALSE) + 
					geom_tiplab(offset=0, hjust=0, 
								align=F, as_ylab = F,
								geom = "text") + 
					xlim_tree(9) +
					xlim_expand(12, 1) +
					vexpand(.04, -1)
					#  + theme(plot.margin = unit(c(l=0,t=0,r=-0.5,b=0), "cm"))
v

colnames(df2) <- c("Susceptible","Susceptible\nvia experiment","With individual\n pathology data")

v_counts <- gheatmap(v, df2, offset=11, color=NULL, width=4, colnames_offset_y=-1) +
				scale_fill_viridis_c(option = "G", begin=1, end=0, 
					na.value = "gray60", breaks=c(0,1,5, 15, 40, 100), 
					name = "Number of Bat and Rodent Species\n", trans="log1p")+
				theme(legend.position = c(0.675,-0.05),legend.direction = "horizontal") + 
				theme(plot.margin = unit(c(l=0,t=0,r=1.9,b=0), "cm"))

v_counts

ggsave("../plots_tables/virus_tree_with_counts.pdf", v_counts, width=7, height=11)
# ggsave("../plots_tables/virus_tree_with_counts.png", v_counts, width=7, height=11)



# Gap analysis via life histories
# COMBINE https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3344

combine <- read.csv("../raw_data/COMBINE/trait_data_imputed.csv")
combine$phylacine_binomial <- gsub(" ","_",combine$phylacine_binomial)

# intersect(dat$Host_Upham, combine$phylacine_binomial)#80
# setdiff(dat$Host_Upham, combine$phylacine_binomial)# Tamias_amoenus, Liomys_salvini
combine$phylacine_binomial[grep("Neotamias_amoenus",combine$phylacine_binomial)] <- "Tamias_amoenus"
combine$phylacine_binomial[grep("Heteromys_salvini",combine$phylacine_binomial)] <- "Liomys_salvini"

combine$Host_Upham <- combine$phylacine_binomial


# Identify key life-history traits

# subset to bats and rodents
combine <- combine[combine$order%in%c("Chiroptera","Rodentia"),]
apply(combine, 2, function(x) sum(!is.na(x))/length(x))

traits <- c("adult_mass_g","max_longevity_d","age_first_reproduction_d",
				"gestation_length_d","litter_size_n","litters_per_year_n",
				"weaning_age_d", "hibernation_torpor", 
				"dphy_invertebrate","dphy_vertebrate", "dphy_plant")


combine <- dplyr::select(combine, c(Host_Upham, order, traits))
combine <- combine %>% filter(!Host_Upham%in%c("Not_recognized")) %>% unique()

# for duplicated host species take average of traits
sum(duplicated(combine$Host_Upham))#104
combine <- combine %>% group_by(Host_Upham, order) %>% summarise_all(funs(mean))
sum(duplicated(combine$Host_Upham))#1
combine$Host_Upham[duplicated(combine$Host_Upham)]# Not_recognised


df1_hosts <- rbind(df1_bat, df1_rodent)
head(df1_hosts)
names(df1_hosts)[1] <- "Host_Upham"

df1_hosts <- left_join(df1_hosts, combine)

df1_hosts <- dplyr::select(df1_hosts, Host_Upham, order, sampled, traits)

# trying to update sampled to be three factors - unrepresented; represented (bats); represented (rodents)
df1_hosts$sampled <- as.character(df1_hosts$sampled)
df1_hosts$sampled[df1_hosts$sampled=="0"] <- "Unrepresented"
df1_hosts$sampled[df1_hosts$sampled=="1" & df1_hosts$order=="Chiroptera"] <- "Represented (bat)"
df1_hosts$sampled[df1_hosts$sampled=="1" & df1_hosts$order=="Rodentia"] <- "Represented (rodent)"

df1_hosts$sampled <- as.factor(df1_hosts$sampled)
df1_hosts <- df1_hosts %>% filter(!is.na(order))

colors <- scale_fill_manual(values=c(colours_BR,"darkgrey"), name="")

dist_theme <- theme_classic() 
scy <- scale_y_continuous(expand = c(0, 0))

mass_p <- ggplot(df1_hosts, aes(x = adult_mass_g, fill = sampled)) +  scale_x_log10() + colors + 
        labs(fill = "Experimentally \n infected hosts") + ylab("") +
        geom_density(alpha = 0.4) + dist_theme + scy + xlab("Adult body mass (g)")

mass_p <- mass_p + facet_wrap(~order, scales="free") + 
                              theme(
                                    strip.background = element_blank(),
                                    strip.text.x = element_text(size=14)
                                  )


options(scipen = 999)
long_p <- ggplot(df1_hosts, aes(x = max_longevity_d, fill = sampled)) + colors + 
        labs(fill = "Experimentally \n infected hosts") + ylab("") +
        geom_density(alpha = 0.4) + dist_theme + scy + xlab("Max longevity (days)")

long_p <- long_p + facet_wrap(~order, scales="free") + 
                              theme(
                                    strip.background = element_blank(),
                                    strip.text.x = element_blank()
                                  )




trait_plot <- plot_grid(
            mass_p + theme(legend.position = c(0.9, 0.8)) + theme(plot.margin=grid::unit(c(0,0,10,0), "mm")) + theme(axis.title.x=element_text(size=14, margin=margin(t=3, unit="mm"))), 
            long_p + theme(legend.position = "none") + theme(plot.margin=grid::unit(c(-4,0,0,0), "mm")) + theme(axis.title.x=element_text(size=14, margin=margin(t=3, unit="mm"))), 
            nrow=2, 
        labels = c('A)','B)'), align="v", vjust=c(1.4,-1.1), label_fontface="plain")

trait_plot

ggsave("../plots_tables/traits_sampled_species.pdf", trait_plot, width=8, height=8)
# ggsave("../plots_tables/traits_sampled_species.png", trait_plot, width=8, height=8, bg = 'white')


