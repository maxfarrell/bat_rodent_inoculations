# data_summaries.R

require(ape);packageVersion("ape")#5.5
require(dplyr);packageVersion("dplyr")#1.1.4
require(tidyr); packageVersion("tidyr")#1.3.0
require(xtable); packageVersion("xtable")#1.8.4
require(ggplot2);packageVersion("ggplot2")#3.4.4
require(ggtree);packageVersion("ggtree")#3.0.4
require(patchwork);packageVersion("patchwork")#1.1.2
require(cowplot);packageVersion("cowplot")#1.1.1


# experimental data
primary <- read.csv("../raw_data/primary_screen.csv")
speciesonly <- read.csv("../raw_data/species_only_data.csv")
indiv <- read.csv("../raw_data/individual_data.csv")

indiv %>% filter(Susceptible_YN=="Y") %>% summarise(n_papers=n_distinct(PaperID))
# 105 papers with individual data on susceptible individuals

speciesonly %>% filter(Susceptible_YN=="Y") %>% summarise(n_papers=n_distinct(PaperID))
# 32 papers from speciesonly dataframe with susceptible individuals

# 105+32

# re-create column for individualID
# Create numbering variable
indiv <- indiv %>%                              
		  group_by(PaperID) %>%
  		  mutate(IndividualID = paste(row_number(),unique(PaperID), sep="_"))

# host data
tree <- read.nexus("../raw_data/upham_tree_666.nex")
tax <- read.csv("../raw_data/taxonomy_mamPhy_5911species.csv")

bats <- tax[tax$ord=="CHIROPTERA",]
rodents <- tax[tax$ord=="RODENTIA",]

bat_tree <- drop.tip(tree, setdiff(tree$tip.label, bats$Species_Name))
rodent_tree <- drop.tip(tree, setdiff(tree$tip.label, rodents$Species_Name))

hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)

# virus names and taxonomy
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- left_join(viruses, vtax)

#subset to only host reported and upham
hosts <- hosts %>% dplyr::select(Host_reported, Host_Upham) %>% unique()

# add names
primary <- dplyr::left_join(primary, viruses)
primary <- dplyr::left_join(primary, hosts)
speciesonly <- dplyr::left_join(speciesonly, viruses)
speciesonly <- dplyr::left_join(speciesonly, hosts)
indiv <- dplyr::left_join(indiv, viruses)
indiv <- dplyr::left_join(indiv, hosts)

# merge speciesonly with individual data
dat <- full_join(indiv, speciesonly)

# removing hosts and viruses with NA names
dat <- dat[!is.na(dat$Virus_ICTV),]
dat <- dat[!is.na(dat$Host_Upham),]

# recoding disease outcomes
dat$ClinicalDisease <- dplyr::recode(dat$ClinicalDisease, "Y"= 1, "N" = 0)
dat$Susceptible_YN <- dplyr::recode(dat$Susceptible_YN, "Y"= 1, "N" = 0)
dat$Mortality_YN <- dplyr::recode(dat$Mortality_YN, "Y"= 1, "N" = 0)

# presence of symptoms
dat_symptoms <- read.csv("../clean_data/symptoms_severity.csv")

dat_symptoms <- dat_symptoms %>% rowwise() %>% mutate(

							Disease_YN = (severity_rank>1)*1,

              				Visible_YN = any(aggression, lethargy, vocalization, sensitivity, 
                          					appetite, skin_issue, weight_loss, movement, 
                          					external_damage, eye_nose, pneumonia)*1,

							Symptoms_sum = sum(aggression, lethargy, vocalization, sensitivity, 
													appetite, skin_issue, weight_loss, movement, 
													internal_damage, external_damage, eye_nose, pneumonia))

dat_symptoms_sm <- dplyr::select(dat_symptoms, PaperID, IndividualID, 
							Virus_ICTV, Virus_strain_isolate, Host_Upham, Disease_YN, Visible_YN, Symptoms_sum, severity_rank, internal_damage)

dat_symptoms_sm <- unique(dat_symptoms_sm)
dim(dat_symptoms_sm)# 2030
dim(dat)#2923

dat <- left_join(dat, dat_symptoms_sm)
dim(dat)# 2923


###############
## summaries ##
###############

# entire database
n_distinct(dat$PaperID)# 137 papers
n_distinct(dat$IndividualID, na.rm=TRUE)# 2870 individuals with individual-level data (including non-susceptibles and those with no pathology dat reported)
n_distinct(dat$IndividualID, na.rm=TRUE) + sum(dat$N_individuals, na.rm=TRUE)# 6994 individuals sampled (including non-susceptible interactions and those with no pathology reported)
n_distinct(dat$Virus_ICTV)# 67 viruses
n_distinct(dat$Host_Upham)# 93 hosts
n_distinct(paste(dat$Virus_ICTV,dat$Host_Upham))# 197 combinations
n_distinct(paste(dat$Virus_ICTV,dat$Host_Upham))/(n_distinct(dat$Virus_ICTV)*n_distinct(dat$Host_Upham))# 0.0316


# subset to only susceptible individuals and studied with pathlogy reported
dat_patho <- left_join(dat, primary)
dat_patho <- dat_patho[dat_patho$Susceptible_YN==1,]
dat_patho <- dat_patho[dat_patho$PathologyReported_YN=="Y",]
dat_patho <- dat_patho[!is.na(dat_patho$Virus_ICTV),]
dat_patho <- dat_patho[!is.na(dat_patho$Host_Upham),]


# summaries (studies with any susceptible individuals and pathology reported!)
n_distinct(dat_patho$PaperID, na.rm=TRUE)#112

# identify which articles these are for bibliography
articles <- read.csv("../raw_data/article_metadata.csv")
included <- articles[articles$PaperID %in% dat_patho$PaperID,]
dim(unique(included))#112
# View(included)
# write.csv(included, "../clean_data/articles_susceptible_pathology.csv", row.names=FALSE)

sum(dat_patho$N_individuals, na.rm=TRUE) + n_distinct(dat_patho$IndividualID, na.rm=TRUE)#5,348
n_distinct(dat_patho$Host_Upham, na.rm=TRUE)#85

dat_patho %>% group_by(Host_Upham) %>% summarise(n_indiv_total = sum(N_individuals, na.rm=T) + n_distinct(IndividualID, na.rm=TRUE)) %>% View()


n_distinct(dat_patho$Host_Upham[dat_patho$Host_order=="Chiroptera"], na.rm=TRUE)#32
n_distinct(dat_patho$Host_Upham[dat_patho$Host_order=="Rodentia"], na.rm=TRUE)#53

n_distinct(dat_patho$Virus_ICTV[dat_patho$Host_order=="Chiroptera"], na.rm=TRUE)#35
n_distinct(dat_patho$Virus_ICTV[dat_patho$Host_order=="Rodentia"], na.rm=TRUE)#28

n_distinct(dat_patho$PaperID[dat_patho$Host_order=="Chiroptera"], na.rm=TRUE)#71
n_distinct(dat_patho$PaperID[dat_patho$Host_order=="Rodentia"], na.rm=TRUE)#41

n_distinct(dat_patho$Virus_ICTV, na.rm=TRUE)#54
n_distinct(paste0(dat_patho$Host_Upham,dat_patho$Virus_ICTV), na.rm=TRUE)#149
n_distinct(paste0(dat_patho$PaperID,dat_patho$Host_Upham,dat_patho$Virus_ICTV), na.rm=TRUE)#192

# studies with individual data only
dat_patho_indiv <- dat_patho[!is.na(dat_patho$IndividualID),]
n_distinct(dat_patho_indiv$PaperID)#84
n_distinct(dat_patho_indiv$IndividualID, na.rm=TRUE)#1620
n_distinct(dat_patho_indiv$Host_Upham)#75
n_distinct(dat_patho_indiv$Virus_ICTV)#44
n_distinct(paste0(dat_patho_indiv$Host_Upham,dat_patho_indiv$Virus_ICTV))#128
n_distinct(paste0(dat_patho_indiv$PaperID,dat_patho_indiv$Host_Upham,dat_patho_indiv$Virus_ICTV))#149



###################
## Viral passage ##
###################
dat_pass <- dat %>% select(PaperID, Host_order, Virus_ICTV, Passage_host) %>% unique()
unique(dat_pass$Passage_host) %>% sort()


# sort(unique(dat_pass$Passage_host))

dat_pass$Passage_host <- factor(dat_pass$Passage_host, 
	levels=c("Artiodactyla","Carnivora", "Chiroptera", "Diptera","Lagomorpha","Primate",
		"Primate/Diptera/Rodentia","Rodentia","Rodentia/Diptera","Rodentia/Primate","Rodentia/Aves",
		"Aves","Aves/Primate","None","Unknown"))

colors <- c(
			"#785EF0",
			"limegreen",
			"#1B85BF",
			"black",
			"#20908C",
			"#D2B48C",
			"#FFB000",
			"#AB1808",
			"#00274C",
			"#F28500",
			"#659A32",
			"#0352fc",
			"#AD72D6",
			"#88CAFC",
			"grey55"
	)

pass_stack <- ggplot(dat_pass, aes(fill=Passage_host,  x=Passage_host)) + 
		geom_bar(position="dodge", stat="count") +
		# geom_point(position = position_jitter(seed = 1, width = 0.2), 
		# 	alpha=0.4) + 
		# scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		# theme(legend.title=element_text("Passage host")) + 
		labs(fill = "Passage host order") + 
		# xlab("Inoculated host order") + 
		xlab("") + 
		ylab("Number of article-virus-host combinations") +
		theme_bw() + theme(panel.grid.major.x = element_blank()) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		facet_wrap(~Host_order)

pass_stack

# ggsave("../plots_tables/virus_passage_host_order.pdf", pass_stack, width=9, height=4.5)

dat_pass %>% filter(Passage_host=="Chiroptera")
# only one study used viruses passaged in bats


######################
# Number of passages #
######################


dat_npass <- dat %>% select(PaperID, Host_order, Virus_ICTV, Passage_host, N_passages) %>% unique()
unique(dat_pass$Passage_host) %>% sort()

dat_npass$Passage_host <- factor(dat_npass$Passage_host, 
	levels=c("Artiodactyla","Carnivora", "Chiroptera", "Diptera","Lagomorpha","Primate",
		"Primate/Diptera/Rodentia","Rodentia","Rodentia/Diptera","Rodentia/Primate","Rodentia/Aves",
		"Aves","Aves/Primate","None","Unknown"))

unique(dat_npass$N_passages)
# Take the bigger of the two numbers for studies with ranges of passages
dat_npass$N_passages[grep("-", dat_npass$N_passages)]
dat_npass$N_passages <- gsub("[0-9]*-", "", dat_npass$N_passages)

# Add passages when multiple different passage souces were used
dat_npass$N_passages[grep("/", dat_npass$N_passages)]
dat_npass$N_passages[dat_npass$N_passages=="24/2"] <- "26"
dat_npass$N_passages[dat_npass$N_passages=="1/1"] <- "2"
dat_npass$N_passages[dat_npass$N_passages=="1/2"] <- "3"
dat_npass$N_passages[dat_npass$N_passages=="1/1/1"] <- "3"
dat_npass$N_passages[dat_npass$N_passages=="1/5"] <- "6"

# set "Multiple" and "Low" to "Unknown"
dat_npass$N_passages[dat_npass$N_passages=="Multiple"] <- "Unknown"
dat_npass$N_passages[dat_npass$N_passages=="Low"] <- "Unknown"

# Set "Unknown" to NA
dat_npass$N_passages[dat_npass$N_passages=="Unknown"] <- NA

# Set to numeric
dat_npass$N_passages <- as.numeric(dat_npass$N_passages)

dat_npass %>% group_by(Host_order) %>% 
			reframe(mean=mean(N_passages, na.rm=TRUE),
					median=median(N_passages, na.rm=TRUE), 
					min=min(N_passages, na.rm=TRUE),
					max=max(N_passages, na.rm=TRUE))

t.test(dat_npass$N_passages ~ dat_npass$Host_order, mu = 0, alternative = "two.sided")

colours_BR <- c("#1B85BF", "#AB1808")

pass_box <- dat_npass %>% ggplot(., aes(y=N_passages, x= Host_order, 
		color=Host_order, fill=Host_order)) + 
		geom_boxplot(outlier.shape=NA, alpha=0.2) +
		geom_jitter(width=0.3, alpha=0.7)+
		scale_color_manual(values=colours_BR) + 
		scale_fill_manual(values=colours_BR) + 
		theme_bw() + theme(panel.grid.major.x = element_blank(), legend.position = "none") +
		scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(-10, 0, 1, 2, 3, 5, 10, 20, 40, 100, 200, 500)) +
		# labs(color = "Inoculated host order")+
		xlab("\nInoculated host order") + ylab("Number of passages prior to inoculation") +
		theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

pass_box

# ggsave("../plots_tables/n_viral_passages.pdf", pass_box, width=5.5, height=5)



##################################################
# Numbers of species showing any sign of disease #
##################################################


dat_patho %>% select(-c(PaperID, IndividualID)) %>% unique() %>% group_by(Virus_ICTV, Host_Upham, Host_order, family) %>%
				summarise(any_patho=(sum(Disease_YN, na.rm=TRUE)>0)*1)


# find interactions with any reported disease
dat_patho <- dat_patho %>% select(-c(PaperID, IndividualID)) %>% unique() %>% group_by(Virus_ICTV, Host_Upham, Host_order, family) %>%
				summarise(any_patho=(sum(Disease_YN, na.rm=TRUE)>0)*1)

dat_patho <- dat_patho[!is.na(dat_patho$Host_order),]

dat_patho %>% group_by(Host_order) %>% summarise(n_hosts=n_distinct(Host_Upham), n_hv=n_distinct(paste(Host_Upham, Virus_ICTV)), n_diseased=sum(any_patho), n_not_diseased=sum(any_patho==0))
# 48 bat-virus interactions showed signs of disease (18 did not) 
# 40 rodent-virus interactions showed signs of disease (43 did not)
# 48/(48+18)# 73%
# 40/(40+43)# 48%

patho_v_fam <- dat_patho %>% group_by(family) %>% summarise(
	n_diseased_bats=sum(any_patho[Host_order=="Chiroptera"]), n_not_diseased_bats=sum(any_patho[Host_order=="Chiroptera"]==0),
	n_diseased_rats=sum(any_patho[Host_order=="Rodentia"]), n_not_diseased_rats=sum(any_patho[Host_order=="Rodentia"]==0))

patho_v_fam <- patho_v_fam %>% group_by(family) %>% summarise(
	bats_diseased=paste0(n_diseased_bats,"/",sum(n_diseased_bats,n_not_diseased_bats)),
	rats_diseased=paste0(n_diseased_rats,"/",sum(n_diseased_rats,n_not_diseased_rats)))

names(patho_v_fam) <- c("Virus family", "Bats", "Rodents")


latex <- print.xtable(xtable(patho_v_fam, caption = "Proportions of host-virus pairs showing any sign of disease across virus families and host order"), 
						include.rownames=FALSE, print.results = FALSE)

# writeLines(
#   c(
#     "\\documentclass[11pt]{article}",
#     "\\begin{document}",
#     "\\thispagestyle{empty}",
#     latex,
#     "\\end{document}"
#   ),
#   "../plots_tables/patho_by_virus_fam.tex"
# )

# tools::texi2pdf("../plots_tables/patho_by_virus_fam.tex", clean = TRUE)



# # REMOVING LYSSAVIRUSES
# dat_nolyss <- dat_patho[grep("Lyssa", dat_patho$Virus_ICTV, invert=TRUE),]
# dat_nolyss %>% group_by(Host_order) %>% summarise(n_hosts=n_distinct(Host_Upham), n_hv=n_distinct(paste(Host_Upham, Virus_ICTV)), n_diseased=sum(any_patho), n_not_diseased=sum(any_patho==0))
# # 31 bat-virus interactions showed signs of disease (17 did not)
# # 39 rodent-virus interactions showed signs of disease (43 did not)

# patho_v_fam <- dat_nolyss %>% group_by(family) %>% summarise(
# 	n_diseased_bats=sum(any_patho[Host_order=="Chiroptera"]), n_not_diseased_bats=sum(any_patho[Host_order=="Chiroptera"]==0),
# 	n_diseased_rats=sum(any_patho[Host_order=="Rodentia"]), n_not_diseased_rats=sum(any_patho[Host_order=="Rodentia"]==0))

# patho_v_fam <- patho_v_fam %>% group_by(family) %>% summarise(
# 	bats_diseased=paste0(n_diseased_bats,"/",sum(n_diseased_bats,n_not_diseased_bats)),
# 	rats_diseased=paste0(n_diseased_rats,"/",sum(n_diseased_rats,n_not_diseased_rats)))

# patho_v_fam

# names(patho_v_fam) <- c("Virus family", "Bats", "Rodents")

# latex <- print.xtable(xtable(patho_v_fam, caption = "Proportions of host-virus pairs showing any sign of disease across virus families and host order. Lyssaviruses removed."), 
# 						include.rownames=FALSE, print.results = FALSE)

# # writeLines(
# #   c(
# #     "\\documentclass[11pt]{article}",
# #     "\\begin{document}",
# #     "\\thispagestyle{empty}",
# #     latex,
# #     "\\end{document}"
# #   ),
# #   "../plots_tables/patho_by_virus_fam_nolyss.tex"
# # )

# # tools::texi2pdf("../plots_tables/patho_by_virus_fam_nolyss.tex", clean = TRUE)



##################################################
# Numbers of species showing any clinical disease #
##################################################


# subset to only susceptible individuals and studied with pathlogy reported
dat_clin <- left_join(dat, primary)
dat_clin <- dat_clin[dat_clin$Susceptible_YN==1,]
dat_clin <- dat_clin[dat_clin$PathologyReported_YN=="Y",]
dat_clin <- dat_clin[!is.na(dat_clin$Virus_ICTV),]
dat_clin <- dat_clin[!is.na(dat_clin$Host_Upham),]


# find interactions with any reported clinical disease
dat_clin <- dat_clin %>% select(-c(PaperID, IndividualID)) %>% unique() %>% group_by(Virus_ICTV, Host_Upham, Host_order, family) %>%
				summarise(any_clin=(sum(ClinicalDisease, na.rm=TRUE)>0)*1)

dat_clin <- dat_clin[!is.na(dat_clin$Host_order),]
dat_clin <- dat_clin[!is.na(dat_clin$Host_order),]

# which paramyxivurs-host associations showed clinical disease and any disease?

dat_clin %>% filter(any_clin==1, family=="Paramyxoviridae") %>% unique()
dat_clin %>% filter(any_clin==1, family=="Filoviridae") %>% unique()

dat_patho %>% filter(any_patho==1, family=="Paramyxoviridae") %>% unique()
dat_patho %>% filter(any_patho==1, family=="Filoviridae") %>% unique()

dat_patho %>% filter(any_patho==0, family=="Paramyxoviridae") %>% unique()
dat_patho %>% filter(any_patho==0, family=="Filoviridae") %>% unique()


dat_clin %>% group_by(Host_order) %>% summarise(n_hosts=n_distinct(Host_Upham), n_hv=n_distinct(paste(Host_Upham, Virus_ICTV)), n_clin_diseased=sum(any_clin), n_not_clin_diseased=sum(any_clin==0))
# 35 bat-virus interactions showed signs of disease (31 did not) 
# 32 rodent-virus interactions showed signs of disease (51 did not)
# 35/(35+31)# 53%
# 32/(32+51)# 38.5%

clin_v_fam <- dat_clin %>% group_by(family) %>% summarise(
	n_clinicallydiseased_bats=sum(any_clin[Host_order=="Chiroptera"]), n_not_clinicallydiseased_bats=sum(any_clin[Host_order=="Chiroptera"]==0),
	n_clinicallydiseased_rats=sum(any_clin[Host_order=="Rodentia"]), n_not_clinicallydiseased_rats=sum(any_clin[Host_order=="Rodentia"]==0))

clin_v_fam <- clin_v_fam %>% group_by(family) %>% summarise(
	bats_clinicallydiseased=paste0(n_clinicallydiseased_bats,"/",sum(n_clinicallydiseased_bats,n_not_clinicallydiseased_bats)),
	rats_clinicallydiseased=paste0(n_clinicallydiseased_rats,"/",sum(n_clinicallydiseased_rats,n_not_clinicallydiseased_rats)))

names(clin_v_fam) <- c("Virus family", "Bats", "Rodents")

latex <- print.xtable(xtable(clin_v_fam, caption = "Proportions of host-virus pairs with author reported clinical disease across virus families and host order"), 
						include.rownames=FALSE, print.results = FALSE)

# writeLines(
#   c(
#     "\\documentclass[11pt]{article}",
#     "\\begin{document}",
#     "\\thispagestyle{empty}",
#     latex,
#     "\\end{document}"
#   ),
#   "../plots_tables/clinicaldisease_by_virus_fam.tex"
# )

# tools::texi2pdf("../plots_tables/clinicaldisease_by_virus_fam.tex", clean = TRUE)


#########################

# Counts of host-virus interactions with sub-clinical disease only

# subset to only susceptible individuals and studied with pathlogy reported
dat_subclin <- left_join(dat, primary)
dat_subclin <- dat_subclin[dat_subclin$Susceptible_YN==1,]
dat_subclin <- dat_subclin[dat_subclin$PathologyReported_YN=="Y",]
dat_subclin <- dat_subclin[!is.na(dat_subclin$Virus_ICTV),]
dat_subclin <- dat_subclin[!is.na(dat_subclin$Host_Upham),]

# find interactions with max severity rank is 2 (i.e. subclinical infection)
dat_subclin <- dat_subclin %>% select(-c(PaperID, IndividualID)) %>% unique() %>% group_by(Virus_ICTV, Host_Upham, Host_order, family) %>%
				summarise(any_subclin=(sum(max(severity_rank)==2, na.rm=TRUE)>0)*1)

dat_subclin <- dat_subclin[!is.na(dat_subclin$Host_order),]
dat_subclin <- dat_subclin[!is.na(dat_subclin$Host_order),]

# which paramyxivurs and filovirus-host associations max severity was subclinical disease?

dat_subclin %>% filter(any_subclin==1, family=="Paramyxoviridae") %>% unique()
dat_subclin %>% filter(any_subclin==1, family=="Filoviridae") %>% unique()

dat_subclin %>% group_by(Host_order) %>% summarise(n_hosts=n_distinct(Host_Upham), n_hv=n_distinct(paste(Host_Upham, Virus_ICTV)), n_subclin_diseased=sum(any_subclin), n_not_subclin_diseased=sum(any_subclin==0))
# 10 bat-virus interactions showed signs of subclinical disease (56 did not) 
# 10/(10+56)# 15%
# 0 rodent-virus interactions showed signs of disease (83 did not)

subclin_v_fam <- dat_subclin %>% group_by(family) %>% summarise(
	n_subclinicallydiseased_bats=sum(any_subclin[Host_order=="Chiroptera"]), n_not_subclinicallydiseased_bats=sum(any_subclin[Host_order=="Chiroptera"]==0),
	n_subclinicallydiseased_rats=sum(any_subclin[Host_order=="Rodentia"]), n_not_subclinicallydiseased_rats=sum(any_subclin[Host_order=="Rodentia"]==0))

subclin_v_fam <- subclin_v_fam %>% group_by(family) %>% summarise(
	bats_subclinicallydiseased=paste0(n_subclinicallydiseased_bats,"/",sum(n_subclinicallydiseased_bats,n_not_subclinicallydiseased_bats)),
	rats_subclinicallydiseased=paste0(n_subclinicallydiseased_rats,"/",sum(n_subclinicallydiseased_rats,n_not_subclinicallydiseased_rats)))

names(subclin_v_fam) <- c("Virus family", "Bats", "Rodents")

latex <- print.xtable(xtable(subclin_v_fam, caption = "Proportions of host-virus pairs with subclinical disease across virus families and host order"), 
						include.rownames=FALSE, print.results = FALSE)

# writeLines(
#   c(
#     "\\documentclass[11pt]{article}",
#     "\\begin{document}",
#     "\\thispagestyle{empty}",
#     latex,
#     "\\end{document}"
#   ),
#   "../plots_tables/subclinicaldisease_by_virus_fam.tex"
# )

# tools::texi2pdf("../plots_tables/subclinicaldisease_by_virus_fam.tex", clean = TRUE)




### COMBINED TABLE OF ANY DISEASE, CLINICAL DISEASE, AND SUB-CLINICAL DISEASE

combo_tab <- cbind(patho_v_fam, clin_v_fam[2:3], subclin_v_fam[2:3])
latex <- print.xtable(xtable(combo_tab, caption = "Proportions of host-virus pairs across disease reporting levels"), 
						include.rownames=FALSE, print.results = FALSE)
# writeLines(
#   c(
#     "\\documentclass[11pt]{article}",
#     "\\begin{document}",
#     "\\thispagestyle{empty}",
#     latex,
#     "\\end{document}"
#   ),
#   "../plots_tables/combo_disease_by_virus_fam.tex"
# )



#########################

# which interactions were not susceptible:
dat_suscep <- dat %>% group_by(Virus_ICTV, Host_Upham, Host_order) %>% summarise(suscep=any(as.logical(Susceptible_YN), na.rm=TRUE))
dat_suscep[dat_suscep$suscep==FALSE,]


# How many viruses infect both a bat and a rodent?
vho <- data.frame(with(dat_patho, table(Virus_ICTV, Host_order)))
vho$Freq[vho$Freq>1] <- 1
vho <- tidyr::pivot_wider(vho, names_from=Host_order, values_from=Freq)
vho$both <- vho$Chiroptera * vho$Rodentia
sum(vho$both)# 9 viruses in both a rodent and a bat (susceptible with pathology reported)

br_viruses <- vho$Virus_ICTV[vho$both==1] %>% unique()
sort(br_viruses)

dat_br <- dat[dat$Virus_ICTV %in% br_viruses,]
dat_br <- dat_br[!is.na(dat_br$severity_rank),]
length(unique(dat_br$Virus_ICTV))#9
sort(unique(dat_br$Virus_ICTV))#9
length(unique(dat_br$PaperID))#55
length(unique(dat_br$Host_Upham))#53

# how many studies inoculate both bats and rodents with the same virus?

dat_bandr <- dat_br %>% group_by(PaperID) %>% summarise(n_orders=n_distinct(Host_order))
n_distinct(dat_bandr$PaperID)#5
dat_bandr %>% filter(n_orders>1) %>% summarise(n_papers=n_distinct(PaperID)) # ZERO


# distribution of sample sizes across studies
dat %>% filter(Susceptible_YN==1 & !is.na(severity_rank)) %>% select(family,  PaperID, IndividualID, N_individuals) %>% unique() %>% 
				summarise(n_individuals = sum(n_distinct(IndividualID, na.rm=TRUE), N_individuals, na.rm=TRUE)) %>%
				arrange(desc(n_individuals)) %>% select(n_individuals) %>% sum() #5622 total individuals



# investigating reservoir data #
virus_traits <- read.csv("../clean_data/Virus_traits.csv")

virus_traits$Reservoir[virus_traits$Reservoir%in%c("Vespbat", "Pterobat")] <- "Chiroptera"
virus_traits$Reservoir[virus_traits$Reservoir%in%c("Rodent")] <- "Rodentia"

vtraits <- dplyr::select(virus_traits, Virus_ICTV, Reservoir)

dat <- left_join(dat, vtraits, relationship="many-to-many") %>% unique()
# dim(dat)# 2445 because some viruses are associated with multiple reservoirs 

# if any known reservoir order match the infected host order, set match to 1
dat$Reservoir_match <- NA

for (i in 1:nrow(dat)) {

	virus <- dat$Virus_ICTV[i]
	reservoirs <- vtraits$Reservoir[vtraits$Virus_ICTV%in%virus]
	match <- any(dat$Host_order[i] %in% reservoirs)
	dat$Reservoir_match[i] <- (match*1) 

}

# orphan viruses 
dat$Virus_ICTV[dat$Reservoir=="Orphan"] %>% sort() %>% unique()
# notable include MERS and ebolaviruses...

# which ones have NA for reservoir
dat$Virus_ICTV[is.na(dat$Reservoir)] %>% sort() %>% unique()

# setting viruses with no clear reservoir to be NA for Reservoir_match
dat$Reservoir_match[dat$Reservoir %in% c("Orphan", "MULTIPLE")] <- NA

dat$Reservoir[dat$Reservoir %in% c("Anseriformes", 
										 	"Charadriiformes",
											"Columbiformes",
											"Galliformes",
											"Passeriformes",
											"Pelecaniformes",
											"Suliformes")] <- "Aves"

dat$Reservoir[dat$Reservoir %in% c("Human","NonHumanPrimate")] <- "Primates"
dat$Reservoir[dat$Reservoir %in% c("Artiodactyl")] <- "Artiodactyla"
dat$Reservoir[dat$Reservoir %in% c("Carnivore")] <- "Carnivora"


# how many viruses with susceptible hosts and pathology data have reservoir information?
dat_patho_res <- left_join(dat_patho, dat %>% select(Virus_ICTV, Reservoir))
dat_patho_res$Reservoir[dat_patho_res$Reservoir %in% c("Orphan", "Multiple")] <- NA

length(unique(dat_patho_res$Virus_ICTV[!is.na(dat_patho_res$Reservoir)]))# 36
length(unique(dat_patho_res$Virus_ICTV[is.na(dat_patho_res$Reservoir)]))# 18
unique(dat_patho_res$Virus_ICTV[is.na(dat_patho_res$Reservoir)])# 18


# are Host order and reservoir match independent?
# host-virus level
dat_sm <- dat %>% filter(Susceptible_YN==1, !is.na(severity_rank))
dat_sm <- dat_sm[!dat_sm$Reservoir%in%c("Orphan"),]
dat_sm <- dat_sm[!is.na(dat_sm$Reservoir),]

dat_sm <- dplyr::select(dat_sm, Virus_ICTV, Host_Upham, Host_order, Reservoir_match) %>% unique()
dat_sm <- ungroup(dat_sm) %>% dplyr::select(-PaperID) %>% unique()

table(dat_sm$Host_order, dat_sm$Reservoir_match)

# do the values of one Reservoir_match depend on the value of Host order? 
# chi square
# Null hypothesis: There are no relationships between the categorical variables. 
# Alternative hypothesis: There are relationships between the categorical variables. 
# Therefore, knowing the value of one variable does help you predict the value of another variable.

chisq.test(dat_sm$Host_order, dat_sm$Reservoir_match, correct=FALSE)
# p=0.02, therefore reject null - there is evidence of a relationship between variables
# 35/(25+35)# 0.58
# 24/(24+40)# 0.375


### Reservoir match versus Host Evolutionary Isolation

# adding host specificity data
host_range <- read.csv("../clean_data/host_range_data.csv")

dat_hr <- left_join(dat, host_range)

dat_hr <- dat_hr %>% filter(Susceptible_YN==1, !is.na(severity_rank)) %>%
			ungroup() %>% select(-PaperID) %>%
			select(Host_order, Virus_ICTV, Host_Upham, Reservoir, Reservoir_match, 
					ses_pd, EvoIso, mean_reservoir_dist) %>% unique() %>%
			filter(!is.na(Reservoir_match)) %>% filter(!is.na(Reservoir))

dat_hr_noRes <- dat_hr %>% select(-Reservoir) %>% unique()
cor(dat_hr_noRes$Reservoir_match, dat_hr_noRes$EvoIso)# -0.55

colours_BR <- c("#1B85BF", "#AB1808")

rhei_bplot <- ggplot(dat_hr_noRes, aes(x=factor(Reservoir_match), y=EvoIso, color=Host_order)) + 
				geom_boxplot(width=0.2, outlier.shape=NA) +
				geom_jitter(width=0.2, alpha=0.3) + facet_grid(~Host_order) +
				scale_color_manual(values=colours_BR) + theme_bw() +
				theme(legend.position="none") +
				xlab("Reservoir match") +
				ylab("Host Evolutionary Isolation")

rhei_bplot
# ggsave("../plots_tables/HEI_by_ResMatch.pdf", width=6, height=3)



rhei_res <- ggplot(dat_hr, aes(x=factor(Reservoir), y=EvoIso, color=Host_order, group=factor(Reservoir_match))) + 
				geom_jitter(aes(shape=factor(Reservoir_match)), size=3, alpha=0.5, width=0.3) +
				scale_color_manual(values=colours_BR) +
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, vjust=1)) +
				labs(shape="Reservoir match", colour="Host order") +
				xlab("Presumed reservoir") + 
				ylab("Host Evolutionary Isolation")

rhei_res
# ggsave("../plots_tables/HEI_by_Reservoir.pdf", width=7, height=4)



colours_YP <-c("#FFC20A","#800080")
dat_hr_sm <- dat_hr %>% ungroup() %>% select(-Reservoir) %>% unique()

rhei_hists <-  ggplot(dat_hr_sm, aes(x=EvoIso, fill=factor(Reservoir_match))) + 
					geom_histogram(alpha=0.5, position="identity") +
					facet_grid(Host_order~.) + 
					scale_fill_manual(values=colours_YP) + 
					theme_bw() + 
					scale_y_continuous(expand = c(0, 0), limits=c(0,14)) + 
					labs(fill="Reservoir match") + 
					ylab("Number of unique host-virus interactions")	


cor(dat_hr_sm$EvoIso, dat_hr_sm$Reservoir_match)
cor(dat_hr_sm$EvoIso[dat_hr_sm$Host_order=="Chiroptera"], dat_hr_sm$Reservoir_match[dat_hr_sm$Host_order=="Chiroptera"])
cor(dat_hr_sm$EvoIso[dat_hr_sm$Host_order=="Rodentia"], dat_hr_sm$Reservoir_match[dat_hr_sm$Host_order=="Rodentia"])




#######################################################
# Plot accumulation of hosts and viruses through time #
#######################################################

dat %>% group_by(Host_order) %>%
			summarise(n_distinct(paste0(Host_Upham, Virus_ICTV)))

dat %>% filter(Susceptible_YN==1, !is.na(severity_rank)) %>% group_by(Host_order) %>%
			summarise(n_distinct(paste0(Host_Upham, Virus_ICTV)))


meta <- read.csv("../raw_data/article_metadata.csv")

id_year <- meta %>% select(PaperID, year)

# hist(id_year$year)
range(id_year$year, na.rm=T)
# looks like we have studies 1919-2022 in "Papers"
# not all of these are included in our study

# dat <- dat %>% filter(Susceptible_YN==1)
# do we want to do this?
# maybe we want to include things people looked at but didnt find successful infection?

dat_year <- dplyr::left_join(dat, id_year)

dat_year <- dat_year %>% filter(Susceptible_YN==1, !is.na(severity_rank)) %>% 
					select(year, Virus_ICTV, Host_Upham, Host_order, PaperID) %>% unique() %>% arrange(year)

range(dat_year$year, na.rm=TRUE)
# 1931 - 2022


dat_year
range(dat_year$year)
# hist(dat_year$year, breaks=20)
# dip in effort 1986-1996
# almost biphasic

# plot of number of studies for bats and rodents
npaper_year <- dat_year %>% select(PaperID, Host_order, year) %>% unique() %>%
				group_by(Host_order, year) %>% mutate(n_papers=n_distinct(PaperID))


ppyear <- ggplot(npaper_year, aes(x=year, y=n_papers, fill=Host_order)) + geom_bar(stat="identity")
# ppyear

dat_year <- dat_year %>%
  group_by(Host_order)%>% # if you have a third variable and you want to achieve the same results for each group
  mutate(uniq_viruses = cumsum(!duplicated(Virus_ICTV)),
  			uniq_hosts = cumsum(!duplicated(Host_Upham)),
  			uniq_hp = cumsum(!duplicated(paste0(Virus_ICTV, Host_Upham)))) %>%
  group_by(year, Host_order) %>% # add group variable for more layers
  summarise(uniq_viruses = last(uniq_viruses),
  			uniq_hosts = last(uniq_hosts),
  			uniq_hp	= last(uniq_hp))


colours_BR <- c("#1B85BF", "#AB1808")
colours_YB <-c("#FFC20A","#0C7BDC")

# ggplot(dat_year) + geom_line(aes(x=year, y=uniq_hosts, color=Host_order))

dat_long <- gather(dat_year, variable, value, uniq_viruses:uniq_hp, factor_key=TRUE)

levels(dat_long$variable) <- c("Viruses", "Hosts", "Host-Virus combinations")

ggplot(dat_long) + geom_line(aes(x=year, y=value, color=Host_order, linetype=Host_order), linewidth=1.2) + facet_wrap(~variable) +
					theme_classic() + scale_color_manual(values=colours_BR) +
					scale_y_continuous(expand = c(0, 0), limits = c(0, 90), breaks=seq(0,90,10)) + 
					scale_x_continuous(limits = c(1930, 2022), breaks=seq(1930,2020,20)) + 
					labs(color = element_blank()) + ylab("Total number investigated") + xlab("Year") + 
					theme(legend.position = c(0.09, 0.88)) +
					guides(colour = guide_legend(override.aes = list(colour = colours_BR, linetype = c(1,2)))) +
					scale_linetype(guide = "none") +
					theme(legend.key.width=unit(1.2,"cm"), legend.key.height=unit(0.4,"cm"))

ggsave("../plots_tables/accumulation_curves.pdf", width=8.0, height=3.2)
# ggsave("../plots_tables/accumulation_curves.png", width=9.2, height=3.2)

# Ten year gap in rodent accumulation between 1986-1996...



###########################################
## Counts of research effort across taxa ##
###########################################

# number of papers, hosts, and individuals by virus
dat_virus <- dat %>% filter(Susceptible_YN==1) %>% 
						select(PaperID, Virus_ICTV, Host_Upham, IndividualID, Host_order) %>%
						unique() %>%
						group_by(Virus_ICTV) %>%
						mutate( n_papers=n_distinct(PaperID),
								n_hosts=n_distinct(Host_Upham),
								n_individuals=n_distinct(IndividualID)) %>%
						select(Virus_ICTV, n_papers, n_hosts, n_individuals) %>%
						unique() %>% arrange(desc(n_hosts), desc(n_individuals))
dat_virus


# by viral family
dat_vf <- dat %>% select(PaperID, Virus_ICTV, Host_Upham, family, IndividualID, N_individuals) %>%
				filter(!is.na(PaperID)) %>% 
				unique() %>% group_by(family) %>% 
				mutate(   n_viruses=sum(n_distinct(Virus_ICTV, na.rm = TRUE)),
							 n_hosts=sum(n_distinct(Host_Upham, na.rm = TRUE)),
							 n_individuals=sum(n_distinct(IndividualID, na.rm=TRUE), N_individuals, na.rm = TRUE),
							 n_papers=sum(n_distinct(PaperID, na.rm = TRUE))) %>%
				select(-c(PaperID:N_individuals, )) %>% unique() %>%
				arrange(desc(n_viruses))

dat_vf <- as.data.frame(dat_vf)
names(dat_vf) <- c("Family","N viruses", "N Hosts", "N individuals","N papers")

latex <- print.xtable(xtable(dat_vf), include.rownames=FALSE, print.results = FALSE)

# writeLines(
#   c(
#     "\\documentclass[11pt]{article}",
#     "\\begin{document}",
#     "\\thispagestyle{empty}",
#     latex,
#     "\\end{document}"
#   ),
#   "../plots_tables/table_sampling_virus_family.tex"
# )

# tools::texi2pdf("../plots_tables/table_sampling_virus_family.tex", clean = TRUE)


# by viral order and family
dat_vtax <- dat %>% select(PaperID, Virus_ICTV, Host_Upham, kingdom, phylum, class, order, family, IndividualID, N_individuals) %>%
				unique %>% group_by(kingdom, phylum, class, order, family) %>% 
				mutate(  n_viruses=sum(n_distinct(Virus_ICTV, na.rm = TRUE)),
						 n_hosts=sum(n_distinct(Host_Upham, na.rm = TRUE)),
						 n_individuals=sum(n_distinct(IndividualID, na.rm = TRUE), N_individuals, na.rm=TRUE),
						 n_papers=sum(n_distinct(PaperID, na.rm = TRUE)),
						 sampling=(n_viruses * n_hosts * log(n_individuals) * n_papers)) %>%
				select(-c(PaperID:N_individuals, )) %>% unique() %>%
				arrange(desc(n_viruses))
# View(dat_vtax)

# Bias for ss vs ds viruses / DNA/RNA viruses


# number of papers, viruses, individuals by host
dat_host <- dat %>% select(PaperID, Virus_ICTV, Host_Upham, Host_order, IndividualID, N_individuals) %>%
						unique() %>%
						group_by(Host_order, Host_Upham) %>%
						mutate( n_papers=n_distinct(PaperID),
								n_viruses=n_distinct(Virus_ICTV),
								n_individuals=sum(n_distinct(IndividualID, na.rm=TRUE), N_individuals, na.rm=TRUE)) %>%
						select(Host_order, Host_Upham, n_papers, n_viruses, n_individuals) %>%
						unique() %>% arrange(desc(n_individuals))
dat_host

# what proportion of bat families are sampled?
#
bat_families <- unique(tax$fam[tax$ord=="CHIROPTERA"])

tax_test <- tax
names(tax_test)[1] <- "Host_Upham"
dat_hosttax <- left_join(dat, tax_test)

dat_hosttax <- dat_hosttax %>% select(Host_Upham, Host_order, fam, Virus_ICTV, IndividualID, PaperID) %>%
						unique() %>%
						group_by(Host_order, fam) %>%
						mutate(  n_papers=sum(n_distinct(PaperID)),
						 n_species=sum(n_distinct(Host_Upham)),
						 n_viruses=sum(n_distinct(Virus_ICTV)),
						 n_individuals=sum(n_distinct(IndividualID))) %>%						
						select(Host_order, fam, n_species, n_viruses, n_individuals, n_papers) %>%
						unique() %>% arrange(Host_order, -n_species)


dat_hosttax

# intersect(dat_hosttax$fam, bat_families)
# only four families sampled out of 18 bat families
# PTEROPODIDAE
# VESPERTILIONIDAE
# PHYLOSTOMIDAE
# MOLOSSIDAE
# 4/18 = 22%


# what proportion of rodent families are sampled?
#
rodent_families <- unique(tax$fam[tax$ord=="RODENTIA"])
intersect(dat_hosttax$fam, rodent_families)
# 11 rodent families sampled out of 34 families
length(intersect(dat_hosttax$fam, rodent_families))/length(rodent_families)
# 11/34 = 32%


# host families by size
tax %>% select(Species_Name, fam, ord) %>%
			filter(ord=="CHIROPTERA") %>%
			group_by(fam) %>% unique() %>%
			mutate(n_species=n_distinct(Species_Name)) %>%
			select(fam, n_species) %>% unique() %>%
			arrange(n_species)

# Sampling is among four most speciose bat families

# We aren't sampling rhinolophids
# looks like there are two studies of Rhinolophids, one at sp. level, the other no individuals

# how many viruses have bats / rodents as reservoirs?
# should do this separately to create vtraits
# make sure to include the primary screen viruses

dat %>% select(PaperID, Virus_ICTV, Host_Upham, Host_order, IndividualID) %>%
				unique %>% group_by(Host_order) %>% 
				mutate(  n_papers=sum(n_distinct(PaperID)),
						 n_hosts=sum(n_distinct(Host_Upham)),
						 n_viruses=sum(n_distinct(Virus_ICTV)),
						 n_individuals=sum(n_distinct(IndividualID))) %>%
				select(-c(PaperID:IndividualID, )) %>% unique()

# sense check
length(sort(unique(dat$Host_Upham[dat$Host_order=="Rodentia"])))# 59
length(sort(unique(dat$IndividualID[dat$Host_order=="Chiroptera"])))# 1995 (one is NA for papers with no individual level data reported)


# Table for export
dat_hosttax$fam <- Hmisc::capitalize(tolower(dat_hosttax$fam))
names(dat_hosttax) <- c("Order", "Family", "N species", "N viruses", "N individuals", "N papers")

latex <- print.xtable(xtable(dat_hosttax), include.rownames=FALSE, print.results = FALSE)

# writeLines(
#   c(
#     "\\documentclass[11pt]{article}",
#     "\\begin{document}",
#     "\\thispagestyle{empty}",
#     latex,
#     "\\end{document}"
#   ),
#   "../plots_tables/table_sampling_host_family.tex"
# )

# tools::texi2pdf("../plots_tables/table_sampling_host_family.tex", clean = TRUE)



##########################
## Host x Virus Heatmap ##
##########################


# Maximum observed severity

colours_BR <- c("#1B85BF", "#AB1808")


# make virus tree
vtax <- as.data.frame(unclass(vtax), stringsAsFactors=TRUE)
frm <- ~superkingdom/realm/kingdom/phylum/class/order/family/genus/Virus_ICTV
vtree <- as.phylo(frm, data = vtax, collapse=FALSE)
vtree$edge.length <- rep(1, nrow(vtree$edge))

# include only viruses in subset data (e.g. after removing mole)
vtree <- drop.tip(vtree, setdiff(vtree$tip.label, dat$Virus_ICTV))
# plot(vtree)


# load host phylo
tree <- read.nexus("../raw_data/upham_tree_666.nex")
host_tree <- drop.tip(tree, setdiff(tree$tip.label, dat$Host_Upham))
host_tree$tip.label <- gsub("_"," ",host_tree$tip.label)


# Heatmap with severity_rank 
# filter out non-susceptible and studies where pathology was not reported

dat_patho <- left_join(dat, primary)
dat_patho <- dat_patho[dat_patho$Susceptible_YN==1,]
dat_patho <- dat_patho[dat_patho$PathologyReported_YN=="Y",]
dat_patho <- dat_patho[!is.na(dat_patho$Virus_ICTV),]
dat_patho <- dat_patho[!is.na(dat_patho$Host_Upham),]

# updating this to filter out non-susceptible and studies where pathology was not reported...
com2 <- dat_patho %>% ungroup() %>%
				filter(Susceptible_YN==1) %>%
				filter(PathologyReported_YN=="Y") %>%
				dplyr::select(PaperID, Virus_ICTV, Host_Upham, IndividualID, N_individuals, severity_rank, Host_order) %>% 
				unique() %>% 
				group_by(Virus_ICTV, Host_Upham, Host_order) %>%
				summarise(n_individuals = as.numeric(sum(n_distinct(IndividualID, na.rm=TRUE) + sum(N_individuals, na.rm=TRUE))),
							max_severity=max(severity_rank, na.rm=TRUE),
							mean_severity=mean(severity_rank, na.rm=TRUE)) 


sum(com2$n_individuals)

com2$Host_Upham <- gsub("_", " ", com2$Host_Upham)
com2 <- data.frame(com2)
com2$Virus_ICTV %>% unique() %>% sort()

# switch to ICTV virus names for plotting (not latin binomial names)
com2$Virus_ICTV[grep("Avian ortho",com2$Virus_ICTV)] <- "Newcastle disease virus"
com2$Virus_ICTV[grep("Middle East",com2$Virus_ICTV)] <- "MERS coronavirus"
com2$Virus_ICTV[grep("Severe acute",com2$Virus_ICTV)] <- "SARS coronavirus"
com2$Virus_ICTV <- gsub("Saint ", "St. ", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("Venezuelan equine encephalitis", "VEE", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("Western equine encephalitis", "WEE", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("Eastern equine encephalitis", "EEE", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("aravan", "Aravan virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("australis", "Australian bat lyssavirus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("lagos", "Lagos bat virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("caucasicus", "West Caucasian bat virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("irkut", "Irkut virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("khujand", "Khujand virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("rabies", "Rabies virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("^Lyssavirus ", "", com2$Virus_ICTV)


com2$Virus_ICTV <- gsub("ebolavirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("orthobunyavirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("vesiculovirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("coltivirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("orthohantavirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("henipavirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("marburgvirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("phlebovirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("mammarenavirus", "virus", com2$Virus_ICTV)
com2$Virus_ICTV <- gsub("pararubulavirus", "virus", com2$Virus_ICTV)

com2$Virus_ICTV <- gsub(" disease$", "", com2$Virus_ICTV)

# remove " virus" from end to save space other than SARS coronavirus
com2$Virus_ICTV <- gsub(" virus$", "" , com2$Virus_ICTV)



# switch to ICTV virus names for plotting (not latin binomial names)
vtree$tip.label[grep("Avian ortho",vtree$tip.label)] <- "Newcastle disease virus"
vtree$tip.label[grep("Middle East",vtree$tip.label)] <- "MERS coronavirus"
vtree$tip.label[grep("Severe acute",vtree$tip.label)] <- "SARS coronavirus"
vtree$tip.label <- gsub("Saint ", "St. ", vtree$tip.label)
vtree$tip.label <- gsub("Venezuelan equine encephalitis", "VEE", vtree$tip.label)
vtree$tip.label <- gsub("Western equine encephalitis", "WEE", vtree$tip.label)
vtree$tip.label <- gsub("Eastern equine encephalitis", "EEE", vtree$tip.label)
vtree$tip.label <- gsub("aravan", "Aravan virus", vtree$tip.label)
vtree$tip.label <- gsub("australis", "Australian bat lyssavirus", vtree$tip.label)
vtree$tip.label <- gsub("lagos", "Lagos bat virus", vtree$tip.label)
vtree$tip.label <- gsub("caucasicus", "West Caucasian bat virus", vtree$tip.label)
vtree$tip.label <- gsub("irkut", "Irkut virus", vtree$tip.label)
vtree$tip.label <- gsub("khujand", "Khujand virus", vtree$tip.label)
vtree$tip.label <- gsub("rabies", "Rabies virus", vtree$tip.label)
vtree$tip.label <- gsub("^Lyssavirus ", "", vtree$tip.label)

vtree$tip.label <- gsub("ebolavirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("orthobunyavirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("vesiculovirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("coltivirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("orthohantavirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("henipavirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("marburgvirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("phlebovirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("mammarenavirus", "virus", vtree$tip.label)
vtree$tip.label <- gsub("pararubulavirus", "virus", vtree$tip.label)

vtree$tip.label <- gsub(" disease$", "", vtree$tip.label)

# remove " virus" from end to save space other than SARS coronavirus
vtree$tip.label <- gsub(" virus$", "" , vtree$tip.label)


# remove viruses not in data
vtree <- drop.tip(vtree, setdiff(vtree$tip.label, com2$Virus_ICTV))

# remove hosts not in data
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, com2$Host_Upham))

phylo_v <- vtree$tip.label
phylo_h <- host_tree$tip.label

com2 <- com2 %>% mutate(
  host = factor(Host_Upham, levels = phylo_h),
  virus = factor(Virus_ICTV, levels = phylo_v))


plot_net <- ggplot(com2, aes(y = host, x = virus, fill = max_severity)) +
  geom_tile() +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(option="B", name="Maximum severity", 
    					na.value = "grey", begin=0.15,
    					end = 0.92, direction=-1,
    					breaks=c(1,5,20,75,250)) +
  labs(fill = "Maximum severity") +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1, size = 11, vjust=0.5),
    # axis.text.y = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")) + geom_hline(yintercept=32.5, linetype=2, color="gray", linewidth=0.5)
 plot_net


# colour host clades
MRCA(host_tree, "Sigmodon hispidus", "Glis glis")#118
MRCA(host_tree, "Antrozous pallidus", "Cynopterus sphinx")#87

tree2 <- groupClade(host_tree, c(87, 118))

ggtree_host <- ggtree(tree2, aes(color=group), ladderize=FALSE) + 
  	scale_color_manual(values=c("black", colours_BR)) + 
  	theme(legend.position="none") + 
	geom_tiplab(offset=200, hjust=1, 
				align=F, as_ylab = F,
				geom = "text") +
   			
		   			 theme(plot.margin = unit(c(l=0,t=0,r=-10.5,b=0), "cm"))

ggtree_host
ggtree_host$data

h_net <- ggtree_host + plot_net + plot_layout(widths=c(1,3))
# h_net


ggtree_virus <- ggtree(vtree, ladderize=FALSE) + 
				geom_tiplab(hjust=1, offset=17, angle=90, align=F) + 
				coord_flip() + 
				xlim_expand(25, 1) +
				theme(plot.margin = unit(c(l=-0.20,t=0.30,r=0,b=0.07), "cm"))
				# xlim(0, 11)
ggtree_virus

legend_b <- cowplot::get_legend(
  plot_net + 
    theme(legend.direction = "horizontal") + 
    guides(colour = guide_legend(title.position = "top"))
    # guides(color = guide_legend(nrow = 2, title.position = "top")) +
)

h_net

# first align the top-row plot (h_net) with the right-most plot of the
# bottom row (ggtree_virus)
plots <- cowplot::align_plots(h_net + theme(legend.position="none"), ggtree_virus +
				xlim(90, 1), 
						align = 'v', axis = 'r')

# then build the bottom row
bottom_row <- cowplot::plot_grid(legend_b, ggtree_virus, axis = 'r', rel_widths=c(1,3))

# then combine with the top row for final plot
a <- cowplot::plot_grid(plots[[1]], bottom_row, ncol = 1, rel_heights = c(4,0.85))
a
# ggsave("../plots_tables/heatmap_severity_pathology_Reported.pdf", a, width=12.5, height=16.5)
# ggsave("../plots_tables/heatmap_severity_pathology_Reported.png", a, width=12.5, height=16.5, bg = "white")



# adding side barplot for n_individuals
cust_breaks <- c(10, 100, 1000)

com3 <- com2 %>% group_by(host, Host_order) %>% summarise(n_individuals=sum(n_individuals))

bp <- 	ggplot(data=com3, aes(y=host, x=log10(n_individuals+0.1), fill=Host_order)) +
  		geom_bar(stat="identity", alpha=0.4)+
  		theme_bw() +  
  		# theme(panel_border(color = "grey85"))+
  		theme( 	axis.title.y=element_blank(),
        		axis.text.y=element_blank(),
         		axis.ticks.y=element_blank(),
         		axis.title.x=element_blank(),
        		# axis.text.x=element_blank(),
         		# axis.ticks.x=element_blank(),
         		legend.position="top",
     		    panel.grid.minor.x = element_blank(),
     		    # panel.grid.major.x = element_blank(),
	            panel.grid.minor.y = element_blank()
	            # panel.grid.major.y = element_blank()
         		) +
  		xlab("# individuals") + 
  		scale_x_continuous(expand = c(0, 0),
  			breaks = c(log10(cust_breaks)), limits=log10(c(1,2500)),
  			labels = c(cust_breaks)) +

  		scale_fill_manual(values=colours_BR) +
  		theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  		# scale_x_log10(breaks=cust_breaks, limits=c(0.1,2300), expand = c(0, 0)) + 
  		annotation_logticks(
					  base = 10,
					  sides = "b",
					  outside = FALSE,
					  scaled = TRUE,
					  short = unit(0.1, "cm"),
					  mid = unit(0.2, "cm"),
					  long = unit(0.3, "cm"),
					  colour = "grey",
					  linetype = 1,
					  alpha = 1,
					  color = NULL)
        # scale_y_log10(breaks = cust_breaks, minor_breaks = minor_breaks) +
        # annotation_logticks()
        # coord_equal() 

  		# scale_y_discrete(limits=rev) # to reverse order and match the tree
bp

# re-aligning
h_net_bp <- ggtree_host + plot_net + theme(legend.position="none") + bp + plot_layout(widths=c(1,3,0.4))
h_net_bp

# first align centre of the top-row plot (h_net) with the middle plot of the
# bottom row (ggtree_virus)
plots2 <- cowplot::align_plots(h_net_bp + theme(legend.position="none"), ggtree_virus +
				xlim(90, 1), 
						align = 'v', axis = 'r')


# then build the bottom row
bottom_row2 <- cowplot::plot_grid(legend_b, 
									ggtree_virus, #+ theme(plot.margin=margin(t=-30)), 
									NULL,  ncol=3, axis = 'l', rel_widths=c(1,3,0.4))
				

# then combine with the top row for final plot
b <- cowplot::plot_grid(plots2[[1]], bottom_row2 + theme(plot.margin = unit(c(-0.37,0.05,1.0,0.00), "cm")), ncol = 1, rel_heights = c(4,0.8))
b


bat <- magick::image_read("../plots_tables/desmodus.png") %>%
  # image_resize("570x380") %>%
  magick::image_colorize(100, colours_BR[1])

rat <- magick::image_read("../plots_tables/rattus.png") %>%
  # image_resize("570x380") %>%
  magick::image_colorize(100, colours_BR[2])

# b_imgs <- ggdraw() + draw_plot(b) + draw_image(bat, width=0.05, hjust = 0, vjust = 0.025) +
# 						  draw_image(rat, width=0.04, hjust = -0.19, vjust = -0.48,)

b_imgs <- ggdraw() + draw_plot(b) + draw_image(bat, width=0.05, hjust = 0.05, vjust = 0.24) +
						  draw_image(rat, width=0.04, hjust = -0.0, vjust = -0.20,)
b_imgs


ggsave("../plots_tables/heatmap_severity_Nindividuals.pdf", b_imgs, width=14, height=18.25)





## Average reported severity per host-virus pair

plot_net <- ggplot(com2, aes(y = host, x = virus, fill = mean_severity)) +
  geom_tile() +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(option="B", name="Mean severity", 
    					na.value = "grey", begin=0.15,
    					end = 0.92, direction=-1,
    					breaks=c(1,5,20,75,250)) +
  labs(fill = "Mean severity") +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1, size = 11, vjust=0.5),
    # axis.text.y = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")) + geom_hline(yintercept=32.5, linetype=2, color="gray", linewidth=0.5)
 plot_net

h_net <- ggtree_host + plot_net + plot_layout(widths=c(1,3))
# h_net

legend_b <- cowplot::get_legend(
  plot_net + 
    theme(legend.direction = "horizontal") + 
    guides(colour = guide_legend(title.position = "top"))
    # guides(color = guide_legend(nrow = 2, title.position = "top")) +
)

h_net

# first align the top-row plot (h_net) with the right-most plot of the
# bottom row (ggtree_virus)
plots <- cowplot::align_plots(h_net + theme(legend.position="none"), ggtree_virus +
				xlim(90, 1), 
						align = 'v', axis = 'r')

# then build the bottom row
bottom_row <- cowplot::plot_grid(legend_b, ggtree_virus, axis = 'r', rel_widths=c(1,3))

# then combine with the top row for final plot
a <- cowplot::plot_grid(plots[[1]], bottom_row, ncol = 1, rel_heights = c(4,0.85))
a


# re-aligning
h_net_bp <- ggtree_host + plot_net + theme(legend.position="none") + bp + plot_layout(widths=c(1,3,0.4))
h_net_bp

# first align centre of the top-row plot (h_net) with the middle plot of the
# bottom row (ggtree_virus)
plots2 <- cowplot::align_plots(h_net_bp + theme(legend.position="none"), ggtree_virus +
				xlim(90, 1), 
						align = 'v', axis = 'r')


# then build the bottom row
bottom_row2 <- cowplot::plot_grid(legend_b, 
									ggtree_virus, #+ theme(plot.margin=margin(t=-30)), 
									NULL,  ncol=3, axis = 'l', rel_widths=c(1,3,0.4))
				

# then combine with the top row for final plot
b <- cowplot::plot_grid(plots2[[1]], bottom_row2 + theme(plot.margin = unit(c(-0.37,0.05,1.0,0.00), "cm")), ncol = 1, rel_heights = c(4,0.8))
b


bat <- magick::image_read("../plots_tables/desmodus.png") %>%
  # image_resize("570x380") %>%
  magick::image_colorize(100, colours_BR[1]
)
rat <- magick::image_read("../plots_tables/rattus.png") %>%
  # image_resize("570x380") %>%
  magick::image_colorize(100, colours_BR[2])

# b_imgs <- ggdraw() + draw_plot(b) + draw_image(bat, width=0.05, hjust = 0, vjust = 0.025) +
# 						  draw_image(rat, width=0.04, hjust = -0.19, vjust = -0.48,)

b_imgs <- ggdraw() + draw_plot(b) + draw_image(bat, width=0.05, hjust = 0.05, vjust = 0.24) +
						  draw_image(rat, width=0.04, hjust = -0.0, vjust = -0.20,)
b_imgs


ggsave("../plots_tables/heatmap_mean_severity_Nindividuals.pdf", b_imgs, width=14, height=18.25)

