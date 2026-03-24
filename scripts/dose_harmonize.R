# dose_harmonize.R

require(dplyr);packageVersion("dplyr")#1.1.4
require(ggplot2);packageVersion("ggplot2")#3.4.4
require(cowplot);packageVersion("cowplot")#1.1.1

require(brms);packageVersion("brms")#2.21.6
require(bayesplot);packageVersion("bayesplot")#1.8.1
require(ape);packageVersion("ape")#5.5


# data
primary <- read.csv("../raw_data/primary_screen_Mar5_2024.csv")
dat <- read.csv("../raw_data/individual_data_Mar20_2024.csv")
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- dplyr::left_join(viruses, vtax)
hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)
#subset to only host reported and upham
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()

# construct column for individualID
dat <- dat %>%                              # Create numbering variable
		  group_by(PaperID) %>%
  		  mutate(IndividualID = paste(row_number(),unique(PaperID), sep="_"))


# add names
primary <- dplyr::left_join(primary, viruses)
primary <- dplyr::left_join(primary, hosts)
dat <- dplyr::left_join(dat, hosts)
dat <- dplyr::left_join(dat, viruses)

# remove NA viruses
primary <- primary[!is.na(primary$Virus_ICTV),]
dat <- dat[!is.na(dat$Virus_ICTV),]

# adding mean body size
# via COMBINE https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3344

combine <- read.csv("../raw_data/COMBINE/trait_data_imputed.csv")
combine$phylacine_binomial <- gsub(" ","_",combine$phylacine_binomial)

# intersect(dat$Host_Upham, combine$phylacine_binomial)#82
# setdiff(dat$Host_Upham, combine$phylacine_binomial)# Tamias_amoenus, Liomys_salvini, NA
combine$phylacine_binomial[grep("Neotamias_amoenus",combine$phylacine_binomial)] <- "Tamias_amoenus"
combine$phylacine_binomial[grep("Heteromys_salvini",combine$phylacine_binomial)] <- "Liomys_salvini"

combine$Host_Upham <- combine$phylacine_binomial
bodymass <- select(combine, Host_Upham, adult_mass_g) %>% unique()

# for duplicated host species take average of mass
bodymass <- bodymass %>% group_by(Host_Upham) %>% mutate(mean_mass=mean(adult_mass_g))
bodymass$adult_mass_g <- bodymass$mean_mass
bodymass <- select(bodymass, -mean_mass)
bodymass <- unique(bodymass)

dat <- left_join(dat, bodymass)


# cleaning dose amount and converting to numeric
dat$Dose_amount[dat$Dose_amount=="Unknown"] <- NA

dat$Dose_amount <- sapply(dat$Dose_amount, function(txt) eval(parse(text=txt)))

# harmonizing dose unit
sort(unique(dat$Dose_unit))
dat$Dose_unit[dat$Dose_unit%in%
				c("Adult mouse LD50","Median MICLD50","Median_MICLD50", "MICLD50",
				"MICLD50_newborns","MICLD50_suckling","MIPLD50","Mouse LD50")] <- "Mouse LD50"

dat$Dose_unit[dat$Dose_unit%in%
				c("Mean_TCID50","Median_TCID50","TCID50")] <- "TCID50"


sort(table(dat$Dose_unit), decreasing=TRUE)
# 1262 Mouse LD50
# 567 PFU
# 540 TCID50
# 116 FFU


# Harmonize Inoculation Route
dat <- dat %>% tidyr::separate(col=Inoculation_route, 
									into=c("Route_type","Route_location"), 
									sep=" [(]", extra="merge")

sort(unique(dat$Route_type))
sort(unique(dat$Route_location))

dat$Route_type[grep("Intranasal", dat$Route_location)] <- "Subcutaneous; Intranasal"
dat$Route_type[grep("Contact", dat$Route_location)] <- "Bite; Contact"

sort(unique(dat$Route_type))

# "Subcutaneous injections are administered in the fat layer, underneath the skin. 
# Intramuscular injections are delivered into the muscle. 
# Intradermal injections are delivered into the dermis, 
# or the skin layer underneath the epidermis (which is the upper skin layer)."

# rough collapsing
sort(table(dat$Route_type), decreasing=TRUE)

# Injection into body other than brain:
# Subcutaneous / Intramuscular / Intraperitoneal / Intradermal / Subdermal / Intravenous

dat$Route_type[dat$Route_type%in%
				c("Subcutaneous","Intramuscular","Intraperitoneal",
					"Intradermal","Subdermal","Intravenous",
					"Footpad","Subcutaneous; Intraperitoneal",
					"Intraperitoneal; Intravenous",
					"Intracardial; Intraperitoneal")] <- "Injection (body)"

# Injection into brain:
dat$Route_type[dat$Route_type%in%
				c("Intracerebral","Intracranial")] <- "Injection (brain)"


# Oral / nasal inoculation:
dat$Route_type[dat$Route_type%in%
				c("Oronasal","Oral","Intranasal", "Nasal")] <- "Oronasal"

# Simluated "Natural" infections
dat$Route_type[dat$Route_type%in%
				c("Transplacental","Contact","Skin","Vector",
					"Bite","Bite; Contact", "Aerosol","Natural")] <- "Natural"

# setting others (multiple routes or vague "injected") to NA
dat$Route_type[dat$Route_type%in%
				c("Injected - no further details","Subcutaneous; Intranasal")] <- NA

sort(table(dat$Route_type), decreasing=TRUE)


colors <- c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000")

dat_rt <- dat[!is.na(dat$Route_type),]

rt_stack <- ggplot(dat_rt, aes(fill=Route_type,  x=Host_order)) + 
		geom_bar(position="stack", stat="count") +
		# geom_point(position = position_jitter(seed = 1, width = 0.2), 
		# 	alpha=0.4) + 
		# scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "none",legend.title=element_blank()) + 
		xlab("") + ylab("Number of individuals")  

# rt_stack

rt_fill <- ggplot(dat_rt, aes(fill=Route_type, x=Host_order)) + 
		geom_bar(position="fill", stat="count") +
		# geom_point(position = position_jitter(seed = 1, width = 0.2), 
		# 	alpha=0.4) + 
		# scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "none",legend.title=element_blank()) + 
		xlab("") + ylab("Proportion of individuals")  

# rt_fill

combo_plot <- cowplot::plot_grid(rt_stack, rt_fill)

# extract legend from plot1
legend <- get_legend(
  rt_stack + theme(legend.position = "bottom")
)
  
# Combine combined plot and legend using plot_grid()
plot_grid(combo_plot, legend, ncol=1, rel_heights = c(1, .05))
ggsave("../plots_tables/Inoculation_routes.pdf", height=6, width=6)
ggsave("../plots_tables/Inoculation_routes.png", height=6, width=6)


# adjusting for body mass
# dat_sm$Dose_mass <- log10(dat_sm$Dose_amount) / log10(dat_sm$adult_mass_g)
# dat_sm$Dose_mass <- log10(dat_sm$Dose_amount) / (dat_sm$adult_mass_g)
dat$Dose_mass <- (dat$Dose_amount) / (dat$adult_mass_g)
hist(dat$Dose_mass)


# Graphs with only major dose units
dat_sm <- filter(dat, Dose_unit%in%c("Mouse LD50","PFU","TCID50","FFU"))

# keeping only sucessful infections  
dat_sm <- filter(dat_sm, Susceptible_YN%in%c("Y"))

dat_sm$Dose_unit <- factor(dat_sm$Dose_unit, levels = c("Mouse LD50", "PFU", "TCID50","FFU"))



colors <- c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000")

p <- ggplot(dat_sm, aes(Host_order, Dose_amount, fill=factor(Dose_unit), color=factor(Dose_unit))) + 
		geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, width = 0.2), 
			alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "none") + xlab("Host order") + ylab("Dose") + 
		facet_grid(cols=vars(Dose_unit)) + 
				stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.7,
               colour = "black")
p

ggsave("../plots_tables/Doses_by_type.pdf", p, width=8, height=5)
ggsave("../plots_tables/Doses_by_type.png", p, width=8, height=5)


# per body size
p2 <- ggplot(dat_sm, aes(Host_order, Dose_mass, fill=Dose_unit, color=Dose_unit)) + 
		geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, width = 0.2), 
			alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		# theme_bw() +
		theme(legend.position = "none") + xlab("Host order") + ylab("Dose/g") + 
		facet_grid(cols=vars(Dose_unit)) + 
				stat_summary(fun = "mean",   
            geom = "crossbar", 
               width = 0.7,
               colour = "black")
p2

ggsave("../plots_tables/Dose_per_g_by_type.pdf", p2, width=8, height=5)
ggsave("../plots_tables/Dose_per_g_by_type.png", p2, width=8, height=5)


# Dose per viral family

p3 <- ggplot(dat_sm, aes(Host_order, Dose_mass, fill=factor(Dose_unit), color=factor(Dose_unit))) + 
		geom_violin(alpha=0.4) +
		# geom_point(position = position_jitter(seed = 1, width = 0.1), 
		# 	alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "bottom", legend.title=element_blank()) + ylab("Dose/g") + xlab("") + 
		facet_grid(cols=vars(family)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
		# + 
		# 		stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              inherit.aes = TRUE)
p3

ggsave("../plots_tables/Dose_by_viral_fam.pdf", p3, height=4, width=16)
ggsave("../plots_tables/Dose_by_viral_fam.png", p3, height=4, width=16)


p4 <- ggplot(dat_sm, aes(Host_order, Dose_mass, fill=factor(Dose_unit), color=factor(Dose_unit))) + 
		geom_violin(alpha=0.4) +
		# geom_point(position = position_jitter(seed = 1, width = 0.1), 
		# 	alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "bottom", legend.title=element_blank()) + ylab("Dose/g") + xlab("") + 
		facet_grid(cols=vars(order)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
		# + 
		# 		stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              inherit.aes = TRUE)
p4

ggsave("../plots_tables/Dose_by_viral_order.pdf", p4, height=4, width=12)
ggsave("../plots_tables/Dose_by_viral_order.png", p4, height=4, width=12)

dose_dat <- select(dat, PaperID, IndividualID, Virus_ICTV, Host_Upham,
					adult_mass_g, Dose_amount, Dose_unit, Dose_mass, Route_type) %>% unique()


dim(dose_dat)#2885

write.csv(dose_dat, "../clean_data/dose_data.csv", row.names=FALSE)


# modelling dose differences between bats and rodents

bayesplot_theme_set(theme_bw())

dat_dose <- dose_dat

# dat_dose <- left_join(dat, dose_dat)
# dat_dose <- dat_dose[dat_dose$Susceptible_YN==1,]
# dat_dose <- dat_dose[!is.na(dat_dose$Susceptible_YN),]

# # log transform dose/g
# dat_dose <- dat_dose %>% group_by(Dose_unit) %>% mutate(Dose_mass = log(Dose_mass))

dat_dose$Host_name <- dat_dose$Host_Upham
dat_dose$Virus_name <- dat_dose$Virus_ICTV


# host data
tree <- read.nexus("../raw_data/upham_tree_666.nex")
hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)

# virus names and taxonomy
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- left_join(viruses, vtax)

#subset to only host reported and upham
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()

# make subset tree to sampled hosts
host_tree <- drop.tip(tree, setdiff(tree$tip.label, dat$Host_Upham))

# make virus tree
vtax <- as.data.frame(unclass(vtax), stringsAsFactors=TRUE)
frm <- ~superkingdom/realm/kingdom/phylum/class/order/family/genus/Virus_ICTV
vtree <- as.phylo(frm, data = vtax, collapse=FALSE)
vtree$edge.length <- rep(1, nrow(vtree$edge))
http://www.banningeyre.com/
# include only viruses in subset data (e.g. after removing mole)
vtree <- drop.tip(vtree, setdiff(vtree$tip.label, dat$Virus_ICTV))
# plot(vtree)

# phylogenetic correlation structure
host_cov <- vcv(host_tree, corr=TRUE)
virus_cov <- vcv(vtree, corr=TRUE)


# adding host_order

dat_dose <- left_join(dat_dose, select(dat, c(Host_Upham, Host_order), ))

names(dat_dose)

hist(log(dat_dose$Dose_mass))

# if (!file.exists("../fit_models/dose_model.rds")) {

  dose_m1 <- brm(log(Dose_mass) ~ Host_order +
  					  (1|Dose_unit) +  (1|Route_type), # + 
                      # (1|Host_name) + (1|Virus_name) +  
                      # (1|Host_Upham) + (1|Virus_ICTV), 
    data=dat_dose, family=gaussian(), 
    iter=4000, thin=1, 
    # cov_ranef = list(Host_Upham = host_cov, Virus_ICTV=virus_cov),
    control=list(adapt_delta=0.95, max_treedepth=10), cores=4, 
    save_pars=save_pars(all=TRUE))

  saveRDS(dose_m1, "../fit_models/dose_model.rds")

# } else { dose_m1 <- readRDS("../fit_models/dose_model.rds")}


summary(dose_m1)







# merging with severity data to plot severity by dose

dat_symptoms <- read.csv("../clean_data/symptoms_severity.csv")

dat <- left_join(dat_symptoms, dat_sm)
dat <- dat[!is.na(dat$Dose_unit),]

colours_BR <- c("#1B85BF", "#AB1808")

p5 <- ggplot(dat, aes(y=severity_rank, x=Dose_mass, fill=Host_order, color=Host_order, shape=Host_order)) + 
		# geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, height = 0.15), 
			alpha=0.4) + 
		scale_x_continuous(trans='log10') +
		scale_fill_manual(values=colours_BR) + scale_color_manual(values=colours_BR) +
		# theme(legend.position = "none") + 
	    labs(color="  Host order", shape="  Host order", fill="  Host order") + ylab("Severity") + xlab("Dose/g") + 
		facet_grid(cols=vars(Dose_unit)) + theme(legend.key=element_blank()) 
				# stat_summary(fun = "mean",
    #            geom = "crossbar", 
    #            width = 0.7,
    #            colour = "black")
p5

ggsave("../plots_tables/severity_by_dose.pdf", p5, width=10, height=4)
ggsave("../plots_tables/severity_by_dose.png", p5, width=10, height=4)
