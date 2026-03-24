# plotting_functions.R
require(brms)
require(ggplot2)


# Code for scale and unscale with mean=0 sd=0.5 (following Gelman 2008)
scale_half <- function(x){
  return( (x-mean(x, na.rm=TRUE)) * 0.5/sd(x, na.rm=TRUE))
}

unscale_half <- function(x,original){
    x*sd(original, na.rm=TRUE)*2 + mean(original, na.rm=TRUE)  
}


cond_eff_plot <- function(model, ...){
    
    draws <- linpred_draws(
                  model,
                  newdata=model$data,
                  value = ".linpred",
                  ndraws = 4000,
                  seed = NULL,
                  re_formula = NULL,
                  category = ".category",
                  dpar = NULL) 
    
    ggplot(draws, aes(x=Host_order, y=.linpred, fill=Host_order)) + stat_eye(.width=c(0.95,0.8), alpha=0.7) + 
          theme_bw() + scale_fill_manual(values=colours_BR) + ylab("Conditional log-odds") + xlab("Host order") + theme(  legend.position="none")
}


# parameters table
param_tab <- function(model, ...){
  parameters::parameters(model, centrality="mean", ci=0.9, digits=2, diagnostic="Rhat", effects="all", test=NULL, verbose = FALSE) %>% insight::export_table(format = "markdown")
}


ar_theme <-   theme_bw() + 
              theme(legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.y=element_text(size=rel(1.1)), axis.title=element_text(size=9))

# intervals plot colors
color_scheme_set("teal")


# joint conditionals plot

joint_conditionals <- function(p1,p2,p3,global_ylim,padding){

      plot_grid(
      
        p1 + 
        global_ylim + 
        padding + 
        ylab("Conditional log-odds of Disease Presence") +
        xlab("") +
        theme_bw() + 
        theme(legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()),
      
        p2 +
        global_ylim + 
        padding + 
        ylab("Conditional log-odds of Mortality") +
        xlab("") +
        theme_bw() + 
        theme(legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()),
      
        p3 +
        global_ylim + 
        padding + 
        ylab("Conditional log-odds of Disease Severity") +
        xlab("") +
        theme_bw() + 
        theme(legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()),
      
                ncol=3, scale=0.90, 
                labels = c('A)', 'B)', "C)"), 
                label_size = 12, hjust=0.02,
                # labels = c('Disease Presence', 'Mortality', "Severity"), 
                align="hv"
      
               )
}


# joint conditionals + intervals combined plot

joint_conditionals_intervals <- function(p1,p2,p3,global_ylim,padding,ar1,ar2,ar3){

      egg::ggarrange(
        (p1 + 
        global_ylim + 
        padding + 
        ggtitle("A) Disease presence")+
        ylab("Conditional log-odds") +
        xlab("") +
        theme_bw() + 
        theme(legend.position = "none",
              axis.title.y=element_text(size=12),
              axis.text.x=element_text(size=12),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())),
        
        (p2 +
        global_ylim + 
        padding + 
        ggtitle("B) Mortality")+
        ylab("") +
        xlab("") +
        theme_bw() + 
        theme(legend.position = "none",
              axis.text.x=element_text(size=12),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())),
      
        (p3 +
        global_ylim + 
        padding + 
        ggtitle("C) Severity")+
        ylab("") +
        xlab("") +
        theme_bw() + 
        theme(legend.position = "none",
              axis.text.x=element_text(size=12),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())),
        
        (ar1 + padding + theme(axis.text.y=element_text(size=12))), 
        (ar2 + padding + theme(axis.text.y=element_blank())), 
        (ar3 + padding + theme(axis.text.y=element_blank())),
         ncol=3, padding=0.5)

}

# Intervals plot with relabelling

labels <- c("b_Host_orderChiroptera" = "Host order (Chiroptera)", 
            "b_n_individuals" = "Number of individuals",
            "b_ses_pd" = "Host phylogenetic diversity", 
            "b_EvoIso" = "Host evolutionary isolation",
            "b_Dose_mass" = "Dose/g",
            "b_Dose_amount" = "Dose",
            "b_Reservoir_match" = "Reservoir match",
            "b_adult_mass_g" = "Body mass",
            "b_Host_orderChiroptera:Reservoir_match" = "Host order x Reservoir Match")


intervals_plot <- function(model, ...) {
  ip <- mcmc_intervals(as.matrix(model), regex_pars = "b_[^I]", prob=0.8, prob_outer=0.95, point_est="mean", point_size = 3.5, outer_size=1, inner_size=3) +  ar_theme + xlab("Estimated effect size") 

  ip <- ip + scale_y_discrete(labels = labels, limits=rev(ip$data$parameter))

  return(ip)
}


get_intervals <- function(draws) {
  bayesplot::mcmc_intervals_data(draws, prob=0.8, prob_outer=0.95, regex_pars = "b_[^I]")
}

compare_posteriors <- function(..., dodge_width = 0.7) {
  dots <- rlang::dots_list(..., .named = TRUE)
  draws <- lapply(dots, function(x) {
    if (class(x)[1] == "stanreg") {
        posterior::subset_draws(posterior::as_draws(x$stanfit),
            variable = names(fixef(x))
        )
    } else if (class(x)[1] == "brmsfit") {
        brm_draws <- posterior::subset_draws(posterior::as_draws(x$fit),
            variable = paste0("b_", rownames(fixef(x)))
        )
    } else {
        stop(paste0(class(x)[1], " objects not supported."))
    }
  })
  intervals <- lapply(draws, get_intervals)
  combined <- dplyr::bind_rows(intervals, .id = "model")

  # Generalized version: count predictors per model dynamically
  n_predictors <- sapply(intervals, nrow)
  model_names <- names(dots)
  combined$model <- rep(model_names, times = n_predictors)
  combined$model <- factor(combined$model, levels = rev(model_names))
  
  ggplot(combined, aes(x = m, y = parameter, color = model, group = model)) +
    geom_vline(xintercept = 0, color="lightgrey") +
    geom_linerange(aes(xmin = l, xmax = h), size = 3, position = position_dodge(dodge_width)) +
    geom_linerange(aes(xmin = ll, xmax = hh), position = position_dodge(dodge_width)) +
    
    geom_point(aes(color = model), fill="white", shape=21, size=3.5, position = position_dodge(dodge_width)) +
    
    geom_point(aes(color = model, fill=model), shape=21, size=3.5, position = position_dodge(dodge_width)) +
    
    scale_y_discrete(labels = labels, limits=rev(unique(as.character(combined$parameter)))) + 
    
    scale_color_manual(values=c("#e08214","#8073ac"), guide = guide_legend(reverse = TRUE)) +
    
    scale_fill_manual(values=alpha(c("#e08214","#8073ac"),0.5),guide = guide_legend(reverse = TRUE))+
    
    ar_theme + 
    
    labs(y="", x="Estimated effect size", col="", fill="")

}


# six panel forest plot

six_panel_forest_plot <- function(s_dis,i_dis,s_mort,i_mort,s_sev,i_sev){

  padding <- theme(plot.margin = unit(c(t=0.1, r=0.1, b=0, l=0), "cm"))
    
  egg::ggarrange(        
      (compare_posteriors(s_dis, i_dis) + padding + theme(axis.text.y=element_text(size=12)) + ggtitle("A) Disease presence")), 
      (compare_posteriors(s_mort, i_mort) + padding + theme(axis.text.y=element_blank()) + ggtitle("B) Mortality")), 
      (compare_posteriors("Species-level"=s_sev, "Individual-level"=i_sev) + padding + theme(axis.text.y=element_blank(), legend.position = "right") + ggtitle("C) Severity")),
      ncol=3, padding=0.5)
}


# taxa-level hierarchical effects plots

virus_hierarchical_effects <- function(re){

      # extract and format hierarchical effects
      re_virus <- re$Virus_ICTV+re$Virus_name
      re_virus <- as.data.frame(re_virus)
      re_virus$Virus <- row.names(re_virus)
      
      # clean up names
      names(re_virus) <- gsub("[.]Intercept","",names(re_virus))
      
      # plot next to virus taxonomic tree
      
      # virus taxonomy
      vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
      vtax <- as.data.frame(unclass(vtax), stringsAsFactors=TRUE)
      
      viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
      viruses <- dplyr::left_join(viruses, vtax)
      
      # make tree
      frm <- ~superkingdom/realm/kingdom/phylum/class/order/family/genus/Virus_ICTV
      vtree <- as.phylo(frm, data = vtax, collapse=FALSE)
      vtree$edge.length <- rep(1, nrow(vtree$edge))
      
      
      # remove viruses with no susceptibe bat or rodent hosts
      vtree <- drop.tip(vtree, setdiff(vtree$tip.label, re_virus$Virus))
      
      # reoder virus factors to match tree
      re_virus$Virus <- factor(re_virus$Virus, levels=rev(vtree$tip.label))

      # shorten names on taxonomy
      vtree$tip.label <- gsub(" virus$", "", vtree$tip.label)
      vtree$tip.label[grep("Avian ortho",vtree$tip.label)] <- "Newcastle Disease"
      vtree$tip.label[grep("Middle East",vtree$tip.label)] <- "MERS-related coronavirus"
      vtree$tip.label[grep("Severe acute",vtree$tip.label)] <- "SARS-related coronavirus"
      vtree$tip.label[grep("ebolavirus", vtree$tip.label, invert=T)] <- gsub(" [a-z]*virus$", "", vtree$tip.label[grep("ebolavirus", vtree$tip.label, invert=T)])
      vtree$tip.label <- gsub("Venezuelan equine encephalitis", "VEE", vtree$tip.label)
      vtree$tip.label <- gsub("Western equine encephalitis", "WEE", vtree$tip.label)
      vtree$tip.label <- gsub("Eastern equine encephalitis", "EEE", vtree$tip.label)
      vtree$tip.label <- gsub(" disease$", "", vtree$tip.label)

      v_re_plot <- ggplot(re_virus, aes(x=Estimate, y=Virus)) + 
                    geom_linerange(aes(xmin=Q5, xmax=Q95), color="gray65", linewidth=0.9) + 
                    geom_linerange(aes(xmin=Q25, xmax=Q75), color="gray45", linewidth=2.0) + 
                    geom_point(aes(x=Estimate), size=2) + xlab("Estimated hierarchical effect") + 
                    ylab(NULL) + theme_bw() + 
                    theme(axis.text.y=element_blank(), axis.title=element_text(size=9),
                        plot.margin = margin(t = 0,  # Top margi
                                             r = 0.1,  # Right mrg
                                             b = 0,  # Bottom magi
                                             l = 0,  # Left marg
                                             unit = "cm"))

    v <- ggtree(ape::rotateConstr(vtree, rev(vtree$tip.label)), ladderize=FALSE) + 
              geom_tiplab(offset=0.2, hjust=0, 
                    align=F, as_ylab = F,
                    geom = "text", size=3) + 
              xlim_tree(8) +
              xlim_expand(15, 1) +
              # vexpand(.04, -1) +                   
              theme( plot.margin = margin(t = 0.1,  # Top margin
                                          r = 0,  # Right margin
                                          b = 0,  # Bottom margin
                                          l = 0,  # Left margin
                                          unit = "cm"))

    plot_grid(v, NULL, v_re_plot, nrow=1, rel_widths = c(0.8, -0.05, 1), align="h", scale=1.0)

}


host_hierarchical_effects <- function(re){

    # extract and format hierarchical effects
    re_host <- re$Host_Upham+re$Host_name
    re_host <- as.data.frame(re_host)
    re_host$Host <- row.names(re_host)
    
    # clean up names
    names(re_host) <- gsub("[.]Intercept","",names(re_host))
    
    # plot next to virus taxonomic tree
    
    # remove hosts with no susceptibe bat or rodent hosts
    htree <- drop.tip(tree, setdiff(tree$tip.label, re_host$Host))
    
    
    # adding host order for coloring points
    re_host <- left_join(re_host, dat %>% ungroup() %>% select(Host_Upham, Host_order) %>% rename(Host = Host_Upham)%>% unique())
    
    # reoder host factors to match tree
    re_host$Host <- factor(re_host$Host, levels=rev(htree$tip.label))
    
    # remobe underscore in host names on tree
    htree$tip.label <- gsub("_"," ", htree$tip.label)
    
    h_re_plot <-  ggplot(re_host, aes(x=Estimate, y=Host, color=Host_order)) + 
                  geom_linerange(aes(xmin=Q5, xmax=Q95, color=Host_order), linewidth=0.9, alpha=0.25) + 
                  geom_linerange(aes(xmin=Q25, xmax=Q75, color=Host_order), linewidth=2.0, alpha=0.50) + 
                  geom_point(aes(x=Estimate), size=2) + 
                  scale_colour_manual(values=colours_BR) + xlab("Estimated hierarchical effect") + ylab(NULL) + theme_bw() + 
                  theme(legend.position = "none",
                        axis.text.y=element_blank(), axis.title=element_text(size=10),
                        plot.margin = margin(t = 0,  # Top margin
                                             r = 0.1,  # Right margin
                                             b = 0,  # Bottom margin
                                             l = 0,  # Left margin
                                             unit = "cm"))

    h <- ggtree(ape::rotateConstr(htree, rev(htree$tip.label)), ladderize=FALSE) + 
              geom_tiplab(offset=2, hjust=0, 
                    align=F, as_ylab = F,
                    geom = "text", size=3) + 
              xlim_tree(10) +
              xlim_expand(160, 1) +
              # vexpand(.04, -1) +                   
              theme( plot.margin = margin(t = 0.1,  # Top margin
                                          r = 0,  # Right margin
                                          b = 0,  # Bottom margin
                                          l = 0,  # Left margin
                                          unit = "cm"))


    plot_grid(h, NULL, h_re_plot, nrow=1, rel_widths = c(0.8, -0.07, 1), align="h", scale=1.0)

}


# general hierarchical effects plot

re_intercepts_plot <- function(re, level_order){

    # convert RE list of arrays into single dataframe
    convert_arrays_to_long <- function(list_of_arrays) {
      
      long_df <- purrr::map_dfr(list_of_arrays, .id = "group_name", .f = function(arr) {
        df <- as.data.frame.table(arr, responseName = "Value")
        return(df)
      })
      long_df$group_name <- as.factor(long_df$group_name)
      return(long_df)
     }
    
    long_df <- convert_arrays_to_long(re)
    
    # convert to wide for plot
    wide_df <- pivot_wider(long_df, names_from=Var2, values_from=Value)
        
    wide_df$group_name <- factor(wide_df$group_name, levels = rev(level_order))
    
    re_labels <- c(
                "Host_Upham" = "Host (phylogenetic)",
                "Host_name" = "Host (non-phylogenetic)", 
                "Virus_ICTV" = "Virus (taxonomic)",
                "Virus_name" = "Virus (non-taxonomic)",
                "PaperID" = "Article",
                "Route_type" = "Inoculation route",
                "Dose_unit" = "Dose unit"
                )
    
    wide_df$group_name <- plyr::revalue(wide_df$group_name, re_labels, warn_missing=FALSE)
    


   ggplot(wide_df, aes(x=Estimate, y=group_name)) + 
                    geom_point(aes(x=Estimate), size=3, alpha = 0.2) + xlab("Estimated hierarchical effect") + 
                    ylab(NULL) + theme_bw() + 
                    theme(axis.title=element_text(size=9),
                        plot.margin = margin(t = 0,  # Top margi
                                             r = 0.1,  # Right mrg
                                             b = 0,  # Bottom magi
                                             l = 0,  # Left marg
                                             unit = "cm"))

}


# conditional effects plots for host order X reservoir match interaction 
# shows 80% and 95% credible intervals

plot_conditional_effects <- function(model, y_label) {
  
  # Extract data for 80% and 95% CIs
  ci_95 <- conditional_effects(
    model, effects = "Host_order:Reservoir_match", int_conditions = list(Reservoir_match = c(0, 1)), prob = 0.95)[[1]]
  
  ci_80 <- conditional_effects(model, effects = "Host_order:Reservoir_match", int_conditions = list(Reservoir_match = c(0, 1)), prob = 0.80 )[[1]]
  
  # visual tweaks
  my_colors <- c("0" = "#355CE8", "1" = "#E8C135")
  pdodge <- position_dodge(width = 0.4)
  
  ggplot(ci_95, aes(x = Host_order, y = estimate__, 
                    color = factor(Reservoir_match), 
                    group = factor(Reservoir_match))) +

    # Outer 95% CI
    geom_linerange(aes(ymin = lower__, ymax = upper__), 
                   position = pdodge, 
                   linewidth = 0.8, 
                   alpha = 0.5) +
    # Inner 80% CI
    geom_linerange(data = ci_80, 
                   aes(ymin = lower__, ymax = upper__), 
                   position = pdodge, 
                   linewidth = 2,
                   alpha = 0.8) +
    # Mean Point Estimate
    geom_point(position = pdodge, 
               size = 3) +

    # Scales and Labels
    scale_color_manual(values = my_colors, 
                       name = "Reservoir Match",
                       labels = c("No Match (0)", "Match (1)")) +

    theme_bw() + labs(y = y_label, x = "Host Order") + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  }