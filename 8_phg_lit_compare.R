#===========================================================================================================#
# 8_phg_lit_compare.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | January 23, 2021
# Background: Compare measured PHg (µg/g) from this study w/ observations in global rivers
#===========================================================================================================#

  # DATA PROCESSING & ANALYSIS ####

    # Data ####

      # Read
        hgdf <- read.csv(paste0(dir, "lit_comparison.csv"), header=T)
        hgdf <- droplevels(subset(hgdf, hgdf$region!="Africa"))
        hgdf$region <- factor(hgdf$region, levels=c("Asia","Europe","America"))
        hgdf$ref <- factor(hgdf$ref, levels=c("other","this_study"))
        
      # Subset data
        this <- droplevels(subset(hgdf, hgdf$ref=="this_study", select=c("region","ocean","river","ref_num","phg_ugg")))
        other <- droplevels(subset(hgdf, hgdf$ref=="other", select=c("region","ocean","river","ref_num","phg_ugg")))
  
      # Summary stats
        # Overall range
          range(na.omit(this$phg_ugg)); hist(na.omit(this$phg_ugg))
          range(na.omit(other$phg_ugg)); hist(na.omit(other$phg_ugg))
        # mean & SE
          tapply(this$phg_ugg, this$region, FUN=mean, na.rm=T); tapply(this$phg_ugg, this$region, FUN=std.error, na.rm=T) # this study
        # n obs
          for(i in unique(this$region)){print(paste0(i," = ",length(na.omit(this$phg_ugg[this$region==i]))," obs"))} # this study
        # mean & SE
          tapply(other$phg_ugg, other$region, FUN=mean, na.rm=T); tapply(other$phg_ugg, other$region, FUN=std.error, na.rm=T) # other studies
        # n obs
          for(i in unique(other$region)){print(paste0(i," = ",length(na.omit(other$phg_ugg[other$region==i]))," obs"))} # other studies
        # mean & range
          tapply(hgdf$phg_ugg, hgdf$ref, FUN=mean, na.rm=T); tapply(hgdf$phg_ugg, hgdf$ref, FUN=range, na.rm=T); for(i in unique(hgdf$ref)){print(paste0(i," = ",length(na.omit(hgdf$phg_ugg[hgdf$ref==i]))," obs"))}
          
        
    # Plotting ####
      
      # Density distribution- export as 4"x5" landscape PDF
        # Parameters
          cols_region <- c("white", "blue", "red")
          cols_ref <- c("grey40","white")
          ltype <- c(2,1)
        # Plot
          ggplot(data=hgdf, aes(x=phg_ugg, fill=ref, linetype=ref)) +
            # Scales
              scale_x_continuous(limits=c(-0.2,5.5), breaks=seq(0,5,1)) +
            # Themes
              # theme_void() + # <<== TURN ON OR OFF FOR EMPTY PLOTTING SPACE
              theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
                #axis.title.x=element_blank() + # <<== TURN ON OR OFF FOR EMPTY PLOTTING SPACE
                #axis.ticks.x=element_blank() + # <<== TURN ON OR OFF FOR EMPTY PLOTTING SPACE
                #axis.text.x=element_blank() + # <<== TURN ON OR OFF FOR EMPTY PLOTTING SPACE
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                panel.background=element_rect(fill='white'),
                panel.border=element_rect(colour="black", fill=NA, size=1),
                plot.background=element_rect(colour="white", size=1),
                text=element_text(size=19),
                axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
                axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0)),
                axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10)),
                axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
                axis.title.x=element_text(margin=margin(t=0, r=0, b=0, l=0)),
                legend.position="none") +
            # Data
              scale_fill_manual(values=cols_ref) +
              scale_linetype_manual(values=ltype) +
              labs(y="density", x=PHg~(µg~g^-1)) +
              geom_density(data=hgdf, aes(x=phg_ugg, y=), alpha=0.8)
          