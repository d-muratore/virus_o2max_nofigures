###############Figure 3###############
  ##Data
  f3_heatmap_data <- read.csv("../virus_o2max_data/figure-3-rpob.csv", row.names = 1)
  f3_barplot_data <- read.csv("../virus_o2max_data/figure-3-sum-rpob.csv")
  
  ##Main heatmap
   f3a<-pheatmap(f3_heatmap_data[c(2:9)],
               cluster_rows = T, 
               cluster_cols = T,
               show_rownames = T, 
               annotation_row = f3_heatmap_data[c(1)],
               annotation_colors=list(taxa_group_annotation=c(Bacteria='#66c2a5',
                                      Eukaryota='#fc8d62',
                                      Archaea='#8da0cb')),
               scale = "row",
               fontsize = 10, 
               show_colnames = T, 
               hclust= "correlation",
               color = magma(10),
               border_color = "black",
               cellwidth = 15, cellheight =10)
    
    ##Barplot
     f3b<-f3_barplot_data %>%
        ggplot(aes(x=rpoB_vst_sum, y = reorder(order, heatmap_order))) +
        geom_col() +
        scale_x_continuous(expand = c(0, 1)) +
        theme_bw() + 
        theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6.7, colour = "black")) 
    
     
    ggsave('figures/Figure_3a.pdf',f3a,scale=2)
    ggsave('../figures/Figure_3b.pdf',f3b)
    
###############Figure 4###############
     ##Data
     f4_data <- read.csv("../data/figure-4-flow.csv")
     
     #boxplot
     f4<-f4_data %>%
       ggplot(aes(y = depth_ID, x  = cells_per_mL)) +
       geom_boxplot() +
       geom_point(pch = 21, fill = "black", position=position_jitter(width=0.1, height=0.1), size = 2.5, alpha = 0.8) +
       scale_y_discrete(limits = c("DCM", "SOM", "BML", "SRF")) +
       facet_wrap(~measurement, scales = "free", ncol = 2) +
       scale_fill_manual(values = c("800" = "white", "2000" = "grey")) +
       theme_Publication() +
       theme(axis.text.x = element_text( hjust = 1, size = 10)) 
     
     ggsave('../figures/Figure_4.pdf',f4)
     
###############Figure 5 a & c###############
    #Data
    f5ac_polony <- read.csv("../data/figure-5ac-polony.csv")
     
    #boxplot
    f5ac<-f5ac_polony  %>%
       ggplot(aes(y = depth_ID, x  = cyanophage_per_mL)) +
       geom_boxplot() +
       geom_point(pch = 21, fill = "black", position=position_jitter(width=0.1, height=0.1), size = 2.5, alpha = 0.8) +
       scale_y_discrete(limits = c("DCM", "SOM", "BML", "SRF")) +
       facet_wrap(~measurement, scales = "free", ncol = 2) +
       scale_fill_manual(values = c("800" = "white", "2000" = "grey")) +
       scale_x_continuous(limits = c(0, 800000)) +
       theme_Publication() +
       theme(axis.text.x = element_text( hjust = 1, size = 10)) 
     
###############Figure 5 b & d ###############
    #Data
     f5bd_ipolony <- read.csv("../data/figure-5bd-ipolony.csv")
     
     #boxplot
     f5bd<-f5bd_ipolony  %>%
      ggplot(aes(y = depth_ID, x = cells_per_mL )) +
       geom_boxplot() +
       geom_point( pch = 21, fill = "black", position=position_jitter(width=0.1, height=0.1), size = 2.5, alpha = 0.8) +
       scale_y_discrete(limits = c("DCM", "SOM", "BML", "SRF")) +
       facet_wrap(~measurement, scales = "free", ncol = 2) +
       scale_x_continuous(limits = c(0, 8000)) +
       theme_Publication() +
       theme(axis.text.x = element_text( hjust = 1, size = 10)) 
     
     ggsave('../figures/Figure_5.pdf',f5ac/f5bd)

  
###############Supplementary Figure 3###############
     #working on it
     
###############Supplementary Figure 4a###############
     #Data
     fs4a_data <- read.csv("../data/supp-figure-4a-Pro-rpob.csv")
     
     fs4<-fs4a_data %>%
       filter(phylogeny  == "Prochlorococcus" ) %>%
       ggplot(aes(fill=clade, x=vst, y=depth_ID)) + 
       geom_bar(position="fill", stat="identity") +
       scale_y_discrete(limits=c("DCM", "SOM", "BML","SRF" )) +
       scale_fill_brewer(palette = "Paired") +
       theme(legend.position = "none") +
       theme_bw() +
       ggtitle("Prochlorococcus SAGs")
     
     ggsave('../figures/Figure_S4.pdf',fs4)
  
###############Supplementary Figure 5###############
     #data
     fs5_data <- read.csv("../data/supp-figure-5-rpob-rpb1.csv")
     
     ##Summarize the mean vst by depth, during all daytime points
     depth_means <- fs5_data  %>%
       filter(day_night == "day") %>%
       group_by(query, depth_ID) %>%
       summarise_if(is.numeric, .funs = mean)
     
     ##convert to wide to make comparisons by depth
     depth_means <- depth_means[c(1,2,14)] %>% #select depth, query (a rpob or RPB1 hit), vst
       pivot_wider(values_from = vst, names_from = depth_ID)
     
     ###Retain a list of only rpob/rpb1 hits that have the highest mean vst at the SOM
     SOM_rpobs <- depth_means %>%
       filter(SOM > SRF & SOM > BML & SOM > DCM) %>%
       select(query) %>%
       inner_join(fs5_data , by = "query") #merge with all vst data for stats
     
     ##Statistically test whether each rpob/rpb1 hit is significantly higher at the SOM compared to all other depths 
     #First, Log() transform
     SOM_rpobs$log.vst <- log(SOM_rpobs$vst)
     
     #Kruskal wallace test of each rpob/rpb1 hit by depth
     SOM_rpobs_kw <- SOM_rpobs %>%
       filter(day_night == "day") %>% ##since we are only interested in daytime samples
       nest(data = -c(query)) %>%
       mutate(kruskal_raw = map(data, ~ kruskal.test(.x$log.vst, .x$depth_ID)),
              kruskal = map(kruskal_raw, broom::tidy)) %>%
       select(-data) %>%
       unnest(kruskal)
     
     ##Get only significant different rpbs, pvalue < 0.05
     SOM_rpobs_kw <- SOM_rpobs_kw %>%
       filter(p.value <= 0.05) %>% #699 rpoBs
       inner_join(fs5_data, by = "query")  #merge with all vst data for stats
     
     
     #now perform multiple testing using the Dunn's test
     ##Dunn's test, method = "bh" is the false discovery rate p-adj method aka Benjamini-Hochberg
     ##see https://community.rstudio.com/t/applying-dunn-test-using-purrr-map/15155/5 
     #need to add log() vst column again
     SOM_rpobs_kw$log.vst <- log(SOM_rpobs_kw$vst)
     
     SOM_rpobs_dunn <- SOM_rpobs_kw %>%
       filter(day_night == "day") %>%
       group_by(query) %>%
       nest() %>% 
       mutate(model = map(data, ~dunn.test(x = .x$log.vst, g = .x$depth_ID, method = "bh") %>% as_tibble())) %>% 
       unnest(model)
     
     ##Keep only those rpob/rpbs that are significantly (padj<=0.1) increased at the SOM relative to the SRF, BML, and DCM
     SOM_rpobs_dunn_df <- SOM_rpobs_dunn[c(1,6,7)] #query, padj and comparison column
     
     SOM_rpobs_final <- SOM_rpobs_dunn_df %>%
       filter(P.adjusted <= 0.1) %>%
       filter(!comparisons == "BML - DCM") %>%
       filter(!comparisons == "BML - SRF") %>%
       filter(!comparisons == "DCM - SRF") %>% #remove comparisons that don't include the SOM
       mutate(count = 1) %>% # add a count column 
       group_by(query) %>% 
       summarise_if(is.numeric, .funs = sum) %>%
       filter(count == 3) %>% ##only retain queries with significant comparisons between SOM - SRF, SOM - BML, and SOM - DCM (so, that would make the total number of occurances in the final table = 3).  #329 rpob/rpb1 hits are significantly increased ONLY at the SOM.
       inner_join(fs5_data , by = "query") #merge with all vsts for plot
     
     
     #plot total vst of each RPOB/RPB1 hit  as a barplot 
     ##Summarize by taxonomy first before plotting
     tax_key <- read.csv("../data/taxonomy_key.csv") #lineage key
     
     SOM_rpobs_final_sums <-  SOM_rpobs_final %>%
       group_by(family) %>%
       summarise_if(is.numeric, .funs = sum) %>%
       inner_join(tax_key, by = "family")
     
     fs5<-SOM_rpobs_final_sums %>%
       filter(!domain == "NA") %>%
       filter(vst > 200) %>% #retain only top-most abundant families
       ggplot(aes(y = vst, x = reorder(family,vst))) +
       coord_flip() +
       geom_bar(aes(fill = class), alpha = 0.7, stat = "identity") +
       facet_wrap( ~domain, scales="free") +
       scale_fill_brewer(palette = "Paired") +
       theme_bw() + 
       theme(legend.position = "bottom",panel.grid = element_blank(), axis.text.x = element_text(size = 6.7,angle = 45, colour = "black"), strip.background = element_blank(), strip.text.x = element_text(size = 20)) +
       ggtitle("Summed vst of sig. dif. rpoB/RPB1 candidates")
     
     ggsave('../figures/Figure_S5.pdf',fs5)
     

###############Supplementary Figure 6###############
     
  #Data
     kegg_curated_pathways <- read.csv("../data/kegg_pathway_assignments_final.csv") #manually curated KEGG pathways
     
     
     rpob_KO_ratios <- read.csv("../data/supp-figure-6-ko-rpoB-ratios.csv")
        #remove first column 
          rpob_KO_ratios$X <- NULL
        #The normalized transcript ("vst") for each gene annotated with a KO was summed by taxonomy, then divided by the rpoB vst value for the corresponding taxonomy to generate an "rpoB ratio". The ratios are represented by "org_KO", meaning "organism _ KEGG Orthology"
     
     #Find significantly increased or decreased KO's (organism, KO) during day time point (8 am)
       #Filter out org_KO's that have all rpob ratio's of 1 across all samples at 8 am 
       rpob_KO_ratios_8AM <- rpob_KO_ratios %>%
         filter(time == "800") %>%
         group_by(org_KO) %>%
         filter(!mean(ratio.rpob) == 1)
         
       ##Perform KW+Dunn multiple comparisons statistical test for each org_KO during the day (8 am) 
        rpob_KO_ratios_dunn <- rpob_KO_ratios_8AM  %>%
         group_by(org_KO) %>%
         nest() %>% 
         mutate(model = map(data, ~dunn.test(x = .x$ratio.rpob, g = .x$depth_ID, method = "bh") %>% as_tibble())) %>% 
         unnest(model)
       
     #Filter significant ones (padj<=0.1) and retain only comparisons that include the SOM
     rpob_KO_ratios_dunn_sig <- rpob_KO_ratios_dunn[c(1,6,7)] %>%
       filter(P.adjusted <= 0.1) %>%
       filter(!comparisons == "BML - DCM") %>%
       filter(!comparisons == "BML - SRF") %>%
       filter(!comparisons == "DCM - SRF") %>%
       mutate(count = 1)
     
     ##Filter to only retain KOs where all three comparisons (SRF - SOM, BML - SOM, and DCM - SOM) are significant
     #count occurance of sig dif comparisons for each org_KO
     rpob_KO_ratios_dunn_sig_count <-  rpob_KO_ratios_dunn_sig %>%
       group_by(org_KO) %>% summarise_if(is.numeric, .funs = sum) %>%
       filter(count == 3 ) #5,198 final KO's 
     
     #Filter the dunn results table to include only those KO's with all 3 sig comparisons
     rpob_KO_ratios_dunn_sig <- merge(rpob_KO_ratios_dunn_sig, rpob_KO_ratios_dunn_sig_count[c(1)])
     
     
     ##Take the average vst by depth for significant KOs to see whether it was UP or DOWN at the SOM, during the day
     rpob_KO_ratios_dunn_sig <- merge(rpob_KO_ratios_dunn_sig, rpob_KO_ratios, by = "org_KO")
     
     
     #filter to retain day time points, and summarise means by depth
     rpob_KO_ratios_dunn_sig <- rpob_KO_ratios_dunn_sig %>%
       filter(time == "800") %>%
       group_by(org_KO, depth_ID) %>%
       summarise_if(is.numeric, .funs = mean)
     
     rpob_KO_ratios_dunn_sig_means <-   rpob_KO_ratios_dunn_sig[c(1,2,5)] %>%
       pivot_wider(names_from = depth_ID, values_from = ratio.rpob)
     
     ##calculate the difference in means between the SOM vs SRF, and the SOM vs BML
     
     rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF <- (rpob_KO_ratios_dunn_sig_means$SOM - rpob_KO_ratios_dunn_sig_means$SRF)
     
     rpob_KO_ratios_dunn_sig_means$SOM_minus_BML <- (rpob_KO_ratios_dunn_sig_means$SOM - rpob_KO_ratios_dunn_sig_means$BML)
     
     rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM <- (rpob_KO_ratios_dunn_sig_means$SOM - rpob_KO_ratios_dunn_sig_means$DCM)
     
     #Postivie values = significantly increased at the SOM, negative values = significantly decreased at the SOM
     #replace the value with increased or decreased to compare between depth comparisons
     rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF <- replace(rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF, which(rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF > 0), "increased")
     rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF <- replace(rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF, which(rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF < 0), "decreased")
     
     rpob_KO_ratios_dunn_sig_means$SOM_minus_BML <- replace(rpob_KO_ratios_dunn_sig_means$SOM_minus_BML, which(rpob_KO_ratios_dunn_sig_means$SOM_minus_BML > 0), "increased")
     rpob_KO_ratios_dunn_sig_means$SOM_minus_BML <- replace(rpob_KO_ratios_dunn_sig_means$SOM_minus_BML, which(rpob_KO_ratios_dunn_sig_means$SOM_minus_BML < 0), "decreased")
     
     
     rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM <- replace(rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM, which(rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM > 0), "increased")
     rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM <- replace(rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM, which(rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM < 0), "decreased")
     
     ##determine consensus. remove org_KOs if the consensus between SOM v SRF and SOM v BML is not the same
     #create a column combining the results of both comparisons
     rpob_KO_ratios_dunn_sig_means$consensus <- paste( rpob_KO_ratios_dunn_sig_means$SOM_minus_SRF, "-",  rpob_KO_ratios_dunn_sig_means$SOM_minus_BML)
     
     rpob_KO_ratios_dunn_sig_means$consensus <- paste( rpob_KO_ratios_dunn_sig_means$consensus, "-",  rpob_KO_ratios_dunn_sig_means$SOM_minus_DCM)
     
     #keep only those that are "increased - increased - increased" or "decreased - decreased - decreased"
     rpob_KO_ratios_dunn_sig_means <- rpob_KO_ratios_dunn_sig_means %>%
       filter(consensus == "increased - increased - increased" | consensus == "decreased - decreased - decreased") #final set of KOs = 1,105
     
     #get KO annotations from the inital dataset. Also keep taxonomic lineage for summary. 
     
     ko_to_orgKo <- read.csv("../data/orgKO_to_lineage.csv")
     
     rpob_KO_ratios_dunn_sig_means <- rpob_KO_ratios_dunn_sig_means %>%
       left_join(ko_to_orgKo, by = "org_KO")
     
     #merge with kegg pathways
     rpob_KO_ratios_dunn_sig_means_withpathway <-  rpob_KO_ratios_dunn_sig_means %>%
       left_join(kegg_curated_pathways, by = "KEGG_KO") %>%
       mutate(count = 1) #add count column for summary
     
     
     ###Increased or decreased KO pathways at the SOM, Bacteria vs EUks
     rpob_KO_ratios_dunn_sig_pathway_count <-   rpob_KO_ratios_dunn_sig_means_withpathway %>%
       group_by(domain, consensus, consensus_pathway) %>%
       summarise_if(is.numeric, .funs = sum)
     
     #plot 
     fs6<-rpob_KO_ratios_dunn_sig_pathway_count %>%
       filter(!consensus_pathway == "NA") %>% #don't include annotations without curated pathways
       filter(!consensus_pathway == "[]") %>% #don't include ambiguous ko pathways
        # filter(count > 1) %>%
       ggplot(aes(x =count, y = reorder(consensus_pathway, count))) +
       geom_bar(aes(fill = consensus), color = "black", stat = "identity", position =  position_dodge2(preserve = "single")) +
       scale_x_continuous(expand = c(0, 1)) +
       theme_Publication() +
       theme(panel.grid.minor.y = element_line(colour = "black"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_blank(), axis.line.y.right = element_blank(), axis.line.x.bottom  = element_blank(), panel.border = element_blank()) +
       facet_wrap(~domain, scales =  "free") 
     
     ggsave('../figures/Figure_S6.pdf',fs6,scale=2)
###############Supplementary Figure 7a & b###############
     #Uses same data from Supp Figure 6
     
     #Add column : averaged vst:rpoB ratio (across depths)
     rpob_KO_ratios_dunn_sig_means_withpathway <- rpob_KO_ratios_dunn_sig_means_withpathway %>%
       mutate(ave_vst_rpob_ratio = (SRF + BML + SOM + DCM)/4)
     
     #for the figure, if the vst:rpob ratio was decreased at the SOM, make the value negative 
     rpob_KO_ratios_dunn_sig_means_withpathway <-  rpob_KO_ratios_dunn_sig_means_withpathway  %>%
       mutate(ave_vst_rpob_ratio_trend = case_when(consensus == "decreased - decreased - decreased" ~  ave_vst_rpob_ratio * -1,
                                                   consensus== "increased - increased - increased" ~  ave_vst_rpob_ratio * 1))
     
     #Supplementary Figure 7a
     fs7a<-rpob_KO_ratios_dunn_sig_means_withpathway %>%
       filter(org == "Prochlorococus") %>%
       ggplot(aes(x = ave_vst_rpob_ratio_trend, y = reorder(description, ave_vst_rpob_ratio_trend))) +
       geom_bar(aes(fill = consensus), col = "black", stat = "identity") +
       geom_vline(xintercept = 0) +
       scale_x_continuous(name = "Averaged ratio", limits = c(-0.5,2.5), breaks = seq(-0.5,2.5, by = 0.5)) +
       theme_Publication() +
       theme(panel.grid.major.x = element_blank(),  axis.title.y = element_blank(), axis.text.x = element_text(size = 10))
  
    #Supplementary Figure 7b
     #plot Prochlorococcus amt rpob ratio transcripts across all depths for day vs night
     fs7b<-rpob_KO_ratios %>%
       filter(org_KO == "Prochlorococus _ ko:K03320") %>% #Pro AMT KO
       filter(time == "800" | time == "2000") %>%
       ggplot(aes(y=depth_ID, x=ratio.rpob, fill = factor(time))) + 
       geom_boxplot() +
       geom_point(position =position_dodge(width=0.75), aes(group = factor(time)), size = 0.2) +
       scale_fill_manual(values = c("800" = "white", "2000" = "grey")) +
       scale_y_discrete(limits=c( "DCM", "SOM",  "BML", "SRF")) +
       theme_bw() + 
       theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 10)) +
       ggtitle(label = "amt, AMT, MEP, ammonium transporter, Amt family")
     
     ggsave('../figures/Figure_S7.pdf',fs7a/fs7b,scale=2)
###############Supplementary Figure 8a & b ###############
     #Data
     fs8_ipolonyperc <- read.csv("../data/supp-figure-8-ipolonyperc.csv")
     
     #boxplot
     fs8<-fs8_ipolonyperc  %>%
       ggplot(aes(y = depth_ID, x = count_or_percent )) +
       geom_boxplot() +
       geom_point( pch = 21, fill = "black", position=position_jitter(width=0.1, height=0.1), size = 2.5, alpha = 0.8) +
       scale_y_discrete(limits = c("DCM", "SOM", "BML", "SRF")) +
       facet_wrap(~measurement, scales = "free", ncol = 2) +
       #scale_x_continuous(limits = c(0, 8000)) +
       theme_Publication() +
       theme(axis.text.x = element_text( hjust = 1, size = 10)) 
     
     ggsave('../figures/Figure_S8.pdf',fs8)
  
###############Supplementary Figure 9###############
     #Data
     fs9_t4t7 <- read.csv("../data/supp-fig-9-data.csv")
     #scaffolds were assigned T4 or T7 taxonomy based on presence of hallmark proteins (see methods)
     
     ###Average vst's by depth for each vOTU scaffold 
       ###Calculate mean and sd
       data_summary <- function(data, varname, groupnames){
         require(plyr)
         summary_func <- function(x, col){
           c(mean = mean(x[[col]], na.rm=TRUE),
             sd = sd(x[[col]], na.rm=TRUE))
         }
         data_sum<-ddply(data, groupnames, .fun=summary_func,
                         varname)
         data_sum <- rename(data_sum, c("mean" = varname))
         return(data_sum)}
       
       fs9_t4t7_means <- data_summary(fs9_t4t7, varname="vst", groupnames=c("depth_ID" ,"vOTU_scaffold", "taxonomic_assignment")) 
     
     #heatmap format
     fs9_t4t7_means  <-  as.data.frame(fs9_t4t7_means[c(1:4)] %>%
       pivot_wider(names_from = depth_ID, values_from = vst))
     rownames( fs9_t4t7_means ) <-  fs9_t4t7_means$vOTU_scaffold
    
     #Heatmap
     colfunc<-colorRampPalette(c("white","lightgray", "black")) 
     
     fs9<-pheatmap(fs9_t4t7_means[c(3:6)],
              annotation_row = fs9_t4t7_means[c(2)],
              scale = "row",
              cluster_rows = T, 
              cluster_cols = T,
              fontsize = 6, 
              border_color = 'darkgray',
              show_rownames = F, 
              show_colnames = T,
              color = colfunc(50)) 

     ggsave('../figures/Figure_S9.pdf',fs9)     

     ## Supplemental Figure 3
      #metadata
     metadata <- read.table("../data/metadata.txt", sep= "\t", header = T)
     
     #data
     fs3_data <- read.table("../data/supp-figure-3-vsts.txt", sep = "\t", header = T)
     
     #PCA analysis
     pca <- prcomp((t(fs3_data)), center = T, scale = F)
     
     percentVar <- round(((pca$sdev) ^ 2 / sum((pca$sdev) ^ 2)* 100), 2) 
     
     pca.dat <- as.data.frame(pca$x)
     pca.dat$sample <- row.names(pca.dat)
     pca.dat<- merge(pca.dat, metadata, by = "sample")
     
    #Plot PCAs 
     makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")
     
     #panel a - color coded by depth pc1pc2
     ggplot(pca.dat, aes(PC1,PC2)) +
       geom_point(aes(size = 8,  color = depth_ID), alpha = 0.8) +
       xlab(makeLab(percentVar[1],1))+ 
       ylab(makeLab(percentVar[2],2))+
       scale_color_manual(values = c("SRF" = "#ac7db7", "BML" = "#dd7e1b", "DCM" = "#6ab244", "SOM" = "#d7549f")) +
       theme_bw()+
       theme(plot.title = element_text(hjust = 0.5),text = element_text(size=15), legend.text = element_text(size = 16, face = "bold"),
             legend.title = element_text(size = 16, colour = "black", face = "bold"), legend.position = "right") +
       guides(colour = guide_legend(override.aes = list(size=8))) 
     
     #panel b - color coded by depth pc2pc3
     ggplot( pca.dat, aes(PC2,PC3)) +
       geom_point(aes(size = 8,  color = depth_ID), alpha = 0.8) +
       xlab(makeLab(percentVar[2],2))+ 
       ylab(makeLab(percentVar[3],3))+
       scale_color_manual(values = c("SRF" = "#ac7db7", "BML" = "#dd7e1b", "DCM" = "#6ab244", "SOM" = "#d7549f")) +
       theme_bw()+
       theme(plot.title = element_text(hjust = 0.5),text = element_text(size=15), legend.text = element_text(size = 16, face = "bold"),
             legend.title = element_text(size = 16, colour = "black", face = "bold"), legend.position = "right") +
       guides(colour = guide_legend(override.aes = list(size=8))) 
     
     
     #panel c - color coded by day/night pc1pc2
     ggplot( pca.dat, aes(PC1,PC2)) +
       geom_point(aes(size = 8,  color = day_night), alpha = 0.8) +
       xlab(makeLab(percentVar[1],1))+ 
       ylab(makeLab(percentVar[2],2))+
       theme_bw()+
       theme(plot.title = element_text(hjust = 0.5),text = element_text(size=15), legend.text = element_text(size = 16, face = "bold"),
             legend.title = element_text(size = 16, colour = "black", face = "bold"), legend.position = "right") +
       guides(colour = guide_legend(override.aes = list(size=8))) 
     
     #panel d - color coded by depth pc2pc3
     ggplot( pca.dat, aes(PC2,PC3)) +
       geom_point(aes(size = 8,  color = day_night), alpha = 0.8) +
       xlab(makeLab(percentVar[2],2))+ 
       ylab(makeLab(percentVar[3],3))+
       theme_bw()+
       theme(plot.title = element_text(hjust = 0.5),text = element_text(size=15), legend.text = element_text(size = 16, face = "bold"),
             legend.title = element_text(size = 16, colour = "black", face = "bold"), legend.position = "right") +
       guides(colour = guide_legend(override.aes = list(size=8))) 
     
     
     #panel e- color coded by depth pc2pc3, facet by depthID
     ggplot(pca.dat, aes(PC2,PC3)) +
       geom_point(aes(size = 8,  color = day_night), alpha = 0.8) +
       xlab(makeLab(percentVar[2],2))+ 
       ylab(makeLab(percentVar[3],3))+
       theme_bw()+
       theme(plot.title = element_text(hjust = 0.5),text = element_text(size=15), legend.text = element_text(size = 16, face = "bold"),
             legend.title = element_text(size = 16, colour = "black", face = "bold"), legend.position = "right") +
       guides(colour = guide_legend(override.aes = list(size=8))) +
       facet_wrap(~depth_ID, nrow = 4)
     

     
