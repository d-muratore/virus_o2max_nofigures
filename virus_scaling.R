## Calculating seasonal emergence of 
## scaling relationship between Prochlorococcus and VLP
## abundances with SOM appearance for Figure S12
## Running bats_reanalysis.R prior to this script
## will have all correct files loaded
bats_virus<-data.table::fread('../data/bats_virus_data.tsv',
                              skip=1,
                              col.names=c('station',
                                          'cruise_ID',
                                          'date_in',
                                          'decyear',
                                          'lat_in',
                                          'lon_in',
                                          'depth',
                                          'depth_nom',
                                          'depth_mixed',
                                          'abund_vir',
                                          'abund_vir_sd'))

bottle_virus<-left_join(bottle_and_ctd,
                        bats_virus %>% mutate(decyear=as.numeric(decyear),
                                              depth_nom=as.numeric(depth_nom)),
                        by=c('dec_year.x'='decyear','rounded_depth'='depth_nom'))

## Removing casts with missing data and formatting for scaling analysis:
nice_bottle_virus<-bottle_virus %>% filter(rounded_depth<=150,
                                           ox>0,
                                           !(cast %in% bad_cal_casts),
                                           pro>0,
                                           as.numeric(abund_vir)>0) %>%
  mutate(month=month(as.POSIXct((dec_year.x-1988)*60*60*24*365.25,origin='1988-01-01')),
         log_vir=log10(as.numeric(abund_vir)*1e6),
         log_pro=log10(pro),
         lin_vir=10^log_vir,
         lin_pro=10^log_pro)

split_nice_virus<-split(nice_bottle_virus,
                        nice_bottle_virus$month)


## Estimating separate models for each month of year

virus_model2_law<-lapply(split_nice_virus,
                         function(x) lmodel2(log_vir~log_pro,
                                             data=x,
                                             range.x='relative',
                                             range.y='relative',
                                             nperm=1e4))

## Collecting model results
do.call(c,lapply(virus_model2_law,function(x) (x$regression.results[1,3]-x$regression.results[2,3])/x$regression.results[1,3]))


## Preparing model results for plotting
virus_model2_out<-do.call(rbind,lapply(virus_model2_law,
                                       function(x) data.frame(exp=x$regression.results[2,3],
                                                              exp_low=x$confidence.intervals[2,4],
                                                              exp_high=x$confidence.intervals[2,5],
                                                              model_signif=x$regression.results[2,5],
                                                              int=x$regression.results[2,2],
                                                              int_low=x$confidence.intervals[2,2],
                                                              int_high=x$confidence.intervals[2,3])))
virus_model2_out$month<-names(virus_model2_law)


exponents_model2<-ggplot(virus_model2_out)+
  geom_errorbar(data=virus_model2_out %>% 
                  filter(month %in% c('8','9','10','11')),
                aes(x=as.numeric(month),ymin=exp_low,
                    ymax=exp_high))+
  geom_point(data=virus_model2_out %>% 
               filter(month %in% c('8','9','10','11')),
             aes(x=as.numeric(month),
                 y=exp,
                 fill=log10(model_signif)),
             size=4,
             shape=21,
             col='black')+
  geom_hline(yintercept=0,linetype='dashed')+
  scale_fill_gradient2(name=expression(paste('Log'[10],' p-Value')),
                       midpoint=log10(0.05),low='yellow3',
                       high='navyblue',mid='white',
                       lim=c(-4,-2))+
  theme_bw()+
  xlab('Month of Year')+
  scale_x_continuous(breaks=seq(2,12))+
  ylab(expression(atop(paste('Power-Law Exponent ',beta),
                       paste('Viruses=',italic(Prochlorococcus)^beta))))+
  theme(text=element_text(size=16),
        legend.position='bottom',
        legend.text=element_text(size=10))

scalings<-ggplot(bottle_virus %>% 
                    filter(rounded_depth<=150,ox>0,!(cast %in% bad_cal_casts)) %>%
                    mutate(month=month(as.POSIXct((dec_year.x-1988)*60*60*24*365.25,origin='1988-01-01'))))+
  geom_point(aes(x=log10(pro),
                 y=log10(as.numeric(abund_vir)*1e6),
                 col=o2_sat),
             size=2.5)+
  facet_wrap(~month,
             labeller=as_labeller(c('1'='Jan',
                                    '2'='Feb',
                                    '3'='Mar',
                                    '4'='Apr',
                                    '5'='May',
                                    '6'='Jun',
                                    '7'='Jul',
                                    '8'='Aug',
                                    '9'='Sep',
                                    '10'='Oct',
                                    '11'='Nov',
                                    '12'='Dec')))+
  geom_abline(data=virus_model2_out %>% filter(month %in% c('8','9','10','11')) %>%
                mutate(month=month(as.Date(paste0('1988/',month,'/01')))),
              aes(intercept=int,slope=exp,group=month))+
  scale_color_viridis_c(option='B',name='Oxygen Saturation (%)',lim=c(75,120))+
  xlab(expression(paste('Log'[10],' ',italic(Prochlorococcus),'/mL')))+
  scale_x_continuous(limits=c(3,6))+
  scale_y_continuous(limits=c(5,7.5))+
  ylab(expression(paste('Log'[10],' Viruses/mL')))+
  theme_bw()+
  theme(strip.text=element_text(size=14),
        strip.background=element_blank(),
        legend.position='bottom')

ggsave('../figures/Figure_S12.pdf',scalings|exponents_model2,
       scale=1.25)
