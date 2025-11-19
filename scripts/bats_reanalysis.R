## Re-analysis of BATS CTD data used for main text Figure 2 and Figure S2

bats_casts<-data.table::fread('../data/all_bats_ctd.txt',
                              fill=TRUE,
                              col.names=c('cast','dec_year','lat','lon','pressure','depth',
                                          'temp','conduct','salinity','ox','beamc',
                                          'fluor','par'))

## Editing casts that had different file format to match other casts
weird_casts<-bats_casts[is.na(bats_casts$dec_year),1]
weird_casts_fixed<-matrix(NA,nrow=nrow(weird_casts),ncol=ncol(bats_casts))
for(i in 1:nrow(weird_casts)){
  all_values<-str_split(weird_casts[i,1],pattern=' ')[[1]]
  all_values<-as.numeric(all_values[-which(all_values=='')])
  weird_casts_fixed[i,]<-all_values
  print(i/nrow(weird_casts))
}
weird_casts_fixed<-data.frame(weird_casts_fixed)
colnames(weird_casts_fixed)<-c('cast','dec_year','lat','lon','pressure','depth',
                               'temp','conduct','salinity','ox','beamc',
                               'fluor','par')
bats_casts[is.na(bats_casts$dec_year),]<-weird_casts_fixed

## Doing CTD data quality checks

broken_up_casts<-split(bats_casts,bats_casts$cast)
busted_conductivity<-do.call(c,lapply(broken_up_casts,function(x) min(x$conduct)<0))
busted_temperature<-do.call(c,lapply(broken_up_casts,function(x) min(x$temp)<0))
busted_pressure<-do.call(c,lapply(broken_up_casts,function(x) min(x$pressure)<0))
has_surface<-do.call(c,lapply(broken_up_casts,function(x) min(x$pressure)>10))
n_scans<-do.call(c,lapply(broken_up_casts,function(x) nrow(x)<100))
all_clear_casts<-Reduce(`|`,list(busted_conductivity,
                                 busted_temperature,
                                 busted_pressure,
                                 has_surface,
                                 n_scans))
okay_casts<-do.call(rbind,broken_up_casts[!all_clear_casts])

## NOTE:: Doing this analysis leverages the R librarty reticulate,
## which allows for utilization of python packages in R
## The package we will be using for our CTD reanalysis is
## the python release of the GSW Gibbs Seawater Toolbox v 3.4.0
## you will need to set this up on your system however works
## best for you

use_condaenv('base')
gsw<-import('gsw')

## Now using GSW v3.4.0 to analyze CTD data
thermo_corrected<-okay_casts %>%
  group_by(cast) %>%
  mutate(rounded_depth=round(depth),
         corrected_depth=gsw$z_from_p(pressure,lat),
         sp=gsw$SP_from_C(conduct*10,temp,pressure),
         sa=gsw$SA_from_SP(sp,pressure,lon,lat),
         ct=gsw$CT_from_t(sa,temp,pressure),
         sigma_theta=gsw$sigma0(sa,ct),
         o2_sol=gsw$O2sol(sa,ct,pressure,lon,lat),
         o2_sat=100*ox/o2_sol,
         mld=corrected_depth[min(which(!(abs(sigma_theta-sigma_theta[min(which(pressure %in% c(10,11)))])<=0.125) & pressure>10))-1],
         mld=ifelse(is.infinite(mld),min(corrected_depth),mld),
         som=ifelse(mld<=-300,
                    NA,corrected_depth[which(ox==max(ox[(corrected_depth>=-300 & corrected_depth<mld)]) & corrected_depth>=-300 & corrected_depth<=mld)]),
         som_intensity=ifelse(is.na(som),
                              NA,
                              ox[corrected_depth==som]),
         dcm=ifelse(max(fluor)==-999,
                    NA,
                    corrected_depth[between(corrected_depth,
                                            unique(mld)-200,
                                            unique(mld))][which.max(fluor[between(corrected_depth,
                                                                                  unique(mld)-200,
                                                                                  unique(mld))])]),
         som=ifelse(som_intensity==-999,NA,som),
         som_intensity=ifelse(som_intensity==-999,NA,som_intensity))

write.csv(thermo_corrected,'../intermediate_files/all_bats_thermo_corrected.csv')

## Using Winkler bottle data to remove CTD casts with poorly calibrated
## optode oxygen data
bottle_data<-data.table::fread('../virus_o2max_data/bats_bottle.txt',
                               skip=57,
                               col.names=c('cast','date_raw','dec_year',
                                           'time','lat','lon','depth','temp',
                                           'conduct','salinity','sigma_theta',
                                           'ox_winkler','oxfixtemp','ox_anom',
                                           'dco2','alkalinity','nitrate','nitrite',
                                           'phosphate','silica','poc','pon','toc',
                                           'tn','bact','pop','tdp','srp','bio_si',
                                           'litho_si','pro','syn','pes','nes')) %>%
  mutate(cast=as.character(cast))

bottle_and_ctd<-bottle_data %>%
  mutate(cast=as.character(cast),
         cast=gsub('[0-9][0-9]$','',cast),
         rounded_depth=round(depth/2)*2) %>%
  left_join(thermo_corrected,by=c('cast'='cast',
                                  'rounded_depth'='rounded_depth')) %>%
  drop_na(!c(dcm))

oxygen_cal<-lm(data=filter(bottle_and_ctd,ox>0,ox_winkler>0),formula=ox~ox_winkler)
summary(oxygen_cal)
standard_residuals<-oxygen_cal$residuals/sqrt(mean(oxygen_cal$residuals^2))
bad_cals<-filter(bottle_and_ctd,ox>0,ox_winkler>0)[which(abs(standard_residuals)>3),]
bad_cal_casts<-unique(bad_cals$cast)

## Constructing calendar of daily average oxygen profiles over
## BATS climatology:

daily_avg<-filter(thermo_corrected,between(ox,195,255),pressure<=300,!(cast %in% bad_cal_casts)) %>%
  mutate(cruise_date=as.POSIXct((dec_year-1988)*60*60*24*365.25,origin='1988-01-01')) %>%
  mutate(julian_day=yday(cruise_date),
         rounded_date=round(julian_day/12)*12) %>%
  group_by(rounded_date,pressure) %>%
  summarize(avg_o2=mean(ox),
            avg_mld=mean(mld*(-1)),
            depth=mean(rounded_depth)*-1) %>%
  group_by(rounded_date) %>%
  mutate(avg_mld=mean(avg_mld,na.rm=TRUE))

figure_2a<-ggplot(daily_avg %>% mutate(super_round_depth=-1*round(depth/5)*5) %>%
         group_by(super_round_depth,rounded_date) %>%
         summarize(avg_mld=mean(avg_mld),
                   avg_o2=mean(avg_o2),
                   avg_pressure=mean(pressure)))+
  metR::geom_contour_fill(aes(x=rounded_date,y=super_round_depth,z=avg_o2))+
  geom_smooth(aes(x=rounded_date,y=avg_mld),
              col='black',size=1.5)+
  scale_y_reverse(lim=c(250,10),expand=c(0,0),name='Depth [m]',
                  breaks=c(10,seq(50,250,by=50)))+
  scale_x_continuous(expand=c(0,0),name='Julian Day',
                     breaks=seq(0,350,by=50),
                     sec.axis=sec_axis(name='',
                                       breaks=seq(0.5,11.5),
                                       trans=~./30,
                                       labels=c('Jan',
                                                'Feb',
                                                'Mar',
                                                'Apr',
                                                'May',
                                                'Jun',
                                                'Jul',
                                                'Aug',
                                                'Sep',
                                                'Oct',
                                                'Nov',
                                                'Dec')))+
  scale_fill_cmocean(name='balance',limits=c(195,250),
                     breaks=c(200,220,240))+
  guides(fill=guide_colorbar(title=expression(paste('O'[2], ' ', mu, 'mol/kg'))))+
  theme_bw()+
  theme(text=element_text(size=18))+
  coord_fixed(ratio=0.5)

figure_2b<-ggplot(filter(thermo_corrected,
                         between(ox,195,255),
                         pressure<=300,
                         !(cast %in% bad_cal_casts)))+
  geom_point(aes(x=as.POSIXct((dec_year-1988)*60*60*24*365.25,origin='1988-01-01'),
                 y=corrected_depth*(-1),col=ox),
             size=1)+
  geom_line(data=(filter(thermo_corrected,
                         between(ox,195,255),pressure<=300,!(cast %in% bad_cal_casts),
                         is.na(som)==FALSE) %>% group_by(cast) %>% slice_head()),
            aes(x=as.POSIXct((dec_year-1988)*60*60*24*365.25,origin='1988-01-01'),
                y=-1*mld),col='black',size=1.5)+
  scale_y_reverse(lim=c(300,0))+
  scale_x_datetime(date_breaks='4 years',
                   date_labels='%Y')+
  scale_color_cmocean(name='balance',limits=c(195,255))+
  facet_wrap(~month(as.POSIXct((dec_year-1988)*60*60*24*365.25,origin='1988-01-01')),
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
  guides(color=guide_colorbar(title=expression(paste('O'[2], ' ', mu, 'mol/kg'))))+
  theme_bw()+
  ylab('Depth [m]')+
  xlab('')+
  theme(strip.background=element_blank(),
        text=element_text(size=18,face='bold'),
        axis.text.x=element_text(angle=90,face='plain'))

ggsave('../figures/Figure_2.pdf',figure_2a/figure_2b,scale=1.5)

## Calculating annual seasonal cycle in SOM emergence

good_o2_casts<-thermo_corrected %>% filter(!(cast %in% bad_cal_casts))

o2_modeling<-good_o2_casts %>%
  group_by(cast) %>%
  dplyr::filter(between(-1*pressure,mean(mld)-10,mean(mld)),
                ox>0,
                between(o2_sat,0,120)) %>%
  summarize(mean_o2s=mean(o2_sat),
            mld=mean(mld),
            dec_year=mean(dec_year),
            mean=MASS::fitdistr(o2_sat,'normal')$estimate[1],
            sd=MASS::fitdistr(o2_sat,'normal')$estimate[2],
            month=factor(ifelse(round(12*(dec_year-floor(dec_year)))>0,
                                round(12*(dec_year-floor(dec_year))),12), 
                         levels=1:12))

o2_cycle<-nls(mean_o2s~a+b*sin(2*pi*(d+(dec_year-floor(dec_year)))),
              data=o2_modeling)

summary(o2_cycle)

cycle_model<-ggplot(o2_modeling)+
  geom_point(aes(x=12*(dec_year-floor(dec_year)),
                 y=mean_o2s,
                 col=-1*mld))+
  geom_line(aes(x=12*(dec_year-floor(dec_year)),
                y=predict(o2_cycle,dec_year-floor(dec_year))),
            col='red',size=1.5)+
  scale_color_viridis_c(name='ML Depth [m]',direction=-1,
                        trans='reverse')+
  xlab('Month of Year')+
  ylab(expression(paste('Mean ','O'[2], ' Saturation in 10m Below Mixed Layer')))+
  theme_bw()
ggsave('../figures/Figure_S2.pdf',cycle_model)

