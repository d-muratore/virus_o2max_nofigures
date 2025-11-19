### This script will do the necessary processing of AE1926
### cruise CTD data in order to generate the panels in main text
## figure 1 and supplemental figure 1

## Read in data
ctd_data<-read.csv('../virus_o2max_data/AE1926_ctd.csv') %>%
  dplyr::select(-X) 

## Prepare for decibar binning and smoothing
ctd_smooth_ready<-ctd_data %>%
  mutate(rounded_depth=round(z/0.5)*0.5) %>%
  group_by(Cast_ID,rounded_depth,ISO_DateTime_UTC) %>%
  dplyr::summarize(across(everything(),
                   list(median),
                   .names='binned_{.col}')) %>%
  mutate(time=as.POSIXct(gsub('Z','',ISO_DateTime_UTC),
                         format='\t%Y-%m-%dT%H:%M:%S',
                         tz='UTC'),
         cast=as.numeric(gsub('^.*InVirTC','',Cast_ID))) %>%
  ungroup() %>%
  group_by(cast) %>%
  mutate(lat=mean(binned_Latitude),
         lon=mean(binned_Longitude)) %>%
  dplyr::select(-c(binned_Latitude,binned_Longitude)) %>%
  mutate(time_local=with_tz(time,tzone='Etc/GMT+4')) %>%
  dplyr::select(c(cast,time,time_local,rounded_depth,starts_with('binned')))

variable_columns<-grep('binned',colnames(ctd_smooth_ready))
ctd_smoothed<-split(ctd_smooth_ready,ctd_smooth_ready$cast)

## Perform KS smoothing with a 5 db bandwidth
sample_smoothing<-function(wt,variable_columns){
  wt[,variable_columns]<-apply(wt[,variable_columns],
                               2,
                               function(x) ksmooth(wt$rounded_depth,
                                                   x,
                                                   bandwidth=5)$y)
  return(wt)
}
ctd_smoothed<-do.call(rbind,lapply(ctd_smoothed,function(x) sample_smoothing(wt=x,variable_columns=variable_columns)))

## Estimate mixed layer depth for each cast
ctd_smoothed<-ctd_smoothed %>%
  group_by(cast) %>%
  mutate(mld=max(rounded_depth[binned_sigma_theta>=(0.125+binned_sigma_theta[(round(rounded_depth)==-10)][1])]))

## Generate Figure S1
beam_ts<-ggplot(ctd_smoothed)+
  geom_point(aes(x=binned_CStarAt0,
                 y=binned_o2_mmkg,
                 col=binned_Fluoresence_2,
                 group=cast),
             size=1)+
  scale_color_cmocean(name='delta',direction=1)+
  guides(col=guide_colorbar(title='Chla Fluorescence'))+
  xlab('Beam Attenuation')+
  ylab(expression(paste('O'[2], ' ', mu, 'mol/kg')))+
  theme_bw()+
  theme(strip.background=element_blank(),
        text=element_text(size=14),
        legend.position='bottom')

ggsave('../figures/Figure_S1.pdf',beam_ts)

## Generate contour plots for main text F1
chl_contours<-ggplot(ctd_smoothed)+
  geom_contour_filled(aes(x=time_local,y=round(rounded_depth),z=binned_Fluoresence_2))+
  geom_line(aes(x=time_local,y=mld),linetype='dashed',col='white',size=1.5)+
  geom_point(aes(x=time_local,y=mld),col='white',size=3)+
  scale_fill_cmocean(name='delta',discrete=T)+
  coord_cartesian(ylim=c(-200,0))+
  scale_y_continuous(expand=c(0,0),name='Depth [m]',
                     breaks=c(seq(-300,0,by=50),c(-5,-40,-60,-110)),
                     labels=c(seq(300,0,by=-50),c('SRF','BML','SOM','DCM')))+
  scale_x_datetime(expand=c(0,0),name='Local Time')+
  theme_bw()+
  guides(fill=guide_colorsteps(title=str_wrap('Chlorophyll Fluorescence [mg/m^3]',width=10),
  ))+
  theme(text=element_text(size=16),
        legend.position='bottom',
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.key.width=unit(0.06,'npc'))

oxygen_contours<-ggplot(ctd_smoothed)+
  geom_contour_filled(aes(x=time_local,y=round(rounded_depth),z=binned_oxygen_sat*100))+
  geom_line(aes(x=time_local,y=mld),linetype='dashed',col='white',size=1.5)+
  geom_point(aes(x=time_local,y=mld),col='white',size=3)+
  scale_fill_cmocean(name='amp',discrete=T)+
  coord_cartesian(ylim=c(-200,0))+
  scale_y_continuous(expand=c(0,0),name='Depth [m]',
                     breaks=c(seq(-300,0,by=50),c(-5,-40,-60,-110)),
                     labels=c(seq(300,0,by=-50),c('SRF','BML','SOM','DCM')))+
  scale_x_datetime(expand=c(0,0),name='Local Time')+
  theme_bw()+
  guides(fill=guide_colorsteps(title='Oxygen Saturation (%)'))+
  theme(text=element_text(size=16),
        legend.position='bottom',
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.key.width=unit(0.06,'npc'))

particle_contours<-ggplot(ctd_smoothed)+
  geom_contour_filled(aes(x=time_local,y=round(rounded_depth),z=binned_CStarAt0),breaks=seq(0,0.05,by=0.005))+
  geom_line(aes(x=time_local,y=mld),linetype='dashed',col='white',size=1.5)+
  geom_point(aes(x=time_local,y=mld),col='white',size=3)+
  scale_fill_cmocean(name='turbid',discrete=T)+
  coord_cartesian(ylim=c(-200,0))+
  scale_y_continuous(expand=c(0,0),name='Depth [m]',
                     breaks=c(seq(-300,0,by=50),c(-5,-40,-60,-110)),
                     labels=c(seq(300,0,by=-50),c('SRF','BML','SOM','DCM')))+
  scale_x_datetime(expand=c(0,0),name='Local Time')+
  theme_bw()+
  guides(fill=guide_colorsteps(title=str_wrap('Beam Attenuation [1/m]',width=10)))+
  theme(text=element_text(size=16),
        legend.position='bottom',
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.key.width=unit(0.06,'npc'))

sigma_theta_contours<-ggplot(ctd_smoothed)+
  geom_contour_filled(aes(x=time_local,y=round(rounded_depth),z=binned_sigma_theta),breaks=seq(23.5,27,by=0.25))+
  scale_fill_cmocean(name='dense',discrete=TRUE)+
  geom_line(aes(x=time_local,y=mld),linetype='dashed',col='white',size=1.5)+
  geom_point(aes(x=time_local,y=mld),col='white',size=3)+
  coord_cartesian(ylim=c(-200,0))+
  scale_y_continuous(expand=c(0,0),name='Depth [m]',
                     breaks=c(seq(-300,0,by=50),c(-5,-40,-60,-110)),
                     labels=c(seq(300,0,by=-50),c('SRF','BML','SOM','DCM')))+
  scale_x_datetime(expand=c(0,0),name='Local Time')+
  theme_bw()+
  guides(fill=guide_colorsteps(title='Sigma Theta [kg/m^3]'))+
  theme(text=element_text(size=16),
        legend.position='bottom',
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.key.width=unit(0.06,'npc'))

f1_bottom_panels<-(sigma_theta_contours/chl_contours)|(oxygen_contours/particle_contours)
ggsave('figures/Figure_1cdef_updated.pdf',f1_bottom_panels,scale=1.5)

## Load in ship navigation data to generate cruise map:
nav_track<-data.table::fread('../data/ship_navigation.csv') %>%
  mutate(local_time=lubridate::with_tz(iso_time,tzone='Etc/GMT+4'))

## Loading in CTD sampling locations
surface_ctd<-ctd_smoothed %>% group_by(cast) %>% slice_head(n=1)

world<-map_data('world')

track_component<-ggplot()+
  geom_polygon(data=filter(world,region=='Bermuda'),
               aes(x=long,y=lat,group=group),
               fill='cornsilk1',
               col='darkolivegreen4')+
  coord_sf(ylim=c(31.5,32.5),
           xlim=c(-65,-64))+
  ggsn::scalebar(data=filter(world,region=='Bermuda'),
                 dist_unit='km',model='WGS84',transform=TRUE,
                 dist=25,location='bottomleft',
                 anchor=c('y'=31.575,'x'=-64.95),
                 height=0.25,st.dist=0.2)+
  geom_path(data=nav_track,
            mapping=aes(x=ship_longitude,y=ship_latitude),
            size=1,
            col='darkgrey')+
  geom_point(data=surface_ctd,
             aes(x=binned_lon,
                 y=binned_lat),
             size=2)+
  ggrepel::geom_label_repel(data=filter(surface_ctd,cast==2),
                            aes(x=binned_lon,
                                y=binned_lat,
                                label='12 Oct 8pm'),
                            nudge_y=0.2)+
  ggrepel::geom_label_repel(data=filter(surface_ctd,cast==29),
                            aes(x=binned_lon,
                                y=binned_lat,
                                label='17 Oct 8am'),
                            nudge_x=-0.2)+
  geom_point(aes(y=31.66,
                 x=-64.16),
             size=4,col='red')+
  ggrepel::geom_label_repel(aes(y=31.66,
                                x=-64.16,
                                label='BATS'),nudge_y=-0.075)+
  ylab('Latitude Degrees North')+
  xlab('Longitude Degrees East')+
  theme_minimal()+
  theme(text=element_text(size=18),
        panel.background=element_rect(fill='#F1F8FF',
                                      color='grey65'),
        panel.grid.major=element_line(color='grey45'))

ggsave('../figures/Figure_1a.pdf',track_component)

## Adding depth layer-averaging mentioned by stats in main text:

srf_stats<-ctd_smoothed %>%
  group_by(cast) %>%
  filter(between(rounded_depth,-40,-5)) %>%
  summarize(avg_o2=mean(binned_o2_mmkg),
            avg_o2s=mean(binned_oxygen_sat),
            avg_cp=mean(binned_CStarAt0),
            avg_chl=mean(binned_Fluoresence_2))

dcm_stats<-ctd_smoothed %>%
  group_by(cast) %>%
  filter(between(rounded_depth,-120,-100)) %>%
  summarize(avg_o2=mean(binned_o2_mmkg),
            avg_o2s=mean(binned_oxygen_sat),
            avg_cp=mean(binned_CStarAt0),
            avg_chl=mean(binned_Fluoresence_2))

som_stats<-ctd_smoothed %>%
  group_by(cast) %>%
  filter(between(rounded_depth,-60,-50)) %>%
  summarize(avg_o2=mean(binned_o2_mmkg),
            avg_o2s=mean(binned_oxygen_sat),
            avg_cp=mean(binned_CStarAt0),
            avg_chl=mean(binned_Fluoresence_2))

