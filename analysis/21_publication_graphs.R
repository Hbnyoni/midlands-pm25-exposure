# Statistical graphics authored in R/ggplot2 at a fixed 180 mm journal width.
# No smoothing, inferred clinical grades, local health fits or invented intervals.
options(encoding='UTF-8',warn=1)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
res <- file.path(root,'results'); out <- file.path(root,'figures')
rd <- function(name) read.csv(file.path(res,name),check.names=FALSE)
ev <- function(name) read.csv(file.path(root,'evidence',name),check.names=FALSE)
blue <- '#176B99'; orange <- '#CF6B3B'; ink <- '#202B33'; teal <- '#257D79'
pal <- c(GHAP=blue,ACAG=orange,'Agache 2024'=blue,'Lin 2024'=orange,
         Admissions=blue,'ED visits'=orange,'Respiratory mortality'=blue,'All-cause mortality'=orange)
shapes <- c(GHAP=16,ACAG=15,'Agache 2024'=16,'Lin 2024'=15,Admissions=16,'ED visits'=15,
            'Respiratory mortality'=16,'All-cause mortality'=15)
C <- rd('cumulative_area_profiles.csv')
area_order <- C$area[C$area!='Midlands Province'][order(-C$cumulative_excess_ug_m3_days[C$area!='Midlands Province'])]
area_y <- function(x) 9-match(x,area_order)
sty <- theme_classic(base_size=10.5,base_family='Arial') + theme(
  text=element_text(colour=ink),axis.text=element_text(size=9.5,colour=ink),
  axis.title=element_text(size=10.5),axis.title.x=element_text(margin=margin(t=7)),
  axis.title.y=element_text(margin=margin(r=7)),axis.line=element_line(linewidth=.55,colour=ink),
  axis.ticks=element_line(linewidth=.5,colour=ink),axis.ticks.length=unit(1.5,'mm'),
  panel.grid.major.x=element_line(colour='#E4E9EC',linewidth=.20),
  legend.position='top',legend.title=element_blank(),legend.text=element_text(size=9.5),
  legend.key.width=unit(7,'mm'),legend.key.height=unit(4,'mm'),legend.margin=margin(0,0,3,0),
  plot.tag=element_text(size=12,face='bold'),plot.tag.position='topleft',
  plot.margin=margin(7,8,6,7),plot.background=element_rect(fill='white',colour=NA))
theme_set(sty)
product_scales <- function() list(scale_colour_manual(values=pal,breaks=names(pal)),scale_shape_manual(values=shapes,breaks=names(shapes)))
dist_y <- function(labels=TRUE) list(scale_y_continuous(breaks=1:8,labels=if(labels) rev(area_order) else NULL,
                                                      limits=c(.5,8.5),expand=c(0,0)),labs(y=NULL),
                                    theme(axis.ticks.y=element_blank(),axis.line.y=element_blank()))
row_guides <- function() geom_hline(yintercept=seq(1.5,7.5,1),colour='#E9EDEF',linewidth=.2)
records <- list()
publish <- function(p,name,height) {
  if(length(args)>1 && !name %in% args[-1]) return(invisible(NULL))
  for(ext in c('png','tif','pdf','svg')) {
    f <- file.path(out,paste0(name,'.new.',ext))
    if(ext=='png') ragg::agg_png(f,width=180,height=height,units='mm',res=600,background='white')
    if(ext=='tif') ragg::agg_tiff(f,width=180,height=height,units='mm',res=600,background='white',compression='lzw')
    if(ext=='pdf') grDevices::cairo_pdf(f,width=180/25.4,height=height/25.4,family='Arial',bg='white')
    if(ext=='svg') svglite::svglite(f,width=180/25.4,height=height/25.4,bg='white',system_fonts=list(sans='Arial'))
    # Construct text metrics only after a device with native font support is open.
    print(p);dev.off()
    target <- file.path(out,paste0(name,'.',ext))
    if(!file.copy(f,target,overwrite=TRUE)) stop('Could not write ',target)
    unlink(f)
  }
  editable <- file.path(out,'editable_R');dir.create(editable,showWarnings=FALSE)
  saveRDS(p,file.path(editable,paste0(name,'.rds')))
  records[[length(records)+1]] <<- data.frame(figure=name,width_mm=180,height_mm=height,
      minimum_font_pt=9,dpi=600,engine='R ggplot2 and patchwork')
  cat('Created',name,'\n');flush.console()
}

# Figure 4: daily chronology and monthly population-time with a common percentage scale.
d <- rd('daily_area_exposure_and_scenarios.csv');d <- d[d$area=='Midlands Province',];d$date <- as.Date(d$date)
m <- rd('monthly_exposure.csv');m <- m[m$area=='Midlands Province',]
m$pct <- 100*m$mean_population_fraction_above_15
m$yr <- factor(m$year,levels=2022:2017);m$mo <- factor(m$month,levels=1:12,labels=month.abb)
a <- ggplot(d,aes(date,pw_pm25)) + geom_line(colour=blue,linewidth=.50,na.rm=FALSE) +
  geom_hline(yintercept=15,colour=orange,linetype='longdash',linewidth=.65) +
  scale_x_date(date_breaks='1 year',date_labels='%Y',expand=c(0,0),limits=as.Date(c('2017-01-01','2022-12-31'))) +
  scale_y_continuous(breaks=seq(0,60,10),limits=c(0,60),expand=c(0,0)) +
  labs(x=NULL,y='PM\u2082.\u2085 (\u00b5g/m\u00b3)',tag='a') +
  theme(panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour='#E4E9EC',linewidth=.2))
b <- ggplot(m,aes(mo,yr,fill=pct)) + geom_tile(colour='white',linewidth=.4) +
  geom_text(aes(label=sprintf('%.0f',pct),colour=pct>=55),size=3.25,family='Arial',show.legend=FALSE) +
  scale_colour_manual(values=c('FALSE'=ink,'TRUE'='white')) +
  scale_fill_gradientn(colours=c('#F1F6F8','#B6DAE4','#61AEC4','#276B94','#183A5B'),limits=c(0,100),
    breaks=c(0,25,50,75,100),name='Exposed population-days (%)',
    guide=guide_colourbar(direction='horizontal',title.position='top',barwidth=unit(75,'mm'),barheight=unit(3,'mm'))) +
  labs(x=NULL,y=NULL,tag='b') + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),legend.position='bottom',legend.title=element_text(size=9.5))
publish(a/b+plot_layout(heights=c(1.15,1)),'Figure_4_temporal_exposure',148)

# Figure 5: sustained concentrations versus population-time, with common district order.
z <- C[C$area!='Midlands Province',];z$y <- area_y(z$area)
pts <- rbind(data.frame(y=z$y,value=z$ghap_six_year_mean,product='GHAP'),
             data.frame(y=z$y,value=z$acag_six_year_mean,product='ACAG'))
a <- ggplot(z) + row_guides() + geom_segment(aes(x=ghap_six_year_mean,xend=acag_six_year_mean,y=y,yend=y),colour='#A7B7C1',linewidth=.75) +
  geom_point(data=pts,aes(value,y,colour=product,shape=product),size=2.65) + product_scales()+dist_y() +
  scale_x_continuous(breaks=seq(10,14,1),limits=c(10,14),expand=c(.015,0)) + labs(x='Six-year mean PM\u2082.\u2085\n(\u00b5g/m\u00b3)',tag='a')
z$millions <- z$population_days_above15_annual_weights/1e6
b <- ggplot(z,aes(millions,y)) + row_guides() + geom_col(fill=blue,width=.48,orientation='y') +
  geom_text(aes(label=sprintf('%.1f',millions)),hjust=-.18,size=3.3,family='Arial') + dist_y(FALSE) +
  scale_x_continuous(limits=c(0,270),breaks=seq(0,250,50),expand=c(0,0)) + labs(x='Population-days >15 \u00b5g/m\u00b3\n(millions)',tag='b')
publish(a+b+plot_layout(guides='collect') & theme(legend.position='top'),'Figure_5_district_exposure',131)

# Figure 6: all six years and the cross-product six-year rank comparison.
comp <- rd('exposure_product_sensitivity.csv');comp <- comp[comp$area=='Midlands Province',]
ann <- rd('annual_exposure_summary.csv');ann <- ann[ann$area=='Midlands Province',]
pt <- rbind(data.frame(year=comp$year,value=comp$ghap_v1_pw_pm25,product='GHAP'),data.frame(year=comp$year,value=comp$comparator_pw_pm25,product='ACAG'))
a <- ggplot(pt,aes(year,value,colour=product,shape=product)) + geom_hline(yintercept=5,colour='#687B87',linetype='longdash',linewidth=.6) +
  geom_line(linewidth=.7)+geom_point(size=2.45)+product_scales()+scale_x_continuous(breaks=2017:2022,labels=17:22) +
  scale_y_continuous(limits=c(4,16),breaks=c(5,10,15))+labs(x='Year (20xx)',y='Annual PM\u2082.\u2085 (\u00b5g/m\u00b3)',tag='a')
pt <- rbind(data.frame(year=ann$year,value=100*ann$population_above_annual_10/ann$population,product='GHAP'),
             data.frame(year=comp$year,value=100*comp$comparator_population_above_10/comp$comparator_population,product='ACAG'))
b <- ggplot(pt,aes(year,value,colour=product,shape=product))+geom_line(linewidth=.7)+geom_point(size=2.45)+product_scales()+
  scale_x_continuous(breaks=2017:2022,labels=17:22)+scale_y_continuous(limits=c(0,105),breaks=seq(0,100,25))+
  labs(x='Year (20xx)',y='Population >10 \u00b5g/m\u00b3 (%)',tag='b')
s <- rd('district_priority_stability.csv');s$y <- area_y(s$area)
pt <- rbind(data.frame(y=s$y+.12,value=s$ghap_rank,product='GHAP'),data.frame(y=s$y-.12,value=s$acag_rank,product='ACAG'))
c <- ggplot(s)+row_guides()+geom_segment(aes(x=ghap_rank,xend=acag_rank,y=y+.12,yend=y-.12),colour='#A7B7C1',linewidth=.75)+
  geom_point(data=pt,aes(value,y,colour=product,shape=product),size=2.6)+product_scales()+dist_y()+
  scale_x_continuous(breaks=1:8,limits=c(.6,8.4),expand=c(0,0))+labs(x='Six-year concentration rank (1 = highest)',tag='c')
b <- b+guides(colour='none',shape='none');c <- c+guides(colour='none',shape='none')
publish((a+b)/c+plot_layout(heights=c(1,1.22),guides='collect') & theme(legend.position='top'),'Figure_6_product_sensitivity',174)

# Figure 8: matched-source asthma, COPD and long-term IHD. Intervals retain their meanings.
H <- rd('matched_asthma_evidence_sensitivity.csv');R <- rd('outcome_specific_health_profiles.csv')
forest <- function(z,letter,xlab,labels=FALSE,prediction=FALSE) {
  z$base <- area_y(z$area)
  z$y <- z$base + ifelse(z$group %in% c('Agache 2024','GHAP'),.14,ifelse(z$group %in% c('Lin 2024','ACAG'),-.14,0))
  z$est <- 100*(z$rr_estimate-1);z$lo <- 100*(z$ci95_lower-1);z$hi <- 100*(z$ci95_upper-1)
  p <- ggplot(z,aes(y=y,colour=group,shape=group))+row_guides()+geom_vline(xintercept=0,colour='#6D7880',linetype='longdash',linewidth=.55)
  if(prediction) p <- p+geom_segment(aes(x=100*(pi80_lower-1),xend=100*(pi80_upper-1),yend=y),linewidth=2.3,alpha=.16)
  p <- p+geom_segment(aes(x=lo,xend=hi,yend=y),linewidth=.65)+
    geom_segment(aes(x=lo,xend=lo,y=y-.065,yend=y+.065),linewidth=.6)+
    geom_segment(aes(x=hi,xend=hi,y=y-.065,yend=y+.065),linewidth=.6)+
    geom_point(aes(x=est),size=2.35)+product_scales()+dist_y(labels)+labs(x=xlab,tag=letter)+
    scale_x_continuous(breaks=scales::breaks_pretty(n=4),expand=expansion(mult=c(.025,.06)))+
    theme(legend.position='top',legend.direction='vertical',legend.key.height=unit(3.6,'mm'),legend.spacing.y=unit(0,'mm'))
  p
}
z <- H[H$year==2022&H$endpoint=='asthma_ED'&H$area!='Midlands Province',];z$group <- z$source
a <- forest(z,'a','Asthma ED\ncontrast (%)',TRUE)
z <- R[R$period=='2022'&R$summary=='p95'&R$function_id=='lin2024_COPD_admissions_w7'&R$area!='Midlands Province',];z$group <- 'GHAP'
b <- forest(z,'b','COPD admissions\ncontrast (%)')+guides(colour='none',shape='none')
z <- R[R$summary=='six_year_mean'&R$function_id=='orellano2024_IHD_mortality'&R$area!='Midlands Province',];z$group <- z$product
c <- forest(z,'c','IHD mortality\ncontrast (%)',prediction=TRUE)
publish(a+b+c+plot_layout(widths=c(1,1,1)),'Figure_8_health_profiles',151)

# Supplement S2: both kinds of coefficient interval remain explicit.
r <- rd('transferred_rr_scenarios.csv');r <- r[r$area=='Midlands Province'&r$window=='annual'&r$reference_cap==5,]
r <- r[r$function_id %in% c('orellano2024_respiratory','orellano2024_allcause'),]
r$endpoint <- ifelse(r$function_id=='orellano2024_respiratory','Respiratory mortality','All-cause mortality')
r$x <- r$year+ifelse(r$endpoint=='Respiratory mortality',-.12,.12)
a <- ggplot(r,aes(x,100*(rr_estimate-1),colour=endpoint,shape=endpoint))+
  geom_linerange(aes(ymin=100*(pi80_lower-1),ymax=100*(pi80_upper-1)),linewidth=2.7,alpha=.16)+
  geom_linerange(aes(ymin=100*(ci95_lower-1),ymax=100*(ci95_upper-1)),linewidth=.7)+geom_point(size=2.4)+
  geom_hline(yintercept=0,linetype='longdash',colour='#687B87',linewidth=.6)+product_scales()+
  scale_x_continuous(breaks=2017:2022)+labs(x=NULL,y='RR contrast above 1 (%)',tag='a')
ds <- d;ds$valid <- is.finite(ds$liu2019_respiratory_rr)&is.finite(ds$liu2019_respiratory_ci_lower)&is.finite(ds$liu2019_respiratory_ci_upper)
ds$segment <- cumsum(!ds$valid);ds <- ds[ds$valid,]
b <- ggplot(ds,aes(date,100*(liu2019_respiratory_rr-1),group=segment))+
  geom_ribbon(aes(ymin=100*(liu2019_respiratory_ci_lower-1),ymax=100*(liu2019_respiratory_ci_upper-1)),fill=blue,alpha=.18)+
  geom_line(colour=blue,linewidth=.5)+scale_x_date(date_breaks='1 year',date_labels='%Y',expand=c(0,0))+
  labs(x=NULL,y='RR contrast above 1 (%)',tag='b')
publish(a/b+plot_layout(heights=c(1,1)),'Figure_S2_supplementary_risk',146)

# S3: separate average-window fits and separate individual lags; never joined or summed.
f <- ev('health_profile_functions.csv');mm <- ev('agache2024_asthma_functions.csv')
winplot <- function(ends,letter) {
  z <- f[f$endpoint %in% ends & f$window_days>0,]
  z$kind <- ifelse(grepl('admissions',z$endpoint),'Admissions','ED visits')
  z$x <- match(z$window_days,c(1,3,7))+ifelse(z$kind=='Admissions',-.075,.075)
  ggplot(z,aes(x,100*(rr_per_10-1),colour=kind,shape=kind))+
    geom_linerange(aes(ymin=100*(ci95_lower_per_10-1),ymax=100*(ci95_upper_per_10-1)),linewidth=.7)+geom_point(size=2.5)+product_scales()+
    scale_x_continuous(breaks=1:3,labels=c('0','0\u20132','0\u20136'),limits=c(.65,3.35))+
    scale_y_continuous(limits=c(0,15),breaks=c(0,5,10,15))+
    labs(x='Averaging window (lag days)',y='Excess rate per 10 \u00b5g/m\u00b3 (%)',tag=letter)
}
lagplot <- function(endpoint,letter) {
  z <- mm[mm$endpoint==endpoint,];z$kind <- ifelse(grepl('admissions',z$endpoint),'Admissions','ED visits')
  ggplot(z,aes(lag_days,100*(rr-1),colour=kind,shape=kind))+
    geom_hline(yintercept=0,linetype='longdash',colour='#687B87',linewidth=.6)+
    geom_linerange(aes(ymin=100*(ci95_lower-1),ymax=100*(ci95_upper-1)),linewidth=.7)+geom_point(size=2.5)+product_scales()+
    scale_x_continuous(breaks=0:4)+scale_y_continuous(limits=c(-1,4),breaks=-1:4)+
    labs(x='Single lag (days)',y='Relative excess per 10 \u00b5g/m\u00b3 (%)',tag=letter)
}
a <- winplot(c('asthma_admissions','asthma_ED'),'a');b <- winplot(c('COPD_admissions','COPD_ED'),'b')
c <- lagplot('asthma_ED','c');dd <- lagplot('asthma_admissions','d')
b <- b+guides(colour='none',shape='none');c <- c+guides(colour='none',shape='none');dd <- dd+guides(colour='none',shape='none')
publish((a+b)/(c+dd)+plot_layout(guides='collect') & theme(legend.position='bottom'),'Figure_S3_published_windows',174)

# S4: complete, unsmoothed ranking matrices.
ar <- rd('annual_district_rank_stability.csv');wr <- rd('acute_window_rank_stability.csv')
rank_scale <- function() scale_fill_gradientn(colours=rev(c('#F0F5F8','#C0DCE7','#7DB8CF','#3C89AE','#205B84','#153A60')),
    limits=c(1,8),breaks=1:8,name='District rank (1 = greatest exposure)',
    guide=guide_colourbar(direction='horizontal',title.position='top',barwidth=unit(85,'mm'),barheight=unit(3,'mm')))
matrix_theme <- theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),
   legend.position='bottom',legend.title=element_text(size=9.5),axis.title.x=element_text(size=9.5))
matrix_annual <- function(product,letter,labels) {
  z <- ar[ar$product==product,];z$y <- factor(z$area,levels=rev(area_order))
  ggplot(z,aes(factor(year),y,fill=rank))+geom_tile(colour='white',linewidth=.4)+
    geom_text(aes(label=rank,colour=rank<=4),size=3.25,family='Arial',show.legend=FALSE)+
    scale_colour_manual(values=c('FALSE'=ink,'TRUE'='white'))+rank_scale()+
    scale_x_discrete(labels=17:22,expand=c(0,0))+scale_y_discrete(labels=if(labels) waiver() else NULL,expand=c(0,0))+
    labs(x='Year (20xx)',y=NULL,tag=letter)+matrix_theme
}
a <- matrix_annual('GHAP','a',TRUE);b <- matrix_annual('ACAG','b',FALSE)
wr$y <- factor(wr$area,levels=rev(area_order));wr$window <- factor(wr$window_days,levels=c(1,3,7),labels=c('1 day','3 days','7 days'))
c <- ggplot(wr,aes(factor(year),y,fill=rank))+geom_tile(colour='white',linewidth=.3)+
  geom_text(aes(label=rank,colour=rank<=4),size=3.15,family='Arial',show.legend=FALSE)+
  scale_colour_manual(values=c('FALSE'=ink,'TRUE'='white'))+rank_scale()+
  facet_grid(.~window,switch='x')+scale_x_discrete(labels=17:22,expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
  labs(x='Year (20xx)',y=NULL,tag='c')+matrix_theme+
  theme(strip.background=element_blank(),strip.placement='outside',strip.text=element_text(size=9.5),panel.spacing=unit(3,'mm'))
publish((a+b)/c+plot_layout(heights=c(1,1),guides='collect') & theme(legend.position='bottom'),'Figure_S4_rank_stability',175)

# S5: population-weighted empirical distributions, not kernel density estimates.
e <- rd('province_exposure_ecdf.csv');e <- e[e$year==2022,]
stopifnot(all(e$population_at_or_below_percent>=0),max(e$population_at_or_below_percent)<=100+1e-8)
e$population_at_or_below_percent <- pmin(100,e$population_at_or_below_percent)
xmax <- ceiling(max(e$concentration))
for(prod in c('GHAP','ACAG')) e <- rbind(e,data.frame(year=2022,product=prod,
    concentration=c(5,xmax),population_at_or_below_percent=c(0,100)))
e <- e[order(e$product,e$concentration),]
a <- ggplot(e,aes(concentration,population_at_or_below_percent,colour=product))+
  geom_vline(xintercept=c(5,10,15),linetype='longdash',colour='#9AA7AF',linewidth=.5)+
  geom_step(linewidth=.8,direction='hv',show.legend=FALSE)+scale_colour_manual(values=pal,guide='none')+
  scale_x_continuous(limits=c(5,ceiling(max(e$concentration))),breaks=c(5,7.5,10,12.5,15))+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=c(0,0))+
  labs(x='Annual PM\u2082.\u2085 (\u00b5g/m\u00b3)',y='Population at or below (%)',tag='a')
z <- rd('population_exposure_distribution_summary.csv');z <- z[z$year==2022 & z$area!='Midlands Province',]
z$y <- area_y(z$area)+ifelse(z$product=='GHAP',.13,-.13)
b <- ggplot(z,aes(y=y,colour=product,shape=product))+row_guides()+
  geom_vline(xintercept=10,linetype='longdash',colour='#9AA7AF',linewidth=.55)+
  geom_segment(aes(x=p10,xend=p90,yend=y),linewidth=.7)+
  geom_segment(aes(x=p25,xend=p75,yend=y),linewidth=2.1,alpha=.3)+geom_point(aes(x=p50),size=2.6)+product_scales()+dist_y()+
  labs(x='Annual PM\u2082.\u2085 (\u00b5g/m\u00b3)',tag='b')+scale_x_continuous(breaks=seq(6,16,2))
publish(a/b+plot_layout(heights=c(.82,1.18),guides='collect') & theme(legend.position='top'),'Figure_S5_exposure_distributions',180)
record_file <- file.path(root,'qa/ggplot_design_record.csv')
record_table <- do.call(rbind,records)
if(length(args)>1 && file.exists(record_file)) {
  old <- read.csv(record_file);record_table <- rbind(old[!old$figure %in% record_table$figure,],record_table)
}
write.csv(record_table,record_file,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'analysis/R_graphics_session.txt'))
