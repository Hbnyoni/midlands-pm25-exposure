# Independent base-R verification of the added sensitivity analyses.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '.'
read <- function(f) read.csv(file.path(root,f),stringsAsFactors=FALSE,check.names=FALSE)
out <- list()
check <- function(label,got,expected) {
 stopifnot(length(got)==length(expected))
 ok <- (is.na(got)&is.na(expected)) | (!is.na(got)&!is.na(expected)&abs(got-expected)<=1e-8*pmax(1,abs(expected)))
 out[[length(out)+1L]] <<- data.frame(check=label,values=length(got),passed=all(ok),maximum_difference=max(c(0,abs(got-expected)),na.rm=TRUE))
}
p <- subset(read('results/exposure_product_sensitivity.csv'),area!='Midlands Province')
r <- read('results/annual_district_rank_stability.csv')
for(product in c('GHAP','ACAG')) for(y in 2017:2022) {
 a <- r[r$product==product & r$year==y,]
 b <- p[p$year==y,]; b <- b[match(a$area,b$area),]
 v <- if(product=='GHAP') b$ghap_v1_pw_pm25 else b$comparator_pw_pm25
 check(paste('annual concentrations',product,y),a$concentration,v)
 check(paste('annual ranks',product,y),a$rank,rank(-v,ties.method='average'))
}
g <- read('results/annual_product_rank_agreement.csv')
for(y in 2017:2022) {
 z <- p[p$year==y,]
 check(paste('Spearman',y),g$spearman[g$year==y],cor(z$ghap_v1_pw_pm25,z$comparator_pw_pm25,method='spearman'))
}
d <- read('results/health_daily_window_exposure.csv')
w <- read('results/acute_window_rank_stability.csv')
for(y in 2017:2022) for(win in c(1,3,7)) {
 a <- w[w$year==y & w$window_days==win,]
 v <- sapply(a$area,function(area) quantile(d$delta_cap15[d$year==y & d$window_days==win & d$area==area],.95,na.rm=TRUE,type=7))
 check(paste('acute p95',y,win),a$delta_cap15,v)
 check(paste('acute ranks',y,win),a$rank,rank(-v))
}
c <- subset(read('results/cumulative_area_profiles.csv'),area!='Midlands Province')
s <- read('results/district_priority_stability.csv');c <- c[match(s$area,c$area),]
check('six-year GHAP ranks',s$ghap_rank,rank(-c$ghap_six_year_mean))
check('six-year ACAG ranks',s$acag_rank,rank(-c$acag_six_year_mean))
for(type in c('annual','acute')) {
 a <- if(type=='annual') r else w
 for(fun in c('min','max','median')) {
  v <- sapply(s$area,function(area) get(fun)(a$rank[a$area==area]))
  check(paste(type,fun),s[[paste0(type,'_rank_',fun)]],v)
 }
}
m <- read('evidence/agache2024_asthma_functions.csv')
check('Agache ED Table 2 original ratios',m$rr[m$endpoint=='asthma_ED'],c(1.012,1.013,1.014,1.023,1.016))
check('Agache admissions original lower limits',m$ci95_lower[m$endpoint=='asthma_admissions'],c(.992,1,.997,1,.997))
f <- read('evidence/health_profile_functions.csv')
h <- read('results/matched_asthma_evidence_sensitivity.csv')
for(endpoint in c('asthma_ED','asthma_admissions')) for(source in c('Agache 2024','Lin 2024')) {
 z <- h[h$endpoint==endpoint & h$source==source,]
 coef <- if(source=='Agache 2024') m[m$endpoint==endpoint & m$lag_days==0,] else f[f$endpoint==endpoint & f$window_days==1,]
 delta <- mapply(function(area,y) quantile(d$delta_cap15[d$area==area & d$year==y & d$window_days==1],.95,na.rm=TRUE),z$area,z$year)
 check(paste('matched p95',endpoint,source),z$delta_pm25,delta)
 for(limits in c('rr_estimate','ci95_lower','ci95_upper')) {
  col <- if(limits=='rr_estimate') 'rr' else limits
  check(paste(endpoint,source,limits),z[[limits]],exp(log(coef[[col]])*delta/coef$increment_ug_m3))
 }
}
dd <- read('results/daily_area_exposure_and_scenarios.csv');dd$month <- as.integer(format(as.Date(dd$date),'%m'))
t <- read('results/seasonal_population_exposure.csv')
v <- mapply(function(area,mo) sum(dd$population_days_above_15[dd$area==area & dd$month==mo],na.rm=TRUE),t$area,t$month)
den <- mapply(function(area,mo) sum(dd$valid_exposure_population[dd$area==area & dd$month==mo],na.rm=TRUE),t$area,t$month)
check('monthly population-days',t$exposed_population_days,v)
check('monthly denominators',t$observed_population_days,den)
check('monthly exposure fraction',t$exposed_fraction,v/den)
share <- v/ave(v,t$area,FUN=sum)
check('monthly distribution of population-days',t$share_of_area_exposed_population_days,share)
result <- do.call(rbind,out)
write.csv(result,file.path(root,'qa/independent_R_inference_checks.csv'),row.names=FALSE)
print(data.frame(groups=nrow(result),values=sum(result$values),passed=all(result$passed),maximum_difference=max(result$maximum_difference)))
stopifnot(all(result$passed))
