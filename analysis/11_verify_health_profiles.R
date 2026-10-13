# Independent base-R checks of cumulative totals, moving windows, IQR scaling,
# coefficient intervals, cross-year population weighting and outcome transforms.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'Quantitative/Midlands_Public_Exposure'
read <- function(f) read.csv(file.path(root,f),stringsAsFactors=FALSE,check.names=FALSE)
d <- read('results/health_daily_window_exposure.csv')
old <- read('results/daily_area_exposure_and_scenarios.csv')
c <- read('results/cumulative_area_profiles.csv')
r <- read('results/outcome_specific_health_profiles.csv')
f <- read('evidence/health_profile_functions.csv')
j <- read('qa/health_january_window_inputs.csv')
ans <- list(); k <- 0L
check <- function(label,got,expected) {
 stopifnot(length(got)==length(expected))
 ok <- (is.na(got)&is.na(expected)) | (!is.na(got)&!is.na(expected)&abs(got-expected)<=1e-8*pmax(1,abs(expected)))
 k <<- k+1L
 dif <- abs(got-expected); md <- if(all(is.na(dif))) 0 else max(dif,na.rm=TRUE)
 ans[[k]] <<- data.frame(check=label,values=length(got),passed=all(ok),maximum_difference=md)
}
check('IQR conversion per 10',f$rr_per_10,exp(log(f$rr)*10/f$increment_ug_m3))
key <- function(date,area) paste(date,area,sep='|')
s <- d[d$window_days==1,]; ix <- match(key(s$date,s$area),key(old$date,old$area))
check('Same-day excess against original independent export',s$delta_cap15,old$pw_daily_excess_above_15[ix])
for(area in unique(d$area)) {
 x <- d[d$area==area & d$window_days==1,]; x <- x[order(x$date),]
 z <- c[c$area==area,]
 check(paste(area,'cumulative fixed-population intensity'),z$cumulative_excess_ug_m3_days,sum(x$delta_cap15_fixed2022,na.rm=TRUE))
 o <- old[old$area==area,]
 check(paste(area,'population-days'),z$population_days_above15_annual_weights,sum(o$population_days_above_15,na.rm=TRUE))
 check(paste(area,'population excess concentration days'),z$population_excess_concentration_days_annual_weights,sum(o$population*o$pw_daily_excess_above_15,na.rm=TRUE))
 for(w in c(3,7)) {
  q <- d[d$area==area & d$window_days==w,];q<-q[order(q$date),]
  calc <- rep(NA_real_,nrow(x))
  for(t in seq_len(nrow(x))) if(t>=w) calc[t]<-mean(x$delta_cap15_fixed2022[(t-w+1):t])
  check(paste(area,w,'fixed 2022 complete windows'),q$delta_cap15_fixed2022,calc)
  current <- rep(NA_real_,nrow(x))
  for(t in seq_len(nrow(x))) {
   if(as.integer(format(as.Date(x$date[t]),'%j'))<=7) {
    vals <- j[j$date==x$date[t]&j$area==area&j$lag<w,]
    stopifnot(nrow(vals)==w);current[t]<-mean(vals$excess_at_current_year_population)
   } else if(t>=w) current[t]<-mean(x$delta_cap15[(t-w+1):t])
  }
  check(paste(area,w,'current year weights including January'),q$delta_cap15,current)
 }
}
for(column in c('population_2022','population_days_above15_annual_weights','population_excess_concentration_days_annual_weights'))
 check(paste('District sum',column),sum(c[c$area!='Midlands Province',column]),c[c$area=='Midlands Province',column])
fi <- match(r$function_id,f$id)
for(pair in list(c('rr_estimate','rr'),c('ci95_lower','ci95_lower'),c('ci95_upper','ci95_upper'),c('pi80_lower','pi80_lower'),c('pi80_upper','pi80_upper'))) {
 calc<-exp(log(f[fi,pair[2]])*r$delta_pm25/f$increment_ug_m3[fi])
 check(paste('All health transformations',pair[1]),r[,pair[1]],calc)
}
for(area in unique(d$area)) for(w in c(1,3,7)) {
 vals<-d$delta_cap15[d$area==area&d$year==2022&d$window_days==w]
 rr<-r[r$area==area&r$period=='2022'&r$window_days==w&r$summary=='p95',]
 check(paste(area,w,'p95 reduction'),rr$delta_pm25,rep(as.numeric(quantile(vals,.95,na.rm=TRUE,type=7)),nrow(rr)))
}
out<-do.call(rbind,ans);write.csv(out,file.path(root,'qa/independent_R_health_checks.csv'),row.names=FALSE)
stopifnot(all(out$passed))
cat(sum(out$values),'health-profile values checked in',nrow(out),'groups; all PASS\n')
