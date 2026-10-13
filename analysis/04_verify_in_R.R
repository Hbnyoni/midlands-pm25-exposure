# Independent arithmetic verification using base R. No fitted health model.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
root <- normalizePath(file.path(dirname(script), '..'), winslash='/')
res <- file.path(root, 'results')
cells <- read.csv(file.path(res, 'native_cell_annual_metrics.csv'))
annual <- read.csv(file.path(res, 'annual_exposure_summary.csv'))
daily <- read.csv(file.path(res, 'daily_area_exposure_and_scenarios.csv'))
scenarios <- read.csv(file.path(res, 'transferred_rr_scenarios.csv'))
functions <- read.csv(file.path(root, 'evidence', 'published_risk_functions.csv'))
bands <- read.csv(file.path(res, 'annual_who_population_bands.csv'))
checks <- list()
check <- function(name, error, tolerance=1e-8) {
  checks[[length(checks)+1]] <<- data.frame(check=name, maximum_absolute_error=error, tolerance=tolerance, passed=error<=tolerance)
  if (!is.finite(error) || error>tolerance) stop(name, ' failed: ', error)
}
for (y in sort(unique(cells$year))) {
  c <- cells[cells$year==y,]
  a <- annual[annual$year==y & annual$area=='Midlands Province',]
  check(paste(y,'weighted annual concentration'), abs(weighted.mean(c$annual_pm25,c$population)-a$population_weighted_pm25))
  check(paste(y,'population total'), abs(sum(c$population)-a$population),1e-6)
  check(paste(y,'population weighted exceedance days'),abs(weighted.mean(c$days_above_15,c$population)-a$population_weighted_exceedance_days))
  cls <- cut(c$annual_pm25,breaks=c(-Inf,5,10,15,25,35,Inf),right=TRUE,labels=FALSE)-1
  check(paste(y,'WHO bands'),max(abs(cls-c$annual_who_band)),0)
  b <- bands[bands$year==y & bands$area=='Midlands Province',]
  sums <- sapply(0:5,function(k)sum(c$population[cls==k]))
  check(paste(y,'WHO population classes'),max(abs(sums-b$population)),1e-6)
  districts <- annual[annual$year==y & annual$area!='Midlands Province',]
  check(paste(y,'eight reporting districts'),abs(nrow(districts)-8),0)
  for (field in c('area_km2','population','population_above_annual_5','population_above_annual_10',
                  'population_above_annual_15','population_days_above_15','total_population_days')) {
    check(paste(y,'district sum',field),abs(sum(districts[[field]])-a[[field]]),max(1e-6,abs(a[[field]])*1e-12))
  }
  check(paste(y,'district weighted mean'),abs(weighted.mean(districts$population_weighted_pm25,districts$population)-a$population_weighted_pm25))
  check(paste(y,'population days from cell exceedance counts'),abs(sum(c$days_above_15*c$population)-a$population_days_above_15),1e-6)
}
dd <- (daily$pw_daily_excess_above_15+daily$pw_previous_day_excess_at_current_population)/2
check('Lag 0 1 mean difference',max(abs(dd-daily$pw_lag01_reduction),na.rm=TRUE))
check('Lag missingness at first day and gaps',sum(is.na(dd)!=is.na(daily$pw_lag01_reduction)),0)
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i,]
  f <- functions[functions$id==s$function_id,]
  for (col in c('rr_per_10','ci95_lower','ci95_upper','pi80_lower','pi80_upper')) {
    if (is.na(f[[col]])) next
    estimate <- exp(log(f[[col]])*s$delta_pm25/10)
    observed <- if (col=='rr_per_10') s[['rr_estimate']] else s[[col]]
    check(paste('RR',i,col),abs(estimate-observed))
  }
}
# Exact threshold equality and immediately higher exposures must classify differently.
v <- c(5,5.001,10,10.001,15,15.001,25,25.001,35,35.001)
check('Upper inclusive threshold boundaries',max(abs((cut(v,c(-Inf,5,10,15,25,35,Inf),right=TRUE,labels=FALSE)-1)-c(0,1,1,2,2,3,3,4,4,5))),0)
out <- do.call(rbind, checks)
write.csv(out,file.path(root,'qa','independent_R_checks.csv'),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),paste('Checks passed:',nrow(out))),file.path(root,'qa','R_session_and_verification.txt'))
cat(nrow(out),'checks passed. Maximum numerical discrepancy:',max(out$maximum_absolute_error),'\n')
