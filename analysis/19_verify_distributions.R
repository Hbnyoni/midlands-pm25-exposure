# Independently recompute weighted empirical distributions from native inputs.
options(warn=1)
a <- commandArgs(trailingOnly=TRUE)
root <- if(length(a)) normalizePath(a[1]) else normalizePath('.')
d <- read.csv(file.path(root,'qa/distribution_native_inputs.csv'),check.names=FALSE)
s <- read.csv(file.path(root,'results/population_exposure_distribution_summary.csv'))
groups <- split(seq_len(nrow(d)),interaction(d$year,d$product,d$area,drop=TRUE))
checks <- list()
for(idx in groups) {
  z <- d[idx,]; j <- order(z$concentration); z <- z[j,]
  w <- z$population; x <- z$concentration; cw <- cumsum(w)/sum(w)
  ref <- s[s$year==z$year[1] & s$product==z$product[1] & s$area==z$area[1],]
  got <- c(population=sum(w),mean_pm25=sum(w*x)/sum(w))
  for(p in c(10,25,50,75,90)) got[paste0('p',p)] <- x[which(cw >= p/100)[1]]
  for(cut in c(5,10,15,25,35)) got[paste0('population_percent_above_',cut)] <- 100*sum(w[x>cut])/sum(w)
  got['p90_minus_p10'] <- got['p90']-got['p10']
  expected <- as.numeric(ref[1,names(got)]); err <- abs(got-expected)
  checks[[length(checks)+1]] <- data.frame(year=z$year[1],product=z$product[1],area=z$area[1],
    values=length(got),maximum_absolute_difference=max(err),passed=all(err <= 1e-9*pmax(1,abs(expected))))
}
checks <- do.call(rbind,checks)
write.csv(checks,file.path(root,'qa/independent_R_distribution_checks.csv'),row.names=FALSE)
stopifnot(nrow(checks)==108, all(checks$passed))
cat('PASS:',nrow(checks),'groups;',sum(checks$values),'values; maximum difference',max(checks$maximum_absolute_difference),'\n')
