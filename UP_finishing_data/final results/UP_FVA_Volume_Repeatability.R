#!/usr/bin/env Rscript
# Under Pressure: FVA repeatability at different shot-volume thresholds
# Base R only. Uses player_seasons.csv from five-season pipeline.
args <- commandArgs(trailingOnly=TRUE)
opt <- list(input='player_seasons.csv',output='fva_volume_repeatability',bootstrap='1000',seed='2026')
if(length(args) %% 2) stop('Use --input FILE --output DIR --bootstrap N')
if(length(args)) for(i in seq(1,length(args),by=2)) {
 key <- sub('^--','',args[i]); if(!key %in% names(opt)) stop('Unknown option: ',key)
 opt[[key]] <- args[i+1]
}
B <- as.integer(opt$bootstrap); if(is.na(B)||B<0) stop('Invalid bootstrap count')
set.seed(as.integer(opt$seed));dir.create(opt$output,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(opt$input,stringsAsFactors=FALSE)
req <- c('season','player_id','player_name','shots','fva_per_shot','xg_per_shot','fva_usable_shots','fva_missing_shots')
if(length(setdiff(req,names(d)))) stop('Missing: ',paste(setdiff(req,names(d)),collapse=', '))
d$year <- as.integer(substr(d$season,1,4));d$player_id <- as.character(d$player_id)
d <- d[!is.na(d$player_id)&d$player_id!='-1'&is.finite(d$year),]
if(anyDuplicated(d[c('player_id','year')])) stop('Duplicate player-season rows')
d <- d[order(d$player_id,d$year),]
# All adjacent seasons with valid xGOT-derived FVA, then apply each threshold to BOTH years.
pairs <- list()
for(id in unique(d$player_id)) {
 z <- d[d$player_id==id,];if(nrow(z)<2) next
 for(i in seq_len(nrow(z)-1)) if(z$year[i+1]==z$year[i]+1) {
  pairs[[length(pairs)+1]] <- data.frame(player_id=id,player_name=z$player_name[i+1],
    from_season=z$season[i],to_season=z$season[i+1],
    shots_prior=z$shots[i],shots_next=z$shots[i+1],
    usable_prior=z$fva_usable_shots[i],usable_next=z$fva_usable_shots[i+1],
    missing_prior=z$fva_missing_shots[i],missing_next=z$fva_missing_shots[i+1],
    fva_prior=z$fva_per_shot[i],fva_next=z$fva_per_shot[i+1],
    xg_prior=z$xg_per_shot[i],xg_next=z$xg_per_shot[i+1])
 }
}
if(!length(pairs)) stop('No adjacent-season pairs found')
p <- do.call(rbind,pairs)
p <- p[is.finite(p$fva_prior)&is.finite(p$fva_next)&is.finite(p$xg_prior)&is.finite(p$xg_next)&p$usable_prior>0&p$usable_next>0,]
write.csv(p,file.path(opt$output,'all_adjacent_pairs.csv'),row.names=FALSE)
thresholds <- c(30,50,75,100)
# Cluster bootstrap resamples players, not individual pairs, as some players contribute several pairs.
cor_safe <- function(a,b,method='pearson') {
 ok <- is.finite(a)&is.finite(b)
 if(sum(ok)<4 || sd(a[ok])==0 || sd(b[ok])==0) return(NA_real_)
 suppressWarnings(cor(a[ok],b[ok],method=method))
}
boot_ci <- function(z,a,b) {
 if(B==0) return(c(NA_real_,NA_real_,0))
 groups <- split(z,z$player_id);ids <- names(groups)
 if(length(ids)<4) return(c(NA_real_,NA_real_,0))
 v <- replicate(B,{
   chosen <- sample(ids,length(ids),replace=TRUE)
   sampled <- do.call(rbind,unname(groups[chosen]))
   cor_safe(sampled[[a]],sampled[[b]])
 })
 v <- v[is.finite(v)]
 if(length(v)<max(30,B*.5)) return(c(NA_real_,NA_real_,length(v)))
 c(unname(quantile(v,c(.025,.975))),length(v))
}
results <- list(); membership <- list()
for(t in thresholds) {
 z <- p[p$shots_prior>=t & p$shots_next>=t,]
 membership[[length(membership)+1]] <- data.frame(threshold=t,pairs=nrow(z),players=length(unique(z$player_id)))
 for(metric in c('fva','xg')) {
  a <- paste0(metric,'_prior'); b <- paste0(metric,'_next')
  r <- cor_safe(z[[a]],z[[b]]); ci <- boot_ci(z,a,b)
  results[[length(results)+1]] <- data.frame(threshold=t,metric=metric,
    pairs=nrow(z),players=length(unique(z$player_id)),pearson=r,
    spearman=cor_safe(z[[a]],z[[b]],'spearman'),ci_low=ci[1],ci_high=ci[2],bootstrap_valid=ci[3],
    mean_prior=mean(z[[a]]),mean_next=mean(z[[b]]))
 }
 write.csv(z,file.path(opt$output,paste0('pairs_min',t,'.csv')),row.names=FALSE)
}
r <- do.call(rbind,results);write.csv(r,file.path(opt$output,'volume_correlations.csv'),row.names=FALSE)
# Complete-case sensitivity: no unexpectedly missing xGOT in either season.
sens <- list()
for(t in thresholds) {
 z <- p[p$shots_prior>=t&p$shots_next>=t&p$missing_prior==0&p$missing_next==0,]
 sens[[length(sens)+1]] <- data.frame(threshold=t,pairs=nrow(z),players=length(unique(z$player_id)),
   fva_pearson=cor_safe(z$fva_prior,z$fva_next),xg_pearson=cor_safe(z$xg_prior,z$xg_next))
}
write.csv(do.call(rbind,sens),file.path(opt$output,'complete_xgot_sensitivity.csv'),row.names=FALSE)
# Plot with N annotations; NA when sample too small.
png(file.path(opt$output,'repeatability_by_shot_volume.png'),width=1400,height=900,res=150)
par(mar=c(5,5,4,2));plot(NA,xlim=c(27,103),ylim=c(-.45,1),xlab='Minimum non-penalty shots in EACH consecutive season',ylab='Pearson correlation (season t vs t+1)',main='Is finishing execution more repeatable for high-volume shooters?')
abline(h=0,col='gray65',lty=2);cols<-c(fva='#c34e5b',xg='#148778')
for(metric in c('fva','xg')) {
 z<-r[r$metric==metric,];lines(z$threshold,z$pearson,type='b',pch=19,lwd=2,col=cols[metric])
 if(metric=='fva') for(i in seq_len(nrow(z))) if(is.finite(z$pearson[i]))
  text(z$threshold[i],z$pearson[i],paste0('n=',z$pairs[i]),pos=if(z$pearson[i]<.65) 3 else 1,cex=.8,col=cols[metric])
}
legend('topright',legend=c('Finishing execution: FVA/shot','Shot quality: xG/shot'),col=cols,lwd=2,pch=19,bty='n')
dev.off()
writeLines(c('Under Pressure: shot-volume repeatability','',
 'Each threshold requires the minimum number of non-penalty shots in BOTH adjacent seasons.',
 'Only truly consecutive seasons are compared. A player may contribute multiple pairs.',
 'Pearson and Spearman are descriptive, not estimates of latent skill reliability.',
 'Confidence intervals use player-cluster bootstrap, accounting for repeated player pairs.',
 'Complete-xGOT sensitivity restricts both seasons to zero missing xGOT records.',
 'Small high-volume cohorts yield uncertain correlations; inspect sample sizes and CIs.',
 'Volume threshold selects high-volume players and can alter player/role mix; not a causal effect.',
 'Shot count is a noisy proxy for precision; models and missing data may affect FVA.',
 'Do not interpret zero correlation as proof that finishing skill does not exist or cannot be taught.'),
 file.path(opt$output,'README.txt'))
print(r);cat('\nSaved to: ',normalizePath(opt$output),'\n',sep='')
