#!/usr/bin/env Rscript
# Under Pressure: longitudinal finishing execution and shot selection
# Input: player_seasons.csv from UP_Finishing_Five_Season_Longitudinal.R
# Base R only. Optional columns: age, position (not required).
# Example:
# Rscript UP_Finishing_Longitudinal_Models.R --input finishing_five_seasons/player_seasons.csv --output finishing_models --min-shots 30 --min-seasons 2 --bootstrap 500

args <- commandArgs(trailingOnly=TRUE)
opt <- list(input='player_seasons.csv', output='finishing_models', `min-shots`='30', `min-seasons`='2', bootstrap='500', seed='2026')
if (length(args) %% 2 != 0) stop('Arguments must be --key value pairs')
if (length(args)) for (i in seq(1, length(args), by=2)) {
  key <- sub('^--', '', args[i]); if (!key %in% names(opt)) stop('Unknown option: ', args[i]); opt[[key]] <- args[i+1]
}
min_shots <- as.integer(opt[['min-shots']]); min_seasons <- as.integer(opt[['min-seasons']]); B <- as.integer(opt$bootstrap)
if (is.na(min_shots) || min_shots < 1 || is.na(min_seasons) || min_seasons < 2 || is.na(B) || B < 0) stop('Invalid numeric argument')
set.seed(as.integer(opt$seed)); dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(opt$input, stringsAsFactors=FALSE, check.names=FALSE)
required <- c('season','player_id','player_name','shots','goals','xg','xg_per_shot','fva_total','fva_usable_shots','fva_missing_shots','fva_per_shot')
missing <- setdiff(required,names(d)); if(length(missing)) stop('Missing columns: ',paste(missing,collapse=', '))
for (nm in c('shots','goals','xg','xg_per_shot','fva_total','fva_usable_shots','fva_missing_shots','fva_per_shot')) d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
d$season_year <- as.integer(sub('^([0-9]{4}).*','\\1',d$season))
d$player_id <- as.character(d$player_id)
d <- d[is.finite(d$season_year) & !is.na(d$player_id) & d$player_id != '' & d$player_id != '-1', ]
# Unknown player_id -1 can merge different players, so exclude it.
if (anyDuplicated(d[c('player_id','season')])) stop('Duplicate player_id/season rows; resolve before modelling.')
d <- d[order(d$player_id,d$season_year), ]
d$missing_rate <- ifelse(d$shots>0,d$fva_missing_shots/d$shots,NA_real_)
d$valid_fva <- is.finite(d$fva_per_shot) & d$fva_usable_shots>0
# Core cohort: qualifying seasons, not all seasons of each qualifying player.
q <- d[d$shots>=min_shots & d$valid_fva, ]
nseasons <- table(q$player_id); ids <- names(nseasons[nseasons>=min_seasons]); q <- q[q$player_id %in% ids, ]
q <- q[order(q$player_id,q$season_year), ]; q$season_index <- q$season_year-min(q$season_year)
write.csv(q,file.path(opt$output,'qualifying_player_seasons.csv'),row.names=FALSE)

# Weighted player fixed-effects slope using within-player demeaned season and outcome.
# Weights = number of usable shots; this is a precision proxy, not a full variance model.
within_slope <- function(z,y,w='fva_usable_shots') {
 z <- z[is.finite(z[[y]]) & is.finite(z[[w]]) & z[[w]]>0, ]; if(nrow(z)<2) return(NA_real_)
 z$x <- z$season_year
 z$xw <- ave(seq_len(nrow(z)),z$player_id,FUN=function(ii) weighted.mean(z$x[ii],z[[w]][ii]))
 z$yw <- ave(seq_len(nrow(z)),z$player_id,FUN=function(ii) weighted.mean(z[[y]][ii],z[[w]][ii]))
 dx <- z$x-z$xw; dy <- z[[y]]-z$yw
 denom <- sum(z[[w]]*dx^2); if(denom<=0) return(NA_real_)
 sum(z[[w]]*dx*dy)/denom
}
# Note: the helper uses ave over row indices, preserving weights within each player.
cluster_boot <- function(z,y,w='fva_usable_shots',B=500) {
 est <- within_slope(z,y,w)
 if(B==0 || !is.finite(est)) return(c(estimate=est,ci_low=NA_real_,ci_high=NA_real_,boot_n=0))
 groups <- split(z,z$player_id); keys <- names(groups); draws <- rep(NA_real_,B)
 for(b in seq_len(B)) {
   sample_ids <- sample(keys,length(keys),replace=TRUE)
   sampled <- lapply(seq_along(sample_ids),function(i){g<-groups[[sample_ids[i]]];g$player_id<-paste0('bootstrap_',i);g})
   draws[b] <- within_slope(do.call(rbind,sampled),y,w)
 }
 valid <- draws[is.finite(draws)]; ci <- if(length(valid)>10) unname(quantile(valid,c(.025,.975))) else c(NA_real_,NA_real_)
 c(estimate=est,ci_low=ci[1],ci_high=ci[2],boot_n=length(valid))
}

# A second outcome is shot selection, xG per shot; weight by all shots.
fit_results <- list(); cohorts <- list(`2plus`=q)
q3 <- d[d$shots>=min_shots & d$valid_fva, ]; n3 <- table(q3$player_id); q3 <- q3[q3$player_id %in% names(n3[n3>=3]), ]
cohorts[['3plus']] <- q3
for(cohort_name in names(cohorts)) {
 z <- cohorts[[cohort_name]]; if(nrow(z)==0) next
 for(metric in c('fva_per_shot','xg_per_shot')) {
   w <- if(metric=='fva_per_shot') 'fva_usable_shots' else 'shots'
   r <- cluster_boot(z,metric,w,B)
   fit_results[[length(fit_results)+1]] <- data.frame(cohort=cohort_name,metric=metric,players=length(unique(z$player_id)),player_seasons=nrow(z),estimate_per_year=r[1],ci_low=r[2],ci_high=r[3],bootstrap_draws=r[4])
 }
}
# Sensitivity: no missing xGOT in a player-season (restricts older seasons substantially).
q_complete <- q[q$fva_missing_shots==0, ]; nc <- table(q_complete$player_id); q_complete <- q_complete[q_complete$player_id %in% names(nc[nc>=min_seasons]), ]
if(nrow(q_complete)>0) {
 r <- cluster_boot(q_complete,'fva_per_shot',B=B)
 fit_results[[length(fit_results)+1]] <- data.frame(cohort='complete_xgot_only',metric='fva_per_shot',players=length(unique(q_complete$player_id)),player_seasons=nrow(q_complete),estimate_per_year=r[1],ci_low=r[2],ci_high=r[3],bootstrap_draws=r[4])
}
fits <- if(length(fit_results)) do.call(rbind,fit_results) else data.frame()
write.csv(fits,file.path(opt$output,'within_player_trends.csv'),row.names=FALSE)

# Consecutive-season correlations and within-player changes; distinguish adjacent years from gaps.
lagrows <- list()
for(id in unique(q$player_id)) {
 z <- q[q$player_id==id, ]; z <- z[order(z$season_year), ]
 if(nrow(z)<2) next
 for(i in seq_len(nrow(z)-1)) if(z$season_year[i+1]-z$season_year[i]==1) {
   lagrows[[length(lagrows)+1]] <- data.frame(player_id=id,player_name=z$player_name[i+1],from_season=z$season[i],to_season=z$season[i+1],fva_prior=z$fva_per_shot[i],fva_next=z$fva_per_shot[i+1],delta_fva=z$fva_per_shot[i+1]-z$fva_per_shot[i],xg_prior=z$xg_per_shot[i],xg_next=z$xg_per_shot[i+1],delta_xg=z$xg_per_shot[i+1]-z$xg_per_shot[i],shots_prior=z$shots[i],shots_next=z$shots[i+1])
 }
}
lags <- if(length(lagrows)) do.call(rbind,lagrows) else data.frame()
write.csv(lags,file.path(opt$output,'adjacent_season_pairs.csv'),row.names=FALSE)
cor_out <- data.frame()
if(nrow(lags)>2) {
 for(metric in c('fva','xg')) {
   a<-lags[[paste0(metric,'_prior')]]; b<-lags[[paste0(metric,'_next')]]
   cor_out<-rbind(cor_out,data.frame(metric=metric,pairs=sum(complete.cases(a,b)),pearson=cor(a,b,use='complete.obs'),spearman=cor(a,b,method='spearman',use='complete.obs')))
 }
}
write.csv(cor_out,file.path(opt$output,'year_to_year_correlations.csv'),row.names=FALSE)

# Player-specific slopes: descriptive only, not individual significance tests.
player_results <- list()
for(id in unique(q$player_id)) {
 z<-q[q$player_id==id, ]; if(nrow(z)<2) next
 player_results[[length(player_results)+1]] <- data.frame(player_id=id,player_name=tail(z$player_name,1),seasons=nrow(z),total_shots=sum(z$shots),first_season=z$season[1],last_season=tail(z$season,1),fva_slope_per_year=within_slope(z,'fva_per_shot'),xg_slope_per_year=within_slope(z,'xg_per_shot','shots'),first_fva=z$fva_per_shot[1],last_fva=tail(z$fva_per_shot,1))
}
players <- if(length(player_results)) do.call(rbind,player_results) else data.frame()
if(nrow(players)) players<-players[order(players$fva_slope_per_year,decreasing=TRUE,na.last=TRUE), ]
write.csv(players,file.path(opt$output,'player_trends_descriptive.csv'),row.names=FALSE)

# Weighted cohort trend by calendar season, explicitly NOT an age curve.
season_stats <- do.call(rbind,lapply(split(q,q$season),function(z) data.frame(season=z$season[1],players=length(unique(z$player_id)),shots=sum(z$shots),weighted_fva=weighted.mean(z$fva_per_shot,z$fva_usable_shots),weighted_xg_per_shot=weighted.mean(z$xg_per_shot,z$shots),missing_xgot=sum(z$fva_missing_shots))))
season_stats <- season_stats[order(season_stats$season), ]; write.csv(season_stats,file.path(opt$output,'cohort_season_summary.csv'),row.names=FALSE)

# Optional age analysis only if age was supplied. This is observational, not causal.
if('age' %in% names(q)) {
 q$age<-suppressWarnings(as.numeric(q$age)); qa<-q[is.finite(q$age), ]
 if(nrow(qa)>5 && length(unique(qa$player_id))>2) {
   agefit<-cluster_boot(transform(qa,season_year=age),'fva_per_shot',B=B)
   write.csv(data.frame(age_slope=agefit[1],ci_low=agefit[2],ci_high=agefit[3],players=length(unique(qa$player_id))),file.path(opt$output,'age_trend_optional.csv'),row.names=FALSE)
 }
}

# Two figures, using base graphics for portability.
if(nrow(season_stats)>0 && any(is.finite(season_stats$weighted_fva))) {
png(file.path(opt$output,'cohort_trends.png'),width=1300,height=800,res=140)
par(mar=c(5,5,3,5)); plot(seq_len(nrow(season_stats)),season_stats$weighted_fva,type='b',pch=19,lwd=3,col='#087a80',xaxt='n',xlab='Season',ylab='Shot-weighted FVA per shot',main='Finishing execution: qualifying players (2+ seasons)'); axis(1,at=seq_len(nrow(season_stats)),labels=season_stats$season);abline(h=0,col='gray60',lty=2);dev.off()
} else warning('No finite cohort values: chart skipped')
finite_players <- players[is.finite(players$fva_slope_per_year), , drop=FALSE]
if(nrow(finite_players)>0) {
 top<-head(finite_players[order(abs(finite_players$fva_slope_per_year),decreasing=TRUE), ],25)

 top<-top[order(top$fva_slope_per_year), ]; png(file.path(opt$output,'player_slope_examples.png'),width=1400,height=1000,res=140)
 par(mar=c(5,13,3,2));barplot(top$fva_slope_per_year,names.arg=top$player_name,horiz=TRUE,las=1,col=ifelse(top$fva_slope_per_year>=0,'#208c79','#c64e5c'),border=NA,xlab='Within-player FVA per shot change per calendar year',main='Largest observed player slopes (descriptive)');abline(v=0,lwd=1.5);dev.off()
} else warning('No finite player slopes: chart skipped')

readme <- c('UNDER PRESSURE: LONGITUDINAL FINISHING MODELS','',
 paste('Input:',normalizePath(opt$input,mustWork=FALSE)),
 paste('Min shots per qualifying season:',min_shots),paste('Min qualifying seasons per player:',min_seasons),
 paste('Main cohort:',length(unique(q$player_id)),'players;',nrow(q),'player-seasons'),
 '', 'DEFINITION: FVA per shot = sum(xGOT - xG) / number of usable non-penalty shots.',
 'This script uses the FVA values already computed by the earlier pipeline; it does not reconstruct missing shot-level values.',
 'Unknown player ID -1 excluded to prevent merging unrelated players.',
 'Within-player slope: weighted fixed-effects trend across calendar seasons.',
 'Uncertainty: player-cluster bootstrap percentile intervals (resample players, preserving each player trajectory).',
 'Shot volume is a precision proxy; heterogeneity in shot context and model calibration remains unadjusted.',
 'Year-to-year correlations pool adjacent qualifying seasons and are descriptive, not reliability-corrected.',
 'Complete-xGOT sensitivity excludes any player-season with missing xGOT; may select non-representative data.',
 'The goals threshold is intentionally not used to define the main cohort (selection on outcome).',
 'No age or position inference is made unless a trustworthy age column is supplied.',
 'Calendar-year trends can reflect aging, tactical changes, team changes, selection and measurement noise.',
 'Confidence intervals excluding zero are not evidence that coaching caused change.',
 '', 'OUTPUTS: within_player_trends.csv; adjacent_season_pairs.csv; year_to_year_correlations.csv;',
 'player_trends_descriptive.csv; qualifying_player_seasons.csv; cohort_season_summary.csv; figures.',
 '', 'NEXT: add position, date of birth / age, and shot-level covariates for a fuller hierarchical model.')
writeLines(readme,file.path(opt$output,'README.txt'))
cat(paste(readme[1:5],collapse='\n'),'\n\n');print(fits);print(cor_out);cat('\nSaved outputs to ',normalizePath(opt$output),'\n',sep='')
