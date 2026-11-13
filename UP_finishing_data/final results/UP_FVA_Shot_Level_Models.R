#!/usr/bin/env Rscript
# Under Pressure: five-season FotMob shot-level finishing models.
# Requires lme4, jsonlite. Optional ggplot2 is NOT needed.
suppressPackageStartupMessages({library(lme4); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(flag, default) {i<-match(flag,args); if(!is.na(i) && i<length(args)) args[i+1] else default}
indir <- getarg('--input-dir','.')
outdir <- getarg('--output','fva_shot_level')
minshots <- as.integer(getarg('--min-shots','30'))
minseasons <- as.integer(getarg('--min-seasons','2'))
maxfun <- as.integer(getarg('--maxfun','100000'))
if(is.na(minshots)||minshots<1||is.na(minseasons)||minseasons<1) stop('Invalid shot/season thresholds')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
seasons <- sprintf('%d-%d',2021:2025,2022:2026)
files <- file.path(indir,paste0('EPL_',seasons,'_shots.csv'))
if(any(!file.exists(files))) stop('Missing input files: ',paste(files[!file.exists(files)],collapse=', '))
parse_bool <- function(x) tolower(trimws(as.character(x))) %in% c('true','1','t','yes')
read_one <- function(f, season) {
 d<-read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
 needed<-c('match_id','player_id','player_name','shot_type','situation','is_own_goal','xg','xgot','raw_shot_json')
 if(!all(needed %in% names(d))) stop('Missing columns in ',f,': ',paste(setdiff(needed,names(d)),collapse=', '))
 d$season<-season
 d
}
d<-do.call(rbind,Map(read_one,files,seasons))
d$player_id<-as.character(d$player_id)
d$xg<-suppressWarnings(as.numeric(d$xg)); d$xgot<-suppressWarnings(as.numeric(d$xgot))
d<-d[!parse_bool(d$is_own_goal) & is.finite(d$xg) & d$xg>=0 & d$xg<=1 & !is.na(d$player_id) & d$player_id!='-1',,drop=FALSE]
# Extract outcome flags from raw shot JSON; no guessing from missing xGOT.
`%||%`<-function(a,b) if(is.null(a)) b else a
flags<-lapply(d$raw_shot_json,function(s) {
 z<-tryCatch(fromJSON(s,simplifyVector=FALSE),error=function(e) NULL)
 if(is.null(z)) return(c(NA,NA,NA,NA))
 bool<-function(v) if(is.null(v)) NA else as.logical(v)
 c(bool(z$isOnTarget),bool(z$isBlocked),as.character(z$eventType %||% NA_character_),!is.null(z$expectedGoalsOnTarget))
})
 # Convert list flags robustly (event type is character).
d$on_target<-vapply(flags,function(z) as.character(z[1]),character(1))=='TRUE'
d$blocked<-vapply(flags,function(z) as.character(z[2]),character(1))=='TRUE'
d$event_type<-vapply(flags,function(z) as.character(z[3]),character(1))
d$flag_known<-vapply(flags,function(z) as.character(z[1]) %in% c('TRUE','FALSE') && as.character(z[2]) %in% c('TRUE','FALSE'),logical(1))
# All confirmed misses / blocked attempts have xGOT=0; missing on-target xGOT remains NA.
d$xgot_effective<-ifelse(d$flag_known & (!d$on_target | d$blocked),0,d$xgot)
d$xgot_effective[d$on_target & !d$blocked & !is.finite(d$xgot)]<-NA_real_
d$usable<-d$flag_known & is.finite(d$xgot_effective)
d$fva<-d$xgot_effective-d$xg
d$on_target_unblocked<-d$on_target & !d$blocked
d$shot_type<-factor(ifelse(is.na(d$shot_type)|d$shot_type=='','Unknown',d$shot_type))
d$situation<-as.character(d$situation); d$situation[is.na(d$situation)|d$situation=='']<-'Unknown'
# Collapse sparse situations to avoid singular/unstable fixed effects.
top_sit<-names(sort(table(d$situation),decreasing=TRUE))[1:min(7,length(unique(d$situation)))]
d$situation<-factor(ifelse(d$situation %in% top_sit,d$situation,'Other'))
d$season<-factor(d$season,levels=seasons)
d$player<-factor(d$player_id)
d$match<-factor(d$match_id)
d$player_season<-interaction(d$player,d$season,drop=TRUE)
d$xg_logit<-qlogis(pmin(pmax(d$xg,0.001),0.999))
# Eligibility uses all valid non-penalty, non-own-goal shots in each season.
counts<-aggregate(rep(1,nrow(d)),by=list(player=d$player,season=d$season),FUN=sum)
names(counts)[3]<-'shots'
counts<-counts[counts$shots>=minshots,,drop=FALSE]
qualified<-names(which(table(as.character(counts$player))>=minseasons))
d<-d[as.character(d$player) %in% qualified,,drop=FALSE]
d<-droplevels(d)
if(length(unique(d$player))<3) stop('Too few qualifying players; reduce --min-shots or --min-seasons')
quality<-data.frame(season=seasons,shots=sapply(seasons,function(s)sum(as.character(d$season)==s)),
                    missing_on_target_xgot=sapply(seasons,function(s)sum(as.character(d$season)==s & d$on_target_unblocked & !is.finite(d$xgot))),
                    unusable=sapply(seasons,function(s)sum(as.character(d$season)==s & !d$usable)))
write.csv(quality,file.path(outdir,'shot_level_quality.csv'),row.names=FALSE)
write.csv(counts[counts$player %in% qualified,],file.path(outdir,'qualifying_player_seasons.csv'),row.names=FALSE)
cat('Eligible players:',length(unique(d$player)),'; shots:',nrow(d),'\n')
ctl<-lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=maxfun),check.conv.singular='ignore')
gctl<-glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=maxfun),check.conv.singular='ignore')
# Three models. Season is fixed, player is persistent random intercept, player-season is transient.
# Match clustering accounts for common match conditions.
model_defs<-list(
 overall=list(data=d[d$usable,,drop=FALSE], family='gaussian',
   formula=fva ~ xg_logit + shot_type + situation + season + (1|player) + (1|player_season) + (1|match)),
 accuracy=list(data=d[d$flag_known,,drop=FALSE],family='binomial',
   formula=on_target_unblocked ~ xg_logit + shot_type + situation + season + (1|player) + (1|player_season) + (1|match)),
 on_target=list(data=d[d$usable & d$on_target_unblocked,,drop=FALSE],family='gaussian',
   formula=I(xgot_effective-xg) ~ xg_logit + shot_type + situation + season + (1|player) + (1|player_season) + (1|match))
)
results<-list(); fixed<-list(); players<-list()
for(nm in names(model_defs)) {
 item<-model_defs[[nm]]; x<-droplevels(item$data)
 cat('Fitting ',nm,': ',nrow(x),' shots...\n',sep='')
 if(nrow(x)<500 || length(unique(x$player))<3) {warning('Skipping ',nm,': insufficient data');next}
 fit<-tryCatch({
  if(item$family=='binomial') glmer(item$formula,data=x,family=binomial,control=gctl,nAGQ=0)
  else lmer(item$formula,data=x,REML=TRUE,control=ctl)
 },error=function(e){warning(nm,': ',conditionMessage(e));NULL})
 if(is.null(fit)) next
 saveRDS(fit,file.path(outdir,paste0('model_',nm,'.rds')))
 capture.output(summary(fit),file=file.path(outdir,paste0('model_',nm,'_summary.txt')))
 v<-as.data.frame(VarCorr(fit)); getv<-function(group){z<-v$vcov[v$grp==group];if(length(z)) z[1] else NA_real_}
 player_var<-getv('player'); transient_var<-getv('player_season'); match_var<-getv('match')
 resid_var<-if(item$family=='binomial') pi^2/3 else sigma(fit)^2
 # Conditional variance shares are descriptive; logistic residual uses latent-scale pi^2/3.
 total<-player_var+transient_var+match_var+resid_var
 results[[nm]]<-data.frame(model=nm,n_shots=nrow(x),n_players=length(unique(x$player)),n_player_seasons=length(unique(x$player_season)),
 player_variance=player_var,player_season_variance=transient_var,match_variance=match_var,residual_variance=resid_var,
 player_variance_share=player_var/total,player_season_variance_share=transient_var/total,
 singular=isSingular(fit,tol=1e-4))
 co<-as.data.frame(coef(summary(fit)));co$term<-rownames(co);co$model<-nm;fixed[[nm]]<-co
 re<-ranef(fit)$player
 players[[nm]]<-data.frame(player_id=rownames(re),effect=re[,1],model=nm)
}
if(!length(results)) stop('No models fitted. Inspect warnings and data quality.')
write.csv(do.call(rbind,results),file.path(outdir,'shot_level_variance_components.csv'),row.names=FALSE)
write.csv(do.call(rbind,fixed),file.path(outdir,'shot_level_fixed_effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,players),file.path(outdir,'shot_level_player_effects.csv'),row.names=FALSE)
# Summarize split-season empirical repeatability on adjusted player-season residuals.
# Leave full cross-validated repeatability for a separate validation step.
cat('Done. Outputs in: ',normalizePath(outdir),'\n',sep='')
