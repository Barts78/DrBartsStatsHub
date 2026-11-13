#!/usr/bin/env Rscript
# Under Pressure: train on 2021/22-2024/25; predict 2025/26.
# Tests whether historical player identity adds out-of-season predictive value.
suppressPackageStartupMessages({library(lme4); library(jsonlite)})
a <- commandArgs(trailingOnly=TRUE)
arg <- function(k,default) {i<-match(k,a);if(!is.na(i)&&i<length(a)) a[i+1] else default}
indir<-arg('--input-dir','.'); outdir<-arg('--output','fva_holdout')
minshots<-as.integer(arg('--min-shots','30')); minseasons<-as.integer(arg('--min-seasons','2'))
boots<-as.integer(arg('--bootstrap','500')); set.seed(20261020)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
seasons<-sprintf('%d-%d',2021:2025,2022:2026)
fs<-file.path(indir,paste0('EPL_',seasons,'_shots.csv'))
if(any(!file.exists(fs))) stop('Missing files: ',paste(fs[!file.exists(fs)],collapse=', '))
d<-do.call(rbind,Map(function(f,s){x<-read.csv(f,stringsAsFactors=FALSE,check.names=FALSE);x$season<-s;x},fs,seasons))
required<-c('player_id','match_id','xg','xgot','situation','shot_type','is_own_goal','raw_shot_json')
if(!all(required %in% names(d))) stop('Missing: ',paste(setdiff(required,names(d)),collapse=', '))
truth<-function(x) tolower(trimws(as.character(x))) %in% c('true','1','t','yes')
d$player_id<-as.character(d$player_id);d$xg<-suppressWarnings(as.numeric(d$xg));d$xgot<-suppressWarnings(as.numeric(d$xgot))
d<-d[!truth(d$is_own_goal)&is.finite(d$xg)&d$xg>=0&d$xg<=1&!is.na(d$player_id)&d$player_id!='-1',,drop=FALSE]
d$situation<-trimws(as.character(d$situation));d$situation[is.na(d$situation)|d$situation=='']<-'Unknown'
# Exclude penalties BEFORE computing eligibility and model features.
d<-d[tolower(d$situation)!='penalty',,drop=FALSE]
getflags<-function(s){z<-tryCatch(fromJSON(s,simplifyVector=FALSE),error=function(e)NULL);if(is.null(z))return(c(NA,NA));b<-function(v)if(is.null(v))NA else as.logical(v);c(b(z$isOnTarget),b(z$isBlocked))}
f<-lapply(d$raw_shot_json,getflags)
d$flag_known<-vapply(f,function(x)length(x)==2&&all(!is.na(x)),logical(1))
d$on_target<-vapply(f,function(x)isTRUE(x[1]),logical(1));d$blocked<-vapply(f,function(x)isTRUE(x[2]),logical(1))
d$on_target_unblocked<-d$on_target & !d$blocked
d$xgot_eff<-ifelse(d$flag_known & !d$on_target_unblocked,0,d$xgot)
d$xgot_eff[d$on_target_unblocked & !is.finite(d$xgot)]<-NA_real_
d$usable<-d$flag_known & is.finite(d$xgot_eff)
d$fva<-d$xgot_eff-d$xg
d$xg_logit<-qlogis(pmin(pmax(d$xg,.001),.999))
d$shot_type<-as.character(d$shot_type);d$shot_type[is.na(d$shot_type)|d$shot_type=='']<-'Unknown'
# Eligibility is determined ONLY from the training seasons (no future leakage).
train0<-d[d$season!='2025-2026',,drop=FALSE]
ct<-aggregate(rep(1,nrow(train0)),list(player=train0$player_id,season=train0$season),sum)
qualified<-names(which(table(ct$player[ct$x>=minshots])>=minseasons))
d<-d[d$player_id %in% qualified,,drop=FALSE]
# Hold out 2025/26, including only historically observed players.
train<-d[d$season!='2025-2026',,drop=FALSE];test<-d[d$season=='2025-2026',,drop=FALSE]
if(nrow(test)<100)stop('Insufficient test shots')
# Collapse situation levels using training frequencies only.
top<-names(sort(table(train$situation),decreasing=TRUE))[seq_len(min(7,length(unique(train$situation))))]
for(nm in c('situation','shot_type')) {
 if(nm=='situation') {train[[nm]][!train[[nm]]%in%top]<-'Other';test[[nm]][!test[[nm]]%in%top]<-'Other'}
 lev<-unique(train[[nm]]);test[[nm]][!test[[nm]]%in%lev]<-'Other'
 if(!'Other'%in%lev)lev<-c(lev,'Other')
 train[[nm]]<-factor(train[[nm]],levels=lev);test[[nm]]<-factor(test[[nm]],levels=lev)
}
# No season fixed effect: unseen future season has no estimable coefficient.
# A training-only year trend can be added in sensitivity work.
train$player<-factor(train$player_id);test$player<-factor(test$player_id,levels=levels(train$player))
train$player_season<-interaction(train$player_id,train$season,drop=TRUE)
test$player_season<-interaction(test$player_id,test$season,drop=TRUE)
train$match<-factor(train$match_id);test$match<-factor(test$match_id)
ctl<-lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000),check.conv.singular='ignore')
gctl<-glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000),check.conv.singular='ignore')
models<-list(
 overall=list(train=train[train$usable,,drop=FALSE],test=test[test$usable,,drop=FALSE],formula=fva~xg_logit+shot_type+situation+(1|player)+(1|player_season)+(1|match),kind='gaussian'),
 accuracy=list(train=train[train$flag_known,,drop=FALSE],test=test[test$flag_known,,drop=FALSE],formula=on_target_unblocked~xg_logit+shot_type+situation+(1|player)+(1|player_season)+(1|match),kind='binomial'),
 on_target=list(train=train[train$usable&train$on_target_unblocked,,drop=FALSE],test=test[test$usable&test$on_target_unblocked,,drop=FALSE],formula=I(xgot_eff-xg)~xg_logit+shot_type+situation+(1|player)+(1|player_season)+(1|match),kind='gaussian')
)
rows<-list(); predictions<-list()
for(nm in names(models)){
 m<-models[[nm]];tr<-droplevels(m$train);te<-m$test
 # Retain factor levels from training for prediction; unseen grouping levels permitted.
 te$player<-factor(te$player_id,levels=levels(tr$player))
 te$shot_type<-factor(as.character(te$shot_type),levels=levels(tr$shot_type))
 te$situation<-factor(as.character(te$situation),levels=levels(tr$situation))
 if(anyNA(te$shot_type)|anyNA(te$situation))stop('Unseen fixed-effect level in test: ',nm)
 cat('Fitting ',nm,': train=',nrow(tr),' test=',nrow(te),'\n',sep='')
 fit<-if(m$kind=='binomial')glmer(m$formula,data=tr,family=binomial,control=gctl,nAGQ=0) else lmer(m$formula,data=tr,REML=TRUE,control=ctl)
 saveRDS(fit,file.path(outdir,paste0('training_model_',nm,'.rds')))
 # re.form=NA: ignore all player, player-season and match effects.
 # re.form=~(1|player): include ONLY persistent player effects learned from prior seasons.
 p0<-predict(fit,newdata=te,re.form=NA,type=if(m$kind=='binomial')'response' else 'response',allow.new.levels=TRUE)
 p1<-predict(fit,newdata=te,re.form=~(1|player),type=if(m$kind=='binomial')'response' else 'response',allow.new.levels=TRUE)
 y<-if(nm=='accuracy')as.numeric(te$on_target_unblocked) else te$fva
 if(nm=='accuracy'){
  eps<-1e-12;p0<-pmin(pmax(p0,eps),1-eps);p1<-pmin(pmax(p1,eps),1-eps)
  loss0<- -(y*log(p0)+(1-y)*log(1-p0));loss1<- -(y*log(p1)+(1-y)*log(1-p1));metric<-'log_loss'
 }else{loss0<-(y-p0)^2;loss1<-(y-p1)^2;metric<-'mean_squared_error'}
 # Positive improvement means player identity improves predictions.
 improvement<-loss0-loss1
 # Cluster bootstrap by player, preserving all shots for each sampled player.
 ids<-unique(te$player_id);ix<-split(seq_len(nrow(te)),te$player_id)
 boot<-rep(NA_real_,boots)
 if(boots>0)for(b in seq_len(boots)){sample_ids<-sample(ids,length(ids),replace=TRUE);sel<-unlist(ix[sample_ids],use.names=FALSE);boot[b]<-mean(improvement[sel])}
 ci<-if(boots>1)unname(quantile(boot,c(.025,.975),na.rm=TRUE)) else c(NA_real_,NA_real_)
 rows[[nm]]<-data.frame(model=nm,metric=metric,train_shots=nrow(tr),test_shots=nrow(te),test_players=length(ids),baseline_loss=mean(loss0),with_player_loss=mean(loss1),improvement=mean(improvement),ci_low=ci[1],ci_high=ci[2],singular=isSingular(fit),stringsAsFactors=FALSE)
 predictions[[nm]]<-data.frame(model=nm,player_id=te$player_id,season=te$season,observed=y,baseline=p0,with_player=p1,loss_improvement=improvement)
}
write.csv(do.call(rbind,rows),file.path(outdir,'heldout_model_comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,predictions),file.path(outdir,'heldout_shot_predictions.csv'),row.names=FALSE)
cat('Saved validation results in ',normalizePath(outdir),'\n')
