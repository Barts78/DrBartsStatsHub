#!/usr/bin/env Rscript
# Five-season finishing analysis; base R + jsonlite.
# Usage: Rscript UP_Finishing_Five_Season_Longitudinal.R --input-dir . --output finishing_five_seasons --min-shots 30 --min-goals 5
args <- commandArgs(trailingOnly=TRUE)
arg <- function(key,default) {i<-match(paste0('--',key),args);if(is.na(i)) default else {if(i==length(args))stop('Missing value for --',key);args[i+1L]}}
indir<-arg('input-dir','.');out<-arg('output','finishing_five_seasons')
minshots<-as.integer(arg('min-shots','30'));mingoals<-as.integer(arg('min-goals','5'))
if(!requireNamespace('jsonlite',quietly=TRUE))stop('Install jsonlite: install.packages("jsonlite")')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
seasons<-sprintf('%d-%d',2021:2025,2022:2026)
read_season<-function(s) {
 f<-file.path(indir,paste0('EPL_',s,'_shots.csv'))
 if(!file.exists(f))stop('Missing: ',f)
 d<-read.csv(f,stringsAsFactors=FALSE,na.strings=c('','NA','null'))
 needed<-c('match_id','player_id','player_name','is_own_goal','situation','xg','xgot','raw_shot_json')
 if(length(setdiff(needed,names(d))))stop('Missing columns in ',f,': ',paste(setdiff(needed,names(d)),collapse=', '))
 raw<-lapply(d$raw_shot_json,function(x)tryCatch(jsonlite::fromJSON(x),error=function(e)list()))
 d$event<-vapply(raw,function(x)as.character(x$eventType %||% ''),character(1))
 d$blocked<-vapply(raw,function(x)isTRUE(x$isBlocked),logical(1))
 d$on_target<-vapply(raw,function(x)isTRUE(x$isOnTarget),logical(1))
 d$goal<-as.integer(d$event=='Goal')
 d$season<-s
 d<-d[!(tolower(as.character(d$is_own_goal)) %in% c('true','1','t')) & !is.na(d$situation) & d$situation!='Penalty' & !is.na(d$xg),,drop=FALSE]
 d$missing_ontarget<-d$on_target & !d$blocked & is.na(d$xgot)
 d$xgot_effective<-ifelse(d$blocked | !d$on_target,0,d$xgot)
 d$fva<-d$xgot_effective-d$xg
 d$fva_valid<-!is.na(d$fva)
 d
}
`%||%`<-function(a,b)if(is.null(a)||length(a)==0)b else a
all<-do.call(rbind,lapply(seasons,read_season));rownames(all)<-NULL
quality<-do.call(rbind,lapply(split(all,all$season),function(d)data.frame(season=d$season[1],matches=length(unique(d$match_id)),nonpenalty_shots=nrow(d),goals=sum(d$goal),xg=sum(d$xg),missing_ontarget=sum(d$missing_ontarget),fva_usable_shots=sum(d$fva_valid),fva_per_usable_shot=mean(d$fva,na.rm=TRUE))))
write.csv(quality,file.path(out,'season_quality.csv'),row.names=FALSE)
keys<-interaction(all$season,all$player_id,drop=TRUE)
players<-do.call(rbind,lapply(split(all,keys),function(d){valid<-sum(d$fva_valid);data.frame(season=d$season[1],player_id=d$player_id[1],player_name=tail(d$player_name,1),shots=nrow(d),goals=sum(d$goal),xg=sum(d$xg),xg_per_shot=mean(d$xg),fva_total=sum(d$fva,na.rm=TRUE),fva_usable_shots=valid,fva_missing_shots=sum(!d$fva_valid),fva_per_shot=if(valid>0)mean(d$fva,na.rm=TRUE)else NA_real_)}))
rownames(players)<-NULL
write.csv(players,file.path(out,'player_seasons.csv'),row.names=FALSE)
make_cohort<-function(label,condition){
 z<-players[condition & !is.na(players$fva_per_shot),,drop=FALSE]
 ids<-names(which(table(z$player_id)==length(seasons)))
 z<-z[as.character(z$player_id) %in% ids,,drop=FALSE]
 write.csv(z,file.path(out,paste0(label,'_five_season_player_seasons.csv')),row.names=FALSE)
 if(!nrow(z)){message(label,': no players qualify');return(invisible(NULL))}
 # Keep one name per player (avoid reshape warnings when names change).
 wide<-reshape(z[,c('player_id','season','fva_per_shot')],idvar='player_id',timevar='season',direction='wide')
 names_map<-tapply(z$player_name,z$player_id,function(x)tail(x,1))
 wide$player_name<-unname(names_map[as.character(wide$player_id)])
 cols<-paste0('fva_per_shot.',seasons)
 wide$change_2021_to_2025<-wide[[cols[5]]]-wide[[cols[1]]]
 wide<-wide[order(wide$change_2021_to_2025),c('player_id','player_name',cols,'change_2021_to_2025')]
 write.csv(wide,file.path(out,paste0(label,'_five_season_changes.csv')),row.names=FALSE)
 png(file.path(out,paste0(label,'_diverging.png')),width=1500,height=max(800,120*nrow(wide)+240),res=150)
 par(mar=c(5,12,4,2));v<-wide$change_2021_to_2025
 barplot(v,names.arg=wide$player_name,horiz=TRUE,las=1,col=ifelse(v<0,'#bf4854','#238a77'),border=NA,main=paste('Five-season FVA change:',label),xlab='FVA per shot change (2025/26 minus 2021/22)');abline(v=0,lwd=1.3);dev.off()
 png(file.path(out,paste0(label,'_trajectories.png')),width=1500,height=900,res=150)
 mat<-as.matrix(wide[,cols]);matplot(t(mat),type='l',lty=1,col=adjustcolor('grey40',alpha.f=.35),xaxt='n',xlab='Season',ylab='FVA per usable shot',main=paste('Finishing trajectories:',label))
 axis(1,at=1:5,labels=gsub('-','/',seasons));abline(h=0,col='grey50');lines(1:5,apply(mat,2,median),col='#126e7a',lwd=4,type='b',pch=16);legend('topright','Cohort median',col='#126e7a',lwd=3,bty='n');dev.off()
 message(label,': ',nrow(wide),' players, ',sum(v>0),' improved, ',sum(v<0),' declined')
 invisible(wide)
}
make_cohort('shots30',players$shots>=minshots)
make_cohort('goals5',players$shots>=minshots & players$goals>=mingoals)
print(quality)
cat('\nImportant: five-season cohorts require the thresholds in EACH season.\n')
cat('FVA uses confirmed misses/blocks as zero xGOT and excludes unexpectedly missing on-target xGOT.\n')
cat('Goal threshold is outcome-selected, not a direct position filter; compare with shot-volume cohort.\n')
