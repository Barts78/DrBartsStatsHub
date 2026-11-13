#!/usr/bin/env Rscript
# Visualise finishing value added (FVA) by player and season.
# Usage: Rscript UP_FVA_Player_Season_Visualisations.R --input player_seasons.csv --output fva_player_visuals --min-shots 30 --min-seasons 2 --top 25
args <- commandArgs(trailingOnly=TRUE)
opt <- list(input='player_seasons.csv',output='fva_player_visuals',`min-shots`='30',`min-seasons`='2',top='25')
if(length(args)>0) {
 for(i in seq(1,length(args),by=2)) {
  if(i==length(args) || !startsWith(args[i],'--')) stop('Arguments must be --name value pairs')
  key <- substring(args[i],3); if(!key %in% names(opt)) stop('Unknown argument: ',key)
  opt[[key]] <- args[i+1]
 }
}
min_shots <- as.integer(opt[['min-shots']]); min_seasons <- as.integer(opt[['min-seasons']]); top <- as.integer(opt$top)
if(any(is.na(c(min_shots,min_seasons,top))) || min_shots<1 || min_seasons<2 || top<1) stop('Invalid thresholds')
dir.create(opt$output,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(opt$input,stringsAsFactors=FALSE,check.names=FALSE)
required <- c('season','player_id','player_name','shots','fva_per_shot')
if(!all(required %in% names(d))) stop('Missing columns: ',paste(setdiff(required,names(d)),collapse=', '))
d$shots <- as.numeric(d$shots);d$fva_per_shot <- as.numeric(d$fva_per_shot)
d$year <- suppressWarnings(as.integer(substr(d$season,1,4)))
d <- d[is.finite(d$fva_per_shot) & is.finite(d$year) & d$shots>=min_shots & !is.na(d$player_id),,drop=FALSE]
d <- d[order(d$player_id,d$year),,drop=FALSE]
# Use FotMob player ID, not player name, to avoid conflating players.
counts <- table(d$player_id)
d <- d[as.character(d$player_id) %in% names(counts[counts>=min_seasons]),,drop=FALSE]
if(nrow(d)==0) stop('No players qualify; lower --min-shots or --min-seasons')
ids <- unique(d$player_id)
first_last <- list(); changes <- list()
for(id in ids) {
 z <- d[d$player_id==id,,drop=FALSE];z<-z[order(z$year),,drop=FALSE]
 a <- z[1,];b <- z[nrow(z),]
 first_last[[length(first_last)+1]] <- data.frame(player_id=id,player_name=b$player_name,first_season=a$season,last_season=b$season,first_fva=a$fva_per_shot,last_fva=b$fva_per_shot,net_change=b$fva_per_shot-a$fva_per_shot,seasons=nrow(z),total_shots=sum(z$shots))
 if(nrow(z)>1) for(j in 2:nrow(z)) {
  # Only adjacent *calendar* seasons; no misleading multi-year leaps.
  if(z$year[j]-z$year[j-1]!=1) next
  changes[[length(changes)+1]] <- data.frame(player_id=id,player_name=z$player_name[j],from_season=z$season[j-1],to_season=z$season[j],from_fva=z$fva_per_shot[j-1],to_fva=z$fva_per_shot[j],change=z$fva_per_shot[j]-z$fva_per_shot[j-1],from_shots=z$shots[j-1],to_shots=z$shots[j])
 }
}
fl <- do.call(rbind,first_last)
ch <- if(length(changes)) do.call(rbind,changes) else data.frame()
write.csv(fl,file.path(opt$output,'first_last_changes.csv'),row.names=FALSE)
write.csv(d,file.path(opt$output,'qualifying_player_seasons.csv'),row.names=FALSE)
if(nrow(ch)) write.csv(ch,file.path(opt$output,'adjacent_season_changes.csv'),row.names=FALSE)
# Identify each player's biggest *absolute* adjacent-season swing; retain direction and labels.
if(nrow(ch)) {
 ch <- ch[order(ch$player_id,-abs(ch$change),ch$from_season),,drop=FALSE]
 biggest <- ch[!duplicated(ch$player_id),,drop=FALSE]
 write.csv(biggest,file.path(opt$output,'biggest_annual_change_by_player.csv'),row.names=FALSE)
} else biggest <- data.frame()
positive <- '#208c79'; negative <- '#c64e5c'; neutral <- '#747c87'
# Select biggest net changes for legible publication chart.
fl_plot <- head(fl[order(-abs(fl$net_change)),,drop=FALSE],top)
fl_plot <- fl_plot[order(fl_plot$net_change),,drop=FALSE]
png(file.path(opt$output,'first_last_dumbbell.png'),width=1900,height=max(850,85*nrow(fl_plot)+240),res=160)
par(mar=c(5,13,4,6));yr <- seq_len(nrow(fl_plot));
xrange <- range(c(fl_plot$first_fva,fl_plot$last_fva),finite=TRUE); pad <- max(diff(xrange)*0.17,0.01)
plot(NA,xlim=xrange+c(-pad,pad),ylim=c(.5,nrow(fl_plot)+.5),yaxt='n',xlab='FVA per usable shot (xGOT − xG)',ylab='',main='Finishing execution: first vs last qualifying season',bty='n')
abline(v=0,col='gray75',lty=2)
segments(fl_plot$first_fva,yr,fl_plot$last_fva,yr,col=ifelse(fl_plot$net_change>=0,positive,negative),lwd=3)
points(fl_plot$first_fva,yr,pch=21,bg='white',col='gray30',cex=1.5,lwd=2)
points(fl_plot$last_fva,yr,pch=19,col=ifelse(fl_plot$net_change>=0,positive,negative),cex=1.5)
axis(2,at=yr,labels=fl_plot$player_name,las=1,tick=FALSE,cex.axis=.85)
text(fl_plot$first_fva,yr+.22,paste0(substr(fl_plot$first_season,3,4),'/',substr(fl_plot$first_season,8,9)),cex=.55,col='gray35')
text(fl_plot$last_fva,yr-.22,paste0(substr(fl_plot$last_season,3,4),'/',substr(fl_plot$last_season,8,9)),cex=.55,col='gray35')
legend('topright',legend=c('First qualifying season','Last qualifying season'),pch=c(21,19),pt.bg=c('white',NA),col=c('gray30',positive),bty='n',cex=.8)
mtext('Players may have gaps between qualifying seasons; points are not necessarily consecutive years.',side=1,line=3.2,cex=.7,col='gray40')
dev.off()
if(nrow(biggest)) {
 b <- head(biggest[order(-abs(biggest$change)),,drop=FALSE],top)
 b <- b[order(b$change),,drop=FALSE]
 png(file.path(opt$output,'largest_adjacent_season_swings.png'),width=1900,height=max(850,80*nrow(b)+240),res=160)
 par(mar=c(5,13,4,4));yr<-seq_len(nrow(b));xmax<-max(abs(b$change))*1.28
 barplot(b$change,names.arg=b$player_name,horiz=TRUE,las=1,col=ifelse(b$change>=0,positive,negative),border=NA,xlim=c(-xmax,xmax),xlab='Change in FVA per usable shot',main='Largest single-season changes in finishing execution',cex.names=.85)
 abline(v=0,col='gray25',lwd=1.5)
 # barplot midpoints must be used for labels
 dev.off()
 # Redraw with exact bar midpoints and labels.
 png(file.path(opt$output,'largest_adjacent_season_swings_labelled.png'),width=1900,height=max(850,80*nrow(b)+240),res=160)
 par(mar=c(5,13,4,4))
 mids <- barplot(b$change,names.arg=b$player_name,horiz=TRUE,las=1,col=ifelse(b$change>=0,positive,negative),border=NA,xlim=c(-xmax,xmax),xlab='Change in FVA per usable shot',main='Biggest consecutive-season FVA swing per player',cex.names=.85)
 abline(v=0,col='gray25',lwd=1.3)
 labels <- paste0(substr(b$from_season,3,4),'/',substr(b$from_season,8,9),' → ',substr(b$to_season,3,4),'/',substr(b$to_season,8,9))
 text(ifelse(b$change>=0,b$change+0.002,b$change-0.002),mids,labels,pos=ifelse(b$change>=0,4,2),cex=.65,xpd=NA)
 dev.off()
}
cat('Qualified players:',length(ids),'| first-last pairs:',nrow(fl),'| adjacent season pairs:',nrow(ch),'\n')
cat('Outputs:',normalizePath(opt$output),'\n')
