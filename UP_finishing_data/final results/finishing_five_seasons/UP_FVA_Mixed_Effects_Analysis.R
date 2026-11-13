#!/usr/bin/env Rscript
# Under Pressure | Persistent finishing skill: mixed-effects analysis
# Input: player_seasons.csv from five-season FotMob analysis.
# Install once: install.packages('lme4')
# Example: Rscript UP_FVA_Mixed_Effects_Analysis.R --input player_seasons.csv --output fva_mixed_models --min-shots 30 --min-seasons 2 --bootstrap 200
# Models player-season FVA/shot, not individual shots. Weights approximate precision.

args <- commandArgs(trailingOnly=TRUE)
opt <- list(input='player_seasons.csv', output='fva_mixed_models', `min-shots`='30', `min-seasons`='2', bootstrap='200', seed='2026')
if(length(args) %% 2L) stop('Use --option value argument pairs.')
if(length(args)) for(i in seq(1L,length(args),by=2L)) {
  k <- sub('^--','',args[i]); if(!k %in% names(opt)) stop('Unknown option: ',args[i]); opt[[k]] <- args[i+1L]
}
min_shots <- as.integer(opt[['min-shots']]); min_seasons <- as.integer(opt[['min-seasons']]); B <- as.integer(opt$bootstrap)
if(anyNA(c(min_shots,min_seasons,B)) || min_shots<1 || min_seasons<2 || B<0) stop('Invalid --min-shots, --min-seasons or --bootstrap.')
if(!requireNamespace('lme4',quietly=TRUE)) stop("Package 'lme4' is required. Run install.packages('lme4') in R.")
set.seed(as.integer(opt$seed))
dir.create(opt$output,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(opt$input,stringsAsFactors=FALSE,check.names=FALSE)
needed <- c('season','player_id','player_name','shots','xg_per_shot','fva_per_shot','fva_usable_shots','fva_missing_shots')
if(length(miss <- setdiff(needed,names(d)))) stop('Missing columns: ',paste(miss,collapse=', '))
for(k in c('shots','xg_per_shot','fva_per_shot','fva_usable_shots','fva_missing_shots')) d[[k]] <- suppressWarnings(as.numeric(d[[k]]))
d$player_id <- as.character(d$player_id)
d$season_year <- suppressWarnings(as.integer(substr(d$season,1,4)))
d <- d[is.finite(d$season_year) & !is.na(d$player_id) & d$player_id!='-1' & d$player_id!='' & d$shots>=min_shots & d$fva_usable_shots>0 & is.finite(d$fva_per_shot) & is.finite(d$xg_per_shot),,drop=FALSE]
if(anyDuplicated(d[c('player_id','season')])) stop('Duplicate player-season rows; resolve first.')
keep <- names(which(table(d$player_id)>=min_seasons))
d <- d[d$player_id %in% keep,,drop=FALSE]
if(nrow(d)<10 || length(unique(d$player_id))<5) stop('Not enough qualifying player-seasons for modelling.')
d$season_index <- d$season_year - min(d$season_year)
d$player_id <- factor(d$player_id)
d$season_factor <- factor(d$season)
write.csv(d,file.path(opt$output,'qualifying_player_seasons.csv'),row.names=FALSE)

# Weighted random-intercept models. Weights ~ shot count: working approximation to
# inverse sampling variance, not a known per-shot variance. ICC is conditional on
# the residual scale for the reference shot count (median weight).
fit_one <- function(z,metric,weight_col,tag) {
  z <- z[is.finite(z[[metric]]) & is.finite(z[[weight_col]]) & z[[weight_col]]>0,,drop=FALSE]
  z$y <- z[[metric]]
  z$w <- z[[weight_col]] / median(z[[weight_col]])
  z$player_id <- droplevels(factor(z$player_id))
  if(length(unique(z$player_id))<5 || nrow(z)<10) return(NULL)
  ctrl <- lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))
  # Season fixed effects absorb common league-wide changes and model calibration drift.
  mod <- tryCatch(lme4::lmer(y ~ season_factor + (1|player_id),data=z,weights=w,REML=TRUE,control=ctrl),error=function(e)e)
  if(inherits(mod,'error')) { warning(tag,': ',conditionMessage(mod)); return(NULL) }
  vc <- as.data.frame(lme4::VarCorr(mod))
  player_var <- vc$vcov[vc$grp=='player_id'][1]
  resid_var <- sigma(mod)^2
  icc <- player_var/(player_var+resid_var)
  singular <- lme4::isSingular(mod,tol=1e-4)
  row <- data.frame(model=tag,players=nlevels(z$player_id),player_seasons=nrow(z),
                    median_shots=median(z[[weight_col]]),player_variance=player_var,
                    residual_variance_reference=resid_var,icc_reference=icc,
                    singular=singular,converged=length(mod@optinfo$conv$lme4$messages)==0,
                    icc_ci_low=NA_real_,icc_ci_high=NA_real_,bootstrap_success=0L)
  # Parametric bootstrap keeps observed player/season design fixed. It measures
  # uncertainty under the fitted Gaussian mixed model, not model misspecification.
  if(B>0) {
    sims <- tryCatch(lme4::bootMer(mod,FUN=function(m){
      v <- as.data.frame(lme4::VarCorr(m)); pv <- v$vcov[v$grp=='player_id'][1]
      pv/(pv+sigma(m)^2)
    },nsim=B,use.u=FALSE,type='parametric',parallel='no'),error=function(e)e)
    if(!inherits(sims,'error')) {
      vals <- as.numeric(sims$t[,1]); vals <- vals[is.finite(vals)]
      if(length(vals)>=20) { ci <- quantile(vals,c(.025,.975),na.rm=TRUE); row$icc_ci_low <- unname(ci[1]); row$icc_ci_high <- unname(ci[2]) }
      row$bootstrap_success <- length(vals)
    } else warning('Bootstrap failed for ',tag,': ',conditionMessage(sims))
  }
  capture.output(summary(mod),file=file.path(opt$output,paste0(tag,'_model_summary.txt')))
  fe <- data.frame(model=tag,term=names(lme4::fixef(mod)),estimate=unname(lme4::fixef(mod)))
  re <- lme4::ranef(mod)$player_id
  re_out <- data.frame(model=tag,player_id=rownames(re),random_intercept=re[[1]],row.names=NULL)
  list(row=row,fe=fe,re=re_out,model=mod)
}

results <- list(); fixed <- list(); random <- list()
run <- function(z,metric,weight_col,tag) {
  ans <- fit_one(z,metric,weight_col,tag)
  if(is.null(ans)) return(invisible(NULL))
  results[[length(results)+1L]] <<- ans$row
  fixed[[length(fixed)+1L]] <<- ans$fe
  random[[length(random)+1L]] <<- ans$re
  invisible(ans)
}
run(d,'fva_per_shot','fva_usable_shots','fva_main')
run(d,'xg_per_shot','shots','shot_quality_comparator')
# Sensitivity A: at least 3 qualifying seasons per player.
d3 <- d[d$player_id %in% names(which(table(d$player_id)>=3)),,drop=FALSE]
if(nrow(d3)>0) run(d3,'fva_per_shot','fva_usable_shots','fva_3plus_seasons')
# Sensitivity B: no missing usable xGOT in any included player-season.
d0 <- d[d$fva_missing_shots==0,,drop=FALSE]
d0 <- d0[d0$player_id %in% names(which(table(d0$player_id)>=min_seasons)),,drop=FALSE]
if(nrow(d0)>0) run(d0,'fva_per_shot','fva_usable_shots','fva_complete_xgot')
# Sensitivity C: higher-volume shooters, 50+ in at least 2 seasons.
d50 <- d[d$shots>=50,,drop=FALSE]
d50 <- d50[d50$player_id %in% names(which(table(d50$player_id)>=min_seasons)),,drop=FALSE]
if(nrow(d50)>0) run(d50,'fva_per_shot','fva_usable_shots','fva_50plus_shots')
if(!length(results)) stop('No models could be fitted. Check data and model warnings.')
summary_tbl <- do.call(rbind,results)
write.csv(summary_tbl,file.path(opt$output,'mixed_model_variance_icc.csv'),row.names=FALSE)
write.csv(do.call(rbind,fixed),file.path(opt$output,'season_fixed_effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,random),file.path(opt$output,'player_random_intercepts.csv'),row.names=FALSE)

# Main ICC comparison plot; confidence limits if bootstrap available.
plot_rows <- summary_tbl[summary_tbl$model %in% c('fva_main','shot_quality_comparator'),,drop=FALSE]
png(file.path(opt$output,'finishing_vs_shot_quality_icc.png'),width=1200,height=700,res=150)
par(mar=c(6,17,3,3))
xx <- plot_rows$icc_reference
labels <- ifelse(plot_rows$model=='fva_main','Finishing execution (FVA)','Shot quality (xG/shot)')
plot(xx,seq_along(xx),xlim=c(0,max(.15,xx,plot_rows$icc_ci_high,na.rm=TRUE)*1.12),
     ylim=c(.5,length(xx)+.5),pch=19,cex=2,col=c('#bc4d59','#168b80')[seq_along(xx)],
     yaxt='n',ylab='',xlab='Estimated player ICC (reference shot volume)',main='Persistent player differences: five Premier League seasons')
axis(2,at=seq_along(xx),labels=labels,las=1)
for(i in seq_along(xx)) if(is.finite(plot_rows$icc_ci_low[i]) && is.finite(plot_rows$icc_ci_high[i]))
  segments(plot_rows$icc_ci_low[i],i,plot_rows$icc_ci_high[i],i,lwd=2)
abline(v=0,col='gray70');dev.off()

readme <- c('Under Pressure | Mixed-effects finishing analysis',
 'Primary model: player-season FVA/shot ~ season fixed effects + player random intercept.',
 'Comparator: xG/shot with same design and shot-count precision proxy.',
 'Weights: usable shots / median usable shots for FVA; shots / median shots for xG.',
 'ICC = between-player variance / (between-player variance + residual variance at median shot count).',
 'Weighted models have shot-count-dependent residual variance; ICC is NOT a universal property.',
 'Bootstrap: parametric, conditional on model assumptions. Small ICC estimates may be boundary/singular.',
 'No ages or positions in source data; cannot infer coaching effects or age curves.',
 'Missing xGOT sensitivity changes the cohort; do not interpret as causal effect of missingness.',
 'A large xG/shot ICC can reflect stable positions, tactics and chance types.',
 'Player random intercepts are shrinkage estimates, not rankings of true finishing talent.',
 'Selection: only player-seasons above threshold; survivor and role bias remain.',
 'Model is an exploratory Gaussian approximation to aggregated per-shot FVA; a shot-level hierarchical model would be stronger.')
writeLines(readme,file.path(opt$output,'README_interpretation.txt'))
print(summary_tbl,row.names=FALSE)
cat('\nResults saved in: ',normalizePath(opt$output),'\n',sep='')
