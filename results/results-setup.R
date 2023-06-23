library(dplyr)
library(magrittr)
library(viridis)
library(ggplot2)
library(tidyr)
library(kableExtra)
library(ggh4x)
library(ggsci)
library(mvtnorm)
library(moments)
library(here)
library(gridExtra)
library(ellipse)
library(rWishart)

## Functions for methods

output.iid <- function(y){
  mu.y <- mean(y)
  var.y <- var(y)

  mode <- y - mu.y
  pear.res <- mode/sqrt(var.y)

  Fx <- pgamma(y, shape = mu.y^2/var.y, scale = var.y/mu.y)
  quant.res <- qnorm(Fx)

  out <- list(pears = pear.res, quant = quant.res)
  return(out)
}

output.mvn <- function(y, M, distribution){
  mu.y <- mean(y)
  var.y <- var(y)

  mode <- y - mu.y
  pear.res <- mode/sqrt(var.y)

  L <- t(chol(M))
  if(distribution == "normal"){
    quant.res <- t(solve(L, mode))
  }
  if(distribution == "gamma"){
    r <- qnorm(pgamma(y, mu.y^2/var.y, scale = var.y/mu.y))
    quant.res <- t(solve(L, r))
  }
  out <- list(pears = pear.res, quant = as.vector(quant.res))
  return(out)
}

qqex.plot <- function(sim.y, out.y, title){
  df <- data.frame(Residual = c(out.y$pears, out.y$quant),
                   Type = c(
                     rep("Pearson", length(sim.y)),
                     rep("Quantile", length(sim.y))
                   ))
  axis_titles <- data.frame(
    Residual = c("Pearson", "Quantile"),
    Type = c("Pearson", "Quantile"),
    axis_title = c(
      paste("KS Test:",
            round(ks.test(df$Residual[df$Type=="Pearson"], "pnorm")$p.value,4)),
      paste("KS Test:",
            round(ks.test(df$Residual[df$Type=="Quantile"], "pnorm")$p.value,4))
    )
  )
  qq <- ggplot(df, mapping = aes(sample = Residual)) +
    stat_qq() +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~Type) + labs(x=NULL) +
    geom_text(
      data=axis_titles,
      aes(label=axis_title), hjust=0.5,
      x=-0.5,
      y=max(df$Residual)-0.25, color='red'
    ) +
    coord_cartesian(clip="off") +
    theme(
      plot.margin= margin(b=30)
    ) +
    theme_bw()

  h <- ggplot(data.frame(y = sim.y), aes(y)) +
    geom_histogram() + theme_bw() +
    ggtitle(title)

  out <- list(h = h, qq = qq)

 # out <- grid.arrange(h, qq, ncol=2)
  return(out)

}

mvn.demo <- function(){
  set.seed(1) # 1 is good
  n <- 3
  N <- 5000
  S <- toeplitz((n:1)/n)
  C <- rWishart::rWishart(1,n,Sigma=S)[,,1]
  L <- t(chol(C))
  mean <- rep(0,n)
  obs <- as.numeric(rmvnorm(1, mean=mean, sigma=C))
  draws <- rmvnorm(n=N, mean=mean,sigma=C)
  ## Correctly rotated space
  draws.rotated <- t(solve(L, t(draws)) )
  obs.rotated <- solve(L, obs)
  cols <- c(rgb(0,0,0,.5), rgb(1,0,0,.5))
  par(mfcol=c(n,n), mar=0*c(.65,.65,.65,.65), oma=c(.5,.5,.5,0),
      mgp=c(.5,.1,0), tck=-.01, cex.axis=.6 )
  for(i in 1:n){
    for(j in 1:n){
      xlim <- range(c(draws[, c(i,j)], draws.rotated[,c(i,j)]))
      ylim <- range(c(draws[, c(i,j)], draws.rotated[,c(i,j)]))
      if(i==j){
        x <- seq(xlim[1], xlim[2], len=500)
        y1 <- dnorm(x, 0, sd=sqrt(C[i,i]))
        y2 <- dnorm(x,0,sd=1)
        plot(x,y1, ylim=c(0, max(c(y1,y2))*1.4), type='l',
             col=cols[1], axes=FALSE)
        lines(x,y2, col=cols[2])
        ## hist(draws[,j], xlab='', main='', ylab='', col=cols[1], border=cols[1],
        ##      ylim=c(0,.5), xlim=xlim, freq=FALSE)
        ## abline(v=obs[j], col=cols[1])
        points(c(obs[j], obs[j])[2],
               c(0,dnorm(obs[j],0,sqrt(C[j,j])))[2], col=cols[1], pch=15)
        lines(c(obs[j], obs[j]), c(0,dnorm(obs[j],0,sqrt(C[j,j]))), col=cols[1])
        mtext(line=-1, col=cols[1], paste("Marginal percentile=", round(mean(obs[j] > draws[,j]),3)), cex=.7)
        ## hist(draws.rotated[,j], xlab='', main='', ylab='',
        ##      add=TRUE, col=cols[2], freq=FALSE, border=cols[2])
        ## abline(v=obs.rotated[j], col=cols[2])
        points(c(obs.rotated[j], obs.rotated[j])[2],
               c(0,dnorm(obs.rotated[j],0,sqrt(C[j,j])))[2], col=cols[2], pch=15)
        lines(c(obs.rotated[j], obs.rotated[j]), c(0,dnorm(obs.rotated[j],0,sqrt(C[j,j]))), col=cols[2])
        mtext(line=-2, col=cols[2], paste("Marginal percentile=", round(mean(obs.rotated[j] > draws.rotated[,j]),3)), cex=.7)
        box()
      }
      if(i<j){
        ## plot(draws[,i], draws[,j], ann=FALSE, pch=16, cex=.25,
        ##      col=cols[1])
        plot(obs[i], obs[j],col=4, cex=1, pch=16, xlim=xlim,
             ylim=ylim, axes=FALSE)
        ## points(draws.rotated[,i], draws.rotated[,j], ann=FALSE,
        ##        pch=16, cex=.25,  col=cols[2])
        ## points(obs[i], obs[j],col=3, cex=2, pch=16)
        arrows(obs[i], obs[j], obs.rotated[i], obs.rotated[j],
               length=.05, lwd=1.5, col=4)
        lines(ellipse(C[c(i,j),c(i,j)], centre=mean[c(i,j)]), col=cols[1])
        lines(ellipse(diag(2), centre=mean[c(i,j)]),  col=cols[2])
        box()
      }
      if(i>j) {plot(1,1, type='n', axes=FALSE, ann=FALSE)}
    }
  }
}


## Functions for results

## Define path to results files
path <- paste0(here(), "/results")

## Read in pvalue results
pvals <- lapply(list.files(path, pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
## Read in MLEs
mles <- lapply(list.files(path, pattern='_mles.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
## Read in stats
Stats <- lapply(list.files(path, pattern='_stats.RDS',
                          full.names=TRUE), readRDS) %>% bind_rows

## Models with non-convergence
nc <- Stats[!is.na(Stats$converge) & (Stats$converge==1|Stats$maxgrad>0.1),]
nc.h0 <- nc[nc$version == "h0",]
nc.h1 <- nc[nc$version == "h1",]

table(nc.h0$model, nc.h0$misp, nc.h0$do.true)
table(nc.h1$model, nc.h1$misp, nc.h1$do.true)


## Data Pre-processing

### p-values
pvals$version <- factor(pvals$version,
                        levels = c('h0', 'h1'),
                        labels = c('Correct',
                                   'Mis-specified'))
pvals$type[which(pvals$method == "pears")] <- "pears"
pvals$type[which(pvals$method == "mcmc")] <- "sim"
pvals$type <- factor(pvals$type,
                     levels = c("osa", "sim", "pears"),
                     labels = c("Analytical Methods",
                                "Simulation Methods",
                                "Pearson"))

pvals <- dplyr::filter(pvals,
                       method %in% c(
                         'fg',
                         'osg',
                         'gen',
                         'cdf',
                         'mcmc',
                         'pears',
                         'uncond',
                         'uncond_nrot',
                         'cond',
                         'cond_nrot'))

pvals$level <- factor(pvals$method,
                      level = c(
                        'fg',
                        'osg',
                        'gen',
                        'cdf',
                        'mcmc',
                        'pears',
                        'uncond',
                        'uncond_nrot',
                        'cond',
                        'cond_nrot'),
                      label = c(
                        "Unconditional",
                        "Unconditional",
                        "Unconditional",
                        "Unconditional",
                        "Unconditional",
                        "Conditional",
                        "Unconditional",
                        "Unconditional",
                        "Conditional",
                        "Conditional"))
pvals$method <- factor(pvals$method,
                       level = c(
                        'pears',
                        'gen',
                        'osg',
                        'fg',
                        'cdf',
                        'mcmc',
                        'uncond',
                        'uncond_nrot',
                        'cond',
                        'cond_nrot'),
                      label = c(
                        'Pearson',
                        'one-step Generic',
                        'one-step Gaussian',
                        'full Gaussian',
                        'cdf',
                        'MCMC',
                        "Unconditional ecdf, Rotated",
                        "Unconditional ecdf, Not Rotated",
                        "Conditional ecdf, Rotated",
                        "Conditional ecdf, Not Rotated"
                      ))

#Only filter out non-converging models if model is correctly specified
nc.id <- dplyr::filter(nc, version == "h0")$id

if(length(nc.id) > 0){
  nc.idx <- which(pvals$id %in% nc.id)
  pvals <- pvals[-nc.idx,]
}

## Functions
filter.true <- function(df, mod, test = "GOF.ks", method.vec){
  dplyr::filter(df, model == mod &
                  test == "GOF.ks" &
                  method %in% method.vec &
                  do.true == TRUE)
}

filter.est <- function(df, mod, test = "GOF.ks", method.vec){
  dplyr::filter(df, model == mod &
                  test == "GOF.ks" &
                  method %in% method.vec &
                  do.true == FALSE)
}

plot.pval.hist <- function(df, doTrue){
  no.misp <- length(unique(df$misp))
  p <- ggplot(df, aes(pvalue, fill = version, color=version))
  if(no.misp == 1){
    p <- p + facet_grid2(method ~ version, labeller = label_wrap_gen(12),
                         scales = "free_y", independent = "y")
  }
  if(no.misp > 1){
    p <- p + facet_nested(method ~ version+misp, labeller = label_wrap_gen(12),
                         scales = "free_y", independent = "y")
  }
  p <- p +  geom_histogram(position ='identity', bins = 50,
                           show.legend = FALSE) +
    scale_fill_manual(values = c("#440154FF", "#35B779FF")) +
    scale_color_manual(values = c("#440154FF", "#35B779FF")) +  theme_bw() +
    scale_x_continuous(breaks=c(0,0.5, 1))
  print(p)
}

plot.err.pow <- function(df.true, df.est){
  pvals.true <- df.true %>% filter(version == "correct") %>%
    group_by(misp, method, type) %>%
    summarize(typeIerror = sum(pvalue <= 0.05)/sum(pvalue >= 0))
  pvals.true$restype <- "Theoretical"
  pvals.true$power <- NA

  pvals.err.est <- df.est %>% filter(version == "correct") %>%
    group_by(misp, method, type) %>%
    summarize(typeIerror = sum(pvalue <= 0.05)/sum(pvalue >= 0))
  pvals.err.est$restype <- "Estimated"

  pvals.power <-df.est %>% filter(version == "mis-specified") %>%
    group_by(misp, method, type) %>%
    summarize(power = sum(pvalue <= 0.05)/sum(pvalue >= 0))
  pvals.power$restype <- "Estimated"
  pvals.est <- left_join(pvals.err.est, pvals.power)

  pvals.comb <- rbind(pvals.true, pvals.est) %>%
    pivot_longer(., c(4,6), names_to = "metric", values_to = "pvalue")

  pvals.comb$metric = factor(pvals.comb$metric,
                             levels = c("typeIerror", "power"),
                             labels = c("Type I Error", "Power"))

  data_vline <- tidyr::expand_grid(misp = unique(pvals.comb$misp),
                                   metric = unique(pvals.comb$metric))
  data_vline$vline <- ifelse(data_vline$metric == "Type I Error", 0.05, 0.95)

  p <- pvals.comb  %>%
    ggplot(., aes(x = pvalue, y = method)) +
    geom_vline(data = data_vline,
               aes(xintercept = vline), color = "red",
               size = 0.75, show.legend = FALSE) +
    geom_point(mapping = aes(color = type)) +
    facet_nested(restype ~ misp + metric, labeller = label_wrap_gen(16) )  +
    scale_x_continuous(breaks=c(0,0.5,1)) +
    ylab("") + xlab("") +
    scale_color_viridis_d() +
    theme_bw() + theme(legend.position="bottom")
  print(p)
}

tbl.err.pow <- function(df, caption = NULL){


  pvals.err <- df %>% filter(version == "correct") %>%
    group_by(misp, method, type) %>%
    summarize('Type I Error' = sum(pvalue <= 0.05)/sum(pvalue >= 0))

  pvals.power <-df %>% filter(version == "mis-specified") %>%
    group_by(misp, method, type) %>%
    summarize(Power = sum(pvalue <= 0.05)/sum(pvalue >= 0))

  pvals.est <- left_join(pvals.err, pvals.power)
  misp.names <- as.character(unique(pvals.est$misp))
  nmisp <- length(misp.names)
  misp.header = c(1, rep(2, nmisp))
  names(misp.header) <- c(" ", misp.names)

  if(nmisp == 1){
    tbl <- pvals.est %>%
      as.data.frame() %>% dplyr::select(method, 'Type I Error', Power)  %>%
      kableExtra::kbl(., format = "latex", caption = caption, booktabs = TRUE) %>%
        kableExtra::kable_styling(., "striped", "HOLD_position") %>%
      kableExtra::add_header_above(., misp.header)
  }

  if(nmisp > 1){

    tbl <- pvals.est %>%
      tidyr::pivot_longer(., cols = 4:5, names_to = "metric", values_to = "pvalue") %>%
      dplyr::select(misp, method, metric, pvalue) %>%
      tidyr::pivot_wider(., names_from = c(misp, metric), values_from = pvalue) %>%
      as.data.frame()
    colnames(tbl) <- c("method", rep(c("Type I Error", "Power"), nmisp))
    tbl <- tbl %>%
      kableExtra::kbl(., format = "latex", caption = caption, booktabs = TRUE) %>%
      kableExtra::kable_styling(., "striped", "HOLD_position") %>%
      kableExtra::add_header_above(., misp.header)
  }

  tbl

}

### Type I error and Power
results.simpleGLMM.grps <- readRDS(paste0(path, '/simpleGLMM_missunifcov_obs_sample_sizes.RDS'))
results.simpleGLMM.grps <- readRDS(paste0(path, '/simpleGLMM_missunifcov_grps_sample_sizes.RDS'))
results.simpleGLMM.grps$pvals$model <-
  results.simpleGLMM.grps$runtimes$model <-
  "simpleGLMM.vary.ngrps"

results.simpleGLMM.obs <- readRDS(paste0(path, '/simpleGLMM_missunifcov_grps_sample_sizes.RDS'))
results.simpleGLMM.obs$pvals$model <-
  results.simpleGLMM.obs$runtimes$model <-
  "simpleGLMM.vary.nobs"

results.randomwalk <- readRDS(paste0(path, '/randomwalk_mu0_sample_sizes.RDS'))
results.spatial <- readRDS(paste0(path, '/spatial_mispomega_sample_sizes.RDS'))
runtimes.all <- rbind(results.simpleGLMM.obs$runtimes,
                     # results.simpleGLMM.grps$runtimes,
                      results.randomwalk$runtimes,
                      results.spatial$runtimes)
runtimes.all <- runtimes.all %>% filter(!is.na(med))
pvals.all <- rbind(results.simpleGLMM.obs$pvals,
                    #  results.simpleGLMM.grps$pvals,
                      results.randomwalk$pvals,
                      results.spatial$pvals)

runtimes.all$level <- factor(runtimes.all$type,
                             level = c(
                               'fg',
                               'osg',
                               'gen',
                               'cdf',
                               'mcmc',
                               'pears',
                               'uncond',
                               'uncond_nrot',
                               'cond',
                               'cond_nrot'),
                             label = c(
                               "Unconditional",
                               "Unconditional",
                               "Unconditional",
                               "Unconditional",
                               "Unconditional",
                               "Conditional",
                               "Unconditional",
                               "Unconditional",
                               "Conditional",
                               "Conditional"))

runtimes.all$method <- factor(runtimes.all$type,
                              level = c(
                                'fg',
                                'osg',
                                'gen',
                                'cdf',
                                'mcmc',
                                'pears',
                                'uncond',
                                'uncond_nrot',
                                'cond',
                                'cond_nrot'),
                              label = c(
                                'full Gaussian',
                                'one-step Gaussian',
                                'one-step Generic',
                                'cdf',
                                'MCMC',
                                'Pearson',
                                "Unconditional ecdf, Rotated",
                                "Unconditional ecdf, Not Rotated",
                                "Conditional ecdf, Rotated",
                                "Conditional ecdf, Not Rotated"
                              ))
pvals.all$method <- factor(pvals.all$method,
                           level = c(
                             'pears',
                             'gen',
                             'osg',
                             'fg',
                             'cdf',
                             'mcmc',
                             'uncond',
                             'uncond_nrot',
                             'cond',
                             'cond_nrot'),
                           label = c(
                             'Pearson',
                             'one-step Generic',
                             'one-step Gaussian',
                             'full Gaussian',
                             'cdf',
                             'MCMC',
                             "Unconditional ecdf, Rotated",
                             "Unconditional ecdf, Not Rotated",
                             "Conditional ecdf, Rotated",
                             "Conditional ecdf, Not Rotated"
                           ))


t1.err.plot <- function(df){
  x.name <- "Number of Observations"
  if("simpleGLMM.vary.nobs" %in% df$model |
     "simpleGLMM.vary.ngroups" %in% df$model ){
    x.name <- "Number of Observations x Number of Groups"
  }
  p <- df %>%
    filter(version == "h0") %>%
    group_by(nobs, method) %>%
    summarise(t1_err = sum(pvalue<0.05)/sum(pvalue>=0)) %>%
    ggplot(aes(x = nobs, y = t1_err)) +
    ylab("Type I Error Rate") +
    xlab(x.name) +
    geom_point() + geom_line() +
    facet_wrap(~method, labeller = label_wrap_gen(14)) +
    theme_bw()

  print(p)
}

pow.plot <- function(df){
  x.name <- "Number of Observations"
  if("simpleGLMM.vary.nobs" %in% df$model |
     "simpleGLMM.vary.ngroups" %in% df$model ){
    x.name <- "Number of Observations x Number of Groups"
  }
  p <- df %>%
    filter(version == "h1") %>%
    group_by(nobs, method) %>%
    summarise(power = sum(pvalue<=0.05)/sum(pvalue>=0)) %>%
    ggplot(aes(x = nobs, y = power)) +
      ylab("Power") +
      xlab(x.name)  +
      geom_point() + geom_line() +
      facet_wrap(~method, labeller = label_wrap_gen(14)) +
      theme_bw()

  print(p)
}

