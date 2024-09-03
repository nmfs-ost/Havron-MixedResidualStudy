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
library(stringr)
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
    geom_vline(xintercept = mean(sim.y), color = "red")
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
nc <- Stats[(Stats$converge.status == 1 |
               Stats$converge.hessian == FALSE |
               Stats$converge.maxgrad > 0.01),]
nc.h0 <- nc[nc$version == "h0",]
nc.h1 <- nc[nc$version != "h0",]
nc.cs <- Stats[Stats$converge.status == 1,]
nc.cs.h0 <- nc.cs[nc.cs$version == "h0",]
nc.cs.h1 <- nc.cs[nc.cs$version != "h0",]

table(nc.h0$model, nc.h0$misp, nc.h0$type, nc.h0$do.true)
table(nc.h1$model, nc.h1$misp, nc.h1$type, nc.h1$do.true)

table(nc.cs.h0$model, nc.cs.h0$misp, nc.cs.h0$type, nc.cs.h0$do.true)
table(nc.cs.h1$model, nc.cs.h1$misp, nc.cs.h1$type, nc.cs.h1$do.true)


## Data Pre-processing

### p-values
# use pvals$misp instead
# pvals$version <- factor(pvals$version,
#                         levels = c('h0', 'h1', 'h1', 'h1'),
#                         labels = c('Correct',
#                                    'Mis-specified', 'Mis-specified', ))
pvals$res.type[which(pvals$method == "pears")] <- "pears"
pvals$res.type[which(pvals$method == "mcmc")] <- "sim"
pvals$res.type <- factor(pvals$res.type,
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
pvals$misp.type <- factor(pvals$misp,
                          level = c(
                            "correct",
                            "overdispersion",
                            "missre",
                            "gamma-normal",
                            "mu0",
                            "mispre",
                            "hsk",
                            "nb-pois",
                            "missunifcov",
                            "pois-zip",
                            "ln-error"
                          ),
                          label = c(
                            "Correct",
                            "Misp. Data Model",
                            "Missing RE",
                            "Misp. Data Model",
                            "Misp. Data Model",
                            "Misp. RE Model",
                            "Misp. Data Model",
                            "Misp. Data Model",
                            "Misp. Data Model",
                            "Misp. Data Model",
                            "Misp. Data Model" ))

pvals$do.true <- factor(pvals$do.true,
                        levels = c(TRUE, FALSE),
                        labels = c("theoretical", "estimated"))

pvals$test <- factor(pvals$test,
                     levels = c("GOF.ad",
                                "GOF.ks",
                                "GOF.lf",
                                "outlier",
                                "AOV",
                                "EqVar",
                                "disp",
                                "Auto",
                                "SAC"
                     ),

                     labels = c("Anderson-Darling",
                                "Kolmogorov-Smirnov",
                                "Lilliefors",
                                "outlier",
                                "Analysis of Variance",
                                "Levene's Equal Variance",
                                "disp",
                                "Autocorrelation",
                                "Spatial Autocorrelation"
                     ))

#Only filter out non-converging models if model is correctly specified
nc.id <- dplyr::filter(nc, version == "h0")$id

if(length(nc.id) > 0){
  nc.pvals.idx <- which(pvals$id %in% nc.id)
  nc.mles.idx <- which(mles$id %in% nc.id)
  pvals <- pvals[-nc.pvals.idx,]
  mles <- mles[-nc.mles.idx,]
}


## Functions
plot.mles <- function(df){
  p <- df %>% ggplot() +
      geom_violin(aes(x = misp, y = bias)) +
      geom_hline(aes(yintercept = 0), linetype = "dashed") +
      facet_grid(par~type, scales = 'free_y',
                 labeller = label_parsed) +
      theme_bw()
  return(p)
}


filter.true <- function(df, mod, type_, test_ = NA, method.vec){
  if(is.na(test_)){
    dplyr::filter(df, model == mod & type == type_ &
                    method %in% method.vec &
                    do.true == "theoretical")
  } else {
    dplyr::filter(df, model == mod & type == type_ &
                    test == test_ &
                    method %in% method.vec &
                    do.true == "theoretical")
  }
}

filter.est <- function(df, mod, type_, test_ = NA, method.vec){
  if(is.na(test_)){
    dplyr::filter(df, model == mod & type == type_ &
                    method %in% method.vec &
                    do.true == "estimated")
  } else {
    dplyr::filter(df, model == mod & type == type_ &
                  test == test_ &
                  method %in% method.vec &
                  do.true == "estimated")
  }
}

filter.all <- function(df, mod, type_, method.vec){

    dplyr::filter(df, model == mod & type == type_ &
                    method %in% method.vec)
}

plot.pval.hist <- function(df, doTrue){
  no.misp <- sum(unique(df$misp) != "A: Correct")
  p <- ggplot(df, aes(pvalue, fill = misp, color=misp))
  if(no.misp == 1){
    p <- p + facet_grid2(method ~ misp, labeller = label_wrap_gen(12),
                         scales = "free_y", independent = "y")
  }
  if(no.misp > 1){
    p <- p + facet_nested(method ~ misp, labeller = label_wrap_gen(12),
                         scales = "free_y", independent = "y")
  }
  p <- p +  geom_histogram(position ='identity', bins = 50,
                           show.legend = FALSE) +
     scale_fill_viridis_d(begin = 0, end = 0.9) +
    scale_color_viridis_d(begin = 0, end = 0.9)  +  theme_bw() +
    scale_x_continuous(breaks=c(0,0.5, 1))
  print(p)
}

plot.ecdf <- function(df, doTrue){
  ggplot(pval.df, aes(pvalue, color = method)) + stat_ecdf(geom = "step") + facet_grid(method~misp)
}

histogram.err.pow <- function(df, Type, ResType){
  df %>% dplyr::filter(type == Type & do.true == ResType & model != "linmod" &
                         test != "outlier" & test != "disp") %>%
    group_by(test, model, misp, method) %>%
    summarize(pval = round(sum(pvalue <= 0.05)/sum(pvalue >= 0),3))%>%
    ggplot(., aes(x = method, y = pval, fill = test, group = test)) + geom_col(position = "dodge") +
    facet_grid(misp ~ model) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

histogram.err <- function(df, Type){
  df %<>% dplyr::filter(type == Type & model != "linmod" &
                         test != "outlier" & test != "disp" & misp == "correct") %>%
    group_by(test, model, method, do.true) %>%
    summarize(pval = round(sum(pvalue <= 0.05)/sum(pvalue >= 0),3))
  if(Type == "GLMM"){
    df %<>% dplyr::filter(method == "Pearson" | method == "cdf" |
                            method == "Conditional ecdf, Not Rotated")
  } else {
    df %<>% dplyr::filter(method == "Pearson" | method == "cdf" |
                            method == "Conditional ecdf, Not Rotated" |
                            method == "Unconditional ecdf, Rotated" |
                            method == "Unconditional ecdf, Not Rotated")
  }
  ggplot(df, aes(x = method, y = pval, color = test, group = test)) +geom_point() +
    facet_grid(do.true ~ model) + theme_bw() + geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

histogram.pow <- function(df, Type){
  df %>% dplyr::filter(type == Type & model != "linmod" &
                         test != "outlier" & test != "disp" & misp != "Correct") %>%
    group_by(test, model, method, misp, do.true) %>%
    summarize(pval = round(sum(pvalue <= 0.05)/sum(pvalue >= 0),3)) %>%
    ggplot(., aes(x = method, y = pval, color = test, group = test)) +geom_point() +
    facet_grid(do.true ~ model) +
    theme_bw() + geom_hline(yintercept = 0.95) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

histogram.pow.bymisp <- function(df, Type){
  df %<>% dplyr::filter(type == Type & model != "linmod" &
                         test != "outlier" & test != "disp" & misp != "correct") %>%
    group_by(test, model, method, misp, misp.type, do.true) %>%
    summarize(pval = round(sum(pvalue <= 0.05)/sum(pvalue >= 0),3))
    if(Type == "GLMM"){
      df %<>% dplyr::filter(method == "Pearson" | method == "cdf" |
                              method == "Conditional ecdf, Not Rotated")
    } else {
      df %<>% dplyr::filter(method == "Pearson" | method == "cdf" |
                              method == "Conditional ecdf, Not Rotated" |
                              method == "Unconditional ecdf, Rotated" |
                              method == "Unconditional ecdf, Not Rotated")
    }
    ggplot(df, aes(x = method, y = pval, color = test, group = test)) +geom_point() +
    facet_grid(misp.type ~ model + do.true) +
    theme_bw() + geom_hline(yintercept = 0.95) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "bottom")
}

plot.err.pow <- function(df){

  results <- df %>% filter(test != "outlier") %>%
    group_by(test, misp, method, res.type, do.true) %>%
    summarize(pval = round(sum(pvalue <= 0.05)/sum(pvalue >= 0),3))

  results$do.true <- factor(results$do.true,
                            levels = c(TRUE, FALSE),
                            labels = c("theoretical", "estimated"))
  results$err_type <- ifelse(results$misp == "A: Correct", "Type I Error", "Power")
  results$err_type <- factor(results$err_type,
                             levels = c("Type I Error", "Power"),
                             labels = c("Type I Error", "Power"))
  results$err_line <- ifelse(results$misp == "A: Correct", 0.05, 0.95)

  p <- results  %>%
    ggplot(., aes(x = pval, y = method))  +
    geom_point(mapping = aes(color = do.true)) +
    scale_color_aaas() + xlab("p-value") +
    facet_grid(test~err_type + misp) + theme_bw() +
    theme(legend.position = "top") + #, legend_title = element_blank())+#,
    #      axis.text.x = element_blank(),
    #      axis.ticks.x = element_blank()) +
    scale_x_continuous(breaks = c(0,0.5,1),
                       labels = c("0", "0.5", "1")) +
    geom_vline(mapping = aes( xintercept = results$err_line)) +
    labs(color = "reisduals")


  return(p)
}

err.table <-  function(df1, df2, sig.level, caption = NULL){

  pvals1 <- df1 %>% filter(test != "outlier" & misp == "Correct") %>%
    group_by(test, method, do.true) %>%
    summarize(pval = round(sum(pvalue <= sig.level)/sum(pvalue >= 0),3))
  pvals2 <- df2 %>% filter(test != "outlier" & misp == "Correct") %>%
    group_by(test, method, do.true) %>%
    summarize(pval = round(sum(pvalue <= sig.level)/sum(pvalue >= 0),3))

  pvals.wider1 <- pivot_wider(pvals1, names_from = method, values_from = pval)
  pvals.wider2 <- pivot_wider(pvals2, names_from = method, values_from = pval)

  pvals.wider <- rbind(pvals.wider1, pvals.wider2)
  pvals.wider$do.true <- factor(pvals.wider$do.true,
                          level = c(TRUE, FALSE),
                          label = c("theoretical", "estimated"))
  colnames(pvals.wider)[which(colnames(pvals.wider) == "do.true")] <- "residual type"

  pvals.wider$test <- factor(pvals.wider$test,
                        level = c(
                          "AOV",
                          "EqVar",
                          "GOF.ad",
                          "GOF.ks",
                          "GOF.lf",
                          "Auto",
                          "SAC"
                        ),
                        label = c(
                          "AOV Equal Means",
                          "Levene's Equal Variances",
                          "Anderson Darling",
                          "Kolmogorov-Smirnov",
                          "Lilliefors",
                          "Autocorrelation",
                          "Moran's I"
                        ))

  out.table <- pvals.wider %>% arrange(test) %>%
    kableExtra::kbl(., format = "latex", caption = caption,
                    booktabs = TRUE, midrule = "", escape = FALSE) %>%
    kableExtra::kable_styling(., "striped", "HOLD_position")
  if(unique(df1$type) == "LMM"){
    out.table <- out.table %>%
      column_spec(1, width = "7em") %>%
      column_spec(2, width = "4em")%>%
      column_spec(3, width = "4em")%>%
      column_spec(4, width = "4em")%>%
      column_spec(5, width = "4em")%>%
      column_spec(6, width = "4em")%>%
      column_spec(8, width = "4em")%>%
      column_spec(9, width = "6em")%>%
      column_spec(10, width = "5em")%>%
      column_spec(11, width = "6em") %>%
      column_spec(12, width = "5em") %>%
      collapse_rows(column = 1, latex_hline = "major")
  }
  if(unique(df1$type) == "GLMM"){
    out.table <- out.table %>%
      column_spec(1, width = "7em") %>%
      column_spec(2, width = "4em")%>%
      column_spec(3, width = "4em")%>%
      column_spec(4, width = "4em")%>%
      column_spec(5, width = "4em")%>%
      column_spec(6, width = "4em")%>%
      column_spec(7, width = "6em")%>%
      column_spec(8, width = "6em") %>%
      column_spec(9, width = "6em") %>%
      column_spec(10, width = "6em") %>%
      collapse_rows(column = 1, latex_hline = "major")
  }

  out.table

}

pow.table <-  function(df1, df2, sig.level, caption = NULL){

  misp.names <- as.character(unique(df1$misp)[-1])
  nmethods <- length(unique(df1$method))
  misp.header = 1
  names(misp.header) <- misp.names[1]

  for(i in seq_along(misp.names)){

    pvals1 <- df1 %>% filter(test != "outlier" & misp == misp.names[i]) %>%
      group_by(test, method, do.true) %>%
      summarize(pval = round(sum(pvalue <= sig.level)/sum(pvalue >= 0),3))
    pvals2 <- df2 %>% filter(test != "outlier" & misp == misp.names[i]) %>%
      group_by(test, method, do.true) %>%
      summarize(pval = round(sum(pvalue <= sig.level)/sum(pvalue >= 0),3))

    pvals.wider1 <- pivot_wider(pvals1, names_from = method, values_from = pval)
    pvals.wider2 <- pivot_wider(pvals2, names_from = method, values_from = pval)

    pvals.wider <- rbind(pvals.wider1, pvals.wider2)

    out <- pvals.wider %>% as.data.frame() %>% arrange(test)
     if(i == 1){
       out.df <- rbind(c(misp.names[i], rep(" ", ncol(out)-1)), out)
     } else {
       out.df <- rbind(out.df, c(misp.names[i], rep(" ", ncol(out)-1)), out)
     }
}
  out.df$do.true <- factor(out.df$do.true,
                           level = c(TRUE, FALSE, " "),
                           label = c("theoretical", "estimated", " "))
  colnames(out.df)[which(colnames(out.df) == "do.true")] <- "residual type"

  out.df$test <- factor(out.df$test,
                        level = c(
                          "AOV",
                          "EqVar",
                          "GOF.ad",
                          "GOF.ks",
                          "GOF.lf",
                          "Auto",
                          "SAC",
                          misp.names
                        ),
                        label = c(
                          "AOV Equal Means",
                          "Levene's Equal Variances",
                          "Anderson Darling",
                          "Kolmogorov-Smirnov",
                          "Lilliefors",
                          "Autocorrelation",
                          "Moran's I",
                          misp.names
                        ))

  out.table <- out.df %>%
      kableExtra::kbl(., format = "latex", caption = caption,
                      booktabs = TRUE, midrule = "", escape = FALSE) %>%
      kableExtra::kable_styling(., "striped", "HOLD_position")
  if(unique(df1$type) == "LMM"){
    out.table <- out.table %>%
      column_spec(1, width = "7em") %>%
      column_spec(2, width = "4em")%>%
      column_spec(3, width = "4em")%>%
      column_spec(4, width = "4em")%>%
      column_spec(5, width = "4em")%>%
      column_spec(6, width = "4em")%>%
      column_spec(8, width = "4em")%>%
      column_spec(9, width = "6em")%>%
      column_spec(10, width = "5em")%>%
      column_spec(11, width = "6em") %>%
      column_spec(12, width = "5em") %>%
      collapse_rows(column = 1, latex_hline = "major")
  }
  if(unique(df1$type) == "GLMM"){
    out.table <- out.table %>%
      column_spec(1, width = "7em") %>%
      column_spec(2, width = "4em")%>%
      column_spec(3, width = "4em")%>%
      column_spec(4, width = "4em")%>%
      column_spec(5, width = "4em")%>%
      column_spec(6, width = "4em")%>%
      column_spec(7, width = "6em")%>%
      column_spec(8, width = "6em") %>%
      column_spec(9, width = "6em") %>%
      column_spec(10, width = "6em") %>%
      collapse_rows(column = 1, latex_hline = "major")
  }

  out.table

}
tbl.err.pow <- function(df, caption = NULL){


  pvals <- df %>% filter(test != "outlier") %>%
    group_by(test, misp, method, res.type) %>%
    summarize(pval = round(sum(pvalue <= 0.05)/sum(pvalue >= 0),3))

  pvals.wider <- pivot_wider(pvals, names_from = misp, values_from = pval)
  misp.names <- as.character(unique(pvals$misp))
  nmisp <- length(misp.names)

  if(nmisp == 2){
    misp.header = c(1, rep(2, nmisp))
    names(misp.header) <- c(" ", misp.names)
    out <- pvals.wider %>%
      as.data.frame() %>% dplyr::select(method, 'Type I Error', Power)  %>%
      kableExtra::kbl(., format = "latex", caption = caption, booktabs = TRUE) %>%
        kableExtra::kable_styling(., "striped", "HOLD_position") %>%
      kableExtra::add_header_above(., misp.header)
  }

  if(nmisp > 2){
    misp.header = c(rep(1,3), 3)
    names(misp.header) <- c(rep(" ",2), 'Type I Error', 'Power')

    out <- pvals.wider %>%
      dplyr::select(-res.type) %>%
      as.data.frame()
    out <- out %>%
      kableExtra::kbl(., format = "latex", caption = caption, booktabs = TRUE) %>%
      kableExtra::kable_styling(., "striped", "HOLD_position") %>%
      kableExtra::add_header_above(., misp.header)
  }

  out

}
#
# ### Type I error and Power
results.simpleGLMM.lmm <- readRDS(paste0(path, '/simpleGLMM_missunifcov_LMM_obs_sample_sizes.RDS'))
results.simpleGLMM.lmm$pvals$type = "LMM"
results.simpleGLMM.lmm$mles$type = "LMM"
results.simpleGLMM.glmm <- readRDS(paste0(path, '/simpleGLMM_mispre_GLMM_obs_sample_sizes.RDS'))
results.simpleGLMM.glmm$pvals$type = "GLMM"
results.simpleGLMM.glmm$mles$type = "GLMM"
results.randomwalk.lmm <- readRDS(paste0(path, '/randomwalk_hsk_LMM_sample_sizes.RDS'))
results.randomwalk.lmm$pvals$type <- "LMM"
results.randomwalk.lmm$mles$type <- "LMM"
results.randomwalk.glmm <- readRDS(paste0(path, '/randomwalk_mispre_GLMM_sample_sizes.RDS'))
results.randomwalk.glmm$pvals$type <- "GLMM"
results.randomwalk.glmm$mles$type <- "GLMM"
results.spatial.lmm <- readRDS(paste0(path, '/spatial_mispre_LMM_sample_sizes.RDS'))
results.spatial.lmm$pvals$type <- "LMM"
results.spatial.lmm$mles$type <- "LMM"
results.spatial.glmm <- readRDS(paste0(path, '/spatial_pois-zip_GLMM_sample_sizes.RDS'))
results.spatial.glmm$pvals$type <- "GLMM"
results.spatial.glmm$mles$type <- "GLMM"
runtimes.all <- rbind(results.simpleGLMM.lmm$runtimes,
                      results.randomwalk.lmm$runtimes,
                      results.spatial.lmm$runtimes)
runtimes.all <- runtimes.all %>% filter(!is.na(med))
pvals.all <- rbind(results.simpleGLMM.lmm$pvals,
                   results.simpleGLMM.glmm$pvals,
                   results.randomwalk.lmm$pvals,
                   results.randomwalk.glmm$pvals,
                   results.spatial.lmm$pvals,
                   results.spatial.glmm$pvals)
mles.all <- rbind(results.simpleGLMM.lmm$mles,
                  results.simpleGLMM.glmm$mles,
                  results.randomwalk.lmm$mles,
                  results.randomwalk.glmm$mles,
                  results.spatial.lmm$mles,
                  results.spatial.glmm$mles)


runtimes.all$method <- factor(runtimes.all$method,
                              level = c(
                                'runtime.fg',
                                'runtime.osg',
                                'runtime.gen',
                                'runtime.cdf',
                                'runtime.mcmc',
                                'runtime.uncond',
                                'runtime.uncond_nrot',
                                'runtime.cond',
                                'runtime.cond_nrot'),
                              label = c(
                                'full Gaussian',
                                'one-step Gaussian',
                                'one-step Generic',
                                'cdf',
                                'MCMC',
                                "Unconditional ecdf, Rotated",
                                "Unconditional ecdf, Not Rotated",
                                "Conditional ecdf, Rotated",
                                "Conditional ecdf, Not Rotated"
                              ))
pvals.all$method <- factor(pvals.all$method,
                           level = c(
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


plot_ss_t1err_by_dim <- function(df){

    p <- df %>%
      filter(version == "h0") %>%
      group_by(nobs, method, model, type) %>%
      summarise(t1_err = sum(pvalue<0.05)/sum(pvalue>=0)) %>%
      ggplot(aes(x = nobs, y = t1_err, color = method)) +
      ylab("Type I Error Rate") +
      xlab("Number of Observations") +
      geom_point() + geom_line() +
      geom_hline(yintercept = 0.05, color = 'red', linetype = "dashed") +
      facet_grid(model~type, labeller = label_wrap_gen(14)) +
      theme_bw() +
      theme(legend.position = "bottom")

    print(p)
}

plot_ss_pow_by_dim <- function(df){

  p <- df %>%
    filter(version == "h1") %>%
    group_by(nobs, method, model, type) %>%
    summarise(power = sum(pvalue<=0.05)/sum(pvalue>=0)) %>%
    ggplot(aes(x = nobs, y = power, color = method)) +
      ylab("Power") +
      xlab("Number of Observations")  +
      geom_point() + geom_line() +
      geom_hline(yintercept = 0.95, color = 'red', linetype = "dashed") +
      facet_grid(model~type, labeller = label_wrap_gen(14)) +
      theme_bw() +
      theme(legend.position = "bottom")

  print(p)
}


plot_ss_runtimes_natural <- function(df){
  p <- ggplot(df,
              aes(nobs, med, ymin=lwr, ymax=upr,  color=method)) +
    facet_wrap(~model, scales = "free_y", ncol = 1) +
    geom_line() +
    labs(y='runtime (s)') +
    scale_color_discrete()  +
    theme_bw()
  print(p)
}
plot_ss_runtimes_log <- function(df){
  p <- ggplot(df,
              aes(nobs, med, ymin=lwr, ymax=upr,  color=method)) +
    facet_wrap(~model, ncol = 1) +
    geom_line() +
    geom_pointrange(fatten=2) +
    scale_y_log10()+labs(y='log(runtime (s))') +
    scale_color_discrete()  +
    theme_bw()
  print(p)
}
plot_ss_mles_by_dim <- function(df){
  p <- ggplot(df, aes(factor(nobs), bias)) +
    geom_violin() +
    geom_hline(yintercept=0, color='red') +
    facet_grid(h~par) + labs(y='Absolute error')
  print(p)
}


