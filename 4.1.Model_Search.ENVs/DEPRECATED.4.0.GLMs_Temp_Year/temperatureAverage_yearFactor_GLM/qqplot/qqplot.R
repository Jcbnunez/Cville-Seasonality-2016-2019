### library
  library(data.table)
  library(ggplot2)

### function
  gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

### load GLM data
  load("/Users/alanbergland/glm_flight/glm.out.VA_ch_0.Rdata")


  glm.out.ag <- glm.out[chr!="X",list(chr=chr[order(p.lrt)], inv=(invName!="none")[order(p.lrt)],
                                      obs=-log10(sort(p.lrt)), exp=-log10(ppoints(length(p.lrt)))),
                          list(mod)]

  p1 <- ggplot(data=glm.out.ag[mod=="aveTemp+year_factor"], aes(x=exp, y=obs, group=inv, linetype=inv)) +
  geom_line() +
  geom_abline(aes(intercept=0, slope=1)) +
  ylim(0,10) + xlim(0,10) +
  facet_grid(inv~chr)

  glm.out.ag[mod=="aveTemp+year_factor"][, list(minp=max(exp)), list(chr, inv)]

  ggsave(p1, height=6, w=16, file="~/qqplot_glm.png")
