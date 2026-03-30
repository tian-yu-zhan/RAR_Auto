

## functions for ten groups

z.test.known.var = function(trt.vec.in, pbo.vec.in, sd.in){
  
  # sd.pool = sqrt(((length(trt.vec.in)-1)*var.in + (length(pbo.vec.in)-1)*var.in)/
  #                  (length(pbo.vec.in)+length(trt.vec.in)-2))
  # t.stat = (mean(trt.vec.in)-mean(pbo.vec.in))/sd.pool/sqrt(1/length(trt.vec.in)+1/length(pbo.vec.in))
  z.stat = (mean(trt.vec.in)-mean(pbo.vec.in))/sd.in/sqrt(1/length(trt.vec.in)+1/length(pbo.vec.in))
  # p.out = pt(t.stat, lower.tail = FALSE, df = (length(trt.vec.in)+length(pbo.vec.in)-2))
  return(z.stat)
}

anova.func = function(pbo.vec.in, l.vec.in, h.vec.in){
  data.anova = data.frame("x" = c(pbo.vec.in, l.vec.in, h.vec.in),
                          "group" = c(rep("pbo", length(pbo.vec.in)), rep("l", length(l.vec.in)),
                                      rep("h", length(h.vec.in))))
  fit.anova = anova(lm(x ~ group, data.anova))
  p.anova = fit.anova$`Pr(>F)`[1]

  # t.anova = qt(1-p.anova/2, df =fit.anova$Df[2])
  # if (mean(c(l.vec.in, h.vec.in))<mean(pbo.vec.in)) t.anova = -t.anova
  #p.anova.one.sided = 1 - pt(t.anova, df =fit.anova$Df[2])

  p.anova.one.sided = p.anova/2
  if (mean(c(l.vec.in, h.vec.in))<mean(pbo.vec.in))  p.anova.one.sided = 1-p.anova/2
  return(qnorm(1-p.anova.one.sided))
  
  
  # data.dunttee <- data.frame("sample" = c(pbo.vec.in, l.vec.in, h.vec.in),
  #                            "group" = c(rep("A", length(pbo.vec.in)), rep("B", length(l.vec.in)),
  #                                        rep("C", length(h.vec.in))))
  # 
  # fit.dunttee <- aov(sample ~ group, data.dunttee)
  # 
  # fit.dunttee.mcp <- glht(fit.dunttee, linfct=mcp(group="Dunnett"),
  #                         alternative = "greater")
  # p.adj.dun = as.numeric(summary(fit.dunttee.mcp)$test$pvalues)
  # return(max(qnorm(1-p.adj.dun)))

  # p.anova.one.sided = p.anova/2
  # if (mean(c(l.vec.in, h.vec.in))<mean(pbo.vec.in))  p.anova.one.sided = 1-p.anova/2
  # return(qnorm(1-p.anova.one.sided))
  
}

SS.func = function(beta.in, delta.in, sd.pool.in){
  SS = 4*(qnorm(1-one.sided.alpha) + qnorm(1-beta.in))^2*sd.pool.in^2/(delta.in^2)
  return(SS)
}




