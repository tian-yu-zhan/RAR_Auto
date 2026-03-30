

## three treatment arms + pbo + SSR

# setwd("C:/Users/zhantx/Documents/collaborations/clinical studies/immunology/RZB HS/Cui2007/RAR_manuscript_Lu/results")

setwd("~/RAR/results/type_I_error/")

library(asd)
library(multcomp)
# library(agricolae)

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

M.SSR.func = function(yi.in, n.in, sd.pool.in, beta.SSR.in, w1.SSR.in){
  w2.SSR.in = sqrt(1-w1.SSR.in^2)
  M = 2*sd.pool.in^2/w2.SSR.in^2 * ((qnorm(1-0.025)+w2.SSR.in*qnorm(1-beta.SSR.in))/
                                      yi.in - n.in*w1.SSR.in/sqrt(2*n.in)/sd.pool.in)^2
  return(M)
}

n.ind = 5
n.arm = 4
output = matrix(NA, nrow = n.ind, ncol = 76)
colnames(output) = c("n", "m", "n.adj.initial", "n.adj.min", "n.adj.max", "t.SSR.check", 
                     "beta",
                     "n.look","RAR.ratio", "true_ES", "true_mean", "RAR_theta",
                     "reject_l_unweighted", 
                     "reject_m_unweighted", 
                     "reject_h_unweighted",
                     # "reject_s_unweighted",
                     # "reject_s1_unweighted",
                     # "reject_s2_unweighted",
                     "reject_lm_unweighted",
                     "reject_lh_unweighted",
                     "reject_mh_unweighted", 
                     "reject_l_unweighted_close", 
                     "reject_m_unweighted_close",
                     "reject_h_unweighted_close", 
                     "reject_s_unweighted_close", 
                     "reject_s_two_unweighted_close", 
                     
                     "reject_global",
                     
                     "reject_l_CHW", 
                     "reject_m_CHW", 
                     "reject_h_CHW",
                     # "reject_s_CHW",
                     # "reject_s1_CHW",
                     # "reject_s2_CHW",
                     "reject_lm_CHW",
                     "reject_lh_CHW",
                     "reject_mh_CHW", 
                     "reject_l_CHW_close", 
                     "reject_m_CHW_close",
                     "reject_h_CHW_close", 
                     "reject_lm_CHW_close", 
                     "reject_lh_CHW_close",
                     "reject_mh_CHW_close", 
                     "reject_s_CHW_close", 
                     "reject_s_two_CHW_close", 
                     
                     "reject_l_bonferroni", "reject_m_bonferroni", 
                     "reject_h_bonferroni", "reject_bonferroni", 
                     "reject_lm_bonferroni", "reject_lh_bonferroni",
                     "reject_mh_bonferroni",
                     # "reject_s_two_bonferroni",
                     "prob.bon.lm", "prob.bon.lh", "prob.bon.mh",
                     
                     "reject_l_dunnett", "reject_m_dunnett", 
                     "reject_h_dunnett", "reject_dunnett", 
                     "reject_lm_dun", "reject_lh_dun",
                     "reject_mh_dun",
                     # "reject_s_two_dun",
                     "prob.dun.lm", "prob.dun.lh", "prob.dun.mh",
                     
                     "n.pbo.out", "n.l.out", "n.m.out", 
                     "n.h.out","n.s.out", 
                     "n.lh.out", "n.lm.out", "n.mh.out", "n.s.two.out",
                     
                     "prob.h", "prob.m", "prob.l",
                     "prob.lm", "prob.lh", "prob.mh",
                     "ASN.total", "ASN.bon", "ASN.dun"
                     )
output = data.frame(output)
n.itt = 10^6
one.sided.alpha = 0.025

########### RAR adjustment threthold

for (ind in c(1)){
# for (ind in 2:5){

  # set.seed(1)
  RAR.equal.vec = c(1/4, 1/4, 1/4, 1/4)  
  RAR.vec = RAR.vec = c(7, 7, 1, 1)/16
  # RAR.vec = RAR.vec = c(3, 3, 1, 1)/8
  
  # delta.max = 0.4
  
  # n.adj.max = 16*4
  n.look = 12
  #########################################################
  
  # if (ind==1){
  #   RAR.vec = RAR.vec = c(8, 6, 1, 1)/16
  #   m = 200
  #   n.adj = 160
  #   n.look = 10
  #   n = m+n.adj*n.look
  #   mean.vec = c(0, 0, 0, 0)/sqrt(n/4/2)
  #   theta.in = 0.5
  # }
  
  if (ind==1){
    n.adj = 32
    m = 24
    n = m+n.adj*n.look
    mean.vec = c(0, 0, 0, 0)
    theta.in = 0.5
    n.adj.min = 16
    t.SSR.check = 6
    n.bon = 520
    n.dun = 508
    beta.SSR = 0.1
  }
  
  if (ind==2){
    n.adj = 32
    m = 24
    n = m+n.adj*n.look
    delta.max = 0.4
    mean.vec = c(0, 0.1, delta.max*3/4, delta.max*9/8)
    print(mean.vec)
    theta.in = 0.5
    n.adj.min = 32
    t.SSR.check = 6
    n.bon = 520
    n.dun = 508
    beta.SSR = 0.1
  }
  
  if (ind==3){
    n.adj = 32
    m = 24
    n = m+n.adj*n.look
    delta.max=0.4
    mean.vec = c(0, 0.5*delta.max, 1*delta.max, 1*delta.max)
    print(mean.vec)
    theta.in = 0.5
    n.adj.min = 16
    t.SSR.check = 6
    n.bon = 520
    n.dun = 508
    beta.SSR = 0.1
  }
  
  if (ind==4){
    n.adj = 32
    m = 24
    n = m+n.adj*n.look
    delta.max=0.37
    mean.vec = c(0, delta.max, delta.max, delta.max)
    print(mean.vec)
    theta.in = 0.5
    n.adj.min = 16
    t.SSR.check = 6
    n.bon = 520
    n.dun = 508
    beta.SSR = 0.1
  }
  
  if (ind==5){
    n.adj = 32
    m = 24
    n = m+n.adj*n.look
    delta.max=0.3
    mean.vec = c(0, 0.5*delta.max, delta.max, 1.5*delta.max)
    print(mean.vec)
    theta.in = 0.5
    n.adj.min = 32
    t.SSR.check = 6
    n.bon = 520
    n.dun = 508
    beta.SSR = 0.1
  }
 
  n.adj.max = floor((n.bon-m-(t.SSR.check-1)*n.adj)/(n.look-t.SSR.check+1))
  

  

  
  
  
  
  
  
  #################################################################################################
  
    p.out.l.unweighted = p.out.m.unweighted = p.out.h.unweighted =
    p.out.lm.unweighted = p.out.lh.unweighted = p.out.mh.unweighted =
    p.out.s.unweighted = p.out.s1.unweighted = p.out.s2.unweighted = 
    p.out.s.two.unweighted = 
    s.arm.out = s.two.arm.out = 
    s.arm.bon.out = s.arm.dun.out = 
    p.out.l.CHW = p.out.m.CHW = p.out.h.CHW =
    p.out.s.CHW = p.out.s1.CHW = p.out.s2.CHW =
    p.out.lm.CHW = p.out.lh.CHW = p.out.mh.CHW =
    p.out.s.two.CHW = 
    p.out.global = 
    prob.l.vec = prob.m.vec = prob.h.vec =
    prob.lm.vec = prob.lh.vec = prob.mh.vec =
    prob.bon.lm.vec = prob.bon.lh.vec = prob.bon.mh.vec =
    prob.dun.lm.vec = prob.dun.lh.vec = prob.dun.mh.vec =
    n.pbo.out = n.l.out = n.m.out = n.h.out = n.s.out = 
    n.lm.out = n.lh.out = n.mh.out = n.s.two.out = 
    dec.bon.vec = dec.bon.l.vec = dec.bon.m.vec = dec.bon.h.vec = 
    dec.dun.vec = dec.dun.l.vec = dec.dun.m.vec = dec.dun.h.vec =
    rep(0, n.itt)
    

    # se.pbo.out = se.l.out = se.h.out = se.s.out = rep(0, n.itt)
  
  output$RAR_theta[ind] = theta.in
  output$beta[ind] = beta.SSR
  
  output$n[ind] = n
  output$m[ind] = m
  output$n.adj.initial[ind] = n.adj
  output$n.adj.max[ind] = n.adj.max
  output$n.adj.min[ind] = n.adj.min
  
  output$t.SSR.check[ind] = t.SSR.check
  # n.look = (n-m)/n.adj
  output$n.look[ind] = n.look
  z.cohort.array = array(NA, dim = c(n.look+1, 2*(n.arm-1)+2, n.itt))
  
  for (itt in 1:n.itt){
    print(itt)
    n.vec.initial = RAR.equal.vec*m
    
    pbo.vec = rnorm(n.vec.initial[1], mean.vec[1], 1)
    l.vec = rnorm(n.vec.initial[2], mean.vec[2], 1)
    m.vec = rnorm(n.vec.initial[3], mean.vec[3], 1)
    h.vec = rnorm(n.vec.initial[4], mean.vec[4], 1)
    
    z.cohort.mat = matrix(NA, nrow = (n.look+1), ncol = 2*(n.arm-1)+2)
    
    # z.cohort.mat[1, ] = c(m,
    #                       qnorm(1-t.test(l.vec, pbo.vec, alternative = "greater", 
    #                                      var.equal = TRUE)$p.value),
    #                       qnorm(1-t.test(m.vec, pbo.vec, alternative = "greater", 
    #                                      var.equal = TRUE)$p.value),
    #                       qnorm(1-t.test(h.vec, pbo.vec, alternative = "greater", 
    #                                      var.equal = TRUE)$p.value),
    #                       qnorm(1-t.test(c(l.vec, m.vec), pbo.vec, alternative = "greater", 
    #                                      var.equal = FALSE)$p.value),
    #                       qnorm(1-t.test(c(l.vec, h.vec), pbo.vec, alternative = "greater", 
    #                                      var.equal = FALSE)$p.value),
    #                       qnorm(1-t.test(c(m.vec, h.vec), pbo.vec, alternative = "greater", 
    #                                      var.equal = FALSE)$p.value),
    #                       qnorm(1-t.test(c(l.vec, m.vec, h.vec), pbo.vec, alternative = "greater", 
    #                                      var.equal = FALSE)$p.value)
    #                       )
    
    z.cohort.mat[1, ] = c(m,
                          z.test.known.var(l.vec, pbo.vec,1),
                          z.test.known.var(m.vec, pbo.vec,1),
                          z.test.known.var(h.vec, pbo.vec,1),
                          z.test.known.var(c(l.vec,m.vec), pbo.vec,1),
                          z.test.known.var(c(l.vec,h.vec), pbo.vec,1),
                          z.test.known.var(c(m.vec,h.vec), pbo.vec,1),
                          z.test.known.var(c(l.vec, m.vec,h.vec), pbo.vec,1)
    )
    
    n.adj.current = n.adj
    
    for (i in 1:n.look){
    
      # if (i==t.SSR.check) break 
      
      RAR.measure = c(z.test.known.var(l.vec, pbo.vec, 1),
                      z.test.known.var(m.vec, pbo.vec, 1),
                      z.test.known.var(h.vec, pbo.vec, 1)
      )
      
        # RAR.measure = c(t.test(l.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$statistic,
        #                 t.test(m.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$statistic,
        #                 t.test(h.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$statistic)
        
        # RAR.measure = c(mean(l.vec)*sqrt(length(l.vec)),
        #                 mean(m.vec)*sqrt(length(m.vec)),
        #                 mean(h.vec)*sqrt(length(h.vec))
        # )
        
        RAR.vec.adj = c(RAR.vec[1], (RAR.vec[-1])[match(RAR.measure, sort(RAR.measure, decreasing = TRUE))])
        
      ## SSR  
      if (i>=t.SSR.check){
        
        if (which.max(RAR.measure)==1) s.interim.vec = l.vec
        if (which.max(RAR.measure)==2) s.interim.vec = m.vec
        if (which.max(RAR.measure)==3) s.interim.vec = h.vec
        
        SS.total = SS.func(beta.in = beta.SSR,
                           sd.pool.in = mean(c(sd(s.interim.vec), sd(pbo.vec))),
                           delta.in = mean(s.interim.vec) - mean(pbo.vec))
        
        # SS.total = length(c(s.interim.vec, pbo.vec))+
        #   2*M.SSR.func(yi.in = mean(s.interim.vec) - mean(pbo.vec),
        #            n.in = round(length(c(s.interim.vec, pbo.vec))/2),
        #            sd.pool.in = 1,
        #            beta.SSR.in = 0.1, 
        #            w1.SSR.in = sqrt((n.adj*(i-1)+m)/n))
        # stop("a")
        
        SS.need = SS.total - length(s.interim.vec) - length(pbo.vec)
        n.adj.ratio = SS.need / ((n.look-i+1)*sum(RAR.vec[1:2]*n.adj.min))
        n.adj.ratio = max(1, ceiling(min(floor(n.adj.max/n.adj.min), n.adj.ratio)))
        
        n.adj.current = n.adj.min*n.adj.ratio
        # n.adj.current = max(n.adj.min, floor(n.adj.ratio*n.adj.min))
      }
        
        ## RAR chunk
        n.add.vec = RAR.vec.adj*n.adj.current
        # print(n.add.vec)
        # n.add.vec = pmax(1, round(RAR.vec.adj*n.adj.current))
        
      
      n.pbo.add = n.add.vec[1]
      n.l.add = n.add.vec[2]
      n.m.add = n.add.vec[3]
      n.h.add = n.add.vec[4]
  
      
      pbo.new.vec = rnorm(n.pbo.add, mean.vec[1], 1)
      l.new.vec = rnorm(n.l.add, mean.vec[2], 1)
      m.new.vec = rnorm(n.m.add, mean.vec[3], 1)
      h.new.vec = rnorm(n.h.add, mean.vec[4], 1)
      
      if (n.adj>1){
        ## for combination function
          # z.value.l.in = qnorm(1-t.test(l.new.vec, pbo.new.vec, alternative = "greater", 
          #                               var.equal = TRUE)$p.value)
          # 
          # z.value.m.in = qnorm(1-t.test(m.new.vec, pbo.new.vec, alternative = "greater", 
          #                               var.equal = TRUE)$p.value)
          # 
          # z.value.h.in = qnorm(1-t.test(h.new.vec, pbo.new.vec, alternative = "greater", 
          #                               var.equal = TRUE)$p.value)
          # 
          # z.value.lm.in = qnorm(1-t.test(c(l.new.vec, m.new.vec), pbo.new.vec, alternative = "greater", 
          #                               var.equal = FALSE)$p.value)
          # 
          # z.value.lh.in = qnorm(1-t.test(c(l.new.vec, h.new.vec), pbo.new.vec, alternative = "greater", 
          #                                var.equal = FALSE)$p.value)
          # 
          # z.value.mh.in = qnorm(1-t.test(c(m.new.vec, h.new.vec), pbo.new.vec, alternative = "greater", 
          #                                var.equal = FALSE)$p.value)
          # 
          # z.value.global.in = qnorm(1-t.test(c(l.new.vec, m.new.vec, h.new.vec), pbo.new.vec, alternative = "greater", 
          #                                var.equal = FALSE)$p.value)
        
        z.value.l.in = z.test.known.var(l.new.vec, pbo.new.vec, 1)
        
        z.value.m.in = z.test.known.var(m.new.vec, pbo.new.vec, 1)
        
        z.value.h.in = z.test.known.var(h.new.vec, pbo.new.vec, 1)
        
        z.value.lm.in = z.test.known.var(c(l.new.vec,m.new.vec), pbo.new.vec, 1)
        
        z.value.lh.in = z.test.known.var(c(l.new.vec,h.new.vec), pbo.new.vec, 1)
        
        z.value.mh.in = z.test.known.var(c(m.new.vec,h.new.vec), pbo.new.vec, 1)
        
        z.value.global.in = z.test.known.var(c(l.new.vec, m.new.vec,h.new.vec), pbo.new.vec, 1)
        
        z.cohort.mat[1+i, ] = c(n.adj,
                               z.value.l.in, z.value.m.in, z.value.h.in, z.value.lm.in, 
                               z.value.lh.in, z.value.mh.in, z.value.global.in)
      }
      
      pbo.vec = c(pbo.vec, pbo.new.vec)
      l.vec = c(l.vec, l.new.vec)
      m.vec = c(m.vec, m.new.vec)
      h.vec = c(h.vec, h.new.vec)
      
    }
    
    RAR.measure = c(z.test.known.var(l.vec, pbo.vec, 1),
                    z.test.known.var(m.vec, pbo.vec, 1),
                    z.test.known.var(h.vec, pbo.vec, 1)
    )
    
    # RAR.measure = c(mean(l.vec)*sqrt(length(l.vec)),
    #                 mean(m.vec)*sqrt(length(m.vec)),
    #                 mean(h.vec)*sqrt(length(h.vec))
    # )
    
    # RAR.measure = c(t.test(l.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$statistic,
    #                 t.test(m.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$statistic,
    #                 t.test(h.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$statistic)
    
    z.cohort.array[,,itt] = z.cohort.mat
    
    comb.func = function(n.look.in, s.z.vec.in, s.n.vec.in){
      if (n.look.in==0){
        z.com = s.z.vec.in[1]
      } else {
        for (i in 1:n.look.in){
          if (i==1){
            z1 = s.z.vec.in[n.look.in]
            z2 = s.z.vec.in[n.look.in+1]
            
          } else {
            z1 = s.z.vec.in[n.look.in+1-i]
            z2 = z.com
          }
          
          w1 = sqrt(s.n.vec.in[n.look.in-i+1]/sum(s.n.vec.in[(n.look.in-i+1):(n.look.in+1)]))
          w2 = sqrt(1-w1^2)
          
          z.com = w1*z1 + w2*z2
        }
      }
      
    
      return(z.com)
    }
    
    if (n.adj>1){
      l.z.vec = z.cohort.mat[, 2]
      m.z.vec = z.cohort.mat[, 3]
      h.z.vec = z.cohort.mat[, 4]
      
      lm.z.vec = z.cohort.mat[, 5]
      lh.z.vec = z.cohort.mat[, 6]
      mh.z.vec = z.cohort.mat[, 7]
      
      global.z.vec = z.cohort.mat[, 8]
      
      s.n.vec = z.cohort.mat[,1]
      
      z.com.l = comb.func(sum(!is.na(l.z.vec))-1, l.z.vec[!is.na(l.z.vec)], s.n.vec[!is.na(l.z.vec)])
      z.com.m = comb.func(sum(!is.na(m.z.vec))-1, m.z.vec[!is.na(m.z.vec)], s.n.vec[!is.na(m.z.vec)])
      z.com.h = comb.func(sum(!is.na(h.z.vec))-1, h.z.vec[!is.na(h.z.vec)], s.n.vec[!is.na(h.z.vec)])
      
      z.com.lm = comb.func(sum(!is.na(lm.z.vec))-1, lm.z.vec[!is.na(lm.z.vec)], s.n.vec[!is.na(lm.z.vec)])
      z.com.lh = comb.func(sum(!is.na(lh.z.vec))-1, lh.z.vec[!is.na(lh.z.vec)], s.n.vec[!is.na(lh.z.vec)])
      z.com.mh = comb.func(sum(!is.na(mh.z.vec))-1, mh.z.vec[!is.na(mh.z.vec)], s.n.vec[!is.na(mh.z.vec)])
      
      z.com.global = comb.func(sum(!is.na(global.z.vec))-1, global.z.vec[!is.na(global.z.vec)], 
                               s.n.vec[!is.na(global.z.vec)])
      
    } else {
      z.com.l = z.com.m = z.com.h = z.com.lm = z.com.lh = z.com.mh = z.com.global = NA
    }
    
    # data to estimate variance
    p.out.l.unweighted[itt] = t.test(l.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$p.value
    p.out.m.unweighted[itt] = t.test(m.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$p.value
    p.out.h.unweighted[itt] = t.test(h.vec, pbo.vec, alternative = "greater", var.equal = TRUE)$p.value
    
    p.out.lm.unweighted[itt] = t.test(c(l.vec, m.vec), pbo.vec, alternative = "greater", 
                                      var.equal = FALSE)$p.value
    p.out.lh.unweighted[itt] = t.test(c(l.vec, h.vec), pbo.vec, alternative = "greater", 
                                      var.equal = FALSE)$p.value
    p.out.mh.unweighted[itt] = t.test(c(m.vec, h.vec), pbo.vec, alternative = "greater", 
                                      var.equal = FALSE)$p.value
    
    p.out.global[itt] = 1-pnorm(z.com.global)
    
    p.out.l.CHW[itt] = 1-pnorm(z.com.l)
    p.out.m.CHW[itt] = 1-pnorm(z.com.m)
    p.out.h.CHW[itt] = 1-pnorm(z.com.h)
    
    p.out.lm.CHW[itt] = 1-pnorm(z.com.lm)
    p.out.lh.CHW[itt] = 1-pnorm(z.com.lh)
    p.out.mh.CHW[itt] = 1-pnorm(z.com.mh)
    
    ###################################################
    ## select one dose
    RAR.measure.re = RAR.measure
    if (sum(duplicated(RAR.measure))>0){
      RAR.measure.re = RAR.measure.re + c(0.0002, 0.0001, 0)
    }
    
    max.ind = as.numeric(which.max(RAR.measure.re))
    if (max.ind==1){
      s.arm = "l"
      s.vec = l.vec
      p.out.s.unweighted[itt] = p.out.l.unweighted[itt]
      p.out.s1.unweighted[itt] = p.out.lm.unweighted[itt]
      p.out.s2.unweighted[itt] = p.out.lh.unweighted[itt]
      
      p.out.s.CHW[itt] = p.out.l.CHW[itt]
      p.out.s1.CHW[itt] = p.out.lm.CHW[itt]
      p.out.s2.CHW[itt] = p.out.lh.CHW[itt]
      
    } else if (max.ind==2){
      s.arm = "m"
      s.vec = m.vec
      p.out.s.unweighted[itt] = p.out.m.unweighted[itt]
      p.out.s1.unweighted[itt] = p.out.lm.unweighted[itt]
      p.out.s2.unweighted[itt] = p.out.mh.unweighted[itt]
      
      p.out.s.CHW[itt] = p.out.m.CHW[itt]
      p.out.s1.CHW[itt] = p.out.lm.CHW[itt]
      p.out.s2.CHW[itt] = p.out.mh.CHW[itt]
      
    } else if (max.ind==3){
      s.arm = "h"
      s.vec = h.vec
      p.out.s.unweighted[itt] = p.out.h.unweighted[itt]
      p.out.s1.unweighted[itt] = p.out.lh.unweighted[itt]
      p.out.s2.unweighted[itt] = p.out.mh.unweighted[itt]
      
      p.out.s.CHW[itt] = p.out.h.CHW[itt]
      p.out.s1.CHW[itt] = p.out.lh.CHW[itt]
      p.out.s2.CHW[itt] = p.out.mh.CHW[itt]
      
    }
    s.arm.out[itt] = s.arm

    ###################################################
    ## select two doses
    min.ind = as.numeric(which.min(RAR.measure.re))
    if (min.ind==1){
      
      s.two.arm = "mh"
      p.out.s.two.unweighted[itt] = p.out.mh.unweighted[itt]
      p.out.s.two.CHW[itt] = p.out.mh.CHW[itt]
      n.s.two.out[itt] = length(c(m.vec, h.vec))
      
    } else if (min.ind==2){
      
      s.two.arm = "lh"
      p.out.s.two.unweighted[itt] = p.out.lh.unweighted[itt]
      p.out.s.two.CHW[itt] = p.out.lh.CHW[itt]
      n.s.two.out[itt] = length(c(l.vec, h.vec))
      
    } else if (min.ind==3){

      s.two.arm = "lm"
      p.out.s.two.unweighted[itt] = p.out.lm.unweighted[itt]
      p.out.s.two.CHW[itt] = p.out.lm.CHW[itt]
      n.s.two.out[itt] = length(c(l.vec, m.vec))
      
    }
    s.two.arm.out[itt] = s.two.arm
    
    n.mh.out[itt] = length(c(m.vec, h.vec))
    n.lh.out[itt] = length(c(l.vec, h.vec))
    n.lm.out[itt] = length(c(l.vec, m.vec))
    #######################################################
    
    prob.l.vec[itt] = as.numeric(s.arm=="l")
    prob.m.vec[itt] = as.numeric(s.arm=="m")
    prob.h.vec[itt] = as.numeric(s.arm=="h")
    
    prob.lm.vec[itt] = as.numeric(s.two.arm=="lm")
    prob.lh.vec[itt] = as.numeric(s.two.arm=="lh")
    prob.mh.vec[itt] = as.numeric(s.two.arm=="mh")
    
    n.pbo.out[itt] = length(pbo.vec)
    n.l.out[itt] = length(l.vec)
    n.m.out[itt] = length(m.vec)
    n.h.out[itt] = length(h.vec)
    n.s.out[itt] = length(s.vec)
    
    
    ## Bonferroni and Dunttee
    n.pbo.equal = round(n.bon/4)
    n.l.equal = round(n.bon/4)
    n.m.equal = round(n.bon/4)
    n.h.equal = n.bon - n.pbo.equal - n.l.equal - n.m.equal
    
    sample.pbo.equal = rnorm(n.pbo.equal, mean = mean.vec[1], sd = 1)
    sample.l.equal = rnorm(n.l.equal, mean = mean.vec[2], sd = 1)
    sample.m.equal = rnorm(n.m.equal, mean = mean.vec[3], sd = 1)
    sample.h.equal = rnorm(n.h.equal, mean = mean.vec[4], sd = 1)
    
    p.unadj.l = pnorm(z.test.known.var(sample.l.equal, sample.pbo.equal, 1),lower.tail = FALSE)
    p.unadj.m = pnorm(z.test.known.var(sample.m.equal, sample.pbo.equal, 1),lower.tail = FALSE)
    p.unadj.h = pnorm(z.test.known.var(sample.h.equal, sample.pbo.equal, 1),lower.tail = FALSE)
    
    # p.unadj.l = t.test(sample.l.equal, sample.pbo.equal, alternative = "greater",
    #                    var.equal = TRUE)$p.value
    # p.unadj.m = t.test(sample.m.equal, sample.pbo.equal, alternative = "greater",
    #                    var.equal = TRUE)$p.value
    # p.unadj.h = t.test(sample.h.equal, sample.pbo.equal, alternative = "greater",
    #                    var.equal = TRUE)$p.value
    
    
    p.adj.bon = p.adjust(c(p.unadj.l, p.unadj.m, p.unadj.h), "bonferroni")
    dec.bon.vec[itt] = as.numeric(sum(p.adj.bon<=one.sided.alpha)>0)
    dec.bon.l.vec[itt] = as.numeric(p.adj.bon[1]<=one.sided.alpha)
    dec.bon.m.vec[itt] = as.numeric(p.adj.bon[2]<=one.sided.alpha)
    dec.bon.h.vec[itt] = as.numeric(p.adj.bon[3]<=one.sided.alpha)
    
    p.adj.bon.re = p.adj.bon
    if (sum(duplicated(p.adj.bon))>0){
      p.adj.bon.re = p.adj.bon.re + c(0, 0.0001, 0.0002)
    }
    
    prob.bon.lm.vec[itt] = as.numeric(which.max(p.adj.bon.re)==3)
    prob.bon.lh.vec[itt] = as.numeric(which.max(p.adj.bon.re)==2)
    prob.bon.mh.vec[itt] = as.numeric(which.max(p.adj.bon.re)==1)
    
    s.arm.bon.out[itt] = which.min(p.adj.bon.re)
    
    # Dunttee
    n.pbo.equal = round(n.dun/4)
    n.l.equal = round(n.dun/4)
    n.m.equal = round(n.dun/4)
    n.h.equal = n.dun - n.pbo.equal - n.l.equal - n.m.equal
    
    sample.pbo.equal = rnorm(n.pbo.equal, mean = mean.vec[1], sd = 1)
    sample.l.equal = rnorm(n.l.equal, mean = mean.vec[2], sd = 1)
    sample.m.equal = rnorm(n.m.equal, mean = mean.vec[3], sd = 1)
    sample.h.equal = rnorm(n.h.equal, mean = mean.vec[4], sd = 1)
    
    data.dunttee <- data.frame("sample" = c(sample.pbo.equal, sample.l.equal, sample.m.equal, sample.h.equal),
                       "group" = c(rep("A", n.pbo.equal), rep("B", n.l.equal), rep("C", n.m.equal),
                                   rep("D", n.h.equal)))
    
    fit.dunttee <- aov(sample ~ group, data.dunttee)
    
    fit.dunttee.mcp <- glht(fit.dunttee, linfct=mcp(group="Dunnett"),
                            alternative = "greater")
    p.adj.dun = as.numeric(summary(fit.dunttee.mcp, test=(adjusted(type = "free")))$test$pvalues)
    dec.dun.vec[itt] = as.numeric(sum(p.adj.dun<=one.sided.alpha)>0)
    dec.dun.l.vec[itt] = as.numeric(p.adj.dun[1]<=one.sided.alpha)
    dec.dun.m.vec[itt] = as.numeric(p.adj.dun[2]<=one.sided.alpha)
    dec.dun.h.vec[itt] = as.numeric(p.adj.dun[3]<=one.sided.alpha)
    
    p.adj.dun.re = p.adj.dun
    if (sum(duplicated(p.adj.dun))>0){
      p.adj.dun.re = p.adj.dun.re + c(0, 0.0001, 0.0002)
    }
    
    s.arm.dun.out[itt] = which.min(p.adj.dun.re)
    
    prob.dun.lm.vec[itt] = as.numeric(which.max(p.adj.dun.re)==3)
    prob.dun.lh.vec[itt] = as.numeric(which.max(p.adj.dun.re)==2)
    prob.dun.mh.vec[itt] = as.numeric(which.max(p.adj.dun.re)==1)

  }
  
  print(max(n.pbo.out+n.l.out+n.m.out+n.h.out))  
  
  output$RAR.ratio[ind] = paste(RAR.vec*n.adj,collapse = ": ")
  output$true_ES[ind] = paste(round(mean.vec*sqrt(n/4/2),2),collapse = ", ")
  output$true_mean[ind] = paste(round(mean.vec,3),collapse = ", ")
  
  ############ unweighted statistics
  output$reject_l_unweighted[ind] = mean(p.out.l.unweighted<=one.sided.alpha)
  output$reject_m_unweighted[ind] = mean(p.out.m.unweighted<=one.sided.alpha)
  output$reject_h_unweighted[ind] = mean(p.out.h.unweighted<=one.sided.alpha)
  # output$reject_s_unweighted[ind] = mean(p.out.s.unweighted<=0.025)
  
  output$reject_lm_unweighted[ind] = mean(p.out.lm.unweighted<=one.sided.alpha)
  output$reject_lh_unweighted[ind] = mean(p.out.lh.unweighted<=one.sided.alpha)
  output$reject_mh_unweighted[ind] = mean(p.out.mh.unweighted<=one.sided.alpha)
  # output$reject_s1_unweighted[ind] = mean(p.out.s1.unweighted<=0.025)
  # output$reject_s2_unweighted[ind] = mean(p.out.s2.unweighted<=0.025)
  
  output$reject_global[ind] = mean(p.out.global<=one.sided.alpha)

  output$reject_l_unweighted_close[ind] = sum(
    (p.out.global<=one.sided.alpha)&(p.out.lm.unweighted<=one.sided.alpha)&
      (p.out.lh.unweighted<=one.sided.alpha)&(p.out.l.unweighted<=one.sided.alpha)&
      (s.arm.out=="l"))/n.itt
  
  output$reject_m_unweighted_close[ind] = sum(
    (p.out.global<=one.sided.alpha)&(p.out.lm.unweighted<=one.sided.alpha)&
      (p.out.mh.unweighted<=one.sided.alpha)&(p.out.m.unweighted<=one.sided.alpha)&
      (s.arm.out=="m"))/n.itt
  
  output$reject_h_unweighted_close[ind] = sum(
    (p.out.global<=one.sided.alpha)&(p.out.mh.unweighted<=one.sided.alpha)&
      (p.out.lh.unweighted<=one.sided.alpha)&(p.out.h.unweighted<=one.sided.alpha)&
      (s.arm.out=="h"))/n.itt
  
  output$reject_s_unweighted_close[ind] = mean(
    (p.out.global<=one.sided.alpha)&(p.out.s1.unweighted<=one.sided.alpha)&
      (p.out.s2.unweighted<=one.sided.alpha)&(p.out.s.unweighted<=one.sided.alpha)
  )
  
  output$reject_s_two_unweighted_close[ind] = mean(
    (p.out.global<=one.sided.alpha)&(p.out.s.two.unweighted<=one.sided.alpha)
  )
  

  
  ############## CHW
  output$reject_l_CHW[ind] = mean(p.out.l.CHW<=one.sided.alpha)
  output$reject_m_CHW[ind] = mean(p.out.m.CHW<=one.sided.alpha)
  output$reject_h_CHW[ind] = mean(p.out.h.CHW<=one.sided.alpha)
  # output$reject_s_CHW[ind] = mean(p.out.s.CHW<=0.025)
  
  output$reject_lm_CHW[ind] = mean(p.out.lm.CHW<=one.sided.alpha)
  output$reject_lh_CHW[ind] = mean(p.out.lh.CHW<=one.sided.alpha)
  output$reject_mh_CHW[ind] = mean(p.out.mh.CHW<=one.sided.alpha)
  # output$reject_s1_CHW[ind] = mean(p.out.s1.CHW<=0.025)
  # output$reject_s2_CHW[ind] = mean(p.out.s2.CHW<=0.025)
  
  select.l.vec =  (p.out.global<=one.sided.alpha)&(p.out.lm.CHW<=one.sided.alpha)&
    (p.out.lh.CHW<=one.sided.alpha)&(p.out.l.CHW<=one.sided.alpha)&
    (s.arm.out=="l")
  
  output$reject_l_CHW_close[ind] = sum(select.l.vec)/n.itt
  
  select.m.vec =  (p.out.global<=one.sided.alpha)&(p.out.lm.CHW<=one.sided.alpha)&
    (p.out.mh.CHW<=one.sided.alpha)&(p.out.m.CHW<=one.sided.alpha)&
    (s.arm.out=="m")
  output$reject_m_CHW_close[ind] = sum(select.m.vec
   )/n.itt
  
  select.h.vec =  (p.out.global<=one.sided.alpha)&(p.out.mh.CHW<=one.sided.alpha)&
    (p.out.lh.CHW<=one.sided.alpha)&(p.out.h.CHW<=one.sided.alpha)&
    (s.arm.out=="h")
  output$reject_h_CHW_close[ind] = sum(select.h.vec
     )/n.itt
  
  select.s.vec = (p.out.global<=one.sided.alpha)&(p.out.s1.CHW<=one.sided.alpha)&
    (p.out.s2.CHW<=one.sided.alpha)&(p.out.s.CHW<=one.sided.alpha)
  
  output$reject_s_CHW_close[ind] = mean(
    select.s.vec
  )
  
  output$reject_s_two_CHW_close[ind] = mean(
    (p.out.global<=one.sided.alpha)&(p.out.s.two.CHW<=one.sided.alpha)
  )
  
  output$reject_lm_CHW_close[ind] = sum(
    (p.out.global<=one.sided.alpha)&(p.out.lm.CHW<=one.sided.alpha)&
      (s.two.arm.out=="lm"))/n.itt
  
  output$reject_lh_CHW_close[ind] = sum(
    (p.out.global<=one.sided.alpha)&(p.out.lh.CHW<=one.sided.alpha)&(s.two.arm.out=="lh"))/n.itt
  
  output$reject_mh_CHW_close[ind] = sum(
    (p.out.global<=one.sided.alpha)&(p.out.mh.CHW<=one.sided.alpha)&(s.two.arm.out=="mh"))/n.itt
  
  ###########################################################
  
 
  output$n.pbo.out[ind] = mean(n.pbo.out)
  output$n.l.out[ind] = mean(n.l.out[select.l.vec==1])
  output$n.m.out[ind] = mean(n.m.out[select.m.vec==1])
  output$n.h.out[ind] = mean(n.h.out[select.h.vec==1])
  output$n.s.out[ind] = mean(n.s.out[select.s.vec==1])
  
  output$ASN.total[ind] = mean(n.pbo.out + n.l.out + n.m.out + n.h.out)
  
  output$n.lh.out[ind] = mean(n.lh.out[s.two.arm.out=="lh"])
  output$n.lm.out[ind] = mean(n.lm.out[s.two.arm.out=="lm"])
  output$n.mh.out[ind] = mean(n.mh.out[s.two.arm.out=="mh"])
  output$n.s.two.out[ind] = mean(n.s.two.out)
  
  ## Bonferroni and Dunttee
  output$reject_bonferroni[ind] = mean(dec.bon.vec)
  output$reject_l_bonferroni[ind] = sum(dec.bon.l.vec[s.arm.bon.out==1])/n.itt
  output$reject_m_bonferroni[ind] = sum(dec.bon.m.vec[s.arm.bon.out==2])/n.itt
  output$reject_h_bonferroni[ind] = sum(dec.bon.h.vec[s.arm.bon.out==3])/n.itt
  
  output$reject_lm_bonferroni[ind] = mean(pmax(dec.bon.l.vec, dec.bon.m.vec))
  output$reject_lh_bonferroni[ind] = mean(pmax(dec.bon.l.vec, dec.bon.h.vec))
  output$reject_mh_bonferroni[ind] = mean(pmax(dec.bon.h.vec, dec.bon.m.vec))
  
  output$ASN.bon[ind] = n.bon
  
  # output$reject_s_two_bonferroni[ind] = mean((dec.bon.l.vec+dec.bon.m.vec+dec.bon.h.vec)>=2)
  
  output$prob.bon.lm[ind] = mean(prob.bon.lm.vec)
  output$prob.bon.lh[ind] = mean(prob.bon.lh.vec)
  output$prob.bon.mh[ind] = mean(prob.bon.mh.vec)
  
  # output$reject_s_bonferroni[ind] = mean(dec.bon.s.vec)
  
  output$reject_dunnett[ind] = mean(dec.dun.vec)
  output$reject_l_dunnett[ind] = sum(dec.dun.l.vec[s.arm.dun.out==1])/n.itt
  output$reject_m_dunnett[ind] = sum(dec.dun.m.vec[s.arm.dun.out==2])/n.itt
  output$reject_h_dunnett[ind] = sum(dec.dun.h.vec[s.arm.dun.out==3])/n.itt
  # output$reject_s_dunnett[ind] = mean(dec.dun.s.vec)
  
  output$ASN.dun[ind] = n.dun
  
  output$reject_lm_dun[ind] = mean(pmax(dec.dun.l.vec, dec.dun.m.vec))
  output$reject_lh_dun[ind] = mean(pmax(dec.dun.l.vec, dec.dun.h.vec))
  output$reject_mh_dun[ind] = mean(pmax(dec.dun.h.vec, dec.dun.m.vec))
  
  # output$reject_s_two_dun[ind] = mean(pmax(pmax(dec.dun.l.vec, dec.dun.m.vec),
  #                                         pmax(dec.dun.l.vec, dec.dun.h.vec),
  #                                         pmax(dec.dun.h.vec, dec.dun.m.vec)
  #                                         )==1)
  
  
  output$prob.dun.lm[ind] = mean(prob.dun.lm.vec)
  output$prob.dun.lh[ind] = mean(prob.dun.lh.vec)
  output$prob.dun.mh[ind] = mean(prob.dun.mh.vec)
  
  output$prob.h[ind] = mean(prob.h.vec)
  output$prob.m[ind] = mean(prob.m.vec)
  output$prob.l[ind] = mean(prob.l.vec)
  
  output$prob.lm[ind] = mean(prob.lm.vec)
  output$prob.lh[ind] = mean(prob.lh.vec)
  output$prob.mh[ind] = mean(prob.mh.vec)
}


write.csv(output, "results_three_doses_SSR.csv")
# print(output)





