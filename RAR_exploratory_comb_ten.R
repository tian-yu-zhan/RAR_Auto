

## five treatment arms + pbo


# setwd("C:/Users/zhantx/Documents/collaborations/clinical studies/immunology/RZB HS/Cui2007/RAR_manuscript_Lu/results")

setwd("~/RAR/results/manu/")

library(asd)
library(multcomp)
# library(agricolae)
library(doParallel)

n.cluster = 3

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

n.arm = 10
n.ind = 4
output = matrix(NA, nrow = n.ind, ncol = 29+4*n.arm)
colnames(output) = c("n", "m", "n.adj.initial", "n.adj.min", "n.adj.max",
                     "t.SSR.check", "beta.SSR",
                     "n.look","RAR.ratio", "true_ES", "true_mean", "RAR_theta",
                     "reject_bonferroni", "reject_dunnett", 
                     "reject_bonferroni_four", "reject_dunnett_four", 
                     "reject_bonferroni_two", "reject_dunnett_two", 
                     "reject_global",
                     as.vector(paste0("reject_", 1:n.arm)),
                     as.vector(paste0("reject_bon_", 1:n.arm)),
                     as.vector(paste0("reject_dun_", 1:n.arm)),
                     "reject_s",
                     "reject_s_two",
                     "reject_s_four",
                     "n_bon", "n_dun", "ASN_total",
                     "ASN_pbo", 
                     as.vector(paste0("ASN_", 1:n.arm)),
                     "ASN_s_two",
                     "ASN_s_four",
                     "ASN_s")
output = data.frame(output)
n.itt = 10^4
one.sided.alpha = 0.025

########### RAR adjustment threthold

for (ind in c(1:4)){
  set.seed(2)
# for (ind in 2){
  
  #########################################################
  
  n.adj.min = 12
  # n.adj.max = 120
  t.SSR.check = 4
  sd.common = 1
  
  
  # RAR.vec = c(13, 14, 1, 1, 1)/30
  # RAR.vec = c(13, 8, 7, 1, 1)/30
  # RAR.equal.vec = c(14, 4, 4, 4, 4)/30
  
  
  
  # RAR.equal.vec = c(25, rep(3, n.arm))/55
  RAR.equal.vec = c(rep(5, n.arm+1))/55
  # RAR.vec = c(20, c(23, 23, 24, 24), rep(1,6))/120
  
  # if (ind==1){
  #   m = 55
  #   n.adj = 120
  #   n.look = 11
  #   # n.adj.last = n.bon - m-n.adj*(n.look-1)
  #   n = m+n.adj*n.look
  #   mean.vec = rep(0, n.arm+1)
  #   theta.in = 0.5
  #   beta.SSR = NA
  #   SSR.ind = FALSE
  #   RAR.vec = c(20, 23, 23, 24, 24, rep(1,6))/120
  # }
  
  if (ind==1){
    m = 55+121
    n.adj = 120
    n.look = 10
    mean.vec = c(0, rep(0, 10))
    RAR.vec = c(c(40, 40, 40), rep(0,8))/120
  }
  
  if (ind==2){
    m = 55+121
    n.adj = 120
    n.look = 10
    mean.vec = c(0, rep(0.32, 4), rep(0.3, 4), rep(0.27, 2))
    RAR.vec = c(c(40, 40, 40), rep(0,8))/120
  }
  
  if (ind==3){
    m = 55+121
    n.adj = 120
    n.look = 10
    mean.vec = c(0, rep(0, 10))
    RAR.vec = c(rep(24, 5), rep(0,6))/120
  }
  
  if (ind==4){
    m = 55+121
    n.adj = 120
    n.look = 10
    mean.vec = c(0, rep(0.32, 4), rep(0.3, 4), rep(0.27, 2))
    RAR.vec = c(rep(24, 5), rep(0,6))/120
  }
  

  n = m+n.adj*n.look
  # mean.vec = c(0, rep(0.27, 6), c(rep(0.3, 4)))
  theta.in = 0.5
  beta.SSR = NA
  SSR.ind = FALSE
  
  n.bon = n
  n.dun = n

  n.adj.max = floor((n.bon-m-(t.SSR.check-1)*n.adj)/(n.look-t.SSR.check+1))
  

  #################################################################################################
   dec.bon = dec.dun  = matrix(0, nrow = n.itt, ncol = n.arm)
    dec.RAR  = matrix(0, nrow = n.itt, ncol = n.arm + 1)
    n.RAR.mat = matrix(0, nrow = n.itt, ncol = n.arm + 2)

    dec.bon.four.vec = dec.dun.four.vec = 
    dec.bon.two.vec = dec.dun.two.vec = 
    p.out.global = dec.bon.vec = dec.dun.vec = dec.RAR.two = ASN.RAR.two = dec.RAR.four = ASN.RAR.four = 
    rep(0, n.itt)
    
  
  output$RAR_theta[ind] = theta.in
  
  output$n[ind] = n
  output$m[ind] = m
  output$n.adj.initial[ind] = n.adj
  output$n.adj.min[ind] = n.adj.min
  output$n.adj.max[ind] = n.adj.max
  output$beta.SSR[ind] = beta.SSR
  
  n.look = (n-m)/n.adj
  output$n.look[ind] = n.look
  z.cohort.array = array(NA, dim = c(n.look+2, 2^n.arm, n.itt))
  
  output$t.SSR.check[ind] = t.SSR.check
  
  for (itt in 1:n.itt){
    print(itt)
    
  # cl <- makeCluster(n.cluster)
  # registerDoParallel(cl)
  # pred = foreach(itt = 1:n.itt) %dopar% {    
    
    n.vec.initial = RAR.equal.vec*m
    
    sample.list = lapply(1:(n.arm+1), function(x){rnorm(n.vec.initial[x], mean.vec[x], sd.common)})
    
    z.cohort.mat = matrix(NA, nrow = (n.look+2), ncol = (2^n.arm))
    z.cohort.mat = data.frame(z.cohort.mat)
    
    col.name.ind = 1
    for (i in 1:n.arm){
      combn.temp = combn(2:(n.arm+1), i)
      col.name.length = dim(combn.temp)[2]
      colnames(z.cohort.mat)[col.name.ind:(col.name.ind+col.name.length-1)] = 
        as.vector(apply(combn.temp, 2, function(x){paste0(x, collapse = ",")}))
      col.name.ind = col.name.ind+col.name.length
    }
    colnames(z.cohort.mat)[(2^n.arm)] = "n"
    
    for (i in 1:(2^n.arm-1)){
      combn.arm = colnames(z.cohort.mat)[i]
      combn.arm.ind = as.numeric(strsplit(combn.arm, ",")[[1]])
      pbo.sample.pool = as.vector(unlist(sample.list[1]))
      trt.sample.pool = as.vector(unlist(sample.list[combn.arm.ind]))
      # z.cohort.mat[1, i] = qnorm(1-t.test(trt.sample.pool, pbo.sample.pool, alternative = "greater", 
      #                                     var.equal = FALSE)$p.value)
      z.cohort.mat[1, i] = z.test.known.var(trt.sample.pool, pbo.sample.pool, 1)
    }
    z.cohort.mat[1, 2^n.arm] = length(unlist(sample.list))
    
    n.adj.current = n.adj
    
    for (i in 1:n.look){
    
        # print(n.adj.current)
        #RAR.measure = sapply(2:(n.arm+1), function(x){t.test(sample.list[[x]], sample.list[[1]])$statistic})
      
        RAR.measure = sapply(2:(n.arm+1), function(x){mean(sample.list[[x]])*sqrt(length(sample.list[[x]]))})
        
        RAR.vec.adj = c(RAR.vec[1], (RAR.vec[-1])[match(RAR.measure, sort(RAR.measure, decreasing = TRUE))])
        
        ## if SSR and >= t.SSR.check
        if ((i>=t.SSR.check)&(SSR.ind)){
          
          select.arm.interim = which.max(RAR.measure)+1
          s.interim.vec = sample.list[[select.arm.interim]]
          pbo.interim.vec = sample.list[[1]]
          
          SS.total = SS.func(beta.in = beta.SSR,
                             sd.pool.in = mean(c(sd(s.interim.vec), sd( pbo.interim.vec))),
                             delta.in = mean(s.interim.vec) - mean( pbo.interim.vec))
          
          SS.need = SS.total - length(s.interim.vec) - length(pbo.interim.vec)
          n.adj.ratio = SS.need / ((n.look-i+1)*sum(RAR.vec[1:2]*n.adj.min))
          n.adj.ratio = min(n.adj.max/n.adj.min, n.adj.ratio)
          n.adj.current = max(n.adj.min, floor(n.adj.ratio*n.adj.min))
        }


          n.add.vec = pmax(0, round(RAR.vec.adj*n.adj.current))
        
        # print(sum(n.add.vec))
        
        sample.new.list = sapply(1:(n.arm+1), function(x){rnorm(n.add.vec[x], mean.vec[x], sd.common)})
      
        sample.list = mapply(c, sample.list, sample.new.list, SIMPLIFY=FALSE)
        
        for (j in 1:(2^n.arm-1)){
          combn.arm = colnames(z.cohort.mat)[j]
          combn.arm.ind = as.numeric(strsplit(combn.arm, ",")[[1]])
          pbo.sample.pool = as.vector(unlist(sample.new.list[1]))
          trt.sample.pool = as.vector(unlist(sample.new.list[combn.arm.ind]))
          # z.cohort.mat[1+i, j] = qnorm(1-t.test(trt.sample.pool, pbo.sample.pool, alternative = "greater", 
          #                                     var.equal = TRUE)$p.value)
          z.cohort.mat[1+i, j] = z.test.known.var(trt.sample.pool, pbo.sample.pool, 1)
        }
        z.cohort.mat[1+i, 2^n.arm] = length(unlist(sample.new.list))
        
    }
    
    RAR.measure = sapply(2:(n.arm+1), function(x){mean(sample.list[[x]])*sqrt(length(sample.list[[x]]))})
    # RAR.measure = sapply(2:(n.arm+1), function(x){t.test(sample.list[[x]], sample.list[[1]])$statistic})
    
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
    
    for (j in 1:(2^n.arm-1)){
      z.vec = z.cohort.mat[1:(n.look+1), j]
      # z.n.vec = z.cohort.mat[1:(n.look+1), 2^n.arm]
      z.n.vec = c(m, rep(n.adj, n.look))
      
      z.cohort.mat[n.look+2, j] = comb.func(sum(!is.na(z.vec))-1, z.vec[!is.na(z.vec)], z.n.vec[!is.na(z.vec)])
    }
    
    z.cohort.array[,,itt] = as.matrix(z.cohort.mat)
    
    # p.global = t.test(as.vector(unlist(sample.list[2:(n.arm+1)])), as.vector(unlist(sample.list[1])), 
    #                   alternative = "greater", 
    #                   var.equal = TRUE)$p.value
    p.global = 1-pnorm(z.cohort.mat[n.look+2, 2^n.arm-1])
    
    p.out.global[itt] = p.global
    #######################################################
    ## decision for RAR selecting one arm
    for (i in 2:(n.arm+1)){
      contain.ind = grep(as.character(i), colnames(z.cohort.mat))
      z.contain.vec = z.cohort.mat[n.look+2, contain.ind]
      if ((sum(z.contain.vec<1.96)==0)) dec.RAR[itt, (i-1)] = 1
    }
    
    select.arm = as.numeric(which.max(RAR.measure))
    n.RAR.mat[itt, 1:(n.arm+1)] = sapply(1:(n.arm+1), function(x){length(sample.list[[x]])})
    n.RAR.mat[itt, (n.arm+2)] = n.RAR.mat[itt, (select.arm+1)]
    
    dec.RAR[itt, n.arm+1] = dec.RAR[itt, select.arm]
    
    #######################################################
    ## decision for RAR selecting two arms
    select.two.arm = order(RAR.measure, decreasing = TRUE)[1:2] + 1
    select.two.arm = sort(select.two.arm)
    
    contain.ind.1 = grep(as.character(select.two.arm[1]), colnames(z.cohort.mat))
    contain.ind.2 = grep(as.character(select.two.arm[2]), colnames(z.cohort.mat))
    contain.ind.two = intersect(contain.ind.1, contain.ind.2)
    
    if (sum(dec.RAR[itt, 1:n.arm])>=2) dec.RAR.two[itt] = 1
    # z.contain.two.vec = z.cohort.mat[n.look+2, contain.ind.two]
    # 
    # if ((sum(z.contain.two.vec<1.96)==0)) dec.RAR.two[itt] = 1
    
    
    ASN.RAR.two[itt] = sum(n.RAR.mat[itt, select.two.arm])
    
    #######################################################
    ## decision for RAR selecting four arms
    select.four.arm = order(RAR.measure, decreasing = TRUE)[1:4] + 1
    select.four.arm = sort(select.four.arm)
    
    contain.ind.1 = grep(as.character(select.four.arm[1]), colnames(z.cohort.mat))
    contain.ind.2 = grep(as.character(select.four.arm[2]), colnames(z.cohort.mat))
    contain.ind.3 = grep(as.character(select.four.arm[3]), colnames(z.cohort.mat))
    contain.ind.4 = grep(as.character(select.four.arm[4]), colnames(z.cohort.mat))
    contain.ind.four = intersect(intersect(intersect(contain.ind.1, contain.ind.2), contain.ind.3), contain.ind.4)
    
    # z.contain.four.vec = z.cohort.mat[n.look+2, contain.ind.four]
    # 
    # if ((sum(z.contain.four.vec<1.96)==0)) dec.RAR.four[itt] = 1
    
    if (sum(dec.RAR[itt, 1:n.arm])>=4) dec.RAR.four[itt] = 1
    ASN.RAR.four[itt] = sum(n.RAR.mat[itt, select.four.arm])
    
    
    #####################################################
    ## Bonferroni and Dunttee
    n.arm.equal = round(n.bon/(n.arm+1))
    n.pbo.equal = n.bon - n.arm.equal*n.arm

    sample.pbo.equal = rnorm(n.pbo.equal, mean = mean.vec[1], sd = sd.common)
    sample.trt.equal = sapply(2:(n.arm+1), function(x){rnorm(n.arm.equal, mean.vec[x], sd.common)})
    
    p.unadj = sapply(1:n.arm, function(x){
      pnorm(z.test.known.var(sample.trt.equal[, x], sample.pbo.equal, sd.common),lower.tail = FALSE)})   
    # p.unadj = sapply(1:n.arm, function(x){t.test(sample.trt.equal[, x], sample.pbo.equal, 
    #                                              alternative = "greater", var.equal = TRUE)$p.value})    
    p.adj.bon = p.adjust(p.unadj, "bonferroni")
   
    dec.bon[itt, ] = as.numeric(p.adj.bon<=one.sided.alpha)
    dec.bon.vec[itt] = as.numeric(sum(p.adj.bon<=one.sided.alpha)>0)
    dec.bon.four.vec[itt] = as.numeric(sum(p.adj.bon<=one.sided.alpha)>=4)
    dec.bon.two.vec[itt] = as.numeric(sum(p.adj.bon<=one.sided.alpha)>=2)
    
    # Dunttee
    n.arm.equal = round(n.dun/(n.arm+1))
    n.pbo.equal = n.dun - n.arm.equal*n.arm
    
    sample.pbo.equal = rnorm(n.pbo.equal, mean = mean.vec[1], sd = sd.common)
    sample.trt.equal = sapply(2:(n.arm+1), function(x){rnorm(n.arm.equal, mean.vec[x], sd.common)})
    
    data.dunttee <- data.frame("sample" = c(sample.pbo.equal, as.vector(sample.trt.equal)),
                       "group" = c(rep(1, n.pbo.equal), rep(2:(n.arm+1), each = n.arm.equal)))
    data.dunttee$group = factor(data.dunttee$group)
    
    fit.dunttee <- aov(sample ~ group, data.dunttee)
    
    fit.dunttee.mcp <- glht(fit.dunttee, linfct=mcp(group="Dunnett"),
                            alternative = "greater")
    p.adj.dun = as.numeric(summary(fit.dunttee.mcp, test=(adjusted(type = "free")))$test$pvalues)
    
    dec.dun[itt, ] = as.numeric(p.adj.dun<=one.sided.alpha)
    dec.dun.vec[itt] = as.numeric(sum(p.adj.dun<=one.sided.alpha)>0)
    dec.dun.four.vec[itt] = as.numeric(sum(p.adj.dun<=one.sided.alpha)>=4)
    dec.dun.two.vec[itt] = as.numeric(sum(p.adj.dun<=one.sided.alpha)>=2)

  }
  
  output$RAR.ratio[ind] = paste(RAR.vec*n.adj,collapse = ": ")
  output$true_ES[ind] = paste(round(mean.vec*sqrt(n/(n.arm+1)/2),2),collapse = ", ")
  output$true_mean[ind] = paste(round(mean.vec,3),collapse = ", ")
  output$reject_global[ind] = mean(p.out.global<=one.sided.alpha)
  
  output$reject_bonferroni[ind] = mean(dec.bon.vec)
  output$reject_dunnett[ind] = mean(dec.dun.vec)
  output$reject_bonferroni_four[ind] = mean(dec.bon.four.vec)
  output$reject_dunnett_four[ind] = mean(dec.dun.four.vec)
  output$reject_bonferroni_two[ind] = mean(dec.bon.two.vec)
  output$reject_dunnett_two[ind] = mean(dec.dun.two.vec)
  
  output[ind, c(as.vector(paste0("reject_", 1:n.arm)), "reject_s")] = apply(dec.RAR, 2, mean)
  output[ind, c(as.vector(paste0("reject_bon_", 1:n.arm)))] = apply(dec.bon, 2, mean)
  output[ind, c(as.vector(paste0("reject_dun_", 1:n.arm)))] = apply(dec.dun, 2, mean)
  
  output[ind, c("ASN_pbo", as.vector(paste0("ASN_", 1:n.arm)), "ASN_s")] = 
    c(mean(n.RAR.mat[,1]), sapply(1:(n.arm+1), function(x){mean(n.RAR.mat[dec.RAR[,x]==1,x+1])}))
  
  output$reject_s_two[ind] = mean(dec.RAR.two)
  output$ASN_s_two[ind] = mean(ASN.RAR.two[dec.RAR.two==1])
  
  output$reject_s_four[ind] = mean(dec.RAR.four)
  output$ASN_s_four[ind] = mean(ASN.RAR.four[dec.RAR.four==1])
  
  output$n_bon[ind] = n.bon
  output$n_dun[ind] = n.dun
  output$ASN_total[ind] = sum(apply(n.RAR.mat, 2, mean)[1:(n.arm+1)])
  
}


write.csv(output, paste0("results_", n.arm, "_doses_SSR.csv"))
# print(output)





