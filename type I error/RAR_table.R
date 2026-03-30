
setwd("C:/Users/zhantx/Documents/collaborations/clinical studies/immunology/RZB HS/Cui2007/RAR_manuscript_Lu/results/for_manuscript/")

library(xtable)

##########################################################################################
## two treatment arms
# two.table = read.csv("results_two.csv",stringsAsFactors = FALSE)
# two.table = data.frame(two.table)
# 
# two.table.final = matrix(NA, nrow = 12, ncol = 10)
# colnames(two.table.final) = c("Scenario", "Method", "reject_l", 
#                               "reject_h", "reject_s",
#                               "reject_global", "ASN_pbo", "ASN_l", "ASN_h", "ASN_s")
# two.table.final = data.frame(two.table.final)
# 
# for (i in 1:4){
#   data.temp = two.table[i, ]
#   two.table.final$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   two.table.final$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "step-down Dunnett")
#   # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
#   # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
#   # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
#   # two.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
#   # two.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, "1: 1: 1", "1: 1: 1")
#   
#   two.table.final$reject_l[(1:3)+(i-1)*3] = paste0(sprintf('%.2f',c(data.temp$reject_l, data.temp$reject_l_bonferroni, 
#                                               data.temp$reject_l_dunnett)*100), "%")
#   two.table.final$reject_h[(1:3)+(i-1)*3] = paste0(sprintf('%.2f',c(data.temp$reject_h, data.temp$reject_h_bonferroni, 
#                                               data.temp$reject_h_dunnett)*100), "%")
#   two.table.final$reject_s[(1:3)+(i-1)*3] = paste0(sprintf('%.2f',c(data.temp$reject_s, data.temp$reject_bonferroni, 
#                                               data.temp$reject_dunnett)*100), "%")
#   two.table.final$reject_global[(1:3)+(i-1)*3] = paste0(sprintf('%.2f',c(data.temp$reject_global, 
#                                                                          data.temp$reject_bonferroni, 
#                                               data.temp$reject_dunnett)*100), "%")
#   
#   two.table.final[(1)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_h", "ASN_s")] = 
#                     c(sprintf('%.2f',data.temp$n.pbo.out), sprintf('%.2f',data.temp$n.l.out),
#                       sprintf('%.2f',data.temp$n.h.out), sprintf('%.2f',data.temp$n.s.out))
#   
#   two.table.final[(2)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_h", "ASN_s")] = 
#     c(sprintf('%.2f',data.temp$ASN_bon_pbo), sprintf('%.2f',data.temp$ASN_bon_l),
#       sprintf('%.2f',data.temp$ASN_bon_h), sprintf('%.2f',data.temp$ASN_bon_s))
#   
#   two.table.final[(3)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_h", "ASN_s")] = 
#     c(sprintf('%.2f',data.temp$ASN_dun_pbo), sprintf('%.2f',data.temp$ASN_dun_l),
#       sprintf('%.2f',data.temp$ASN_dun_h), sprintf('%.2f',data.temp$ASN_dun_s))
#   
# }
# 
# print(xtable(two.table.final), include.rownames = FALSE, align="c")

##########################################################################################
## 1 out of 3 CHW
# three.table = read.csv("results_one_out_of_three_doses.csv",stringsAsFactors = FALSE)
# three.table = data.frame(three.table)
# 
# three.table.one = matrix(NA, nrow = 15, ncol = 12)
# colnames(three.table.one) = c("Scenario", "Method", "reject_l", "reject_m",
#                                 "reject_h", "reject_s",
#                                 "reject_global", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h",
#                                 "ASN_s")
# three.table.one = data.frame(three.table.one)
# 
# for (i in 1:5){
#   data.temp = three.table[i, ]
# 
#   three.table.one$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   three.table.one$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")
#   # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
#   # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
#   # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
#   # three.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
#   # three.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, "1: 1: 1: 1", "1: 1: 1: 1")
# 
#   three.table.one$reject_l[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_l_CHW_close,
#                             data.temp$reject_l_bonferroni,
#                             data.temp$reject_l_dunnett)*100), "%")
# 
#   three.table.one$reject_m[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_m_CHW_close,
#                             data.temp$reject_m_bonferroni,
#                             data.temp$reject_m_dunnett)*100), "%")
# 
#   three.table.one$reject_h[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_h_CHW_close,
#                             data.temp$reject_h_bonferroni,
#                             data.temp$reject_h_dunnett)*100), "%")
# 
#   three.table.one$reject_s[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_s_CHW_close,
#                             data.temp$reject_bonferroni,
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.one$reject_global[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_global,
#                             data.temp$reject_bonferroni,
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.one[(1)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.0f',data.temp$n.pbo.out), sprintf('%.0f',data.temp$n.l.out),
#       sprintf('%.0f',data.temp$n.m.out), sprintf('%.0f',data.temp$n.h.out),
#       sprintf('%.0f',data.temp$n.s.out))
# 
#   three.table.one[(2)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
#       sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
#       sprintf('%.0f',data.temp$n/4))
# 
#   three.table.one[(3)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
#       sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
#       sprintf('%.0f',data.temp$n/4))
# 
# 
# 
# }
# 
# print(xtable(three.table.one), include.rownames = FALSE, align="c")

# ##########################################################################################
## 1 out of 3 unweighted
# three.table = read.csv("results_one_out_of_three_doses.csv",stringsAsFactors = FALSE)
# three.table = data.frame(three.table)
# 
# three.table.one.unweighted = matrix(NA, nrow = 15, ncol = 12)
# colnames(three.table.one.unweighted) = c("Scenario", "Method", "reject_l", "reject_m",
#                               "reject_h", "reject_s",
#                               "reject_global", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h",
#                               "ASN_s")
# three.table.one.unweighted = data.frame(three.table.one.unweighted)
# 
# for (i in 1:5){
#   data.temp = three.table[i, ]
# 
#   three.table.one.unweighted$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   three.table.one.unweighted$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")
#   # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
#   # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
#   # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
#   # three.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
#   # three.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, "1: 1: 1: 1", "1: 1: 1: 1")
# 
#   three.table.one.unweighted$reject_l[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.2f',c(data.temp$reject_l_unweighted_close,
#                             data.temp$reject_l_bonferroni,
#                             data.temp$reject_l_dunnett)*100), "%")
# 
#   three.table.one.unweighted$reject_m[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.2f',c(data.temp$reject_m_unweighted_close,
#                             data.temp$reject_m_bonferroni,
#                             data.temp$reject_m_dunnett)*100), "%")
# 
#   three.table.one.unweighted$reject_h[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.2f',c(data.temp$reject_h_unweighted_close,
#                             data.temp$reject_h_bonferroni,
#                             data.temp$reject_h_dunnett)*100), "%")
# 
#   three.table.one.unweighted$reject_s[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.2f',c(data.temp$reject_s_unweighted_close,
#                             data.temp$reject_bonferroni,
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.one.unweighted$reject_global[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.2f',c(data.temp$reject_global,
#                             data.temp$reject_bonferroni,
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.one.unweighted[(1)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.2f',data.temp$n.pbo.out), sprintf('%.2f',data.temp$n.l.out),
#       sprintf('%.2f',data.temp$n.m.out), sprintf('%.2f',data.temp$n.h.out),
#       sprintf('%.2f',data.temp$n.s.out))
# 
#   three.table.one.unweighted[(2)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.2f',data.temp$n/4), sprintf('%.2f',data.temp$n/4),
#       sprintf('%.2f',data.temp$n/4), sprintf('%.2f',data.temp$n/4),
#       sprintf('%.2f',data.temp$n/4))
# 
#   three.table.one.unweighted[(3)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.2f',data.temp$n/4), sprintf('%.2f',data.temp$n/4),
#       sprintf('%.2f',data.temp$n/4), sprintf('%.2f',data.temp$n/4),
#       sprintf('%.2f',data.temp$n/4))
# 
# 
# 
# }
# 
# print(xtable(three.table.one.unweighted), include.rownames = FALSE, align="c")

# ##########################################################################################
# ## 2 out of 3
# three.table = read.csv("results_two_out_of_three_doses.csv",stringsAsFactors = FALSE)
# three.table = data.frame(three.table)
# 
# three.table.two = matrix(NA, nrow = 15, ncol = 12)
# colnames(three.table.two) = c("Scenario", "Method", "reject_lm", "reject_lh",
#                               "reject_mh", "reject_s_two",
#                               "reject_global", "ASN_pbo", "ASN_lm", "ASN_lh", "ASN_mh",
#                               "ASN_s_two")
# three.table.two = data.frame(three.table.two)
# 
# for (i in 1:5){
#   data.temp = three.table[i, ]
#   
#   three.table.two$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   three.table.two$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")
#   # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
#   # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
#   # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
#   # three.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
#   # three.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, "1: 1: 1: 1", "1: 1: 1: 1")
#   
#   three.table.two$reject_lm[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_lm_CHW_close,
#                           data.temp$reject_lm_bonferroni, 
#                         data.temp$reject_lm_dun)*100), "%")
#   
#   three.table.two$reject_lh[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_lh_CHW_close,
#                             data.temp$reject_lh_bonferroni, 
#                             data.temp$reject_lh_dun)*100), "%")
#   
#   three.table.two$reject_mh[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_mh_CHW_close,
#                             data.temp$reject_mh_bonferroni, 
#                             data.temp$reject_mh_dun)*100), "%")
#  
#   three.table.two$reject_s_two[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_s_two_CHW_close,
#                             data.temp$reject_s_two_bonferroni, 
#                             data.temp$reject_s_two_dun)*100), "%")
#   
#   
#   three.table.two$reject_global[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_global,
#                             data.temp$reject_bonferroni, 
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.two[(1)+(i-1)*3, c("ASN_pbo", "ASN_lm", "ASN_lh", "ASN_mh", "ASN_s_two")] = 
#     c(sprintf('%.2f',data.temp$n.pbo.out),
#       sprintf('%.2f',data.temp$n.lm.out), sprintf('%.2f',data.temp$n.lh.out),
#       sprintf('%.2f',data.temp$n.mh.out), sprintf('%.2f',data.temp$n.s.two.out))
#   
#   
#   
#   three.table.two[(2)+(i-1)*3, c("ASN_pbo", "ASN_lm", "ASN_lh", "ASN_mh", "ASN_s_two")] = 
#     c(sprintf('%.2f',data.temp$n/4), sprintf('%.2f',data.temp$n*2/4),
#       sprintf('%.2f',data.temp$n*2/4), sprintf('%.2f',data.temp$n*2/4),
#       sprintf('%.2f',data.temp$n*2/4))
#   
#   three.table.two[(3)+(i-1)*3, c("ASN_pbo", "ASN_lm", "ASN_lh", "ASN_mh", "ASN_s_two")] = 
#     c(sprintf('%.2f',data.temp$n/4), sprintf('%.2f',data.temp$n*2/4),
#       sprintf('%.2f',data.temp$n*2/4), sprintf('%.2f',data.temp$n*2/4),
#       sprintf('%.2f',data.temp$n*2/4))
#   
#  
#   
# }
# 
# print(xtable(three.table.two), include.rownames = FALSE, align="c")


###################################################################################################
## SSR table for selecting 1 out of 3
# three.SSR.table = read.csv("results_three_doses_SSR.csv",stringsAsFactors = FALSE)
# three.SSR.table = data.frame(three.SSR.table)
# 
# three.table.one.SSR = matrix(NA, nrow = 15, ncol = 13)
# colnames(three.table.one.SSR) = c("Scenario", "Method", "reject_l", "reject_m",
#                               "reject_h", "reject_s",
#                               "reject_global", "ASN_total", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h",
#                               "ASN_s")
# three.table.one.SSR = data.frame(three.table.one.SSR)
# 
# for (i in 1:5){
#   data.temp = three.SSR.table[i, ]
# 
#   three.table.one.SSR$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   three.table.one.SSR$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")
#   # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
#   # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
#   # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
#   # three.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
#   # three.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, "1: 1: 1: 1", "1: 1: 1: 1")
# 
#   three.table.one.SSR$reject_l[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_l_CHW_close,
#                             data.temp$reject_l_bonferroni,
#                             data.temp$reject_l_dunnett)*100), "%")
# 
#   three.table.one.SSR$reject_m[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_m_CHW_close,
#                             data.temp$reject_m_bonferroni,
#                             data.temp$reject_m_dunnett)*100), "%")
# 
#   three.table.one.SSR$reject_h[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_h_CHW_close,
#                             data.temp$reject_h_bonferroni,
#                             data.temp$reject_h_dunnett)*100), "%")
# 
#   three.table.one.SSR$reject_s[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_s_CHW_close,
#                             data.temp$reject_bonferroni,
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.one.SSR$reject_global[(1:3)+(i-1)*3] =
#     paste0(sprintf('%.1f',c(data.temp$reject_global,
#                             data.temp$reject_bonferroni,
#                             data.temp$reject_dunnett)*100), "%")
# 
#   three.table.one.SSR[(1)+(i-1)*3, c("ASN_total", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.0f',data.temp$ASN.total),
#       sprintf('%.0f',data.temp$n.pbo.out), sprintf('%.0f',data.temp$n.l.out),
#       sprintf('%.0f',data.temp$n.m.out), sprintf('%.0f',data.temp$n.h.out),
#       sprintf('%.0f',data.temp$n.s.out))
# 
#   three.table.one.SSR[(2)+(i-1)*3, c("ASN_total", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.0f',data.temp$ASN.bon),
#       sprintf('%.0f',data.temp$ASN.bon/4), sprintf('%.0f',data.temp$ASN.bon/4),
#       sprintf('%.0f',data.temp$ASN.bon/4), sprintf('%.0f',data.temp$ASN.bon/4),
#       sprintf('%.0f',data.temp$ASN.bon/4))
# 
#   three.table.one.SSR[(3)+(i-1)*3, c("ASN_total", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
#     c(sprintf('%.0f',data.temp$ASN.dun),
#       sprintf('%.0f',data.temp$ASN.dun/4), sprintf('%.0f',data.temp$ASN.dun/4),
#       sprintf('%.0f',data.temp$ASN.dun/4), sprintf('%.0f',data.temp$ASN.dun/4),
#       sprintf('%.0f',data.temp$ASN.dun/4))
# 
# 
# 
# }
# 
# print(xtable(three.table.one.SSR), include.rownames = FALSE, align="c")

###################################################################################################
## table for 4 out of 10
ten.SSR.table = read.csv("results_10_doses_SSR.csv",stringsAsFactors = FALSE)
ten.SSR.table = data.frame(ten.SSR.table)

ten.table.four.SSR = matrix(NA, nrow = 3*4, ncol = 8)
colnames(ten.table.four.SSR) = c("k", "Scenario", "Method", "reject_s_k",
                                  "reject_global", "ASN_total", "ASN_pbo",
                                 "ASN_s_k")
ten.table.four.SSR = data.frame(ten.table.four.SSR)

ten.table.four.SSR$k = c("2", rep("", 5), "4", rep("", 5))

for (i in 1:4){
  data.temp = ten.SSR.table[i, ]

  ten.table.four.SSR$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")

  if (i<=2){
    ten.table.four.SSR$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
    ten.table.four.SSR$reject_s_k[(1:3)+(i-1)*3] =
      paste0(sprintf('%.1f',c(data.temp$reject_s_two,
                              data.temp$reject_bonferroni_two,
                              data.temp$reject_dunnett_two)*100), "%")
    ten.table.four.SSR$ASN_s_k[(1:3)+(i-1)*3] =
      c(sprintf('%.0f',data.temp$ASN_s_two), sprintf('%.0f',round(data.temp$n/11*2)),
        sprintf('%.0f',round(data.temp$n/11*2)))
  } else {
    ten.table.four.SSR$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i-2), "", "")
    ten.table.four.SSR$reject_s_k[(1:3)+(i-1)*3] =
      paste0(sprintf('%.1f',c(data.temp$reject_s_four,
                              data.temp$reject_bonferroni_four,
                              data.temp$reject_dunnett_four)*100), "%")
    ten.table.four.SSR$ASN_s_k[(1:3)+(i-1)*3] =
      c(sprintf('%.0f',data.temp$ASN_s_four), sprintf('%.0f',round(data.temp$n/11*4)),
        sprintf('%.0f',round(data.temp$n/11*4)))
  }

  ten.table.four.SSR$reject_global[(1:3)+(i-1)*3] =
    paste0(sprintf('%.1f',c(data.temp$reject_global,
                            data.temp$reject_bonferroni,
                            data.temp$reject_dunnett)*100), "%")

  ten.table.four.SSR$ASN_total[(1:3)+(i-1)*3] =
    c(sprintf('%.0f',data.temp$ASN_total),sprintf('%.0f',data.temp$n),sprintf('%.0f',data.temp$n))

  ten.table.four.SSR$ASN_pbo[(1:3)+(i-1)*3] =
    c(sprintf('%.0f',data.temp$ASN_pbo), sprintf('%.0f',round(data.temp$n/11)),
      sprintf('%.0f',round(data.temp$n/11)))


}

print(xtable(ten.table.four.SSR), include.rownames = FALSE, align="c")

####################################################################################
## example SSR + 2 out of 4

# example.table = read.csv("results_4_doses_SSR.csv",stringsAsFactors = FALSE)
# example.table = data.frame(example.table)
# 
# example.table.out = matrix(NA, nrow = 3, ncol = 6)
# colnames(example.table.out) = c("Method", 
#                                 "reject_s_two",
#                                   "reject_global", "ASN_total", "ASN_pbo", 
#                                   "ASN_s_two")
# example.table.out = data.frame(example.table.out)
# 
# for (i in 1){
#   data.temp = example.table[i, ]
# 
#   # example.table.out$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   example.table.out$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")
# 
#   example.table.out$reject_s_two[(1:3)+(i-1)*3] =
#     c(paste0(sprintf('%.1f',c(data.temp$reject_s_two, data.temp$reject_bon_two, data.temp$reject_dun_two)*100), "%"))
# 
#   example.table.out$reject_global[(1:3)+(i-1)*3] =
#     c(paste0(sprintf('%.1f',c(data.temp$reject_global, data.temp$reject_bon,
#                               data.temp$reject_dun)*100), "%"))
# 
# 
#   example.table.out[(1)+(i-1)*3, c("ASN_total", "ASN_pbo", "ASN_s_two")] =
#     c(sprintf('%.0f',data.temp$ASN_total),
#       sprintf('%.0f',data.temp$ASN_pbo),
#       sprintf('%.0f',data.temp$ASN_s_two))
# 
#   example.table.out[(2)+(i-1)*3, c("ASN_total", "ASN_pbo", "ASN_s_two")] =
#     c(sprintf('%.0f',data.temp$n_bon),
#       sprintf('%.0f',data.temp$n_bon/5),
#       sprintf('%.0f',data.temp$n_bon*2/5))
# 
#   example.table.out[(3)+(i-1)*3, c("ASN_total", "ASN_pbo", "ASN_s_two")] =
#     c(sprintf('%.0f',data.temp$n_dun),
#       sprintf('%.0f',data.temp$n_dun/5),
#      sprintf('%.0f',data.temp$n_dun*2/5))
# 
# 
# }
# 
# print(xtable(example.table.out), include.rownames = FALSE, align="c")







##########################################################################################
## five treatment arms
# five.table = read.csv("results_5_doses.csv",stringsAsFactors = FALSE)
# five.table = data.frame(five.table)
# 
# five.table.final = matrix(NA, nrow = 21, ncol = 18)
# colnames(five.table.final) = c("Scenario", "Method", "reject_D1", "reject_D2",
#                                 "reject_D3", "reject_D4", "reject_D5", "reject_each", "reject_s",
#                                "reject_global", "ASN_pbo", "ASN_each", 
#                                 "ASN_D1", "ASN_D2", "ASN_D3", "ASN_D4", "ASN_D5", "ASN_s")
# five.table.final = data.frame(five.table.final)
# 
# for (i in 1:7){
#   data.temp = five.table[i, ]
#   
#   five.table.final$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
#   five.table.final$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "step-down Dunnett")
#   # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
#   # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
#   # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
#   # five.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
#   # five.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, 
#   #                                               rep("1:1:1:1:1:1", 2))
#   
#   five.table.final$reject_D1[(1:3)+(i-1)*3] = 
#     c(paste0(sprintf('%.2f',c(data.temp$reject_1)*100), "%"), "-", "-")
#   
#   five.table.final$reject_D2[(1:3)+(i-1)*3] = 
#     c(paste0(sprintf('%.2f',c(data.temp$reject_2)*100), "%"), "-", "-")
#   
#   five.table.final$reject_D3[(1:3)+(i-1)*3] = 
#     c(paste0(sprintf('%.2f',c(data.temp$reject_3)*100), "%"), "-", "-")
#   
#   five.table.final$reject_D4[(1:3)+(i-1)*3] = 
#     c(paste0(sprintf('%.2f',c(data.temp$reject_4)*100), "%"), "-", "-")
#   
#   five.table.final$reject_D5[(1:3)+(i-1)*3] = 
#     c(paste0(sprintf('%.2f',c(data.temp$reject_5)*100), "%"), "-", "-")
#   
#   five.table.final$reject_each[(1:3)+(i-1)*3] = 
#     sapply(1:3, function(x){paste0("(", 
#                                   c(paste0(sprintf('%.2f',c(data.temp$reject_1)*100), "%"), "-", "-")[x], ", ",
#                                   c(paste0(sprintf('%.2f',c(data.temp$reject_2)*100), "%"), "-", "-")[x], ", ",
#                                   c(paste0(sprintf('%.2f',c(data.temp$reject_3)*100), "%"), "-", "-")[x], ", ",
#                                   c(paste0(sprintf('%.2f',c(data.temp$reject_4)*100), "%"), "-", "-")[x], ", ",
#                                   c(paste0(sprintf('%.2f',c(data.temp$reject_5)*100), "%"), "-", "-")[x], 
#                                   ")")})
#   
#   five.table.final$reject_s[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_s,
#                             data.temp$reject_bonferroni, 
#                             data.temp$reject_dunnett)*100), "%")
#   
#   
#   five.table.final$reject_global[(1:3)+(i-1)*3] = 
#     paste0(sprintf('%.2f',c(data.temp$reject_global,
#                             data.temp$reject_bonferroni, 
#                             data.temp$reject_dunnett)*100), "%")
#   
#   five.table.final[(1)+(i-1)*3, c("ASN_pbo", 
#                                   "ASN_D1", "ASN_D2", "ASN_D3", "ASN_D4", "ASN_D5", "ASN_s")] = 
#     c(sprintf('%.2f',data.temp$ASN_pbo), sprintf('%.2f',data.temp$ASN_1),
#       sprintf('%.2f',data.temp$ASN_2), sprintf('%.2f',data.temp$ASN_3),
#       sprintf('%.2f',data.temp$ASN_4), sprintf('%.2f',data.temp$ASN_5),
#       sprintf('%.2f',data.temp$ASN_s))
#   
#   five.table.final[(2)+(i-1)*3, c("ASN_pbo", 
#                                   "ASN_D1", "ASN_D2", "ASN_D3", "ASN_D4", "ASN_D5", "ASN_s")] = 
#     c(sprintf('%.2f',data.temp$n/6), sprintf('%.2f',data.temp$n/6),
#       sprintf('%.2f',data.temp$n/6), sprintf('%.2f',data.temp$n/6),
#       sprintf('%.2f',data.temp$n/6), sprintf('%.2f',data.temp$n/6),
#       sprintf('%.2f',data.temp$n/6))
#   
#   five.table.final[(3)+(i-1)*3, c("ASN_pbo", 
#                                   "ASN_D1", "ASN_D2", "ASN_D3", "ASN_D4", "ASN_D5", "ASN_s")] = 
#     c(sprintf('%.2f',data.temp$n/6), sprintf('%.2f',data.temp$n/6),
#       sprintf('%.2f',data.temp$n/6), sprintf('%.2f',data.temp$n/6),
#       sprintf('%.2f',data.temp$n/6), sprintf('%.2f',data.temp$n/6),
#       sprintf('%.2f',data.temp$n/6))
#   
#   five.table.final$reject_each[(1:3)+(i-1)*3] = 
#     sapply(1:3, function(x){paste0("(", 
#                                    c(paste0(sprintf('%.2f',c(data.temp$reject_1)*100), "%"), "-", "-")[x], ", ",
#                                    c(paste0(sprintf('%.2f',c(data.temp$reject_2)*100), "%"), "-", "-")[x], ", ",
#                                    c(paste0(sprintf('%.2f',c(data.temp$reject_3)*100), "%"), "-", "-")[x], ", ",
#                                    c(paste0(sprintf('%.2f',c(data.temp$reject_4)*100), "%"), "-", "-")[x], ", ",
#                                    c(paste0(sprintf('%.2f',c(data.temp$reject_5)*100), "%"), "-", "-")[x], 
#                                    ")")})
#   
#   
#   
# }
# 
# print(xtable(five.table.final[, -(3:7)]), 
#       include.rownames = FALSE, align="c")






