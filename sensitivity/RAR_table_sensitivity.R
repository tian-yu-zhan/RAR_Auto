
setwd("C:/TZ/collaborations/clinical studies/immunology/RZB HS/Cui2007/RAR_manuscript_Lu/R code/v2/sensitivity/")

library(xtable)

##########################################################################################
## 1 out of 3 CHW
# three.table = read.csv("results_small_one_out_of_three_doses.csv",stringsAsFactors = FALSE)
# three.table = read.csv("results_unknown_var_one_out_of_three_doses.csv",stringsAsFactors = FALSE)
three.table = read.csv("results_bin_small_one_out_of_three_doses.csv",stringsAsFactors = FALSE)
three.table = data.frame(three.table)

three.table.one = matrix(NA, nrow = 15, ncol = 12)
colnames(three.table.one) = c("Scenario", "Method", "reject_l", "reject_m",
                                "reject_h", "reject_s",
                                "reject_global", "ASN_pbo", "ASN_l", "ASN_m", "ASN_h",
                                "ASN_s")
three.table.one = data.frame(three.table.one)

for (i in 1:5){
  data.temp = three.table[i, ]

  three.table.one$Scenario[(1:3)+(i-1)*3] = c(paste0("S", i), "", "")
  three.table.one$Method[(1:3)+(i-1)*3] = c("SSACTP", "Bonferroni", "Dunnett")
  # two.table.final$n[(1:3)+(i-1)*3] = data.temp$n
  # two.table.final$m[(1:3)+(i-1)*3] = c(data.temp$m, "-", "-")
  # two.table.final$n.adj[(1:3)+(i-1)*3] = c(data.temp$n.adj, "-", "-")
  # three.table.final$true_mean[(1:3)+(i-1)*3] = data.temp$true_mean
  # three.table.final$RAR.ratio[(1:3)+(i-1)*3] = c(data.temp$RAR.ratio, "1: 1: 1: 1", "1: 1: 1: 1")

  three.table.one$reject_l[(1:3)+(i-1)*3] =
    paste0(sprintf('%.2f',c(data.temp$reject_l_CHW_close,
                            data.temp$reject_l_bonferroni,
                            data.temp$reject_l_dunnett)*100), "%")

  three.table.one$reject_m[(1:3)+(i-1)*3] =
    paste0(sprintf('%.2f',c(data.temp$reject_m_CHW_close,
                            data.temp$reject_m_bonferroni,
                            data.temp$reject_m_dunnett)*100), "%")

  three.table.one$reject_h[(1:3)+(i-1)*3] =
    paste0(sprintf('%.2f',c(data.temp$reject_h_CHW_close,
                            data.temp$reject_h_bonferroni,
                            data.temp$reject_h_dunnett)*100), "%")

  three.table.one$reject_s[(1:3)+(i-1)*3] =
    paste0(sprintf('%.2f',c(data.temp$reject_s_CHW_close,
                            data.temp$reject_bonferroni,
                            data.temp$reject_dunnett)*100), "%")

  three.table.one$reject_global[(1:3)+(i-1)*3] =
    paste0(sprintf('%.2f',c(data.temp$reject_global,
                            data.temp$reject_bonferroni,
                            data.temp$reject_dunnett)*100), "%")

  three.table.one[(1)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
    c(sprintf('%.0f',data.temp$n.pbo.out), sprintf('%.0f',data.temp$n.l.out),
      sprintf('%.0f',data.temp$n.m.out), sprintf('%.0f',data.temp$n.h.out),
      sprintf('%.0f',data.temp$n.s.out))

  three.table.one[(2)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
    c(sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
      sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
      sprintf('%.0f',data.temp$n/4))

  three.table.one[(3)+(i-1)*3, c("ASN_pbo", "ASN_l", "ASN_m", "ASN_h","ASN_s")] =
    c(sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
      sprintf('%.0f',data.temp$n/4), sprintf('%.0f',data.temp$n/4),
      sprintf('%.0f',data.temp$n/4))



}

print(xtable(three.table.one), include.rownames = FALSE, align="c")




