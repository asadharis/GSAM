# Simulation results.

process.results <- function(name = "Liver_GSE14520_U133A") {
  #name  = "Lung_GSE19804"
  files <- list.files(name)
  files <- paste0(name,"/",files)
  nsim <- length(files)

  res <- data.frame()
  for(i in 1:nsim) {
    load(files[i])
    res <- rbind(res, results)
  }

  require(ggplot2)

  plot1 <- ggplot(res, aes(x=method_name, y=val_auc,
                           fill = method_name)) +
    theme(text = element_text(size=15)) +
    ylab("Validation AUC") + xlab("") + guides(fill = FALSE) +
    theme(axis.text=element_text(size=14, face = "bold")) +
    geom_boxplot()+scale_y_log10()

  plot2 <- ggplot(res, aes(x=method_name, y=val_mse,
                           fill = method_name)) +
    theme(text = element_text(size=15)) +
    ylab("Test Error") + xlab("") + guides(fill = FALSE) +
    theme(axis.text=element_text(size=14, face = "bold")) +
    geom_boxplot()+scale_y_log10()

  require(dplyr)
  res <- as_tibble(res)
  mu <- res %>%
    group_by(method_name) %>%
    summarise(mean_auc = mean(val_auc), se_auc = sd(val_auc)/sqrt(n()),
              mean_err = mean(val_mse), se_err = sd(val_mse)/sqrt(n()))

  sparse <- res %>%
    group_by(method_name) %>%
    summarise(sp_auc = mean(sparse_auc),
              sp_err = mean(sparse_mse) )
  return(list(plot1, plot2, mu, sparse))
}


plot.results <- function(obj,
                         exclude = "") {
  require(ggplot2)
  obj1 <- data.frame(obj[[3]])
  ind <- !(obj1$method_name %in% exclude)
  obj1 <- obj1[ind,]
  plot1 <- ggplot(obj1,
                  aes(x = method_name,
                      y = mean_auc, color = method_name))+
    geom_point(size = 3) + theme(legend.position = "none") +
    geom_errorbar(aes(ymin=mean_auc-se_auc,
                      ymax=mean_auc+se_auc), width=.2)

  plot1



}

##########################
#### Finalize Figures ####
##########################
# Lung disease
lung <- process.results("Lung_GSE19804")
prostate <- process.results("Prostate_GSE6919_U95B")
breast <- process.results("Breast_GSE70947")
throat <- process.results("Throat_GSE42743")

# Aggregate the data.
dat <- rbind(lung[[3]], prostate[[3]],
             breast[[3]], throat[[3]])
j <- nrow(lung[[3]])
dat$study <- rep(c("Lung", "Prostate", "Breast", "Throat"),
                 each = j)
# Exclude TF-2, spam4 and spam5
exclude <- c("tf2", "spam4", "spam5")
ind <- !(dat$method_name %in% exclude)
dat <- dat[ind,]

# Change order of SpAM appearance.
dat$method_name <- factor(dat$method_name,
                          levels = c("lasso", "spam2", "spam3",
                                     "spam10", "ssp", "tf0",
                                     "tf1"),
                          labels = c("Lasso", "SpAM, M=2", "SpAM, M=3",
                                     "SpAM, M=10", "SSP", "TF, k=0",
                                     "TF, k=1"))

##### Generate Table #####

library(dplyr)
library(tidyr)
mytab <- dat %>%
  select(mean_auc,study,method_name) %>%
  pivot_wider(names_from = study, values_from = mean_auc)

mytab[2:4, ] <- mytab[c(3,4,2),]
mytab

mytabSE <- dat %>%
  select(se_auc, study, method_name) %>%
  pivot_wider(names_from = study, values_from = se_auc)

mytabSE[2:4, ] <- mytabSE[c(3,4,2),]
mytabSE1000 <- mytabSE
mytabSE1000[,-1] <- mytabSE[,-1] * 1000
mytabSE1000

round(as.matrix(mytabSE1000[,-1]),2)

fin_res <- paste0(round(as.matrix(mytab[,-1]),3)," (",
       round(as.matrix(mytabSE1000[,-1]),2), ")")
fin_res <- matrix(fin_res, ncol = 4)

colnames(fin_res) <- colnames(mytab)[-1]
fin_tab <- data.frame("Method" = mytab$method_name, fin_res)

xtable(fin_tab)

# Generate plot
ggplot(data = dat,
       mapping = aes(x = method_name, y = mean_auc,
                     color = method_name)) +
  geom_point(size = 3) + theme(legend.position = "none") +
  geom_errorbar(aes(ymin=mean_auc-se_auc,
                    ymax=mean_auc+se_auc), width=.2)+
  facet_wrap(vars(study), scales = "free_y") +

  theme(text = element_text(size = 16)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text())+
  ylab("Mean AUC")


# plot.results(lung, exclude = c("tf2", "spam4", "spam5"))
# # View Sparsity results
# lung[[4]]
# cbind(lung[[4]], (1 - lung[[4]][,2])*100)
#
# plot.results(prostate, exclude = c("tf2", "spam10"))
# prostate[[4]]
# cbind(prostate[[4]][,1], (1 - prostate[[4]][,2])*100)
#
# plot.results(throat, exclude = "tf2")
# cbind(prostate[[4]][,1], (1 - throat[[4]][,2])*100)
#
# plot.results(breast, exclude = "tf2")
#
# cbind(prostate[[4]][,1], (1 - breast[[4]][,2])*100)
