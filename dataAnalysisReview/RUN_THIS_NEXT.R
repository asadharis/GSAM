# Simulation results.

process.results <- function(name = "Liver_GSE14520_U133A") {
  name  = "Colorectal_GSE44076"
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

  return(list(plot1, plot2, mu))
}

### Gives me good results
prostate <- process.results(name = "Prostate_GSE6919_U95B")
liver <- process.results(name = "Liver_GSE14520_U133A")

#breast2 <- process.results("Breast_GSE22820")

colorectal <- process.results("Colorectal_GSE44076")
liver2 <- process.results("Liver_GSE76427")
lung <- process.results("Lung_GSE19804")
