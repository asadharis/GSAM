process.dat <- function(nvars = 100) {
  filename <- paste0("riboflavin_nvars", nvars)
  files <- paste0(filename, "/", list.files(filename))
  nfiles <- length(files)
  err <- matrix(NA, ncol = 9, nrow = nfiles)
  colnames(err) <- c(paste0("spam",3:7), "ssp", paste0("tf",0:2))
  sparse <- err
  for(i in 1:nfiles) {
    load(files[i])
    err[i, ] <- c(spam3$err["min"], spam4$err["min"], spam5$err["min"],
                  spam6$err["min"], spam7$err["min"], ssp$err["min"],
                  tf0$err["min"], tf1$err["min"], tf2$err["min"])
  }
}