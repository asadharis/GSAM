source("Process_Results.R")
source("Process_Results2.R")

#SPAM 3,6,10,20,30,50

mse.plot <- plot.results()

spam3 <- plot.object(obj.name = "mod.spam3", scen.num = 4, p = 100, n = 500,
            title = "SPAM, M = 3")

spam10 <- plot.object(obj.name = "mod.spam10", scen.num = 4, p = 100, n = 500,
            title = "SPAM, M = 10")
spam50 <- plot.object(obj.name = "mod.spam50", scen.num = 4, p = 100, n = 500,
            title = "SPAM, M = 50")

ssp <- plot.object(obj.name = "mod.ssp", scen.num = 4, p = 100, n = 500,
            title = "Sobolev Smoothness Penalty")

tf0 <- plot.object(obj.name = "mod.tf.k0", scen.num = 4, p = 100, n = 500,
            title = "Trend Filter, k = 0")

tf1 <- plot.object(obj.name = "mod.tf.k1", scen.num = 4, p = 100, n = 500,
            title = "Trend Filter, k = 1")
tf2 <- plot.object(obj.name = "mod.tf.k2", scen.num = 4, p = 100, n = 500,
            title = "Trend Filter, k = 2")

