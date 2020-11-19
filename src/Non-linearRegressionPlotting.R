library(ggplot2)
library(gridExtra)
library(reshape)

setwd("~/Desktop/JiayiQu/UCD_PhD/LDproj")

# Ne_table.csv is a file saving Er2 for different combinations of parameters (shown in data folder)
df_Ne = read.csv("Ne_table.csv")
head(df_Ne)

df_Ne_u0 = df_Ne[df_Ne[, "u"] == 0, ] # mutation = 0
df_Ne_u1 = df_Ne[df_Ne[, "u"] == 1e-9, ] # mutation = 1e-9

myfun = function(a, b, c, Ne){
  res = 1/(a+b*c*Ne)
}

Hillfunc = function(c, N, u){
  rho = 4*N*c
  theta = 4*N*u
  (10 + rho + 4*theta)/(22 + 13*rho + 32*theta + rho^2 + 6*rho*theta + 8*theta^2)
}

recomb = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 6.25e-5)

df_u1 = data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_u1) = c("c", "a", "b", "mse_trans", "mse_sved", "mse_hill")
j = 1
for(i in recomb){
  dat = df_Ne_u1[df_Ne_u1[, "r"] == i, ]
  fit <- nls(E_rsqr ~ 1/(a + b*Ne*r), data = dat,
             start = list(a = 0.5, b = 3))
  a_coef = coef(fit)[1]
  b_coef = coef(fit)[2]
  df_u1[j, c(1:3)] = c(i, a_coef, b_coef)
  df_u1[j, "mse_trans"] = mean((dat$E_rsqr -  1/(a_coef + b_coef*dat$Ne*dat$r))^2)
  df_u1[j, "mse_sved"] = mean((dat$E_rsqr -  1/(1 + 4*dat$Ne*dat$r))^2)
  df_u1[j, "mse_hill"] = mean((dat$E_rsqr - Hillfunc(c = dat$r, N = dat$Ne, u = 1e-9))^2)
  
  j = j + 1
}

df_u1

df_u1[, 4:6] = signif(df_u1[, 4:6], 2)

k = c(5, 10, 14, 15)
recomb_u1 = c(0.1, 0.05, 0.01, 6.25e-5)

df_u1= df_u1[k, ]

mse_description = list()
mse_description[[1]] = c(expression(paste("Calibrated (", "6.5 x 10"^-7, ")")), expression(paste("Sved (", "2.8 x 10"^-3, ")")), expression(paste("Hill (", "6.2 x 10"^-4, ")")))
mse_description[[2]] = c(expression(paste("Calibrated (", "1.6 x 10"^-6, ")")), expression(paste("Sved (", "1.4 x 10"^-2, ")")), expression(paste("Hill (", "3.0 x 10"^-3, ")")))
mse_description[[3]] = c(expression(paste("Calibrated (", "3.3 x 10"^-5, ")")), expression(paste("Sved (", "1.3 x 10"^-1, ")")), expression(paste("Hill (", "2.2 x 10"^-2, ")")))
mse_description[[4]] = c(expression(paste("Calibrated (", "2.8 x 10"^-4, ")")), expression(paste("Sved (", "6.0 x 10"^-1, ")")), expression(paste("Hill (", "5.7 x 10"^-2, ")")))

c_list = list("c = 0.10", "c = 0.05", "c = 0.01",
              expression("c = 6.25 x 10"^-5))


j = 1
#plotlist = list()
for (i in recomb_u1){
  dat = df_Ne_u1[df_Ne_u1[, "r"] == i, ]
  p = ggplot(data = dat, aes(x = Ne, y = E_rsqr)) + geom_point() +
    stat_function(fun = myfun, args = list(a = df_u1$a[j], b = df_u1$b[j], c = i), aes(colour = "Trans")) +
    # Sved
    stat_function(fun = myfun, args = list(a = 1, b = 4, c = i), aes(colour = "Sved" )) +
    # Hill
    stat_function(fun = Hillfunc, args = list(u = 1e-9, c = i), aes(colour ="Hill" )) +  
    xlab("Ne") + ylab(expression("r"[E]^2)) +
    scale_x_continuous(expand = c(0, 5), limits = c(5, 52)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
    scale_color_manual("", breaks = c("Trans", "Sved", "Hill"), 
                       labels = c(mse_description[[j]][1], mse_description[[j]][2], mse_description[[j]][3]),  values = c("#E69F00", "red", "#56B4E9"))  + 
    ggtitle(c_list[[j]]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) + 
    theme(legend.justification=c(0.9,0.9), legend.position=c(.9,.9), legend.text=element_text(size=12), legend.title=element_blank())  
  #plotlist[[j]] = p
  pname <- paste0("compPlot_u1_r=",i, "_calibrate")
  ggsave(paste0(pname,".png"),p, width = 5, height = 5)
  j = j+1
}

#grid.arrange(grobs=plotlist,ncol=4)

###### u0
df_u0 = data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_u0) = c("c", "a", "b", "mse_trans", "mse_sved", "mse_hill")
j = 1
for(i in recomb){
  dat = df_Ne_u0[df_Ne_u0[, "r"] == i, ]
  fit <- nls(E_rsqr ~ 1/(a + b*Ne*r), data = dat,
             start = list(a = 0.5, b = 3))
  a_coef = coef(fit)[1]
  b_coef = coef(fit)[2]
  df_u0[j, c(1:3)] = c(i, a_coef, b_coef)
  df_u0[j, "mse_trans"] = mean((dat$E_rsqr -  1/(a_coef + b_coef*dat$Ne*dat$r))^2)
  df_u0[j, "mse_sved"] = mean((dat$E_rsqr -  1/(1 + 4*dat$Ne*dat$r))^2)
  df_u0[j, "mse_hill"] = mean((dat$E_rsqr - Hillfunc(c = dat$r, N = dat$Ne, u = 1e-9))^2)
  
  j = j + 1
}

df_u0

recomb_u0 = c(0.1, 0.05, 0.01, 6.25e-5)
l = c(5, 10, 14, 15)
df_u0[, 4:6] = signif(df_u0[, 4:6], 2)
df_u0= df_u0[l, ]

mse_description_u0 = list()
mse_description_u0[[1]] = c(expression(paste("Calibrated (", "6.5 x 10"^-7, ")")), expression(paste("Sved (", "7.2 x 10"^-4, ")")), expression(paste("Hill (", "3.9 x 10"^-3, ")")))
mse_description_u0[[2]] = c(expression(paste("Calibrated (", "4.8 x 10"^-5, ")")), expression(paste("Sved (", "1.9 x 10"^-3, ")")), expression(paste("Hill (", "1.4 x 10"^-2, ")")))
mse_description_u0[[3]] = c(expression(paste("Calibrated (", "2.0 x 10"^-4, ")")), expression(paste("Sved (", "5.4 x 10"^-3, ")")), expression(paste("Hill (", "9.2 x 10"^-2, ")")))
mse_description_u0[[4]] = c(expression(paste("Calibrated (", "4.8 x 10"^-11, ")")), expression(paste("Sved (", "1.1 x 10"^-5, ")")), expression(paste("Hill (", "3.0 x 10"^-1, ")")))

c_list_u0 = list("c = 0.10", "c = 0.05", 
              "c = 0.01", expression("c = 6.25 x 10"^-5))


j = 1
#plotlist_u0 = list()
for (i in recomb_u0){
  dat = df_Ne_u0[df_Ne_u0[, "r"] == i, ]
  p = ggplot(data = dat, aes(x = Ne, y = E_rsqr)) + geom_point() + 
    stat_function(fun = myfun, args = list(a = df_u0$a[j], b = df_u0$b[j], c = i), aes(colour = "Trans")) +
    # Sved
    stat_function(fun = myfun, args = list(a = 1, b = 4, c = i), aes(colour = "Sved" )) +
    # Hill
    stat_function(fun = Hillfunc, args = list(u = 1e-9, c = i), aes(colour ="Hill" )) +  
    xlab("Ne") + ylab(expression("r"[E]^2)) + 
    scale_x_continuous(expand = c(0, 5), limits = c(5, 52)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
    scale_color_manual("", breaks = c("Trans", "Sved", "Hill"), 
                       labels = c(mse_description_u0[[j]][1], mse_description_u0[[j]][2], mse_description_u0[[j]][3]),  values = c("#E69F00", "red", "#56B4E9"))  + 
    ggtitle(c_list_u0[[j]]) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) + 
    theme(legend.justification=c(0.9,0.9), legend.position=c(.9,.9), legend.text=element_text(size=12), legend.title=element_blank()) 
  #plotlist_u0[[j]] = p
  pname <- paste0("compPlot_u0_r=",i,"_calibrate")
  ggsave(paste0(pname,".png"),p, width = 5, height = 5)
  j = j+1
}

#grid.arrange(grobs=plotlist_u0,ncol=2)

