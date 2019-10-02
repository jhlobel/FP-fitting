#uses the following package
library(minpack.lm)

#sample data. All concentrations are in micromolar.
#x is the protein concentration
#y is the millipolarization (mP) signal. This is a representative single experiment.
x = c(12.5, 6.25 ,3.125 ,1.5625 ,0.78125, 0.390625, 0.1953125, 0.09765625 ,0.048828125, 0.024414063, 0.012207031, 0.006103516,0.003051758,0.001525879,0.000762939,0.00038147,0.000190735,9.53674E-05,4.76837E-05)
y = c(344, 352, 350, 348, 340, 331, 319, 295, 267, 227, 194, 164,140,131,124,115,113,115,115)

#y is the millipolarization (mP) signal, x is protein concentration, k is the dissocation constant (kd). 
#this is used to fit when kd > probe concentration.
lin = nlsLM(y~(max(y)-min(y))*x/(x+k)+min(y),start=list(k=1))

#prints fits and other statistical information
summary(lin)

#y is the mP signal, x is protein concentration, probe is the probe concentration, k is the dissocation constant (kd), and n is the hill coefficient. 
#this is used to fit when kd and probe concentrations are similar and the binding isotherm exhibits coopertivity.
probe = 0.0005
lin = nlsLM(y~(max(y)-min(y))*(((x^n+probe^n+k^n)-sqrt((x^n+probe^n+k^n)^2-(4*x^n*probe^n)))/(2*probe^n))+min(y),start=list(k=1,n=1))

#prints fits and other statistical information
summary(lin)