BD3 <- read.table(file='https://raw.githubusercontent.com/kornball/darinconb/master/Experimet3.csv',
                    header=TRUE, sep='\t')
attach(BD3)
head(BD3)

#GRAFICO DE BOXPLOTS DE VOLTAJE EN RELACION A LOS FACTORES
windows(50,50)
par(mfrow=c(1,2))
boxplot(BD3$Voltaje)
boxplot(Voltaje ~ Onda, data=BD3, xlab = "Longitud de onda", 
        ylab = "Voltaje", cex.lab = 2)
boxplot(Voltaje ~ Distancia, data=BD3, xlab = "Distancia", 
        ylab = "Voltaje", cex.lab = 2)
boxplot(Voltaje ~ PWM, data=BD3, xlab = "PWM", 
        ylab = "Voltaje", cex.lab = 2)

#GraFICO DE DENSIDAD
windows(50,75)
plot(density(BD3$Voltaje), xlab = "Voltaje", ylab = "Densidad", 
     main="", cex.lab=1.5, lwd = 2)

#GRAFICO DE CORRELACION
require(MPV)
windows(50,50)
pairs(BD3,
      upper.panel = panel.reg,
      diag.panel = panel.hist,
      lower.panel = panel.cor)

#OTROS GRAFICOS
with(BD3, symbols(x= BD3$PWM, y= BD3$Voltaje, circles = BD3$Onda, las=1,
                    inches= 0.1, fg="black", xlab = "Longitud de onda", ylab = "Voltaje"))

with(BD3, symbols(x= BD3$Distancia, y= BD3$Voltaje, circles = BD3$PWM, las=1,
                 inches= 0.1, fg="black", xlab = "Distancia", ylab = "Voltaje"))

#GRAFICO 3D DE DISPERSION
require(rgl)
with(BD3, plot3d(x=BD3$PWM, y=BD3$Distancia, z=BD3$Voltaje,
                xlab = "Longitud de onda", ylab = "Distancia", zlab = "Voltaje"))

#AALISIS DE LINEAL
require(gamlss)
mod.BD <- lm(Voltaje ~ Onda + Distancia + PWM, data = BD3)
summary(mod.BD)

#Grafico de interaccion
windows(50,50)
interaction.plot(x.factor = BD3$PWM, trace.factor = BD3$Distancia, response = BD3$Voltaje)
interaction.plot(x.factor = BD3$PWM, trace.factor = BD3$Onda, response = BD3$Voltaje)
library(scatterplot3d)
scatterplot3d(x=PWM, y=Distancia, z=Voltaje, pch=16, cex.lab=1.5,
              highlight.3d=TRUE, type="h")

# Construyendo la superficie de respuesta
install.packages("rsm")
library(rsm)
image(mod, Onda ~ PWM)
contour(mod, Onda ~ PWM)
persp(mod, Onda ~ PWM, zlab = "Voltaje")

persp(mod, Onda ~ PWM, col = "blue",
      bounds = list(temp=c(150, 300), conc=c(10, 30)),
      zlab = "Rendimiento", 
      contours = list(z="bottom", col="blue"),
      theta = -145, phi = 35)

#ENCONTRANDO LA MEJOR DISTRIBUCION
install.packages("galmss")
require(gamlss)
fit <- fitDist(Voltaje, data=BD3, type="realplus")
fit$fits
fit
par(mfrow=c(2, 2))
windows(50,50)
fitBCPEo<- histDist(y=Voltaje, family=BCPEo, main="BCPEo", cex.lab = 3)
fitBCCGo <- histDist(y=Voltaje, family=BCCGo, main="BCCGo",  cex.lab = 3)
fitBCTo <- histDist(y=Voltaje, family=BCTo, main="BCTo",  cex.lab = 3)
fitexGAUS <- histDist(y=Voltaje, family=exGAUS, main="exGAUS", cex.lab = 3)
fitWEI3 <- histDist(y=Voltaje, family=WEI3, main="WEI3", cex.lab = 3)
fitGA <- histDist(y=Voltaje, family=GA, main="GA", cex.lab = 3)
fitGIG <- histDist(y=Voltaje, family=GIG, main="GIG", cex.lab = 3)
fitLOGNO <- histDist(y=Voltaje, family=LOGNO, main="LOGNO", cex.lab = 3)

horizonte <- formula(~ PWM + Onda + Distancia + 
                       PWM * Distancia + PWM*Onda + Onda*Distancia +
                       I(PWM ^ 2) + 
                       I(Distancia ^ 2) +
                       I(Onda ^ 2))
## Ajuste del modelo inicial con la distribucion BCPEo
con1 <- gamlss.control(c.crit=0.001, n.cyc=5)
bcpe0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=BCPEo(),
                control=con1)
bcpe1 <- stepGAICAll.A(bcpe0, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte),
                       nu.scope=list(lower= ~ 1, upper=horizonte),
                       tau.scope=list(lower= ~ 1, upper=horizonte))
bcpe1 <- refit(bcpe1)

## Ajuste del modelo inicial con la distribucion BCCGo
bccg0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=BCCGo())
bccg1 <- stepGAICAll.A(bccg0, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte),
                       nu.scope=list(lower= ~ 1, upper=horizonte))
bccg1 <- refit(bccg1)

## Ajuste del modelo inicial con la distribucion BCTo
con1 <- gamlss.control(c.crit=0.001, n.cyc=100)
bct0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=BCTo(),
               control=con1)
bct1 <- stepGAICAll.A(bct0, trace=F,
                      scope=list(lower= ~ 1, upper=horizonte),
                      sigma.scope=list(lower= ~ 1, upper=horizonte),
                      nu.scope=list(lower= ~ 1, upper=horizonte),
                      tau.scope=list(lower= ~ 1, upper=horizonte))
bct1 <- refit(bct1)
summary(bct1)
## Ajuste del modelo inicial con la distribucion WEI3
wei30 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=WEI3())
wei31 <- stepGAICAll.A(wei30, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion GA
ga0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=GA())
ga1 <- stepGAICAll.A(ga0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion GIG
gig0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=GIG())
gig1 <- stepGAICAll.A(gig0, trace=F,
                      scope=list(lower= ~ 1, upper=horizonte),
                      sigma.scope=list(lower= ~ 1, upper=horizonte),
                      nu.scope=list(lower= ~ 1, upper=horizonte))
gig1 <- refit(gig1)

## Ajuste del modelo inicial con la distribucion LOGNO
logno0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, 
                 family=LOGNO(mu.link='log'))
logno1 <- stepGAICAll.A(logno0, trace=F,
                        scope=list(lower= ~ 1, upper=horizonte),
                        sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion IG
ig0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=IG())
ig1 <- stepGAICAll.A(ig0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion exGauss
Exga <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD5, family=exGAUS())
Exga1 <- stepGAICAll.A(Exga, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte),
                       nu.scope=list(lower= ~ 1, upper=horizonte))
Exga1 <- refit(Exga1)


## Ajuste del modelo inicial con la distribucion EXP
exp0 <- gamlss(produccion ~ 1, data=datos, family=EXP())
exp1 <- stepGAICAll.A(exp0, trace=F,
                      scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion PARETO2
con1 <- gamlss.control(c.crit=0.001, n.cyc=10000)
pare0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=PARETO2(),
                control=con1)
pare1 <- stepGAICAll.A(pare0, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte))



## Ajuste del modelo inicial con la distribucion NO
no0 <- gamlss(Voltaje ~ 1, sigma.fo= ~ 1, data=BD3, family=NO)
no1 <- stepGAICAll.A(no0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte))
no1 <- refit(no1)

#PLOTs y AIC Todos
plot(bct1) 
plot(Exga1)
plot(bccg1)
plot(wei31, which=1)
plot(gg1)

AIC(bccg1, bct1, wei31,
    Exga1, k=log(nrow(BD5)))

#Worm Plots
windows(50,50)
par(mfrow=c(2, 2), bg='white')
wp(bct1)
title("bct1")
wp(Exga1)
title("Exga1")
wp(bccg1, main = "bccg1", cex = 1, cex.lab = 2)
title("bccg1")
wp(wei31)
title("wei31")

#Diferenciando por colores
Amarillo <- BD3[c(1:1272), c(1,2,4)]
Rojo <- BD3[c(1, 1273:2551), c(1,2,4)]
Azul <- BD3[c(1, 2552:3826), c(1,2,4)]
Verde <- BD3[c(1, 2827:5101), c(1,2,4)]
windows(50,50)
par(mfrow=c(2,2))
plot(Amarillo$Voltaje, Amarillo$PWM, col="yellow", 
     main = "Datos para color amarillo", xlab = "PWM", ylab = "Voltaje", cex.lab = 1.5, cex.main = 1.5)
plot(Rojo$Voltaje, Rojo$PWM, col="red",
     main = "Datos para color rojo",xlab = "PWM", ylab = "Voltaje", cex.lab = 1.5, cex.main = 1.5)
plot(Azul$Voltaje, Azul$PWM, col="blue",
     main = "Datos para color azul",xlab = "PWM", ylab = "Voltaje", cex.lab = 1.5, cex.main = 1.5)
plot(Verde$Voltaje, Verde$PWM, col="green",
     main = "Datos para color verde",xlab = "PWM", ylab = "Voltaje", cex.lab = 1.5, cex.main = 1.5)
 