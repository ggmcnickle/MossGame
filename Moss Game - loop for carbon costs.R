#Gordon G. McNickle 2014
#Biology department
#Wilfrid Laurier Unversity
#gmcnickle@wlu.ca
rm(list=ls(all=TRUE))
library(scatterplot3d)
library(rootSolve)
graphics.off()

###############################################
#Parameters
###############################################
## Rn    = N available, e.g. N concentration uM / m-2 / yr
## Rc    = C available, e.g. max photosynthesis with no competition / m-2/yr
## ut    = tissue strategy (i.e. moss "leaf")
## uh    = height strategy
## ctn   = tissue cost in N units 
## ctc	 = tissue cost in C units (includes respiration)
## chn	 = height cost in N units
## chc 	 = height cost in C units (includes respiration)
## zt	 = size asymmetry constant for tissue carbon competition
## zh	 = size asymmetry constant for height competition
##		NOTE: Larger numbers of z make the size-asymmetry more severe
## alph	 = represents C:N ratio of whole plant. Or plant need for C relative to N. 
## beta  = 1-alph (i.e. this is the N:C ratio)	
################################################

moss.game<-function(u) {
		with (as.list(params),	{
	#THIS IS CURRENTLY COMPETITION BETWEEEN TWO MOSS INDIVIDUALS OF THE SAME SPECIES
	#NOTE: u[1] is tissue & u[2] is height for Moss individual 1, 
	#  and u[3] is tissue & u[4] is height for moss individual 2.

		#NITROGEN GAME		
	    	t   <- u[1]+u[3]	#Total tissue (i.e. total moss "leaves")
		Hn  <- Rn*(1-exp(-t))	#harvest of N by both mosses
		f1  <- u[1]/(t) 	#Moss 1 share of N harvest
		f2  <- u[3]/(t)		#Moss 2 share of N harvest
		
		#N DERIVATIVES	
		dHn <- Rn*exp(-t) 	#Derivative of Hn with respect to u (any)
		df1 <- 1/t-u[1]/(t^2)	#Derivative of f1 with respect to u1
		df2 <- 1/t-u[3]/(t^2)	#Derivative of f2 with respect to u3

		#CARBON GAME
	    	t   <- u[1]+u[3]			#Total tissue (i.e. total moss "leaves"
		tz  <- (u[1]^zt)*(u[2]^zh) + (u[3]^zt)*(u[4]^zh) 	#simplifies coding (denomenator of z(u)) 
		Hc  <- Rc*(1-exp(-t))			#harvest of C by both plants
		z1  <- ((u[1]^zt)*(u[2]^zh))/tz		#P1 share of C harvest
		z2  <- ((u[3]^zt)*(u[4]^zh))/tz		#P2 share of C harvest
		
		#C DERIVATIVES	
		dHc   <- Rc*exp(-t) 	#Derivative of Hc with respect to u (any)
		dz1.t <- (zt*(u[2]^(zh))*(u[1]^(zt-1))*(u[4]^zh)*(u[3]^zt))/(tz^2) 	#Derivative of z1 with respect to u1
		dz2.t <- (zt*(u[4]^(zh))*(u[3]^(zt-1))*(u[2]^zh)*(u[1]^zt))/(tz^2)	#Derivative of z2 with respect to u3
		dz1.h <- (zh*(u[2]^(zh-1))*(u[1]^(zt))*(u[4]^zh)*(u[3]^zt))/(tz^2)	#Derivative of z1 with respect to u2
		dz2.h <- (zh*(u[4]^(zh-1))*(u[3]^(zt))*(u[2]^zh)*(u[1]^zt))/(tz^2)	#Derivative of z2 with respect to u4

		#PROFIT FUNCTIONS
		P1n <- f1*Hn - ctn*u[1] - chn*u[2]   #Net profit N moss1
		P2n <- f2*Hn - ctn*u[3] - chn*u[4]   #Net profit N moss2
		P1c <- z1*Hc - ctc*u[1] - chc*u[2]   #Net profit c moss1
		P2c <- z2*Hc - ctc*u[3] - chc*u[4]   #Net profit c moss2

		#PROFIT DERIVATIVES
		dP1n.dut <- df1*Hn + f1*dHn - ctn
		dP1n.duh <- -chn

		dP2n.dut <- df2*Hn + f2*dHn - ctn
		dP2n.duh <- -chn

		dP1c.dut <- dz1.t*Hc + z1*dHc -ctc
		dP1c.duh <- dz1.h*Hc -chc

		dP2c.dut <- dz2.t*Hc + z2*dHc -ctc
		dP2c.duh <- dz2.h*Hc -chc
				
	#DERIVATIVES G-function where, G=(Ps^alph)*(Pr^beta)
		dG1.dut <- alph*(P1c^(alph-1)) * (P1n^beta)*dP1c.dut +
			  beta*(P1c^alph) * (P1n^(beta-1))*dP1n.dut 	#Moss 1 derivative of G with respect to u1
		dG1.duh <- alph*(P1c^(alph-1)) * (P1n^beta)*dP1c.duh +
			  beta*(P1c^alph) * (P1n^(beta-1))*dP1n.duh	#Moss 1 derivative of G with respect to u2

		dG2.dut <- alph*(P2c^(alph-1)) * (P2n^beta)*dP2c.dut +
			  beta*(P2c^alph) * (P2n^(beta-1))*dP2n.dut	#Moss 2 derivative of G with respect to u3
		dG2.duh <- alph*(P2c^(alph-1)) * (P2n^beta)*dP2c.duh +
			  beta*(P2c^alph) * (P2n^(beta-1))*dP2n.duh	#Moss 2 derivative of G with respect to u4


		return(c(dG1.dut = dG1.dut, dG1.duh = dG1.duh, 
			dG2.dut = dG2.dut, dG2.duh = dG2.duh))
		})} 

#Loop to iterate ESS sover over different amounts of N and C available. 

#Create empty vectors to store output
ct.cost<-as.numeric() 	#empty C vector
ch.cost<-as.numeric() 	#empty N vector
Tissue<-as.numeric()	#empty Tissue vector
Height<-as.numeric()	#empty Height vector

#loop parameters
ct.int<-0.05		#Interval to increase C in loop
ch.int<-0.05		#Interval to increase N in loop
ct.max<-3		#Max C in loop,
ct.min<-.9		#Min C in loop

ch.max<-3		#Max N in loop,
ch.min<-.9		#Min N in loop
ct<-ct.min		#start Rc @ Rc.min		
ch<-ch.min		#start Rn @ Rn.min

alpha<-0.983		#alpha value where, beta=(1-alpha)


#Loop
while (ct<ct.max) {
	ct<-ct+ct.int
	while (ch<ch.max) {
		ch<-ch+ch.int
		params<- c(Rc=250, Rn=50, zt=1.5, zh=1.5, ctc=ct, chc=ch, ctn=.03, chn=0.03, alph=alpha, beta=(1-alpha)) 
		ct.cost<-c(ct.cost, ct)
		ch.cost<-c(ch.cost,ch)
		ESS<-multiroot(f = moss.game, start = c(.1,.1,.1,.1), maxiter=100000, positive = TRUE)$root
		Tissue<-c(Tissue,ESS[1])
		Height<-c(Height,ESS[2])
		}
		ch<-ch.min
	}

out<-data.frame(ct.cost, ch.cost, Tissue, Height)
write.table(out, "d:/out.moss.txt", sep="\t")

Total<-Tissue+Height
fTissue<-Tissue/Total
fHeight<-Height/Total

#plot results
dev.new()
scatterplot3d(x = ct.cost, y = ch.cost, z = Total, highlight.3d=TRUE)
dev.new()
scatterplot3d(x = ct.cost, y = ch.cost, z = Tissue, highlight.3d=TRUE)
dev.new()
scatterplot3d(x = ct.cost, y = ch.cost, z = Height, highlight.3d=TRUE)



