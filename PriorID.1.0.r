
# PriorID
# version 1.0

# functions
# rankings
# rankings.antro
# rankings.cv
# rankings.cv.ran

# Functions
###############################################################################
rankings<-function(V,S){
N<-length(V[,1])
m<-max(V[,2])
I<-length(S)
VV<-matrix(ncol=2,nrow=N)
VV[,1]<-V[,1]
VV[,2]<-V[,2]
V<-VV
V_id<-vector(length=N,mode='numeric')
for(i in 1:I){V_id[S[i]]<-1}
# counting locker [i], the total of potential victims associated with locker i and locker_i [i], total of already solved cases associated with locker i.
# para el locker i contamos locker[i], total de posibles victimas asociadas al locker i y locker_i[i], total de casos resueltos asociados al locker i.
locker<-vector(length=m,mode='numeric')
locker_i<-vector(length=m,mode='numeric')
for(i in 1:N){a<-V[i,2];locker[a]<-locker[a]+1}
for(i in 1:I){b<-V[S[i],2];locker_i[b]<-locker_i[b]+1}
nk<-max(locker)
alpha<-vector(length=m,mode='numeric')
cte<-N/nk-2
alpha<-locker/N*cte
alpha0<-sum(alpha)
P_prior<-alpha/alpha0
P_post<-(alpha+locker_i)/(alpha0+I)

# defining odds_prior and odds_posterior for each individual
odds_prior_uni<-vector(mode="numeric",length=N)
post<-vector(mode="numeric",length=N)
prior<-vector(mode="numeric",length=N)
odds_post<-vector(mode="numeric",length=N)
cte<-locker-locker_i

for(i in 1:N){
if(V_id[i]==0){
loc<-V[i,2]
ifelse(cte[loc]!=0,prior[i]<-P_prior[loc]/cte[loc],0)
ifelse(cte[loc]!=0,post[i]<-P_post[loc]/cte[loc],0)}
}
odds_post<-round(post/(1-post),5)
for(i in 1:N){odds_prior_uni[i]<-ifelse(V_id[i]==0,round(1/(N-I-1),5),0)}
a<-data.frame(V,V_id,odds_prior_uni,odds_post)
b<-subset(a,V_id==0)
names(b)<-c("ID","locker","identified","odds_prior","odds_posterior")
c<-order(b$odds_post,decreasing=T)
b<-b[c,]

# plotting:
pdf("rankings.pdf")
plot(b$odds_post,xlab="Possible Victims",ylab="ODDS Posterior",cex=0.5,pch=19,main="Updated ODDS after learning about the solved cases")
abline(h=odds_prior_uni,col="red")
legend("topright",col="red","Odds Prior",lty=1)
dev.off()

# printing result:
b
}



################################################################################
# Function that generates a prioritized ranking when  anthropological information of the skeleton is used (gender and age estimated at the time of death) to select compatible cases.

rankings.antro<-function(Va,S,gender_s,age_s){
V<-Va
N<-length(V[,1])
m<-max(V[,2])
I<-length(S)
V_id<-vector(length=N,mode='numeric')
for(i in 1:I){V_id[S[i]]<-1}
# counting locker [i], the total of potential victims associated with locker i and locker_i [i], total of already solved cases associated with locker i.
locker<-vector(length=m,mode='numeric')
locker_i<-vector(length=m,mode='numeric')
for(i in 1:N){a<-V[i,2];locker[a]<-locker[a]+1}
for(i in 1:I){b<-V[S[i],2];locker_i[b]<-locker_i[b]+1}

nk<-max(locker)

alpha<-vector(length=m,mode='numeric')
cte<-N/nk-2
alpha<-locker/N*cte
alpha0<-sum(alpha)
P_prior<-alpha/alpha0

P_post<-(alpha+locker_i)/(alpha0+I)

# defining odds_prior and odds_posterior for each individual
odds_prior_uni<-vector(mode="numeric",length=N)
post<-vector(mode="numeric",length=N)
prior<-vector(mode="numeric",length=N)
odds_post<-vector(mode="numeric",length=N)

cte<-locker-locker_i

for(i in 1:N){
if(V_id[i]==0){
loc<-V[i,2]
ifelse(cte[loc]!=0,prior[i]<-P_prior[loc]/cte[loc],0)
ifelse(cte[loc]!=0,post[i]<-P_post[loc]/cte[loc],0)}
}

# selecting compatible cases with gender  and age of the skeletal remains (by extending the range in two years in upper and lower limits)
a<-data.frame(V,V_id,post)
b<-subset(a,V_id==0)
names(b)<-c("ID","locker","age","gender","identified","post")
c<-subset(b,b$gender==gender_s & b$age<=age_s[2]+2 & b$age>=age_s[1]-2)
d<-order(c$post,decreasing=T)
e<-c[d,]

B<-sum(e$post)
e$post<-e$post/B
Nef<-length(e$post)

e$odds_post<-round(e$post/(1-e$post),5)
for(i in 1:Nef){e$odds_prior_uni[i]<-round(1/(Nef-1),5)}

# plotting:
pdf("rankings.pdf")
plot(e$odds_post,xlab="Possible Victims",ylab="ODDS Posterior",cex=0.5,pch=19,main="Updated ODDS after learning about the solved cases")
abline(h=e$odds_prior_uni,col="red")
legend("topright",col="red","Odds Prior",lty=1)
dev.off()

# printing result:
e
}

###############################################################################
rounds<-10
rankings.cv<-function(V,S,rounds){
N<-length(V[,1])
m<-max(V[,2])
I<-length(S)
V_id<-vector(length=N,mode='numeric')
VV<-matrix(ncol=2,nrow=N)
VV[,1]<-V[,1]
VV[,2]<-V[,2]
V<-VV
I_s<-round(I/4)
sc<-0.0
ef<-0.0
sc2<-0.0
ef2<-0.0
for(r in 1:rounds){
S_r<-sample(S,I_s)
S_l<-S[(S %in% S_r)==FALSE]

b<-rankings(V,S_l)

up<-b$ID[b$odds_posterior>b$odds_prior]
up_r<-length(S_r[S_r %in% up])
up<-length(up)

sc<-sc+up_r/I_s
ef<-ef+up_r/up
sc2<-sc2+(up_r/I_s)**2
ef2<-ef2+(up_r/up)**2
}
sc<-sc/rounds
ef<-ef/rounds
s_sc<-sqrt((sc2-rounds*sc**2)/(rounds-1))
s_ef<-sqrt((ef2-rounds*ef**2)/(rounds-1))

result<-c(sc,s_sc,ef,s_ef)
names(result)<-c("success","deviation success","effectiveness","deviation effectiveness")

# printing results
result
}

# cross-validation y compara con elegir al azar los casos que sobresalen:
rankings.cv.random<-function(V,S,rounds){
  N<-length(V[,1])
  m<-max(V[,2])
  I<-length(S)
  VV<-matrix(ncol=2,nrow=N)
  VV[,1]<-V[,1]
  VV[,2]<-V[,2]
  V<-VV
  V_id<-vector(length=N,mode='numeric')
  I_s<-round(I/4)
  sc<-0.0
  ef<-0.0
  sc2<-0.0
  ef2<-0.0
  rsc<-0.0
  ref<-0.0
  rsc2<-0.0
  ref2<-0.0
  for(r in 1:rounds){
    S_r<-sample(S,I_s)
    S_l<-S[(S %in% S_r)==FALSE]
    
    b<-rankings(V,S_l)
    
    up<-b$ID[b$odds_posterior>b$odds_prior]
    up_r<-length(S_r[S_r %in% up])
    up<-length(up)
    
    sc<-sc+up_r/I_s
    ef<-ef+up_r/up
    sc2<-sc2+(up_r/I_s)**2
    ef2<-ef2+(up_r/up)**2
    
    c<-b
    c$odds_posterior<-sample(b$odds_posterior)
    
    rup<-c$ID[c$odds_posterior>c$odds_prior]
    rup_r<-length(S_r[S_r %in% rup])
    rup<-length(rup)
    
    rsc<-rsc+rup_r/I_s
    ref<-ref+rup_r/rup
    rsc2<-rsc2+(rup_r/I_s)**2
    ref2<-ref2+(rup_r/rup)**2
  }
  sc<-sc/rounds
  ef<-ef/rounds
  s_sc<-sqrt((sc2-rounds*sc**2)/(rounds-1))
  s_ef<-sqrt((ef2-rounds*ef**2)/(rounds-1))
  
  rsc<-rsc/rounds
  ref<-ref/rounds
  rs_sc<-sqrt((rsc2-rounds*rsc**2)/(rounds-1))
  rs_ef<-sqrt((ref2-rounds*ref**2)/(rounds-1))
  
  result<-c(sc,s_sc,ef,s_ef,rsc,rs_sc,ref,rs_ef)
  names(result)<-c("success","deviation success","effectiveness","deviation effectiveness","random success","random dev success","random effectiveness","random dev effectiveness")
  
  # printing results
  result
}



# examples
gender_s<-"F"
age_s<-c(20,28)
N<-100
total_lockers<-10
total_id<-10
Va<-matrix(ncol=3,nrow=N)
Va[,2]<-floor(runif(N)*total_lockers+1)
Va[,1]<-seq(1,N)
Va_gender<-ifelse(rbinom(N,1,0.5)==0,"F","M")
Va[,3]<-rpois(100,23)
Va<-data.frame(Va,Va_gender)
names(Va)<-c("ID","locker","age","gender")
S<-c(90,1,25,94,61,87,11,12,44,91)
# rankings()
# total de posibles victimas
N<-100
total_lockers<-10
total_id<-10
# definimos V
V<-matrix(ncol=2,nrow=N)
V[,2]<-floor(runif(N)*total_lockers+1)
V[,1]<-seq(1,N)
# definimos S (puede no estar ordenado)
#S<-floor(runif(total_id)*N+1)
S<-c(90,11,25,94,61,87,31,13,4,100)


# examples rankings:
rankings(V,S)
rankings.antro(Va,S,"F",c(15,20))
rankings.cv(V,S,10)
# example improving the succes by learning from a larger set of already-solved cases S (80 per cent)
S<-sample(V[,1],80,replace=FALSE)
rankings.cv(V,S,10)
rankings.cv.random(V,S,10)
###################################################################################
###################################################################################
