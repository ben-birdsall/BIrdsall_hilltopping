elv=c(0, 400)           
lndscp=150                 
nindvs=50                 
nsteps=500                 
mvup=0.8

prod.land = function(elv){
  land  = matrix(ncol=lndscp, nrow=lndscp)
  ypk = round(runif(1,0,lndscp))
  xpk = round(runif(1,0,lndscp))
  land[xpk, ypk] = elv[2]
land[xpk, 1:(ypk-1)] = round(seq(elv[1], elv[2], (elv[2]-elv[1])/(ypk-2)) + rgamma((ypk-1),0,1), 0)
land[xpk, (ypk+1):lndscp] = round(rev(seq(elv[1], elv[2], (elv[2]-elv[1])/(lndscp-ypk-1)) + rgamma((lndscp-ypk),5,2)), 0)
for(r in (xpk-1):1){
land[r,] = land[(r+1),] - round(rgamma(lndscp, 10, 1), 0)
  }
for(r in (xpk+1):lndscp){
land[r,] = land[(r-1),] - round(rgamma(lndscp, 10, 1), 0)
  }
return(land)
}
land = prod.land(elv)
windows()
image(land)  

pop = function(nindvs, lndscp){
  bflies = matrix(nrow=nindvs, ncol=2)
  variance=25
  x = round(runif(1,min = 0,max = lndscp-variance))
  y = round(runif(1,min = 0,max = lndscp-variance))
  bflies[,1]  = sample(x:(x+variance),nindvs,replace = T)
  bflies[,2]  = sample(y:(y+variance),nindvs,replace = T)
  return(bflies)
}
bflies = pop(nindvs, lndscp)
points(bflies[,1]/150, bflies[,2]/150, pch=12,cex=2.0)

moving=function(x){
  mv=matrix(NA,nsteps,2)
  mv[1,]=bflies[x,]
  possmv=NA
for (k in 2:nsteps) {
  possmv[k]=rbinom(1,1,mvup)
  near=land[ifelse((mv[k-1,1])==1,1,(mv[k-1,1]-1)):ifelse((mv[k-1,1])==150,150,(mv[k-1,1]+1)),ifelse((mv[k-1,2])==1,1,(mv[k-1,2]-1)):ifelse((mv[k-1,2])==150,150,(mv[k-1,2]+1))]
  trans=max(which(near==max(near)))
  mv[k,2]=if(possmv[k]==1){ifelse(trans %in% c(1,2,3),mv[(k-1),2]-1,ifelse(trans %in% c(7,8,9),mv[(k-1),2]+1,mv[(k-1),2]))
    } else mv[(k-1),2]+sample(c(1,0,-1),1)
      mv[k,1]=if(possmv[k]==1) {ifelse(trans %in% c(1,4,7), mv[(k-1),1]-1,ifelse(trans %in% c(3,6,9), mv[(k-1),1]+1,mv[(k-1),1]))
        } else mv[(k-1),1]+sample(c(1,0,-1),1)
}  
  return(mv)
}

windows()
image(land)
points(bflies[,1]/150, bflies[,2]/150, pch=12,cex=2.0)
movements=lapply(1:nindvs,moving)
corridors=lapply(movements,function(x)x/150)
lapply(corridors, lines)


