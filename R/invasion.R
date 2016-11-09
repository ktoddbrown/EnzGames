invasion <- function(parm, tradeoff, rel_invader=0.1, timeFlag=c('all', 'last')[2],
                     timearr = c(0:7, (1:12)*30, 365*(1:10), 365*10*(2:10))){
  
  parm$v2 <- tradeoff(b=parm$b, vmax=parm$vmax, cue=parm$cue2)
  parm$v1 <- tradeoff(b=parm$b, vmax=parm$vmax, cue=parm$cue1)
  
  parm$C.ss <- with(parm, h*k/(cue1*v1-h))
  parm$B1.ss <- with(parm, cue1/h*(I-m*C.ss))
  
  dC <- function(t, y, parms){
    B1 <- y[1]; B2 <- y[2]; C <- y[3]
    ans <- with(parms,{
      c(B1=B1*(cue1*v1*C/(k+C) - h), 
        B2=B2*(cue2*v2*C/(k+C) - h),  
        C=I-C*(m+(v1*B1+v2*B2)/(k+C)))
    })
    names(ans) <- c('B1', 'B2', 'C')
    return(list(ans))
  }
  
  
  evolution <- lsoda(y=c(B1=parm$B1.ss, B2=parm$B1.ss*rel_invader, C=parm$C.ss), 
                     times=timearr, func=dC, parms=parm)
  
  temp <- as.data.frame(evolution)
  temp[abs(temp) < 1e-8] <- 0 #cut off absolute tol
  temp$time <- timearr
  
  if(identical(timeFlag, 'last')){
    return(data.frame(parm, temp[dim(temp)[1],]))
  }else{
    return(list(parm, temp))
  }
}