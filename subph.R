
subph<-function(t,data,start,param,mod,v,TypeNEW,knotsNEW,variablesNEW,mb,pb,ncol1,marki){#needs the param and mod parameter from the main function, the data here should be only data 
  #for each individual ID, also need parameters ncol, needs the ncol1 parameter from main function, 
  
  ##only produce the integration of variables having TD effect, all the other effects do not need the integration to save computer time 
  
  
  subbeta<-rep(NA,v)          
  subgx<-rep(NA,v)
  j2<-0
  j5<-0
  j4<-0
  j6<-0
  j7<-0
  j8<-0
  for(i in 1:v){
    # if no WCE+NL and TD effect
    if(mod$NL[i]==0 & mod$TD[i]==0 & mod$WCE[i]==0){
      
      subbeta[i]<-mod$coefficients[i]
      subgx[i]<-data[data[,TypeNEW[1]]==start,variablesNEW[i]]
      
    }else if(mod$NL[i]==1 & mod$TD[i]==0 & mod$WCE[i]==0){ ## only NL
      j2<-j2+1
      subbeta[i]<-1
      subgx[i]<-data[data[,TypeNEW[1]]==start, ncol1+j2]
      
    }else if(mod$NL[i]==0 & mod$TD[i]==1 & mod$WCE[i]==0){ ## only TD
      subbeta[i]<-tdest(t,mod,i)
      subgx[i]<-data[data[,TypeNEW[1]]==start,variablesNEW[i]]
      
    }else if(mod$NL[i]==0 & mod$TD[i]==0 & mod$WCE[i]==1){ ## only WCE
      ## calculated the weighted sum first, then this value can be used the same way as the TVC
      j5<-j5+1
      subbeta[i]<-mod$coefficients[i]
      subgx[i]<-data[data[,TypeNEW[1]]==start, ncol1+marki[1]+j5]
      
    }else if(mod$NL[i]==1 & mod$TD[i]==1 & mod$WCE[i]==0){ ##NL and TD
      j4<-j4+1
      subbeta[i]<-tdest(t,mod,i)
      subgx[i]<-data[data[,TypeNEW[1]]==start, ncol1+marki[2]+j4]
      
    }else if(mod$NL[i]==1 & mod$TD[i]==0 & mod$WCE[i]==1){##NL and WCE
      j6<-j6+1
      subbeta[i]<-1
      subgx[i]<-data[data[,TypeNEW[1]]==start, ncol1+marki[3]+j6]
      
    }else if(mod$NL[i]==0 & mod$TD[i]==1 & mod$WCE[i]==1){ ## TD+WCE
      j7<-j7+1  
      
      subbeta[i]<-tdest(t,mod,i)
      subgx[i]<-data[data[,TypeNEW[1]]==start, ncol1+marki[4]+j7]
      
    }else{ #NL+tD+WCE
      j8<-j8+1
      
      subbeta[i]<-tdest(t,mod,i)
      subgx[i]<-data[data[,TypeNEW[1]]==start, ncol1+marki[5]+j8]
      
    }
  }
  
  subbase<-splineDesign(knotsNEW,t,ord=pb+1)
  
  subpht<-exp(param%*%t(subbase)+subbeta%*%subgx)
  return(subpht)
}


