ants.ec<-function(stat,fieldtype,df){
  ec<-c(0,0,0,0)
  t<-stat
  if(fieldtype=="T"){
    a<-gamma((df+1)/2)/(((df*pi)^(1/2))*gamma(df/2))
    tintegral<-function(t){a*((1+(t^2)/df)^(-1/2*(df+1)))}
    iec<-integrate(tintegral,lower=t,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/(2*pi))*((1+((t^2)/df))^(-1/2*(df-1)))
    ec[3]<-((((4*log(2))*lgamma((df+1)))/(((2*pi)^(3/2))*((df/2)^(1/2))*lgamma(df/2)))*((1+((t^2)/2))^(-1/2*(df-1))))*t
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^2))*((1+((t^2)/df))^(-1/2*(df-1)))*((((df-1)/df)*(t^2))-1)
  }else if(fieldtype=="F"){
    a<-(gamma((df+k-2)/2)/(gamma(df/2)*gamma(k/2)))
    fintegral<-function(t){a*(((k*t/df)^(1/2*(k-2)))*(1+(k*t/2)^(-1/2*(df+k))))}
    iec<-integrate(fintegral,lower=t,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/((2*pi)^(1/2)))*(((lgamma(df+k-1)/2)*2^(1/2))/(lgamma(df/2)*lgamma(k/2)))*((k*t/df)^(1/2(k-2)))*((1+(k*t/v))^(-1/2*(df+k-2)))
    ec[3]<-((4*log(2))/(2*pi))*((lgamma((v+k-2)))/(lgamma(df/2)*lgamma(k/2)))*((k*t/df)^(1/2*(k-2)))*(1+(k*t/df))^(-1/2(df+k-2))*((v-1)*(k*t/v)-(k-1))
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^(3/2)))*(((lgamma((v+k-2)))*2^(-1/2))/(lgamma(df/2)*lgamma(k/2)))*((k*t/df)^(1/2*(k-3)))*((1+(k*t/df))^(-1/2*(df+k-2)))*((df-1)*(df-2)*((k*t/v)^(2))-(2*df*k-df-k-1)*(k*t/df)+(k-1)*(k-2))
  }else if(fieldtype=="X"){
    xintegral<-function(t){((t^(1/2))*e^(-1/2*t))/((2^(df/2))*(gamma(df/2)))}
    iec<-integrate(xintegral,lower=t,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/(2*pi))*(((t^(1/2*(df-1)))*(e^(-1/2*t)))/((2^(1/2(df-2)))*gamma(df/2)))
    ec[3]<-((4*log(2))/(2*pi))*(((t^(1/2*(df-2)))*(e^(-1/2*t)))/((2^(1/2*(df-2)))*gamma(df/2)))*(t-(v-1))
    ec[4]<-((4*log(2)^(3/2))/((2*pi)^(3/2)))*(((t^(1/2*(df-3)))*(e^(-1/2*t)))/((2^(1/2*(df-2)))*(gamma(df/2))))*((t^2)-((2*df-1)*t)+(df-1)*(df-2))
  }else if(fieldtype=="G"){  
    gintegral<-function(t){((1/((2*pi)^(1/2)))*e^(-1/2*(t^2)))}
    iec<-integrate(gintegral,lower=t,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-((((4*log(2))^(1/2))/(2*pi))*e^(-1/2*((t^2))))
    ec[3]<-((4*log(2))/(2*pi))*t*(e^(-1/2*(t^2)))
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^2))*((t^2)-1)*(e^(-1/2*(t^2)))
  }else{
      stop("Incorrect input for FieldType")
      }
  return(ec)
}
