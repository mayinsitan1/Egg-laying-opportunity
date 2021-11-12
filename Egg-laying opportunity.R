#####-----Egg-laying opportunity-----#####
#####-----Liang Ma-----#####
#####-----Nov 2021-----#####
#overwrite files Rblas.dll and Rlapack.dll in directory Program Files/R/R-3.3.3/bin/x64 with Intel MKL counterpart (for hardware acceleration)
#see https://stackoverflow.com/questions/38090206/linking-intels-math-kernel-library-mkl-to-r-on-windows

library(ncdf4)
library(parallel)
mclapply.hack <- function(...) {
  ## Create a cluster
  size.of.list <- length(list(...)[[1]])
  cl <- makeCluster( min(size.of.list, detectCores()) )
  
  ## Find out the names of the loaded packages 
  loaded.package.names <- c(
    ## Base packages
    sessionInfo()$basePkgs,
    ## Additional packages
    names( sessionInfo()$otherPkgs ))
  
  tryCatch( {
    
    ## Copy over all of the objects within scope to
    ## all clusters. 
    this.env <- environment()
    while( identical( this.env, globalenv() ) == FALSE ) {
      clusterExport(cl,
                    ls(all.names=TRUE, env=this.env),
                    envir=this.env)
      this.env <- parent.env(environment())
    }
    clusterExport(cl,
                  ls(all.names=TRUE, env=globalenv()),
                  envir=globalenv())
    
    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        require(yy , character.only=TRUE)})
    })
    
    ## Run the lapply in parallel 
    return( parLapply( cl, ...) )
  }, finally = {        
    ## Stop the cluster
    stopCluster(cl)
  })
}
library(popbio)
library(data.table)
library(pracma)
setwd("E:/")

####----parameters----####
## sun-earth geometry
rd=180/pi # factor to convert radians into degrees
S_p0=1360 # solar constant (W/m2), p159 in Campbell and Norman 1998
sigma=5.67*10^-8 # W/m2/K4 stefan-boltzmann constant, p281 in Campbell and Norman 1998
alpha_S=0.90  # solar absorptivity, for lizard, Table 11.4 in Campbell and Norman 1998, from Gates 1980
alpha_L=0.97 # longwave absoptivity, for S. occidentalis, Barlett & Gates, 1967
F_d=0.5 # view factor for diffuse solar radiation, p181 in Campbell and Norman 1998
F_r=0.5  # view factor for reflected radiation, p181 in Campbell and Norman 1998
F_a=0.5  # view factor for atmospheric thermal radiation, p181 in Campbell and Norman 1998
F_g=0.5  # view factor for ground thermal radation, p181 in Campbell and Norman 1998
rho_S=0.48 # reflectance of desert sand, Tetzlaff, G. 1983, ground measure of Albedo of the Sahara. 
epsilon_s=0.97 # assumed surface emisivity, p163 in Campbell and Norman 1998
c_p=29.3 # specific heat of air, J/mol/K, p44 in Campbell and Norman 1998
#physiological parameters
svlj=30.8 #snout-vent length (mm) Sceloporus undulatus 30 days after releasing to the field; Parker and Andrews 2007 Oecologia
Mj= 3.55*10^(-5)*svlj^3 #juvenile body mass
Topt=32 #mean body temperatures regulated by gravid Sceloporus jarrovi, Beuchat 1986 (gravid female); 
#Selected body temp of juvenile Sceloporus jarrovi; Patterson et al. 2017 (juvenile)
Tvlow<-29.4 
Tvhigh<-36.3 #central 80% of field body temperature; Angilletta, 2001 Ecology
Cmax=exp(-2.23+0.05*svlj) 
Ej0= 0.1*0.25*39.3 #(KJ) initial fat energy storage of juvenile; Levy et al. 2016
Rfat=0.79 #Efficiency of converting assimilated energy to fat; Barker et al. 1998
vj= 10^(0.044+0.2*log10(Mj)) #...juvenile...
dj= 0.7*vj #...juvenile...
Ingest_max<-30.12*0.005*0.5*dj*0.76*3600/1000 #Hourly energy intake (without considering gut space) (KJ); Buckley 2008
D0 <- -1*(-0.024545)/0.001545 #Estimated developmental 0 is 15.88673 (jives with Christian et al 1986 showing developmental stasis at 15C)
CTmax<-41.5 # oC
svl_hib<-33 # mm svl_before_hibernation, Kearney 2011 Functional Ecology
Mj_hib<-3.55*10^(-5)*svl_hib^3
Energy_content<-Ej0*Mj_hib/Mj
Devsuccess_thre<-0.8
Energy_thre<-Energy_content #KJ
Offstress_thre<-0


####----physiological functions----####
cal.Devrate<-function(temps) { # function for translating temperature to Rate 
  sapply(1:8760,FUN=function(x){
    if (temps[x] > D0){
      Rate <- ((-0.024545+0.001545*temps[x])/24) #numbers are coefficients for regressing 1/duration ~ temperature, and divided by 24 to make hourly
      return (Rate)}
    else {return(0)}
  })
}

# convert hourly temps to hourly developmental time (return a vector)
cal.Devtime<-function(rates){
  sapply(1:8760,FUN=function(x){
    Dev_acc<-cumsum(rates[x:8760]) #calculate accumulated incubation from oviposition to the end of the year
    DT_x<-findInterval(1,Dev_acc) #calculate incubation time (complete 2/3 of development=1 of incubation)
    time<-as.numeric(ifelse(Dev_acc[DT_x+1] >=1 & DT_x < 100*24, DT_x, NA))
    return(time)
  })
}

# function for hatching success, from Rory script from levy 2015 (how to get that from the figure?)
cal.Devsuccess<-function(Devtime,temps){ # Devtime and temps are vectors
  sapply(1:8760,function(x){
    if(!is.na(Devtime[x])){
      temps_x<-temps[x:(x+Devtime[x])] #temps_x is the developmental temperatures experienced by egg laid at hour x
      Tmin<-mean(sapply(split(temps_x[1:(floor(length(temps_x)/24)*24)], rep(1:floor(length(temps_x)/24), each=24)), min))
      Tmax<-mean(sapply(split(temps_x[1:(floor(length(temps_x)/24)*24)], rep(1:floor(length(temps_x)/24), each=24)), max))
      if(Tmin < 25 & Tmax <= 42) {  #42 or 44? levy use 44
        logit_p <- -2.18733 + 0.14268*Tmin   
        survival <- exp(logit_p)/(1+exp(logit_p))
      } 
      if(Tmin >= 25 & Tmax <= 42) survival<-0.8
      if(Tmax > 42) survival<-0
    }
    else survival<-0
    return(survival)
  })
}

cal.Tb<-function(T1,T2,T3,T4){ # function for finding Tmin and Tmax
  Tmin<-min(T1,T2,T3,T4)
  Tmax<-max(T1,T2,T3,T4)
  if(Tmin<Topt && Tmax>Topt) return(Topt)
  if(Tmin>=Topt) return(Tmin)
  if(Tmax<=Topt) return(Tmax)
}


# calculate the sum of budget energy
cal.Energy<-function(Devtime,Budget){ 
  sapply(1:8760,function(x){
    if(!is.na(Devtime[x])) { 
      Budget_x<-ifelse(Budget[(x+Devtime[x]):8760]>0,
                       Budget[(x+Devtime[x]):8760]*Rfat,
                       Budget[(x+Devtime[x]):8760])
      Budget_sum_x<-ifelse(all((cumsum(Budget_x)+Ej0)>0), cumsum(Budget_x)[length(Budget_x)]+Ej0, 0)
    }
    else Budget_sum_x<-0
    return(Budget_sum_x)
  })
}



####----modules----####
##1. Read microclimate data and calculate operative temperatures
path<-"E:/Levy et al. 2016 (extracted)/"
read.micro.data<- function(loc_id,model_time,year_id) { #function for reading microclimate data
  range_hour<-((year_id-1)*8760+1):(year_id*8760)
  read.micro<-function(x){ #reading certain microclimate variable
    climate.nc=list.files(path=paste(path,x,sep=""),pattern=paste(x,"_",loc_id,"_MIC_CLIM_36_",model_time,"_*",sep=""))
    nc_geo=nc_open(paste(path,x,"/",climate.nc,sep=""))
    if(x=="SWDOWN") result<-ncvar_get(nc_geo,varid=x)[range_hour]
    if(x%in%c("Tsurface","WIND10")) result<-ncvar_get(nc_geo,varid=x)[range_hour,]
    if(x%in%c("Tair","Tsoil")) result<-ncvar_get(nc_geo,varid=x)[range_hour,,]
    nc_close(nc_geo)
    result
  }
  variable_names = c("SWDOWN","Tsurface","WIND10","Tair","Tsoil")
  microclimate<-lapply(variable_names,FUN=read.micro)
  x="SWDOWN"
  climate.nc=list.files(path=paste(path,x,sep=""),pattern=paste(x,"_",loc_id,"_MIC_CLIM_36_",model_time,"_*",sep=""))
  nc_geo=nc_open(paste(path,x,"/",climate.nc,sep=""))
  time<-ncvar_get(nc_geo,varid="time")[range_hour]
  years = time %/% 1000000
  months = (time %% 1000000) %/% 10000
  days = (time %% 10000) %/% 100
  hours =time %% 100
  date<-paste(months,days,years)
  date<-as.Date(date, format=c("%m %d %Y"))
  j_day<-as.numeric(format(date, "%j"))
  microclimate[["lat"]] = as.numeric(ncvar_get(nc_geo, varid="lat"))
  microclimate[["lon"]] = as.numeric(ncvar_get(nc_geo, varid="lon"))
  Tsurface_sun<-microclimate[[2]][,1]-273 #2 dimension: time, shade; shade=0%
  Tsurface_shade<-microclimate[[2]][,5]-273 #2 dimension: time, shade; shade=100%
  WIND<-microclimate[[3]][,1]#2 dimension: time, layers_mic; layer=0.03m
  Tair_sun<-microclimate[[4]][,1,1]-273#3 dimension: time, shade, layers_mic; shade=0%;layer=0.03m
  Tair_shade<-microclimate[[4]][,5,1]-273#3 dimension: time, shade, layers_mic; shade=100%;layer=0.03m
  Tair_120cm<-microclimate[[4]][,1,15]-273#3 dimension: time, shade, layers_mic; shade=0%;layer=1.20m
  microclimate[["Tsoil_0"]]<-microclimate[[5]][,1,2]-273#3 dimension: time, shade, layers_mic; shade=0%;layer=0.06m
  #microclimate[["Tsoil_50"]]<-microclimate[[5]][,3,2]-273#3 dimension: time, shade, layers_mic; shade=50%;layer=0.06m
  microclimate[["Tsoil_100"]]<-microclimate[[5]][,5,2]-273#3 dimension: time, shade, layers_mic; shade=100%;layer=0.06m
 
  ###---biophysical model---###
  ## sun-earth geometry
  RevAng = 0.21631 + 2 * atan (0.967 * tan (0.0086 * (-186 + j_day))) # Revolution angle in radians
  DecAng = asin (0.39795 * cos (RevAng))                         # Declination angle in radians    From Buckley 2008       
  Daylength<-2*acos(-tan(abs(as.vector(microclimate[["lat"]]))/rd)*tan(DecAng))/(360/rd/24) # Function 6.8 in Gates 1980, use abs value of lat
  f=(279.575+0.9856*j_day)/rd  # f in radians,function 11.4 in Campbell and Norman
  ET= (-104.7*sin (f)+596.2*sin (2*f)+4.3*sin (3*f)-12.7*sin (4*f)-429.3*cos (f)-2.0*cos (2*f)+19.3*cos (3*f))/3600   # (11.4) Equation of time, fuction 11.4 in Campbell and Norman 1998
  lon_convert<-as.vector(microclimate[["lon"]])+180 # convert longitude to 0~360
  LC<-ifelse(lon_convert%%15<7.5,(lon_convert%%15/15),((lon_convert%%15-15)/15))
  t_0=12-LC-ET # solar noon
  ZEN= sin(DecAng)*sin (as.vector(microclimate[["lat"]])/rd) + cos (DecAng)*cos (as.vector(microclimate[["lat"]])/rd)*cos (pi/12*(hours-t_0)) # zenith angle in radians
  psi=ZEN/rd # zenith angle in radians
  F_p=0.5*cos(psi) # view factor for solar beam, p181 in Campbell and Norman 1998
  # short wave radiation
  S_p=microclimate[[1]]/cos(psi) # direct irradience on a surface perpendicular to the beam
  S_d=0.3*(1-S_p/S_p0)* S_p0*cos (psi)   # sky diffuse radiation, function 11.11 and 11.13 in Campbell and Norman 1998
  S_t=S_p*cos (psi)+S_d # global irradiance (total irradiance, beam plus diffuse), function 11.8 and 11.9 in Campbell and Norman 1998 
  S_r=rho_S*S_t # reflected radiation, function 11.10 in Campbell and Norman 1998
  # long wave radiation
  epsilon_ac= 9.2*10^-6*(Tair_120cm+273)^2 # clear sky emissivity, function 10.11 in Campbell and Norman 1998
  L_a=epsilon_ac*sigma*(Tair_120cm+273)^4  # long wave flux densities from atmosphere, fuction 10.7 in Campbell and Norman 1998
  L_g_sun=epsilon_s*sigma*(Tsurface_sun+273)^4  # long wave flux densities from ground, function 10.7 in Campbell and Norman 1998
  L_g_shade=epsilon_s*sigma*(Tsurface_shade+273)^4  # long wave flux densities from ground, function 10.7 in Campbell and Norman 1998
  ## total radiation
  R_abs_sun= alpha_S*(F_p*S_p+ F_d*S_d + F_r*S_r)+alpha_L*(F_a*L_a+F_g*L_g_sun) # Absorbed radiation, function 11.14 in Campbell and Norman 1998
  R_abs_shade= alpha_L*(F_a*L_a+F_g*L_g_shade) # Absorbed radiation, function 11.14 in Campbell and Norman 1998
  #remove short wave
  ## conductance
  g_r_sun=4*sigma*(Tair_sun+273)^3/c_p
  g_r_shade=4*sigma*(Tair_shade+273)^3/c_p # radiative conductance, function 12.7 in Campbell and Norman 1998
  
  g_Ha_j=1.4*0.135*sqrt(WIND/svlj*1000) # boundary layer conductance for heat, mol/m2/s, p211 in Campbell and Norman 1998
  ## operative temperature
  microclimate[["Te_sun"]]<-Tair_sun+(R_abs_sun-epsilon_s*sigma*(Tair_120cm+273)^4)/(c_p*(g_r_sun+g_Ha_j)) # Operative temperature, function 12.19 in Campbell and Norman 1998
  microclimate[["Te_shade"]]<-Tair_shade+(R_abs_shade-epsilon_s*sigma*(Tair_120cm+273)^4)/(c_p*(g_r_shade+g_Ha_j)) # Operative temperature, function 12.19 in Campbell and Norman 1998
  microclimate[["Daylight"]]<-ifelse(hours>(t_0-0.5*Daylength)&hours<(t_0+0.5*Daylength),1,0)
  
  microclimate
}


##2. Calculate development success, offspring stress and growth energy

cal.growth<-function(loc_id,model_time,year_id) {

  growth<-read.micro.data(loc_id,model_time,year_id)[6:12]
  
  # Calculate body temperature, determine thermal opportunity and heat stress
  # Shuttle among Te_sun, Te_shade, Tsoil_sun, Tsoil_shade, aiming at Topt (as close as possible); Ma et al. 2018
  # discribed in the main text 2.3 Activity and body temperatures
  Tb<-mcmapply(cal.Tb,growth[["Te_sun"]],growth[["Te_shade"]],growth[["Tsoil_0"]],growth[["Tsoil_100"]])
  Tb[which(growth[["Daylight"]]==0)]<-growth[["Tsoil_100"]][which(growth[["Daylight"]]==0)]
  Fopp<-ifelse(Tb>Tvlow & Tb<Tvhigh & growth[["Daylight"]]==1, 1, 0)
  Dopp<-ifelse(Tb>Tvlow & Tb<Tvhigh, 1, 0) #digestion opportunity
  Dopp_daily<-.colSums(Dopp,24,length(Dopp)/24)
  Fopp_daily<-.colSums(Fopp,24,length(Fopp)/24)
 

  # Calculate gut content, foraging, digestion and metabolic rate without considering gestation
  # Calculate reproduction, embryonic development, and energy budget
  # function for calculating hourly "incubating rate" (percentage of incubation period completed per hour)
  # very important: incubation period==2/3 whole developmental time
  # Rory made function for "incubating rate" calculation by collecting S. undulatus incubation length data from 6 studies
  # estimate D0 from the coefficients (set y = 0 and solve for x)
  Digest_daymax<-0.115*Cmax*log(Dopp_daily+1) #???
  Gut_content<-c(0) # A vector of hourly gut content; empty gut at first hour
  Digest_daily<-c() # A vector of daily digestion
  Ingest<-c() # A vector of hourly ingestion
  MR<-c() # A vector of hourly metabolic rate
  
  for(i in 1:8760){
    Ingest[i]<-ifelse(Fopp[i]==1,min(Ingest_max,(Cmax-Gut_content[i])),0)
    Gut_content[i+1]<-Gut_content[i]+Ingest[i]
    if(i%%24==0){
      Digest_daily[i/24]<-min(Digest_daymax[i/24],Gut_content[i+1])
      Gut_content[i+1]<-Gut_content[i+1]-Digest_daily[i/24]
    }
    if(Dopp[i]==0) MR[i]<-exp(-10.0+0.51*log(Mj)+0.12*Tb[i])
    if(Dopp[i]!=0 & Ingest[i]==0 & Gut_content[i+1]==0) MR[i]<-exp(-10.0+0.51*log(Mj)+0.12*Tb[i])
    if(Dopp[i]!=0 & Ingest[i]==0 & Gut_content[i+1]!=0) MR[i]<-1.5*exp(-10.0+0.51*log(Mj)+0.12*Tb[i])
    if(Ingest[i]!=0) MR[i]<-2*exp(-10.0+0.51*log(Mj)+0.12*Tb[i])
  }
  
  Budget<-rep(Digest_daily/24,each=24)-MR # A vector of hourly energy budget
  
  
  Devrate<-as.numeric(cal.Devrate(growth[["Tsoil_0"]]))# Devrate is a vector of Rates of Tsoils
  growth[["Devtime_0"]]<-as.numeric(cal.Devtime(Devrate))
  growth[["Act_daily"]]<-rep(Fopp_daily, each=24)
  growth[["Act_daily"]][which(is.na(growth[['Devtime_0']]))]<-NA
  growth[["Devsuccess_0"]]<-as.numeric(cal.Devsuccess(growth[["Devtime_0"]],growth[["Tsoil_0"]]))
  growth[['Off_stress']]<-sapply(1:8760,FUN=function(x){
    if(!is.na(growth[['Devtime_0']][x])){
      i<-(x+growth[['Devtime_0']][x])
      length<-length(which((CTmax-Tb[i:(i+2*7*24)])<3))
      return(length)}
    else return(NA)
  })
  
  growth[["Energy_0"]]<-as.numeric(cal.Energy(growth[["Devtime_0"]],Budget))
  growth[["Energy_0"]][which(growth[['Devsuccess_0']]==0)]<-0
  
  ####calculate mismatch length
  #growth[['Devsuccess_mis']]<-length(which(growth[["Act_daily"]]>2&growth[["Devsuccess_0"]]<Devsuccess_thre))
  #growth[['Offstress_mis']]<-length(which(growth[["Act_daily"]]>2&growth[["Off_stress"]]>Offstress_thre))
  #growth[['Energy_mis']]<-length(which(growth[["Act_daily"]]>2&growth[["Energy_0"]]<Energy_thre))
  
  ####tanslate the mismatch data to Binomial Distribution
  growth[['Devsuccess_mis']]<-ifelse(growth[["Act_daily"]]>2&growth[["Devsuccess_0"]]<Devsuccess_thre,1,0)
  growth[['Devsuccess_mis']][which(is.na(growth[['Devsuccess_mis']]))]<-0
  growth[['Offstress_mis']]<-ifelse(growth[["Act_daily"]]>2&growth[["Off_stress"]]>Offstress_thre,1,0)
  growth[['Offstress_mis']][which(is.na(growth[['Offstress_mis']]))]<-0
  growth[['Energy_mis']]<-ifelse(growth[["Act_daily"]]>2&growth[["Energy_0"]]<Energy_thre,1,0)
  growth[['Energy_mis']][which(is.na(growth[['Energy_mis']]))]<-0
  growth[c(1:2, 10:15)]
  
  
  ####calculate ON length in 50% and 100% shade
  #for(shade in c(50,100)){
    #Devrate<-as.numeric(cal.Devrate(growth[[paste("Tsoil_",shade,sep="")]]))# Devrate is a vector of Rates of Tsoils
    #Devtime<-as.numeric(cal.Devtime(Devrate))
    #growth[[paste("Devsuccess_",shade,sep="")]]<-as.numeric(cal.Devsuccess(Devtime,growth[[paste("Tsoil_",shade,sep="")]]))
    #growth[[paste("Offstress_",shade,sep="")]]<-sapply(1:8760,FUN=function(x){
     # if(!is.na(Devtime[x])){
      #  i<-(x+Devtime[x])
       # length<-length(which((CTmax-Tb[i:(i+2*7*24)])<3))
      #  return(length)}
      #else return(NA)
    #})
    #growth[[paste("Energy_",shade,sep="")]]<-as.numeric(cal.Energy(Devtime,Budget))
    #growth[[paste("Energy_",shade,sep="")]][which(growth[[paste("Devsuccess_",shade,sep="")]]==0)]<-0
  #}
  #growth[c(1:2, 9:14)]
}  


cal.niche<-function(loc_id,model_time,year_id){
  niche<-cal.growth(loc_id,model_time,year_id)
  
  ####calculate mismatch length
  #niche[["Devsuccess_0"]]<-length(which(niche[["Devsuccess_0"]]>Devsuccess_thre))
  #niche[['Off_stress']]<-length(which(niche[["Off_stress"]]==Offstress_thre))
  #niche[["Energy_0"]]<-length(which(niche[["Energy_0"]]>Energy_thre))
  
  ###calculate ON length in 50% and 100% shade
  #for(shade in c(50,100)){
  #  niche[[paste("Devsuccess_",shade,sep="")]]<-length(which(niche[[paste("Devsuccess_",shade,sep="")]]>Devsuccess_thre))
    #niche[[paste("Offstress_",shade,sep="")]]<-length(which(niche[[paste("Offstress_",shade,sep="")]]==Offstress_thre))
    #niche[[paste("Energy_",shade,sep="")]]<-length(which(niche[[paste("Energy_",shade,sep="")]]>Energy_thre))
  #}

  ####tanslate the ON data to Binomial Distribution
  niche[["Devsuccess_0"]]<-ifelse(niche[["Devsuccess_0"]]>Devsuccess_thre,1,0)
  niche[['Off_stress']]<-ifelse(niche[["Off_stress"]]==Offstress_thre,1,0)
  niche[['Off_stress']][which(is.na(niche[['Off_stress']]))]<-0
  niche[["Energy_0"]]<-ifelse(niche[["Energy_0"]]>Energy_thre,1,0)
  names<-c('Devsuccess_0','Off_stress','Energy_0','Devsuccess_mis','Offstress_mis','Energy_mis')
  for(i in 1:length(names)){
    if(all(niche[[names[i]]]==0)||sum(niche[[names[i]]]==1)){ # only one 1 
      niche[[paste(names[i],'_peaknum', sep='')]]<-NA #number of peaks
      niche[[paste(names[i],'_locave', sep='')]]<-NA  #average peak location 
      niche[[paste(names[i],'_locmax', sep='')]]<-NA  #max peak location
      niche[[paste(names[i],'_locfirst', sep='')]]<-NA  #the first peak location
      niche[[paste(names[i],'_loclast', sep='')]]<-NA  #the last peak location
      niche[[paste(names[i],'_loc50', sep='')]]<-NA   #location at the 50% 1
    }
    else {
      x<-as.factor(niche[[names[i]]])
      y<-1:8760
      result<-cdplot(x~y) #calculate the density distribution
      prob<-1-result[["0"]](1:8760)
      if(all(diff(prob)<2.23e-16)){ #length(prob[which(prob>0.9999999)])>1 || 
        niche[[paste(names[i],'_peaknum', sep='')]]<-1  
        niche[[paste(names[i],'_locave', sep='')]]<-min(which(cumsum(niche[[names[i]]])==floor(sum(niche[[names[i]]])/2)))
        niche[[paste(names[i],'_locmax', sep='')]]<-min(which(cumsum(niche[[names[i]]])==floor(sum(niche[[names[i]]])/2)))
        niche[[paste(names[i],'_locfirst', sep='')]]<-min(which(cumsum(niche[[names[i]]])==floor(sum(niche[[names[i]]])/2)))
        niche[[paste(names[i],'_loclast', sep='')]]<-min(which(cumsum(niche[[names[i]]])==floor(sum(niche[[names[i]]])/2)))
        niche[[paste(names[i],'_loc50', sep='')]]<-min(which(cumsum(niche[[names[i]]])==floor(sum(niche[[names[i]]])/2)))
      }
      if(length(prob[which(prob==max(prob))])>1){
        niche[[paste(names[i],'_peaknum', sep='')]]<-1
        niche[[paste(names[i],'_locave', sep='')]]<-mean(which(prob==max(prob)))
        niche[[paste(names[i],'_locmax', sep='')]]<-mean(which(prob==max(prob)))
        niche[[paste(names[i],'_locfirst', sep='')]]<-mean(which(prob==max(prob)))
        niche[[paste(names[i],'_loclast', sep='')]]<-mean(which(prob==max(prob)))
        niche[[paste(names[i],'_loc50', sep='')]]<-mean(which(prob==max(prob)))
      }
      if(any(diff(prob)>2.23e-16)&length(prob[which(prob==max(prob))])<=1){
        peaks<-findpeaks(prob)
        niche[[paste(names[i],'_peaknum', sep='')]]<-nrow(peaks)
        niche[[paste(names[i],'_locave', sep='')]]<-mean(peaks[,2])
        niche[[paste(names[i],'_locmax', sep='')]]<-peaks[which(peaks==max(peaks[,1])),2] #which(g==max(g[,1]), arr.ind=TRUE) return row and col
        niche[[paste(names[i],'_locfir st', sep='')]]<-peaks[1,2]
        niche[[paste(names[i],'_loclast', sep='')]]<-peaks[nrow(peaks),2]
        niche[[paste(names[i],'_loc50', sep='')]]<-min(which(cumsum(niche[[names[i]]])==floor(sum(niche[[names[i]]])/2)))
      }
    }
  }
  niche[c(1:2, 9:44)]
}


cal.niche.years<-function(loc_id,model_time){
  growth.years<-mclapply.hack(1:20,cal.niche,loc_id=loc_id,model_time=model_time)
  growth.years<-matrix(unlist(growth.years,use.names=T),ncol=38,byrow=TRUE)
  growth.years<-colMeans(growth.years,na.rm=TRUE) #NA omit, mean of the real number,only 1 number will have a mean
  return(growth.years)
}


result.past<-matrix(NA,ncol=38,nrow=0)
result.future<-matrix(NA,ncol=38,nrow=0)


start.time <- Sys.time()
for(i in 1:11407){
  result.past<-rbind(result.past,cal.niche.years(loc_id=i,model_time="past"))
  result.future<-rbind(result.future,cal.niche.years(loc_id=i,model_time="future"))
  print(i)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

write.csv(result.past,"result.past.csv",row.names=F)
write.csv(result.future,"result.future.csv",row.names=F)



##3. Calculate Tsoil mean and IQR
path<-"E:/Levy et al. 2016 (extracted)/Tsoil/"
read.micro.data<-function(loc_id,model_time,year_id){ #function for reading microclimate data
  range_hour<-((year_id-1)*8760+1):(year_id*8760)
  climate.nc=list.files(path,pattern=paste("Tsoil_",loc_id,"_MIC_CLIM_36_",model_time,"_*",sep=""))
  nc_geo=nc_open(paste(path,climate.nc,sep=""))
  microclimate<-list()
  microclimate[["lat"]] = as.numeric(ncvar_get(nc_geo, varid="lat"))
  microclimate[["lon"]] = as.numeric(ncvar_get(nc_geo, varid="lon"))
  Tsoil_0<-as.numeric(ncvar_get(nc_geo, varid="Tsoil")[range_hour,1,2]-273)#3 dimension: time, shade, layers_mic; shade=0%;layer=0.06m
  nc_close(nc_geo)
  microclimate[["Tsoil"]]<-mean(Tsoil_0)
  microclimate[["IQR"]]<-IQR(Tsoil_0)
  microclimate
}

cal.niche.years<-function(loc_id,model_time){
  growth.years<-mclapply.hack(1:20,read.micro.data,loc_id=loc_id,model_time=model_time)
  growth.years<-matrix(unlist(growth.years,use.names=T),ncol=4,byrow=TRUE)
  growth.years<-colMeans(growth.years,na.rm=TRUE)
  return(growth.years)
}

result.past<-matrix(NA,ncol=4,nrow=0)
result.future<-matrix(NA,ncol=4,nrow=0)

start.time <- Sys.time()
for(i in 1:11407){
  result.past<-rbind(result.past,cal.niche.years(loc_id=i,model_time="past"))
  result.future<-rbind(result.future,cal.niche.years(loc_id=i,model_time="future"))
  print(i)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

write.csv(result.past,"Tsoil.past.csv",row.names=F)
write.csv(result.future,"Tsoil.future.csv",row.names=F)