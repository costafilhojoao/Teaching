%Porto 2012, Replication of Greenwood and Wright 'Frontiers of Business Cycle Research - chapter 6', 1995
%code by João Ricardo M. G. costa Filho
close all;

%----------------------------------------------------------------------------------------------------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------------------------------------------------------------------------------------------------
var c cm ch y km hm kh hh w l r xm xh x zm zh;
    
varexo em eh;

parameters 
beta th tk deltam deltah theta eta hmbar hhbar rhom rhoh sigmam sigmah e
rbar kmbar ybar xmbar wbar lbar khbar xhbar xbar cmbar chbar phi kappa a cbar b chi g q;

%----------------------------------------------------------------------------------------------------------------------------------------------------------
% 2. calibration/Steady State
%----------------------------------------------------------------------------------------------------------------------------------------------------------
beta=0.9898;                                                                                                 %discount factor
th=0.35;                                                                                                     %
tk=0.70;                                                                                                     %tax rate on capital income
deltam=0.0235;                                                                                                %market capital depreciation rate
deltah=0.0235;                                                                                                %home capital depreciation rate
theta=0.2944;                                                                                                % 
eta=0.3245;                                                                                                  %exponent of home consumption function
hmbar=0.33;                                                                                                  %market hours
hhbar=0.25;                                                                                                  %home hours
rhom=0.95;                                                                                                   %ar coefficient of market technology
rhoh=0.95;                                                                                                   %ar coefficient of home technology
sigmam=0.007;                                                                                                %market innovation standard deviation
sigmah=0.007;                                                                                                %home innovation standard deviation
e=2/3; gamma=2/3;                                                                                            %model 2
%e=0.4; gamma=0;                                                                                             %model 3
rbar=((1/beta)+deltam-1-tk*deltam)/(1-tk);                                                                   %equilibrium interest rate
kmbar=hmbar*(theta/rbar)^(1/(1-theta));                                                                      %equilibrium market capital stock
ybar=kmbar*(rbar/theta);                                                                                     %equilibrium output
xmbar=deltam*kmbar;                                                                                          %equilibrium market investment
wbar=(ybar/hmbar)*(1-theta);                                                                                 %equilibrium wage
%wbar=(ybar/hmbar);
lbar=1-hmbar-hmbar;                                                                                          %equilibrium leisure
khbar=(beta)/(1-beta*(1-deltah))*hhbar*wbar*(1-th);                                                          %equilibrium home capital stock
xhbar=deltah*khbar;                                                                                          %equilibrium home investment
xbar=xmbar+xhbar;                                                                                            %equilibrium total investment
cmbar=ybar-xbar;                                                                                             %equilibrium market consumption
chbar=(khbar^eta)*(hhbar^(1-eta));                                                                           %equilibrium home-produced-goods consumption
phi=(lbar*(cmbar^(e-1))*wbar*(1-th))^(-1);
kappa=hhbar*(lbar*(chbar^e)*(1-eta))^(-1);
a=phi/(kappa+phi);                                                                                           %market consumption weight on total consumption
cbar=(a*cmbar^e+(1-a)*chbar^e)^(1/e);                                                                        %equilibrium total consumption
b=phi/(phi+a*cbar^(-e));                                                                           
chi=(1-tk)/(rbar*(1-tk)+tk*deltam+(1-deltam));
g=(a*(e-1)*cmbar^(e-1)*(1-deltah))/(a*cmbar^(e-1)*(1-deltah)+(1-a)*chbar^e*eta/khbar);
q=((1-a)*e*(1-eta)*(1/khbar)*chbar^e)/(a*cmbar^(e-1)*(1-deltah)+(1-a)*chbar^e*eta/khbar); 
%----------------------------------------------------------------------------------------------------------------------------------------------------------
% 3. Model 
%----------------------------------------------------------------------------------------------------------------------------------------------------------
model(linear); 

c=(cbar^-e)*(a*(cmbar^e)*cm+(1-a)*(chbar^e)*ch);                                                                               
y=theta*km(-1)+(1-theta)*zm+(1-theta)*hm;                                                                                                  
ch=eta*kh(-1)+(1-eta)*zh+(1-eta)*hh;                                                                                               
zm=rhom*zm(-1)-em;                                                                                                              
zh=rhoh*zh(-1)+eh;                                                                                                             
ybar*y=cmbar*cm+xbar*x;                                                                                                          
xmbar*xm=kmbar*(km-(1-deltam)*km(-1));                                                                                         
xhbar*xh=khbar*(kh-(1-deltah)*kh(-1));                                                                                          
xbar*x=xmbar*xm+xhbar*xh;                                                                                                       
e*c+(1-e)*cm-w=l;
e*ch-e*c=hh-l;                                              
(e-1)*cm(+1)-e*c(+1)+chi*rbar*r(+1)=(e-1)*cm-e*c;
e*c+(1-e)*cm=e*c(+1)-g*cm(+1)-e*q*ch(+1)+q*kh;
y-km(-1)=r;                   
y-hm=w;                   
lbar*l=-(hmbar*hm+hhbar*hh);                        

end;
%----------------------------------------------------------------------------------------------------------------------------------------------------------
shocks;
var em=sigmam;
%var eh=sigmah;
corr em, eh=gamma;
end;

steady;

stoch_simul(irf=100, nograph)  c cm ch y km hm kh hh w l r xm xh x;