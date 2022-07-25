% Adaptado de https://github.com/JohannesPfeifer/DSGE_mod/blob/master/Gali_2008/Gali_2008_chapter_3.mod


%--------------------------------------------------------------------------------------------------------------------------------------
% 0. Condução da política monetária
%--------------------------------------------------------------------------------------------------------------------------------------

% Aqui deve ser definido se o banco central utiliza uma regra de crescimento da quantidade de moeda (1) ou 
% uma regra de política monetária (0). 

@#define money_growth_rule=0


%--------------------------------------------------------------------------------------------------------------------------------------
% 1. Definição das Variáveis
%--------------------------------------------------------------------------------------------------------------------------------------


var pi           ${\pi}$ (long_name='inflation')
    y_gap        ${\tilde y}$ (long_name='output gap')
    y_nat        ${y^{nat}}$ (long_name='natural output')      //(in contrast to the textbook defined in deviation from steady state)
    y            ${y}$ (long_name='output')
    r_nat        ${r^{nat}}$ (long_name='natural interest rate')
    r_real       ${r^r}$ (long_name='//real interest rate')     
    i            ${i}$ (long_name='nominal interrst rate')
    n            ${n}$ (long_name='hours worked')
    m_real       ${m-p}$ (long_name='real money stock')
    m_growth_ann ${\Delta m}$ (long_name='money growth annualized')
    
    @#if money_growth_rule==0
        nu ${\nu}$ (long_name='AR(1) monetary policy shock process')    
    @#else
        money_growth  ${\Delta m_q}$ (long_name='money growth')
    @#endif
    
    a            ${a}$ (long_name='AR(1) technology shock process')
    r_real_ann   ${r^{r,ann}}$ (long_name='annualized real interest rate')
    i_ann        ${i^{ann}}$ (long_name='annualized nominal interest rate')
    r_nat_ann    ${r^{nat,ann}}$ (long_name='annualized natural interest rate')
    pi_ann       ${\pi^{ann}}$ (long_name='annualized inflation rate')
    ;     


varexo eps_a ${\varepsilon_a}$   (long_name='technology shock')

        eps_supply ${\varepsilon_s}$   (long_name='supply shock')
          
        eps_demand ${\varepsilon_d}$   (long_name='demand shock')
  
       @#if money_growth_rule==0
            eps_nu ${\varepsilon_\nu}$   (long_name='monetary policy shock')
       @#else   
            eps_m ${\varepsilon_\m}$   (long_name='money growth rate shock')
       @#endif
       ;

%--------------------------------------------------------------------------------------------------------------------------------------
% 2. Calibração
%--------------------------------------------------------------------------------------------------------------------------------------

parameters alppha ${\alppha}$ (long_name='capital share')
           betta ${\beta}$ (long_name='discount factor')
           rho_a ${\rho_a}$ (long_name='autocorrelation technology shock')
    
       @#if money_growth_rule==0
           rho_nu ${\rho_{\nu}}$ (long_name='autocorrelation monetary policy shock')
       @#else   
           rho_m ${\rho_{m}}$ (long_name='autocorrelation monetary growth rate shock')
       @#endif
           siggma ${\sigma}$ (long_name='log utility')
           phi ${\phi}$ (long_name='unitary Frisch elasticity')
           phi_pi ${\phi_{\pi}}$ (long_name='inflation feedback Taylor Rule')
           phi_y ${\phi_{y}}$ (long_name='output feedback Taylor Rule')
           eta ${\eta}$ (long_name='semi-elasticity of money demand')
           epsilon ${\epsilon}$ (long_name='demand elasticity')
           theta ${\theta}$ (long_name='Calvo parameter')
    ;  

siggma = 1;
phi=1;
phi_pi = 1.5;
phi_y  = .5/4;
theta=2/3;

@#if money_growth_rule==0
    rho_nu =0.5;
@#else   
    rho_m=0.5;
@#endif

rho_a  = 0.9;
betta = 0.99;
eta  =4;
alppha=1/3;
epsilon=6;

%--------------------------------------------------------------------------------------------------------------------------------------
% 3. Modelo
%--------------------------------------------------------------------------------------------------------------------------------------

model(linear);

//Composite parameters
#Omega=(1-alppha)/(1-alppha+alppha*epsilon);  
#psi_n_ya=(1+phi)/(siggma*(1-alppha)+phi+alppha); 
#lambda=(1-theta)*(1-betta*theta)/theta*Omega; 
#kappa=lambda*(siggma+(phi+alppha)/(1-alppha));  

//1. New Keynesian Phillips Curve
pi=betta*pi(+1)+kappa*y_gap + eps_supply ;

//2. Dynamic IS Curve 
y_gap=-1/siggma*(i-pi(+1)-r_nat)+y_gap(+1) + eps_demand;

//3. Interest Rate Rule
@#if money_growth_rule==0
i=phi_pi*pi+phi_y*y_gap+nu;
@#endif

//4. Definition natural rate of interest
r_nat=siggma*psi_n_ya*(a(+1)-a);

//5. Definition real interest rate
r_real=i-pi(+1);

//6. Definition natural output
y_nat=psi_n_ya*a;

//7. Definition output gap
y_gap=y-y_nat;

//8. Monetary policy shock
@#if money_growth_rule==0
    nu=rho_nu*nu(-1)+eps_nu;
@#endif

//9. TFP shock
a=rho_a*a(-1)+eps_a;

//10. Production function
y=a+(1-alppha)*n;

//11. Money growth (derived from
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi);

//12. Real money demand (eq. 4)
m_real=y-eta*i;
@#if money_growth_rule==1
//definition nominal money growth
money_growth=m_real-m_real(-1)+pi;
//exogenous process for money growth
money_growth=rho_m*(money_growth(-1))+eps_m;
@#endif

//13. Annualized nominal interest rate
i_ann=4*i;

//14. Annualized real interest rate
r_real_ann=4*r_real;

//15. Annualized natural interest rate
r_nat_ann=4*r_nat;

//16. Annualized inflation
pi_ann=4*pi;

end;

%--------------------------------------------------------------------------------------------------------------------------------------
% 4. Equilíbrio
%--------------------------------------------------------------------------------------------------------------------------------------

steady;
check;
model_diagnostics;
model_info;

%--------------------------------------------------------------------------------------------------------------------------------------
% 5. Simulação
%--------------------------------------------------------------------------------------------------------------------------------------

shocks;
    @#if money_growth_rule==0
        var eps_nu = 0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    @#else   
        var eps_m = 0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    @#endif
end;

@#if money_growth_rule==0
stoch_simul(order = 1, nograph);
@#else
stoch_simul(order = 1, nograph);
@#endif