%%%%% Introdução aos modelos DSGE
%%%%% Produção Domiciliar e Ciclos de Negócio Reais (RBC)
%%%%% Modelo completo log-linear
%%%%% João Ricardo Costa Filho
%%%%% joaocostafilho.com

close all;

@#define antecipacao=1

%--------------------------------------------------------------------------------------------------------------------------------------
% 1. Definição das Variáveis
%--------------------------------------------------------------------------------------------------------------------------------------

%%% Variáveis endógenas %%%

var c  ${c}$    (long_name='Consumo')
    cm ${c_M}$  (long_name='Consumo no mercado')
    ch ${c_H}$  (long_name='Consumo no domicílio')
    y  ${y}$    (long_name='PIB')
    km ${k_M}$  (long_name='Estoque de capital do mercado')
    kh ${k_H}$  (long_name='Estoque de capital do domicílio')
    hm ${h_M}$  (long_name='Horas trabalhadas no mercado')
    hh ${h_H}$  (long_name='Horas trabalhadas no domicílio')
    w  ${w}$    (long_name='Salário real')
    l  ${\ell}$ (long_name='Lazer')
    r  ${r}$    (long_name='Taxa de juros real')
    xm ${x_M}$  (long_name='Investimento no mercado')
    xh ${x_H}$  (long_name='Investimento no domicílio')
    x  ${x}$    (long_name='Investimento agregado')
    zm ${x_M}$  (long_name='Produtividade no mercado')
    zh ${x_H}$  (long_name='Produtividade no domicílio')
    G  ${G}$    (long_name='Gastos do governo')
;

%%% Variáveis exógenas %%%
    
varexo em ${\varepsilon_M}$   (long_name='Choque de produtividade no mercado')
       eh ${\varepsilon_H}$   (long_name='Choque de produtividade no domicílio')
       eg ${\varepsilon_g}$   (long_name='Choque nos gastos do governo')
;

%--------------------------------------------------------------------------------------------------------------------------------------
% 2. Calibração
%--------------------------------------------------------------------------------------------------------------------------------------

parameters  beta     ${\beta}$     (long_name='fator de desconto')
            th       ${\tau_H}$    (long_name='alíquota de imposto sobre a renda trabalho')
            tk       ${\tau_K}$    (long_name='alíquota de imposto sobre a renda do capital')
            deltam   ${\delta_M}$  (long_name='depreciação do capital do mercado')
            deltah   ${\delta_H}$  (long_name='depreciação do capital do domicílio')
            alpha    ${\alpha}$    (long_name='participação do capital na função de produção do mercado')
            eta      ${\eta}$      (long_name='participação do capital na função de produção do domicílio')
            hmbar    ${\bar{H}_M}$ (long_name='horas trabalhadas no mercado em equilíbrio')
            hhbar    ${\bar{H}_H}$ (long_name='horas trabalhadas no domicílio em equilíbrio')
            rhom     ${\rho_M}$    (long_name='autocorrelação da produtividade no mercado')
            rhoh     ${\rho_H}$    (long_name='autocorrelação da produtividade no domicílio')
            rhog     ${\rho_G}$    (long_name='autocorrelação dos gastos do governo')
            e        ${e}$         (long_name='expoente da função CES do consumo total')
            gs       ${\tau_K}$    (long_name='percentual dos gastos do governo em relação ao PIB em equilíbrio')
;

beta   = 0.97;                                                                                                
th     = 0.34;                                                                                                     
tk     = 0.25;                                                                                                     
deltam = 0.05;                                                                                               
deltah = 0.05;                                                                                               
alpha  = 0.44;                                                                                                
eta    = 0.3245;                                                                                                 
hmbar  = 0.33;                                                                                                  
hhbar  = 0.25;                                                                                                  
rhom   = 0.9;                                                                                                   
rhoh   = 0.9;
rhog   = 0.8;                                                                                                                                                                  
e      = 2 / 3;
gs     = 0.2;                                                                              

%--------------------------------------------------------------------------------------------------------------------------------------
% 3. Modelo
%--------------------------------------------------------------------------------------------------------------------------------------

model(linear);

%% Equilíbrio estacionário

# rbar  = ( ( 1 / beta ) + deltam - 1 - tk * deltam ) / ( 1 - tk );
# kmbar = hmbar * ( rbar / alpha )^( 1 / ( alpha - 1 ) );                                                     
# ybar  = kmbar^alpha * hmbar^( 1 - alpha );
# xmbar = deltam * kmbar;
# wbar  = ( 1 - alpha ) * ( ybar / hmbar );
# lbar  = 1 - hmbar - hhbar;
# khbar = beta * hhbar * wbar * ( 1 - th) / ( 1 - beta * ( 1 - deltah ) * ( 1 - eta ) );
# chbar = khbar^eta * hhbar^( 1 - eta );
# xhbar = deltah * khbar;
# xbar  = xmbar + xhbar;
# Gbar  = gs * ybar;
# cmbar = ybar - xbar - Gbar;                                                                                         
#   phi = ( lbar * ( cmbar^(e-1) ) * wbar * ( 1 - th ) )^(-1); % variável auxiliar
# kappa = hhbar * ( lbar * ( chbar^e ) * ( 1 - eta ) )^(-1);   % variável auxiliar
#     a = phi / ( kappa + phi );
#  cbar = ( a * cmbar^e + ( 1 - a ) * chbar^e )^( 1 / e );                                                                          
#     b = phi / ( phi + a * cbar^(-e) );
#   chi = ( 1 - tk ) / ( rbar * ( 1 - tk ) + tk * deltam + ( 1 - deltam ) );  % variável auxiliar
#     g = ( a * ( e - 1 ) * cmbar^( e - 1 ) * ( 1 - deltah ) ) / ( a *cmbar^( e - 1 ) *( 1 - deltah ) + ( 1 - a ) * chbar^e * eta / khbar ); % variável auxiliar
#     q = ( ( 1 - a ) * e * ( 1 - eta ) * ( 1 / khbar ) * chbar^e ) / ( a * cmbar^( e - 1 )*( 1 - deltah ) + ( 1 - a ) * chbar^e * eta / khbar );  % variável auxiliar

%%%%%%%%%%%%% Famílias %%%%%%%%%%%%% 

[name = 'Oferta de trabalho no mercado']
e * c + ( 1 - e ) * cm - w = l;

[name = 'Oferta de trabalho no domicílio']
e * ch - e * c = hh - l;

[name = 'Equação de Euler no mercado']
( e - 1 ) * cm(+1) - e * c(+1) + chi * rbar * r(+1) = (e-1) * cm - e * c;

[name = 'Equação de Euler no domicílio']
e * c + ( 1 - e ) * cm = e * c(+1) - g * cm(+1) - e * q * ch(+1) + q * kh;

[name = 'Consumo']
c = ( cbar^-e ) * ( a * ( cmbar^e ) * cm + ( 1 - a ) * ( chbar^e ) * ch );                                                                               

[name = 'Consumo no domicílio']
ch = eta * kh(-1) + ( 1 - eta ) * zh + ( 1 - eta ) * hh;                                                                                             

[name = 'Lei de Movimento do capital no mercado']
xmbar * xm = kmbar * ( km - ( 1 - deltam ) * km(-1) );                                                                                         

[name = 'Lei de Movimento do capital no domicílio']
xhbar * xh = khbar * ( kh - ( 1 - deltah ) * kh(-1) );                                                                                          

%%%%%%%%%%%%% Empresas %%%%%%%%%%%%%

[name = 'Função de produção']
y = alpha * km(-1) + ( 1 - alpha ) * zm + ( 1 - alpha) * hm;                                                                                                  

[name = 'Demanda por capital']
y - km(-1) = r; 

[name = 'Demanda por trabalho']
y - hm = w;

%%%%%%%%%%%%% Governo %%%%%%%%%%%%% 
 
[name = 'Gastos do governo'] 
@#if antecipacao==0
    G = rhog * G(-1) + eg;
@#else   
    G = rhog * G(-1) + eg(-3);
@#endif

%%%%%%%%%%%%% Agregação %%%%%%%%%%%%% 

[name = 'Restrição de recursos']
ybar * y = cmbar * cm + xbar * x + Gbar * G;

[name = 'Investimento agregado']
xbar * x = xmbar * xm + xhbar * xh; 

[name = 'Restrição do tempo']
lbar * l = - ( hmbar * hm + hhbar * hh ); 

[name = 'Produtividade no mercado']
zm = rhom * zm(-1) + em;                                                                                                              

[name = 'Produtividade no domicílio']
zh = rhoh * zh(-1) + eh;                                                                                                      

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

var eg;
stderr 1;
end;

stoch_simul(order=1, irf=20);