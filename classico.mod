% Adaptado de https://github.com/JohannesPfeifer/DSGE_mod/blob/master/Gali_2008/Gali_2008_chapter_2.mod



%--------------------------------------------------------------------------------------------------------------------------------------
% 1. Definição das Variáveis
%--------------------------------------------------------------------------------------------------------------------------------------


var C            ${C}$           (long_name='Consumo')
    W_real       ${\frac{W}{P}}$ (long_name='Salário Real')
    Pi           ${\Pi}$         (long_name='Taxa de inflaçºao')
    A            ${A}$           (long_name='Produvidade')
    N            ${N}$           (long_name='Horas trabalhadas')
    R            ${R^n}$         (long_name='Taxa de juros nominal') 
    r            ${r}$           (long_name='Real Interest Rate')
    Y            ${Y}$           (long_name='PIB') 
    m            ${\Delta M}$ (long_name='Taxa de crescimento da moeda')
    ;

varexo 
    ea           ${\varepsilon_A}$   (long_name='Choque de produtividade')
    em           ${\varepsilon_m}$   (long_name='Choque monetário')
    ;   

%--------------------------------------------------------------------------------------------------------------------------------------
% 2. Calibração
%--------------------------------------------------------------------------------------------------------------------------------------

parameters alpha   ${\alpha}$ (long_name='parâmetro da função de produção')
           beta    ${\beta}$ (long_name='fator de desconto')
           rho     ${\rho}$ (long_name='autocorrelação da produtividade')
           sigma   ${\sigma}$ (long_name='curvatura da função utilidade')
           phi     ${\phi}$ (long_name='elasticidade "Frisch"')
           phi_pi  ${\phi_{\pi}}$ (long_name='inflation feedback Taylor Rule')
           eta     ${\eta}$ (long_name='semi-elasticidade da demanda por moeda')
    ;

alpha  = 0.33; 
beta   = 0.99;
rho    = 0.9;
sigma  = 1;
phi    = 1;
phi_pi = 1.5;
eta    = 4;

%--------------------------------------------------------------------------------------------------------------------------------------
% 3. Modelo
%--------------------------------------------------------------------------------------------------------------------------------------

model;

%%%%%%%%%%%%% Famílias %%%%%%%%%%%%% 

[name = 'Oferta de Trabalho']
W_real = C ^ sigma * N ^ phi;

[name = 'Equação de Euler']
1/ R = beta * ( C(+1) / C ) ^ ( -sigma ) / Pi(+1);


%%%%%%%%%%%%% Empresas %%%%%%%%%%%%% 

[name = 'Função de Produção']
Y = A * N^( 1 - alpha );

[name = 'Demanda por Trabalho']
W_real = ( 1 - alpha ) * A * N ^( - alpha);


%%%%%%%%%%%%% Política monetária %%%%%%%%%%%%% 

R = 1 / beta * Pi ^ phi_pi + em;

%%%%%%%%%%%%% Agragação %%%%%%%%%%%%% 

[name = 'Taxa de juros real']
r = R / Pi(+1);

[name = 'Condição de Equilíbrio']
Y = C;

[name = 'Produtividade']
log(A)=rho*log(A(-1)) + ea;

[name = 'Condição de Equilíbrio']
m = 4 * ( log( Y ) - log( Y(-1) ) - eta * ( log(R) - log( R(-1) ) ) + log(Pi) );


end;

%--------------------------------------------------------------------------------------------------------------------------------------
% 4. Equilíbrio
%--------------------------------------------------------------------------------------------------------------------------------------

steady_state_model;

A  = 1;
R  = 1 / beta;
Pi = 1;
r  = R;
N = ( 1 - alpha )^( 1 / ( ( 1 - sigma ) * alpha + phi + sigma ) );
Y = A * N ^(1-alpha);
W_real = ( 1 - alpha ) * A * N ^ ( -alpha);
C = Y;
m = 0;
end;

steady;
check;
model_diagnostics;
model_info;

%--------------------------------------------------------------------------------------------------------------------------------------
% 5. Simulação
%--------------------------------------------------------------------------------------------------------------------------------------

shocks;
var ea; stderr 1;
var em; stderr 1;
end;

stoch_simul(ar=1, order=1, irf=20) Y C W_real Pi R r m;