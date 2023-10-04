%%%%% Introdução aos modelos DSGE
%%%%% Modelo de Ciclos de Negócio Reais (RBC)
%%%%% Modelo básico log-linearizado
%%%%% João Ricardo Costa Filho
%%%%% joaocostafilho.com

%--------------------------------------------------------------------------------------------------------------------------------------
% 1. Definição das Variáveis
%--------------------------------------------------------------------------------------------------------------------------------------

%%% Variáveis endógenas %%%

var c  ${c}$  (long_name='Consumo')
    h  ${h}$  (long_name='Horas trabalhadas')
    A  ${A}$  (long_name='Produtividade Total dos Fatores')
    k  ${k}$  (long_name='Estoque de capital')
;

%%% Variáveis exógenas %%%

varexo e ${\varepsilon_A}$   (long_name='Choque de produtividade')
;

%--------------------------------------------------------------------------------------------------------------------------------------
% 2. Calibração
%--------------------------------------------------------------------------------------------------------------------------------------

parameters phi     ${\phi}$ (long_name='Curvatura da função utilidade em relação às horas trabalhadas')
           psi     ${\psi}$ (long_name='peso da desutilidade do trabalho na função utilidade')
           sigma   ${\sigma}$ (long_name='curvatura da função utilidade')
           alpha   ${\alpha}$ (long_name='parâmetro da função de produção')
           beta    ${\beta}$  (long_name='fator de desconto')
           delta   ${\delta}$ (long_name='taxa de depreciação')
           rho     ${\rho}$   (long_name='autocorrelação da produtividade')

;


phi   = 1;
psi   = 1;
sigma = 2;
alpha = 0.44;
beta  = 0.97;
delta = 0.05;
rho   = 0.9;

%--------------------------------------------------------------------------------------------------------------------------------------
% 3. Modelo
%--------------------------------------------------------------------------------------------------------------------------------------

model(linear);

# Abar = 1;
# rbar = 1 / beta - 1 + delta;
# koh = ( rbar / alpha )^( 1 / ( alpha - 1 ) );
# coh = koh^alpha - delta * koh;
# hbar = ( ( 1 - alpha ) / psi * koh^alpha * coh^(-sigma) )^( 1 / ( phi + sigma - 1 ) );
# kbar = koh * hbar;
# ybar = Abar * kbar^alpha * hbar^(1-alpha);
# cbar = coh * hbar; 

%%%%%%%%%%%%% Famílias %%%%%%%%%%%%% 

[name = 'Oferta de Trabalho']
phi * h + sigma * c  = A + alpha * ( k(-1) - h);

[name = 'Equação de Euler']
c(+1) - c = ( 1 - beta * ( 1 - delta ) ) / sigma * ( A(+1) + ( 1 - alpha ) * ( h(+1) - k ) );

[name = 'Lei de Movimento do Capital']
k = ( 1 - delta ) * k(-1) + ybar / kbar * ( A + alpha * k(-1) + ( 1 - alpha ) * h ) - cbar / kbar * c;

%%%%%%%%%%%%% Agragação %%%%%%%%%%%%% 

[name = 'Produtividade']
A = rho * A(-1) + e;

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

var e;
stderr 1;
end;

stoch_simul(ar=1, order=1, irf=20);


