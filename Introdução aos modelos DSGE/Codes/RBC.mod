%%%%% Introdução aos modelos DSGE
%%%%% Modelo de Ciclos de Negócio Reais (RBC)
%%%%% Modelo básico não-linear
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

parameters phi     ${\phi}$ (long_name='curvatura da função utilidade em relação às horas trabalhadas')
           psi     ${\psi}$ (long_name='peso da desutilidade do trabalho na função utilidade')
           sigma   ${\sigma}$ (long_name='curvatura da função utilidade')
           alpha   ${\alpha}$ (long_name='parâmetro da função de produção')
           beta    ${\beta}$  (long_name='fator de desconto')
           delta   ${\delta}$ (long_name='taxa de depreciação')
           rho     ${\rho}$   (long_name='autocorrelação da produtividade')

;


phi   = 1;
psi   = 2.29;
sigma = 2;
alpha = 0.44;
beta  = 0.97;
delta = 0.05;
rho   = 0.9;

%--------------------------------------------------------------------------------------------------------------------------------------
% 3. Modelo
%--------------------------------------------------------------------------------------------------------------------------------------

model;

# Abar = 1;

%%%%%%%%%%%%% Famílias %%%%%%%%%%%%% 

[name = 'Oferta de Trabalho']
psi * exp(h)^phi * exp(c)^sigma = ( 1 - alpha ) * exp(A) * ( exp(k(-1)) / exp(h) )^alpha;

[name = 'Equação de Euler']
exp(c)^(-sigma) = beta * ( exp(c(+1)) )^(-sigma) * ( 1 + alpha * exp(A(+1)) * ( exp(k) / exp(h(+1)) )^(alpha - 1) - delta );

[name = 'Lei de Movimento do Capital']
exp(k) = ( 1 - delta ) * exp(k(-1)) + exp(A) * exp(k(-1))^alpha * exp(h)^( 1 - alpha ) - exp(c);

%%%%%%%%%%%%% Agragação %%%%%%%%%%%%% 

[name = 'Produtividade']
A = ( 1 - rho ) * Abar + rho * A(-1) + e;

end;

%--------------------------------------------------------------------------------------------------------------------------------------
% 4. Equilíbrio
%--------------------------------------------------------------------------------------------------------------------------------------

initval;
A = 1;
h = 0.35;
c = 1.01;
k = 9.32;
end;

%--------------------------------------------------------------------------------------------------------------------------------------
% 5. Simulação
%--------------------------------------------------------------------------------------------------------------------------------------

shocks;

var e;
stderr 1;
end;

stoch_simul(ar=1, order=1, irf=20);


