%% Retrieving the state-space representation from Dynare
%  João Ricardo Costa Filho, 2022

%{
Correspondence: 

joao.filho@novasbe.pt
joaocostafilho.com

This is a free software: you can redistribute it and/or modify it under                                                                
the terms of the GNU General Public License as published by the Free                                                                   
Software Foundation, either version 3 of the License, or (at your option)                                                              
any later version.  See <http://www.gnu.org/licenses/> for more information.                                                           

* set the path to Dynare via Home -> Set Path -> Add Folder -> chose the matlab-subfolder of Dynare
* set the folder where the .mod-file is saved to yout Matlab-path
*/
%}

%% Housekeeping

close all;
clear;

%% solve the model

dynare NK;

clearvars -except oo_ M_ model;

%% State Space

% state space representation: 
% S_t = A * S_{t-1} + B * e_{t}, 
% X_t = C * S_{t-1} + D * e_{t};

%initiate the matrices

% baseline model
A = [];
B = [];
C = [];
D = [];

%state variables
state = M_.state_var';

%control variables
control = [1:1:size(oo_.dr.ghu,1)]';

%remove state variables from control variables vector
for j = 1:size(state,1)
    
    control( control == state(j) ) = [];

end

A = [oo_.dr.ghx( oo_.dr.inv_order_var( state ), : ) ];

B = [oo_.dr.ghu( oo_.dr.inv_order_var( state ), : ) ];

C = [oo_.dr.ghx( oo_.dr.inv_order_var( control ), : ) ];

D = [oo_.dr.ghu( oo_.dr.inv_order_var( control ), : ) ];

%exogenous variables names
exogenous = M_.exo_names;  

%endogenous variables names
endogenous = M_.endo_names;

%variables names
S_variables_names = M_.endo_names(state);
X_variables_names = M_.endo_names(control);

%% Simulations - Sequence of random shocks

e1 = [ 0 -3 0 0 0 2 0.5 0.5 0 0 0 0 ]'; %TFP shock;
e2 = [ 0 0 0 0 0 0 0 0 0 2 0.5 0  ]';   %Supply shock;
e3 = [ 0 -8 0 0 0 4 2 2 0 0 0 0  ]';   %demand shock;

horizon = size(e1,1)+1;

shocks_names = M_.exo_names;

shocks = zeros(M_.exo_nbr, horizon);

shocks( find( shocks_names == "eps_a"),2:horizon) = e1;
shocks( find( shocks_names == "eps_supply"),2:horizon) = e2;
shocks( find( shocks_names == "eps_demand"),2:horizon) = e3;

Ssim = zeros( size(state,1), horizon);
Xsim = zeros( size(control,1), horizon);

for j = 2:horizon
    
    Ssim(:,j) = A * Ssim(:,j-1) + B * shocks(:,j);
    Xsim(:,j) = C * Ssim(:,j-1) + D * shocks(:,j);
    
end

%% Graphs


% TFP shock;

figure('Position', get(0, 'Screensize'))
subplot(2,2,1)
plot(1:horizon,Ssim( find( S_variables_names == "y") ,:), 'Color',[0.6350 0.0780 0.1840], 'Linewidth', 2 ) 
axis("square")
xlim([1 horizon])
title( 'Output (%)' )
hold on

subplot(2,2,2)
plot(1:horizon, Xsim( find( X_variables_names == "n") ,:), 'Color',[0.6350 0.0780 0.1840], 'Linewidth', 2 ) 
axis("square")
xlim([1 horizon])
title('Hours of work (%)')
hold on

subplot(2,2,3)
plot(1:horizon, Xsim( find( X_variables_names == "i_ann") ,:), 'Color',[0.6350 0.0780 0.1840], 'Linewidth', 2 )
axis("square")
xlim([1 horizon])
title('Nominal interest rate (annualized, p.p.)')
hold on

subplot(2,2,4)
plot(1:horizon, Xsim( find( X_variables_names == "pi_ann") ,:), 'Color',[0.6350 0.0780 0.1840], 'Linewidth', 2 )
axis("square")
xlim([1 horizon])
title('Inflation rate (annualized, p.p.)')
hold off
exportgraphics(gcf,'figure2.jpg','Resolution',300)
