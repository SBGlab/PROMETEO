function dxdt = TS_ODE(t,x,parstr)

% states of the model
mG = x(1); % mRNA Repressor G (concentration)
G = x(2);  % Repressor G (concentration)
mR = x(3); % mRNA Repressor R (concentration)
R = x(4);  % Repressor R (concentration

% parameters and inputs of the model
par = parstr.par;
IG = parstr.IG; % Inducer of G expression
IR = parstr.IR; % Inducer of R expression

kLeakG = par(1);    % Leakage Promoter G
kLeakR = par(2);    % Leakage Promoter R
thetaG = par(3);    % Theta Repressor G
thetaR = par(4);    % Theta Repressor R
thetaIG = par(5);   % Theta Inducer G
thetaIR = par(6);   % Theta Inducer R
nIG = par(7);       % n Inducer G
nIR = par(8);       % n Inducer R 
degmG = par(9);     % degradation rate constant mRNA G
degmR = par(10);    % degradation rate constant mRNA R
degG = par(11);     % degradation Repressor G
degR = par(12);     % degradation Repressor R

nG = par(13);      % n Binding Repressor G to Promoter R
nR = par(14);      % n Binding Repressor R to Promoter G
kG = par(15);      % transcription rate constant  promoter G
kR = par(16);      % transcription rate constant promoter R
kpG = par(17);     % translation rate constant RBS G
kpR = par(18);     % translation rate constant RBS R

% model ODE equations
denG = 1+(R/thetaR*(1/(1+(IG/thetaIG)^nIG)))^nR;  
denR = 1+(G/thetaG*(1/(1+(IR/thetaIR)^nIR)))^nG; 

dmGdt = kLeakG + kG/denG - degmG*mG;  % mass balance mRNA G
dGdt =  kpG*mG - degG*G;              % mass balance G
dmRdt = kLeakR + kR/denR - degmR*mR;  % mass balance mRNA G
dRdt =  kpR*mR - degR*R;              % mass balance R

dxdt = [dmGdt; dGdt; dmRdt; dRdt];