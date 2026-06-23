close all
clear all
clc

% PROMOTER PARAMETERS
% Promoter pBetl
kLeak_pBetl =   0.0141;     % Leakage of the pBetl  (nM/min)
k_pBtel = 14.09;            % transcription rate pBetl  (nM/min)

% Promoter pTtg
kLeak_pTtg = 0.0025;        % Leakage of pTtg
k_pPtg = 2.49;              % transcription rate pPtg


% Promoter pTet
kLeak_pTet =  0.0130;       % Leakage of pTet
k_pTet = 13.01;             % transcription rate pTet 


% REPRESSOR PARAMETERS
% Repressor TtgR
theta_TtgR = 400;           % Theta-Threshold parameter TtgR (in arbitrary units)   


degm_TtgR = 1.386e-1;          % degradation mRNA TtgR 
deg_TtgR = 1.65e-2;            % degradation protein TtgR 
n_TtgR = 3;                    % Cooperativity of the binding TtgR
 
% Repressor BetlR
theta_BetlR = 178;         % Theta-Threshold parameter BetlR (in arbitrary units) 
degm_BetlR = 1.386e-1;          % degradation mRNA BetlR
deg_BetlR = 1.65e-2;            % degradation protein BetlR
n_BetlR = 2.5;                  % Cooperativity of the binding BetlR

% Repressor TetR
theta_TetR = 8;                % Theta-Threshold parameter TetR(in arbitrary units)   
degm_TetR = 1.386e-1;          % degradation mRNA TetR
deg_TetR = 1.65e-2;            % degradation protein TetR
n_TetR = 2.15;                 % Cooperativity of the binding TetR

% INDUCER PARAMETERS
% Inducer Cho
thetaI_Cho = 2.926e-1;          % Theta-Threshold parameter for Cho
nI_Cho = 2.7;                   % Cooperativity Cho 

% Inducer Nar
thetaI_Nar =  278.0613;         % Theta-Threshold parameter for  Nar
nI_Nar = 1.9;                   % Cooperativity Nar

% Inducer aTc
thetaI_aTc = 15;                % Theta-Threshold parameter for  aTc
nI_aTc = 3.8; 

% RBS PARAMETERS
kp_B0034 = 0.39;              % translation rate B0034
kp_RBS_st = 0.09;               % translation rate RBS_st
kp_BCD8 = 0.56;                 % translation rate BCD8


color1 =[51 75 233]/255;
color2 =[235 120 44]/255;

% CONFIG pTS1 (M-M) 

% PROMOTER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Promoter G: pBetl
kLeakG =  kLeak_pBetl;         % Leakage of the pBetl
kG = k_pBtel;               % transcription rate pBetl  

% Promoter R: pTtg
kLeakR = kLeak_pTtg;          % Leakage of pTtg
kR = k_pPtg;                % transcription rate pPtg
 
 % REPRESSOR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Repressor G
thetaG = theta_TtgR;         % Theta-Threshold parameter TtgR (in arbitrary units)   
degmG = degm_TtgR;          % degradation mRNA TtgR 
degG = deg_TtgR;            % degradation protein TtgR 
nG = n_TtgR;                    % Cooperativity of the binding TtgR
 
% Repressor R
thetaR = theta_BetlR;         % Theta-Threshold parameter BetlR (in arbitrary units) 
degmR = degm_BetlR;          % degradation mRNA BetlR
degR = deg_BetlR;            % degradation protein BetlR
nR = n_BetlR;                  % Cooperativity of the binding BetlR
 
% INDUCER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inducer G
thetaIG = thetaI_Cho;          % Theta-Threshold parameter for Cho
nIG = nI_Cho;                   % Cooperativity Cho 

% Inducer R
thetaIR =  thetaI_Nar;         % Theta-Threshold parameter for  Nar
nIR = nI_Nar;                   % Cooperativity Nar

% RBS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters associated to the RBS
kpG = kp_B0034;            % translation rate R (B0034)
kpR = kp_B0034;            % translation rate R (B0034)


par = [kLeakG, kLeakR, thetaG, thetaR, thetaIG, thetaIR, nIG, nIR, degmG, degmR, degG, degR, nG, nR, kG, kR, kpG, kpR];

% EXPERIMENTAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TS21_data = [0	0.038	0.000	0.045	0	   0.000	0.012	0.020	
	        60	0.046	0.040	0.049	60	   0.128	0.162	0.178	
	        120	0.059	0.062	0.065	120	   0.435	0.413	0.463	
	        200	0.075	0.075	0.082	200	   0.509	0.538	0.598	
	        340	0.265	0.197	0.254	340	   0.673	0.685	0.701	
	        600	0.757	0.681	0.796	600	   0.926	0.871	0.956	
	        800	1.000	0.878	0.968	800	   0.988	1.046	0.984	
	       1000	0.867	1.788	1.000	1000    1.090	1.000	1.000];
 
Nar = TS21_data(:,1); 
TS21_F_Mean = mean(TS21_data(:,2:4),2);
TS21_F_SD =std(TS21_data(:,2:4)')';
TS21_B_Mean = mean(TS21_data(:,6:8),2);
TS21_B_SD =std(TS21_data(:,6:8)')';

% DYNAMICS 
IG =  0;    % Induces G
IRv_1 = Nar;  % Induces R

parstr.par = par;
parstr.IG = IG;

% TIME SPAN
tspan = 0:0.1:1400;

% IC
% forward dynamics
x0_f = [0; 0; 0; 0];

% backward dynamics
x0_b = [0; 0; 10; 10]; 

% Simulation forward
for ii=1:1:8
parstr.IR = IRv_1(ii);
[t,x_f{ii}] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
 x_f_end{ii} = x_f{ii}(end,:);
end

% Simulation backward
for ii=1:1:8
parstr.IR = IRv_1(ii);
[t,x_b{ii}] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
x_b_end{ii} = x_b{ii}(end,:);
end

tt = 0:30:1380;

[idx,loc]=ismember(tt,t);

for ii=1:1:8
eval(sprintf('GFP_%d_Sim_f = x_f{%d}(loc,2);',ii,ii))
eval(sprintf('GFP_%d_Sim_b = x_b{%d}(loc,2);',ii,ii))
eval(sprintf('RFP_%d_Sim_f = x_f{%d}(loc,4);',ii,ii))
eval(sprintf('RFP_%d_Sim_b = x_b{%d}(loc,4);',ii,ii))
end


% Dose Response
for ii=1:1:8
 PG_ss_f_1(ii)= x_f_end{ii}(2);
 PG_ss_b_1(ii)= x_b_end{ii}(2);
 PR_ss_f_1(ii)= x_f_end{ii}(4);
 PR_ss_b_1(ii)= x_b_end{ii}(4);
end

PG_ss_f_sim_1 = PG_ss_f_1/max(PG_ss_b_1);
PG_ss_b_sim_1 = PG_ss_b_1/max(PG_ss_b_1);
PR_ss_f_sim_1 = PR_ss_f_1/max(PR_ss_b_1);
PR_ss_b_sim_1 = PR_ss_b_1/max(PR_ss_b_1);

PR_ss_f_data_1 = TS21_F_Mean;
PR_ss_b_data_1 = TS21_B_Mean;

Nar = TS21_data(:,1); 
TS21_F_Mean = mean(TS21_data(:,2:4),2);
TS21_F_SD =std(TS21_data(:,2:4)')';
TS21_B_Mean = mean(TS21_data(:,6:8),2);
TS21_B_SD =std(TS21_data(:,6:8)')';


%%% CONFIG pTS2 S-M

% PROMOTER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Promoter G: pBetl
kLeakG =   kLeak_pBetl;        % Leakage of the pBetl
kG = k_pBtel;               % transcription rate pBetl  
 
% Promoter R: pTet
kLeakR = kLeak_pTet;           % Leakage of pTet
kR = k_pTet;                % transcription rate pTet

% REPRESSOR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repressor G
thetaG = theta_TetR;                % Theta-Threshold parameter TetR(in arbitrary units)   
degmG = degm_TetR;          % degradation mRNA TetR
degG = deg_TetR;            % degradation protein TetR
nG = n_TetR;                 % Cooperativity of the binding TetR

% Repressor R
thetaR = 550;      
degmR = degm_BetlR;          % degradation mRNA BetlR
degR = deg_BetlR;            % degradation protein BetlR 
nR = n_BetlR;                  % Cooperativity of the binding BetlR

% INDUCER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inducer G
thetaIG = thetaI_Cho;        % Theta-Threshold parameter for Cho
nIG = nI_Cho;                 % Cooperativity of the Inducer Cho 

%Inducer R
thetaIR = thetaI_aTc;              % Theta-Threshold parameter for  aTc
nIR = nI_aTc; 

% RBS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters associated to the RBS
kpG = kp_B0034;          % translation rate G (BC0034)
kpR = kp_BCD8;            % translation rate R (BCD8)



par = [kLeakG, kLeakR, thetaG, thetaR, thetaIG, thetaIR, nIG, nIR, degmG, degmR, degG, degR, nG, nR, kG, kR, kpG, kpR];

TS9_data = [	0	0.054	0.035	0.036	0	0.153	0.082	0.023	
                20	0.168	0.106	0.164	20	0.022	0.017	0.065	
                40	0.195	0.102	0.165	40	0.831	0.817	0.801	
                48	0.754	0.722	0.823	48	0.821	0.871	0.897	
                62	0.803	0.900	0.876	62	0.897	0.953	0.985	
                72	0.943	0.941	0.987	72	0.914	0.968	0.974	
                88	0.925	0.993	0.968	88	0.995	1.000	1.000	
                100	1.000	1.000	1.000	100	1.000	1.000	1.00];
            
           
             
aTc_2 = TS9_data(:,1); 
TS9_F_Mean = mean(TS9_data(:,2:4),2);
TS9_F_SD =std(TS9_data(:,2:4)')';
TS9_B_Mean = mean(TS9_data(:,6:8),2);
TS9_B_SD =std(TS9_data(:,6:8)')';

% DYNAMICS 
IG =  0;  % induces G
IRv_2 = aTc_2;%  nominal: hysteresis range from 10 to 22


parstr.par = par;
parstr.IG = IG;


% TIME SPAN
tspan = 0:0.1:1400;

% forward dynamics
x0_f = [0; 0; 0; 0];

% backward dynamics
x0_b = [0; 0; 80; 80]; 

% simulation forward
for ii=1:1:8
parstr.IR = IRv_2(ii);
[t,x_f{ii}] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
 x_f_end{ii} = x_f{ii}(end,:);
end

% simulation backward

for ii=1:1:8
parstr.IR = IRv_2(ii);
[t,x_b{ii}] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
x_b_end{ii} = x_b{ii}(end,:);
end


tt = 0:30:1380;

[idx,loc]=ismember(tt,t);


for ii=1:1:8
eval(sprintf('GFP_%d_Sim_f = x_f{%d}(loc,2);',ii,ii))
eval(sprintf('GFP_%d_Sim_b = x_b{%d}(loc,2);',ii,ii))
eval(sprintf('RFP_%d_Sim_f = x_f{%d}(loc,4);',ii,ii))
eval(sprintf('RFP_%d_Sim_b = x_b{%d}(loc,4);',ii,ii))
end


% dose response
for ii=1:1:8
 PG_ss_f_2(ii)= x_f_end{ii}(2);
 PG_ss_b_2(ii)= x_b_end{ii}(2);
 PR_ss_f_2(ii)= x_f_end{ii}(4);
 PR_ss_b_2(ii)= x_b_end{ii}(4);
end

PG_ss_f_sim_2 = PG_ss_f_2/max(PG_ss_b_2);
PG_ss_b_sim_2 = PG_ss_b_2/max(PG_ss_b_2);
PR_ss_f_sim_2 = PR_ss_f_2/max(PR_ss_b_2);
PR_ss_b_sim_2 = PR_ss_b_2/max(PR_ss_b_2);


PR_ss_f_data_2 = TS9_F_Mean;
PR_ss_b_data_2 = TS9_B_Mean;



%%% CONFIG pTS3  (W-S)
% PROMOTER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Promoter G: pTtg
kLeakG =  kLeak_pTtg;         % Leakage of the pTtg
kG = k_pPtg;                % transcription rate pTtg  
 
% Promoter R: pTet
kLeakR = kLeak_pTet;  
kR = k_pTet;                 % transcription rate pTet


% REPRESSOR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repressor G
thetaG = theta_TetR;                % Theta-Threshold parameter TetR(in arbitrary units) 
degmG = degm_TetR;          % degradation mRNA TetR
degG = deg_TetR;            % degradation protein TetR
nG = n_TetR;                 % Cooperativity of the binding TetR

% Repressor R
thetaR = theta_TtgR ;        % Theta-Threshold parameter TtgR (in arbitrary units)   
degmR = degm_TtgR;          % degradation mRNA TtgR
degR = deg_TtgR;            % degradation protein TtgR
nR = n_TtgR;                    % Cooperativity of the binding TtgR


% % INDUCER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inducer G
thetaIG = thetaI_Nar;        % Theta-Threshold parameter for Nar
nIG = nI_Nar;                % Cooperativity of the Inducer Nar

% 
% % Inducer R
thetaIR = thetaI_aTc;           % Theta-Threshold parameter for  aTc
nIR = nI_aTc; 
% 

% RBS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters associated to the RBS
%kpG = 0.1337;           % translation rate R (RBS_st)

kpG = kp_RBS_st;           % translation rate R (RBS_st)
kpR = kp_BCD8;            % translation rate R (BCD8)



par = [kLeakG, kLeakR, thetaG, thetaR, thetaIG, thetaIR, nIG, nIR, degmG, degmR, degG, degR, nG, nR, kG, kR, kpG, kpR];

TS31_data = [	0	0.022	0.049	0.025	0	0.897	0.865	0.869	
                10	0.121	0.165	0.128	10	0.897	0.864	0.812	
                15	0.421	0.465	0.461	15	0.812	0.845	0.874	
                20	0.828	0.804	0.813	20	0.813	0.894	0.845	
                40	0.833	0.894	0.842	40	0.821	0.902	0.856	
	            48	0.829	0.860	0.876	48	0.834	0.945	0.902	
	            62	0.899	0.937	0.901	62	0.897	0.945	0.918	
	            72	0.925	0.960	0.946	72	0.923	0.968	0.968	
	            88	0.940	0.984	1.000	88	0.964	0.987	0.987
	            100	0.987	0.993	1.000	100	0.972	0.948	1.000];
aTc_3 = TS31_data(:,1); 
TS31_F_Mean = mean(TS31_data(:,2:4),2);
TS31_F_SD =std(TS31_data(:,2:4)')';
TS31_B_Mean = mean(TS31_data(:,6:8),2);
TS31_B_SD =std(TS31_data(:,6:8)')';



% DYNAMICS 
IG =  0;  % induces G
IRv = aTc_3;%  


parstr.par = par;
parstr.IG = IG;


% TIME SPAN
tspan = 0:0.1:1400;

% forward dynamics
x0_f = [0; 0; 0; 0];

% backward dynamics
x0_b = [0; 0; 80; 80]; 


% simulation forward
for ii=1:1:10
parstr.IR = IRv(ii);
[t,x_f{ii}] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
 x_f_end{ii} = x_f{ii}(end,:);
end

% simulation backward

for ii=1:1:10
parstr.IR = IRv(ii);
[t,x_b{ii}] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
x_b_end{ii} = x_b{ii}(end,:);
end


tt = 0:30:1380;

[idx,loc]=ismember(tt,t);


for ii=1:1:10
eval(sprintf('GFP_%d_Sim_f = x_f{%d}(loc,2);',ii,ii))
eval(sprintf('GFP_%d_Sim_b = x_b{%d}(loc,2);',ii,ii))
eval(sprintf('RFP_%d_Sim_f = x_f{%d}(loc,4);',ii,ii))
eval(sprintf('RFP_%d_Sim_b = x_b{%d}(loc,4);',ii,ii))
end


% DOSE RESPONSE
for ii=1:1:10
 PG_ss_f(ii)= x_f_end{ii}(2);
 PG_ss_b(ii)= x_b_end{ii}(2);
 PR_ss_f(ii)= x_f_end{ii}(4);
 PR_ss_b(ii)= x_b_end{ii}(4);
end

PG_ss_f_sim = PG_ss_f/max(PG_ss_b);
PG_ss_b_sim = PG_ss_b/max(PG_ss_b);
PR_ss_f_sim = PR_ss_f/max(PR_ss_b);
PR_ss_b_sim = PR_ss_b/max(PR_ss_b);


PR_ss_f_data = TS31_F_Mean;
PR_ss_b_data = TS31_B_Mean;


% FIGURES 
%--------------------------------------------------------------------------

figure(1)
subplot(2,3,1)
plot(IRv_1, PR_ss_f_sim_1, 'color', color1, 'linewidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
hold on
plot(IRv_1, PR_ss_b_sim_1, 'color', color2, 'linewidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('Nar (\muM)')
ylabel('RFP_{norm}')
axis([0 1000 0 2])
title('pTS1 (M-M)')

subplot(2,3,4)
errorbar(Nar, TS21_F_Mean, TS21_F_SD, 'color', color1, 'linewidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
hold on
errorbar(Nar, TS21_B_Mean, TS21_B_SD, 'color', color2, 'linewidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('Nar (\muM)')
ylabel('RFP_{norm}')
axis([0 1000 0 2])

subplot(2,3,2)
plot(IRv_2, PR_ss_f_sim_2, 'color', color1, 'LineWidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
hold on
plot(IRv_2, PR_ss_b_sim_2, 'color', color2,'LineWidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('aTc (ng/mL)')
ylabel('RFP_{norm}')
axis([0 100 0 1.5])
title('pTS2 S-M')

subplot(2,3,5)
errorbar(aTc_2, TS9_F_Mean, TS9_F_SD, 'color', color1,'LineWidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
hold on
errorbar(aTc_2, TS9_B_Mean, TS9_B_SD, 'color', color2,'LineWidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('aTc (ng/mL)')
ylabel('RFP_{norm}')
axis([0 100 0 1.5])

subplot(2,3,3)
plot(IRv', PR_ss_f_sim, 'color', color1, 'LineWidth', 2)
hold on
plot(IRv', PR_ss_b_sim, 'color', color2, 'LineWidth', 2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('aTc (ng/mL)')
ylabel('RFP_{norm}')
legend('Forward (sim)','Backward (sim)')
title('pTS3 W-S')
axis([0 100 0 1.5])

subplot(2,3,6)
errorbar(aTc_3, TS31_F_Mean, TS31_F_SD, 'color', color1, 'LineWidth', 2)
hold on
errorbar(aTc_3, TS31_B_Mean, TS31_B_SD, 'color', color2, 'LineWidth', 2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('aTc (ng/mL)')
ylabel('RFP_{norm}')
legend('Forward (exp)','Backward (exp)')
axis([0 100 0 1.5])

% Define OBJECTIVE FUNCTION
%-----------------------------------------------------------------------

figure(1)
subplot(2,3,1)
plot(IRv_1, PR_ss_f_sim_1, 'color', color1, 'linewidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
hold on
plot(IRv_1, PR_ss_b_sim_1, 'color', color2, 'linewidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('Nar (\muM)')
ylabel('RFP_{norm}')
axis([0 1000 0 2])
title('pTS1(M-M)')



subplot(2,3,2)
plot(IRv_2, PR_ss_f_sim_2, 'color', color1, 'LineWidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
hold on
plot(IRv_2, PR_ss_b_sim_2, 'color', color2,'LineWidth',2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('aTc (ng/mL)')
ylabel('RFP_{norm}')
axis([0 100 0 1.5])
title('pTS2 (S-M)')


subplot(2,3,3)
plot(IRv', PR_ss_f_sim, 'color', color1, 'LineWidth', 2)
hold on
plot(IRv', PR_ss_b_sim, 'color', color2, 'LineWidth', 2)
set(gca,'FontSize',14, 'FontName', 'Helvetica')
xlabel('aTc (ng/mL)')
ylabel('RFP_{norm}')
legend('Forward (sim)','Backward (sim)')
title('pTS3 (W-S)')
axis([0 100 0 1.5])



