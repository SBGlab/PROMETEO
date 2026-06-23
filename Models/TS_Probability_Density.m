clear all
clc
close all

% install IDESS https://doi.org/10.5281/zenodo.7788692

dimension = 2;

% Hill coefficients
nu_TtgR  = 3;
nu_BetlR = 2.5;
nu_TetR  = 2.15;

% mRNA decay rates
dTtgRm = 0.1386;
dBetlRm = 0.1386;
dTetRm = 0.1386;

% Protein decay rates
dBetlR = 1.65e-2;
dTtgR  = 1.65e-2;
dTetR  = 1.65e-2;

% Promoter strengths
kPBetl = 14.09;
kPTet  = 13.01;
kPTtg  = 2.49;

% Dissociation constants
KBCD8 = 0.56;
KRBST = 0.09;

% Leakage constants
K_leakage_pBetl = 0.0141;
K_leakage_pPtg  = 0.0025;
K_leakage_pTet  = 0.01301;

% Epsilon 
epsilon_pBetl = K_leakage_pBetl / kPBetl;
epsilon_pPtg  = K_leakage_pPtg  / kPTtg;
epsilon_pTet  = K_leakage_pTet  / kPTet;


% ================================================================
TS_label = {'1A','1C','2A','2C','3A','3C'};
labels = TS_label;

problemMesh.spatial_mesh = [0 800 800; ...
                             0 600 600];
problemMesh.time_mesh = [0 700 100 5];  % t0  tf  steps  snapshots

no_of_realizations = 1e5;


% --------------------------------------------------
% CONFIG 1: BetlR – TtgR
% --------------------------------------------------
disp('Running CONFIG 1: BetlR–TtgR');

K_BetlRA = 8;  K_BetlRC = 3;
K_TtgRA  = 3;  K_TtgRC  = 3;

H = {0 nu_TtgR; nu_BetlR 0};
epsilon = { {epsilon_pBetl,1}; {epsilon_pPtg,1} };

disp('pTS1 (1A, W–S)');
rCoef_1A = { kPBetl KBCD8 dTtgRm dTtgR; ...
             kPTtg  KRBST dBetlRm dBetlR };

K = {0 K_TtgRA; K_BetlRA 0};

task = problemData(2,H,K,epsilon,rCoef_1A);
sim_1A = task.SSAsimulationGPU(problemMesh,no_of_realizations);
X_1A = double(sim_1A{end}(2,:));
Y_1A = double(sim_1A{end}(4,:));


disp('pTS1 (1C, S–W)');
rCoef_1C = { kPBetl KRBST dTtgRm dTtgR; ...
             kPTtg  KBCD8 dBetlRm dBetlR };

K = {0 K_TtgRC; K_BetlRC 0};

task = problemData(2,H,K,epsilon,rCoef_1C);
sim_1C = task.SSAsimulationGPU(problemMesh,no_of_realizations);
X_1C = double(sim_1C{end}(2,:));
Y_1C = double(sim_1C{end}(4,:));


% --------------------------------------------------
% CONFIG 2: BetlR – TetR
% --------------------------------------------------
disp('Running CONFIG 2: BetlR–TetR');

K_TetRA = 3;  K_TetRC = 3;
K_BetlRA = 9; K_BetlRC = 3;

H = {0 nu_TetR; nu_BetlR 0};
epsilon = { {epsilon_pBetl,1}; {epsilon_pTet,1} };

disp('pTS2 (2A, W–S)');
rCoef_2A = { kPBetl KBCD8 dTetRm dTetR; ...
             kPTet  KRBST dBetlRm dBetlR };

K = {0 K_TetRA; K_BetlRA 0};

task = problemData(2,H,K,epsilon,rCoef_2A);
sim_2A = task.SSAsimulationGPU(problemMesh,no_of_realizations);

X_2A = double(sim_2A{end}(2,:));
Y_2A = double(sim_2A{end}(4,:));


disp('pTS2 (2C, S–W)');
rCoef_2C = { kPBetl KRBST dTetRm dTetR; ...
             kPTet  KBCD8 dBetlRm dBetlR };

K = {0 K_TetRC; K_BetlRC 0};

task = problemData(2,H,K,epsilon,rCoef_2C);
sim_2C = task.SSAsimulationGPU(problemMesh,no_of_realizations);

X_2C = double(sim_2C{end}(2,:));
Y_2C = double(sim_2C{end}(4,:));


% --------------------------------------------------
% CONFIG 3: TtgR – TetR
% --------------------------------------------------
disp('Running CONFIG 3: TtgR–TetR');

K_TtgRA= 4.5;    K_TetRA =3;  
K_TtgRC = 1.5;   K_TetRC =3;

H = {0 nu_TetR; nu_TtgR 0};
epsilon = { {epsilon_pPtg,1}; {epsilon_pTet,1} };


disp('pTS3 (3A, W–S)');
rCoef_3A = { kPTtg KBCD8 dTetRm dTetR; ...
             kPTet KRBST dTtgRm dTtgR };

K = {0 K_TetRA; K_TtgRA 0};

task = problemData(dimension,H,K,epsilon,rCoef_3A);
sim_3A = task.SSAsimulationGPU(problemMesh,no_of_realizations);

X_3A = double(sim_3A{end}(2,:));
Y_3A = double(sim_3A{end}(4,:));


disp('pTS3 (3C, S–W)');
rCoef_3C = { kPTtg KRBST dTetRm dTetR; ...
             kPTet KBCD8 dTtgRm dTtgR };

K = {0 K_TetRC; K_TtgRC 0};

task = problemData(dimension,H,K,epsilon,rCoef_3C);
sim_3C = task.SSAsimulationGPU(problemMesh,no_of_realizations);

X_3C = double(sim_3C{end}(2,:));
Y_3C = double(sim_3C{end}(4,:));


% ================================================================
disp('--- Analyzing GFP distributions ---');

sim_map = containers.Map();
sim_map('1A') = sim_1A;
sim_map('1C') = sim_1C;
sim_map('2A') = sim_2A;
sim_map('2C') = sim_2C;
sim_map('3A') = sim_3A;
sim_map('3C') = sim_3C;

percentage_green = zeros(1,6);
pop1_all = cell(1,6);
pop2_all = cell(1,6);

green_cutoff = 100;

for ii = 1:6
    label = TS_label{ii};
    sim_data = sim_map(label);

    % Extract GFP
    X = double(sim_data{end}(2,:));

    % KDE estimation
    [f, xi] = ksdensity(X);
    f_smooth = smoothdata(f,'gaussian',15);
    f_norm = f_smooth ./ max(f_smooth);

    % Detect valleys
    [~, locs_min] = findpeaks(-f_norm,'MinPeakProminence',0.01);

    if isempty(locs_min)
        [~, locs_min] = findpeaks(-f_norm,'MinPeakProminence',0.00001);
    end
    if numel(locs_min) > 1
        [~, locs_min] = findpeaks(-f_norm,'MinPeakProminence',0.1);
    end


    % --------------------------------------------------------------
    % Case 1: No valley → unimodal → fallback threshold
    % --------------------------------------------------------------
    if isempty(locs_min)
        warning('No valley detected for %s → using fixed cutoff.', label);
        pop1 = X(X <= green_cutoff);
        pop2 = X(X > green_cutoff);

        percentage_green(ii) = numel(pop2) / numel(X);

        % Store populations
        pop1_all{ii} = pop1;
        pop2_all{ii} = pop2;

    else
        % ----------------------------------------------------------
        % Case 2: Bimodal → use valley threshold
        % ----------------------------------------------------------
        threshold = xi(locs_min(1));

        pop1 = X(X <= threshold);
        pop2 = X(X > threshold);

        percentage_green(ii) = numel(pop2) / numel(X);

        pop1_all{ii} = pop1;
        pop2_all{ii} = pop2;
    end
end


% ================================================================
disp('--- % GREEN CELLS ---');
for ii = 1:6
    fprintf('%s: %.2f %%\n', TS_label{ii}, percentage_green(ii)*100);
end

disp('--- Plotting results ---');

figure(1);

plot_titles = { ...
    'BetlR–TtgR : pTS1(W-S)', ...
    'BetlR–TtgR : pTS1(S-W)', ...
    'BetlR–TetR : pTS2(W-S)', ...
    'BetlR–TetR : pTS2(S-W)', ...
    'TtgR–TetR : pTS3(W-S)', ...
    'TtgR–TetR : pTS3(S-W)' };

X_data = {X_1A, X_2A,  X_3A, X_1C, X_2C, X_3C};
Y_data = {Y_1A, Y_2A,  Y_3A, Y_1C, Y_2C, Y_3C};


plot_to_TS = [1 3 5 2 4 6];   % mapping from plot index → TS_label index

for i = 1:6
    subplot(2,3,i)
    histogram2(X_data{i}, Y_data{i}, 'Numbins',30,'Normalization','pdf','FaceColor','flat');
    xlabel('GFP'); ylabel('RFP'); zlabel('PDF');
    colorbar;

    title(sprintf('%s — %.2f %% green', ...
        plot_titles{plot_to_TS(i)}, percentage_green(plot_to_TS(i))*100));

    set(gca,'FontSize',12);
end



