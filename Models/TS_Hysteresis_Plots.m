close all
clear all
clc


kLeak_pBetl = 0.0141; k_pBtel = 14.09;
kLeak_pTtg  = 0.0025; k_pTtg = 2.49;

theta_TtgR  = 400; degm_TtgR  = 1.386e-1; deg_TtgR  = 1.65e-2; n_TtgR  = 3;
theta_BetlR = 178; degm_BetlR = 1.386e-1; deg_BetlR = 1.65e-2; n_BetlR = 2.5;

thetaI_Cho = 0.2926; nI_Cho = 2.7;
thetaI_Nar = 278.0613; nI_Nar = 1.9;

kpG = 0.39; kpR = 0.39;

k_pBtel_vec = [6 8 10 13 15 17];                  
kpGv        = [0.2 0.4 0.6 0.8 1 1.2 1.4 1.6];   
k_pTtg_vec  = [2.49    3.1233    3.7567    4.39];    
kpRv        =[0.37 0.5 0.63 0.85];

nCurvesG = numel(kpGv);
nCurvesR = numel(kpRv);
nCurvesB = numel(k_pBtel_vec);
nCurvesT = numel(k_pTtg_vec);

%% ---------------- Master Palette ----------------
Nmaster = 256;

h_blue = linspace(0.58, 0.50, Nmaster);
s_blue = linspace(0.85, 0.95, Nmaster);
v_blue = linspace(0.80, 0.70, Nmaster);
cmap_blue = hsv2rgb([h_blue(:), s_blue(:), v_blue(:)]);

h_orange = linspace(0.02, 0.10, Nmaster);
s_orange = linspace(0.90, 0.75, Nmaster);
v_orange = linspace(0.75, 0.85, Nmaster);
cmap_orange = hsv2rgb([h_orange(:), s_orange(:), v_orange(:)]);

%% ---------------- Simulation Setup ----------------
Nar   = 0:5:1000;
IG    = 0;
tspan = 0:0.1:1400;
x0_f  = [0;0;0;0];
x0_b  = [0;0;10;10];

figure(1); clf


subplot(2,2,1); hold on; box on
xlabel('Nar (\muM)'); ylabel('GFP_{norm}');
title('Varying k_G');
set(gca,'FontSize',14,'FontName','Helvetica')
axis([0 600 0 1.2])

for jj = 1:nCurvesB

    color_f = pickColor(cmap_blue, nCurvesB, jj);
    color_b = pickColor(cmap_orange, nCurvesB, jj);

    par = [kLeak_pBetl, kLeak_pTtg, theta_TtgR, theta_BetlR, ...
           thetaI_Cho, thetaI_Nar, nI_Cho, nI_Nar, ...
           degm_TtgR, degm_BetlR, deg_TtgR, deg_BetlR, ...
           n_TtgR, n_BetlR, k_pBtel_vec(jj), k_pTtg, kpG, kpR];

    x_f_end = zeros(length(Nar),4);
    x_b_end = zeros(length(Nar),4);

    for ii = 1:length(Nar)
        parstr.par = par;
        parstr.IG  = IG;
        parstr.IR  = Nar(ii);

        [~,x] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
        x_f_end(ii,:) = x(end,:);

        [~,x] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
        x_b_end(ii,:) = x(end,:);
    end

    PG_f = x_f_end(:,2);
    PG_b = x_b_end(:,2);

    PG_f = PG_f / max(PG_b);
    PG_b = PG_b / max(PG_b);

    plot(Nar, PG_f, '-', 'Color', color_f, 'LineWidth', 2)
    plot(Nar, PG_b, '--', 'Color', color_b, 'LineWidth', 2)

end

subplot(2,2,2); hold on; box on
xlabel('Nar (\muM)'); ylabel('GFP_{norm}');
title('Varying kx_G');
set(gca,'FontSize',14,'FontName','Helvetica')

for jj = 1:nCurvesG

    color_f = pickColor(cmap_blue, nCurvesG, jj);
    color_b = pickColor(cmap_orange, nCurvesG, jj);

    par = [kLeak_pBetl, kLeak_pTtg, theta_TtgR, theta_BetlR, ...
           thetaI_Cho, thetaI_Nar, nI_Cho, nI_Nar, ...
           degm_TtgR, degm_BetlR, deg_TtgR, deg_BetlR, ...
           n_TtgR, n_BetlR, k_pBtel, k_pTtg, kpGv(jj), kpR];

    x_f_end = zeros(length(Nar),4);
    x_b_end = zeros(length(Nar),4);

    for ii = 1:length(Nar)
        parstr.par = par;
        parstr.IG  = IG;
        parstr.IR  = Nar(ii);

        [~,x] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
        x_f_end(ii,:) = x(end,:);

        [~,x] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
        x_b_end(ii,:) = x(end,:);
    end

    PG_f = x_f_end(:,2);
    PG_b = x_b_end(:,2);

    PG_f = PG_f / max(PG_b);
    PG_b = PG_b / max(PG_b);

    plot(Nar, PG_f, '-', 'Color', color_f, 'LineWidth', 2)
    plot(Nar, PG_b, '--', 'Color', color_b, 'LineWidth', 2)
end


subplot(2,2,3); hold on; box on
xlabel('Nar (\muM)'); ylabel('GFP_{norm}');
title('Varying k_R');
set(gca,'FontSize',14,'FontName','Helvetica')
axis([0 600 0 1.2])

for jj = 1:nCurvesT

    color_f = pickColor(cmap_blue, nCurvesT, jj);
    color_b = pickColor(cmap_orange, nCurvesT, jj);

    par = [kLeak_pBetl, kLeak_pTtg, theta_TtgR, theta_BetlR, ...
           thetaI_Cho, thetaI_Nar, nI_Cho, nI_Nar, ...
           degm_TtgR, degm_BetlR, deg_TtgR, deg_BetlR, ...
           n_TtgR, n_BetlR, k_pBtel, k_pTtg_vec(jj), kpG, kpR];

    x_f_end = zeros(length(Nar),4);
    x_b_end = zeros(length(Nar),4);

    for ii = 1:length(Nar)
        parstr.par = par;
        parstr.IG  = IG;
        parstr.IR  = Nar(ii);

        [~,x] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
        x_f_end(ii,:) = x(end,:);

        [~,x] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
        x_b_end(ii,:) = x(end,:);
    end

    PG_f = x_f_end(:,2);
    PG_b = x_b_end(:,2);

    PG_f = PG_f / max(PG_b);
    PG_b = PG_b / max(PG_b);

    plot(Nar, PG_f, '-', 'Color', color_f, 'LineWidth', 2)
    plot(Nar, PG_b, '--', 'Color', color_b, 'LineWidth', 2)
end

subplot(2,2,4); hold on; box on
xlabel('Nar (\muM)'); ylabel('GFP_{norm}');
title('Varying kx_R');
set(gca,'FontSize',14,'FontName','Helvetica')

for jj = 1:nCurvesR

    color_f = pickColor(cmap_blue, nCurvesR, jj);
    color_b = pickColor(cmap_orange, nCurvesR, jj);

    par = [kLeak_pBetl, kLeak_pTtg, theta_TtgR, theta_BetlR, ...
           thetaI_Cho, thetaI_Nar, nI_Cho, nI_Nar, ...
           degm_TtgR, degm_BetlR, deg_TtgR, deg_BetlR, ...
           n_TtgR, n_BetlR, k_pBtel, k_pTtg, kpG, kpRv(jj)];

    x_f_end = zeros(length(Nar),4);
    x_b_end = zeros(length(Nar),4);

    for ii = 1:length(Nar)
        parstr.par = par;
        parstr.IG  = IG;
        parstr.IR  = Nar(ii);

        [~,x] = ode45(@TS_ODE, tspan, x0_f, [], parstr);
        x_f_end(ii,:) = x(end,:);

        [~,x] = ode45(@TS_ODE, tspan, x0_b, [], parstr);
        x_b_end(ii,:) = x(end,:);
    end

    PG_f = x_f_end(:,2);
    PG_b = x_b_end(:,2);

    PG_f = PG_f / max(PG_b);
    PG_b = PG_b / max(PG_b);

    plot(Nar, PG_f, '-', 'Color', color_f, 'LineWidth', 2)
    plot(Nar, PG_b, '--', 'Color', color_b, 'LineWidth', 2)
end

%% ---------------- Local Function ----------------
function c = pickColor(cmap, n, idx)
    inds = round(linspace(1, size(cmap,1), n));
    c = cmap(inds(idx), :);
end
