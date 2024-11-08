clc;
clear all;
close all;
tic;

filename_in = 'Data_electricity1.xlsx'; % Put the input file name 
filename_load = 'Load_.csv'; % Put the input file name 
filename_PV = 'PV.xlsx'; % Put the input file name 
filename_tariff = 'tariffs.xlsx'; % Put the input file name
filename_out = 'Case4_.xlsx'; % Put the output filename 

% Load data and aggregate to hourly
dem1 = xlsread(filename_load, 'Load_', 'A1:BC1440'); % Original demand data in 1-minute intervals
dem = reshape(sum(reshape(dem1, 60, [])), 24, []); % Convert to hourly data by summing every 60-minute block

Pg = xlsread(filename_tariff, 'tariff_minutes', 'A1:A1440'); % Original tariff data in 1-minute intervals
Pg = reshape(mean(reshape(Pg, 60, [])), 24, 1); % Convert to hourly data by averaging every 60-minute block

FiT = xlsread(filename_tariff, 'tariff_minutes', 'B1:B1440'); % Original FiT data in 1-minute intervals
FiT = reshape(mean(reshape(FiT, 60, [])), 24, 1); % Convert to hourly data by averaging every 60-minute block

Ip2p_price = xlsread(filename_tariff, 'tariff_minutes', 'C1:C1440'); % Original P2P buy price data in 1-minute intervals
Ip2p_price = reshape(mean(reshape(Ip2p_price, 60, [])), 24, 1); % Convert to hourly data by averaging every 60-minute block

Xp2p_price = xlsread(filename_tariff, 'tariff_minutes', 'D1:D1440'); % Original P2P sell price data in 1-minute intervals
Xp2p_price = reshape(mean(reshape(Xp2p_price, 60, [])), 24, 1); % Convert to hourly data by averaging every 60-minute block

% Load PV data and aggregate to hourly
Solar1 = xlsread(filename_PV, 'PV_100', 'A1:BC1440'); % Original PV data in 1-minute intervals
Solar = reshape(sum(reshape(Solar1, 60, [])), 24, []); % Convert to hourly data by summing every 60-minute block
Solar = Solar ./ 2; % Convert into kWhr

res = Solar;

% Battery and other parameters remain unchanged
alpha = 7.3; % Maximum charge rate of battery 3.3 kW or 5kW
beta = 7.3; % Maximum discharge rate of battery 3.3 kW or 5kW
Su = 13.4; % Upper bounds of storage level in battery 5kWh or 13.4 kWh battery size
Sl = 1; % Lower bounds of storage level in battery 
eta_ch = 0.95; % Battery charging efficiency 
eta_dis = 0.95; % Battery discharging efficiency

Nh = length(dem(1,:)); % No of houses
psip2p = 1 - 0.0001; % Distribution network losses and conversion of DG for P2P sale

BatteryPlace = xlsread(filename_in, 'Battery Placement', 'B5:BD5'); % Battery location
Sp0 = 2 * ones(1, Nh); % Initial battery state
Spp0 = BatteryPlace .* Sp0; % Initial storage level in battery
S0 = Spp0(1, BatteryPlace == 1);

Nt_end = 24;
Nt_begin = 1; % First timestep of no of timesteps selected for the plot
Nt = 24; % 24 intervals (1 hour each)
% dt = 1; % Time step in hours (1 hour)

% Creating decision variable sets
yalmip('clear')
G = sdpvar(Nt, Nh, 'full'); % Grid import decision variable
Gfeed = sdpvar(Nt, Nh, 'full'); % Grid feed-in decision variable
D = sdpvar(Nt, sum(BatteryPlace), 'full'); % Discharging decision variable
C = sdpvar(Nt, sum(BatteryPlace), 'full'); % Charging decision variable
S = sdpvar(Nt, sum(BatteryPlace), 'full'); % SoC decision variable
O_ch = binvar(Nt, sum(BatteryPlace)); % BESS operating mode decision variable
Ip2p = sdpvar(Nt, Nh, Nh-1, 'full'); % Import P2P
Xp2p = sdpvar(Nt, Nh, Nh-1, 'full'); % Export P2P

% Summing over the third dimension to get total P2P import/export per prosumer
Total_Ip2p = sum(Ip2p, 3); % Total P2P imports for each prosumer at each time step
Total_Xp2p = sum(Xp2p, 3); % Total P2P exports for each prosumer at each time step

% Objective function including P2P trading costs and revenues
Objective = Pg(Nt_begin:Nt_end)' * sum(G, 2) - FiT(Nt_begin:Nt_end)' * sum(Gfeed, 2) ...
    + sum(Ip2p_price(Nt_begin:Nt_end)' * sum(Total_Ip2p, 2)) ...
   - sum(Xp2p_price(Nt_begin:Nt_end)' * sum(Total_Xp2p, 2));

Constraints = S(1,:) == S0 + eta_ch * C(1,:) - (1 / eta_dis) * D(1,:) ;
for i = 2:Nt
    Constraints = [Constraints; S(i,:) == S(i-1,:) + eta_ch * C(i,:)  - (1 / eta_dis) * D(i,:)];
end

Constraints = [Constraints, 0 <= G(:) <= inf, 0 <= Gfeed(:) <= inf];
Constraints = [Constraints, 0 <= D(:) <= beta*(1 - O_ch(:)), 0 <= C(:) <= alpha* O_ch(:)];
Constraints = [Constraints, Sl <= S(:) <= Su, 0 <= Ip2p(:) <= inf, 0 <= Xp2p(:) <= inf];

for t1 = 1:Nt
    if (t1 == 1)
        Constraints = [Constraints, S(t1,:) == S0];     
    elseif (t1 < Nt)
        Constraints = [Constraints, Sl <= S(t1,:) <= Su];  
    else
        Constraints = [Constraints, S(t1,:) == S0];  
    end
end

I = sdpvar(Nt, Nh, 'full');
X = sdpvar(Nt, Nh, 'full');

I = sum(Ip2p, 3);
X = sum(Xp2p, 3);

eq1 = [];
k = 0;
for i = 1:Nh
    if BatteryPlace(i) == 1
        k = k + 1;
        eq1 = [eq1, res(:,i) + G(:,i) + I(:,i) + D(:,k) == dem(:,i) + X(:,i) + C(:,k) + Gfeed(:,i)];
    else
        eq1 = [eq1, res(:,i) + G(:,i) + I(:,i) == dem(:,i) + X(:,i) + Gfeed(:,i)];
    end
end
Constraints = [Constraints, eq1];

IXpind = zeros(Nh-1, Nh);
for i = 1:length(IXpind(1,:))
    for j = 1:length(IXpind(:,1))
        if i <= j
            IXpind(j,i) = j + 1;
        else
            IXpind(j,i) = j;
        end
    end
end

Constraints = [Constraints, psip2p * sum(X,2) == sum(I,2)];

options = sdpsettings('verbose', 1, 'solver', 'gurobi');
diagnostics = optimize(Constraints, Objective, options);

Gv = value(G);
Gfeedv = value(Gfeed);
Cv = value(C);
Dv = value(D);
Sv = value(S);
Ochv = value(O_ch);
Ip2p = value(Ip2p);
Xp2p = value(Xp2p);
Pinj = Gv - Gfeedv;

Tot_Grid_import = sum(Gv,'all');
Tot_Grid_export = sum(Gfeedv,'all');
Peak_demand = max(sum(Gv,2));
Avg_demand = mean(sum(Gv,2));
Cost_Grid_import = Pg(Nt_begin:Nt_end)' * sum(Gv, 2);
Revenue_Grid_export = FiT(Nt_begin:Nt_end)' * sum(Gfeedv, 2);
KPI = [value(Objective(1,1)) Tot_Grid_import Tot_Grid_export 0 Peak_demand Cost_Grid_import Revenue_Grid_export Cost_Grid_import - Revenue_Grid_export];

P2P_Import = sum(Ip2p, 3);
P2P_export = sum(Xp2p, 3);

xlswrite('Case4_without_KPI.xlsx', KPI, 'KPIs_X_All_Jan20', 'D25')
xlswrite(filename_out, Gv, 'G', 'C2')
xlswrite(filename_out, Gfeedv, 'Gfeed', 'C2')
xlswrite(filename_out, Sv, 'SoC', 'C2')
xlswrite(filename_out, Cv, 'Battery Charge', 'C2')
xlswrite(filename_out, Dv, 'Battery Discharge', 'C2')
xlswrite(filename_out, Pinj, 'Net Pinj', 'C2')
xlswrite(filename_out, P2P_Import, 'P2P_Ixport', 'C2')
xlswrite(filename_out, P2P_export, 'P2P_export', 'C2')

toc;
fprintf('Electricity cost for community per day: %f\n', value(Objective(1,1))/100);
fprintf('Average Electricity cost per House per day: %f\n', value(Objective(1,1))/(100*Nh));
fprintf('Peak Demand: %f\n', max(sum(Gv, 2)));
fprintf('Average Demand: %f\n', mean(sum(Gv, 2)));
fprintf('Cost of Grid Import per day: %f\n', Pg(Nt_begin:Nt_end)'*sum(Gv, 2)/100);
fprintf('Revenue from Grid Export day: %f\n', FiT(Nt_begin:Nt_end)'*sum(Gfeedv, 2)/100);
fprintf('P2P buying cost %f\n', value(sum(Ip2p_price(Nt_begin:Nt_end)' * sum(Total_Ip2p, 2)))/100);
fprintf('P2P selling cost %f\n', value(sum((Xp2p_price(Nt_begin:Nt_end)' * sum(Total_Xp2p, 2))))/100);

plot(sum(res', 2), 'LineWidth', 1.15, 'Color', 'r')
hold on
plot(sum(dem', 2), 'LineWidth', 1.15, 'Color', 'g')
hold on
plot(sum(Cv', 2), 'LineWidth', 1.15, 'Color', 'b')
hold on
plot(sum(Dv', 2), 'LineWidth', 1.15, 'Color', 'c')
hold on
plot(sum(value(Gfeed'), 2), 'LineWidth', 1.15, 'Color', 'm')
hold on
plot(sum(value(G'), 2), 'LineWidth', 1.15, 'Color', 'k')
hold on
plot(sum(P2P_Import', 2), 'LineWidth', 1.15, 'Color', 'y')
xlabel('Time (hours)')
ylabel('Power(kW)')
legend({'Solar','Demand','Charging','Discharging','Export','Import','P2P Power'},'Location','northwest','NumColumns',3)

C_Gv = Pg .* Gv;
C_Gfeed = FiT .* Gfeedv;
C_Ip2p = Ip2p_price .* value(Total_Ip2p);
C_Xp2p = Xp2p_price .* value(Total_Xp2p);
Indi_Cost = (sum(C_Gv) - sum(C_Gfeed) + sum(C_Ip2p) - sum(C_Xp2p))/100;
Total_community_cost = sum(Indi_Cost);

w = value(Objective(1,1))/100; % Total community cost in euro
x = value(Objective(1,1))/(100*Nh); % Individual cost of prosumers
y = max(sum(Gv, 2)); % Peak demand
y1 = Pg(Nt_begin:Nt_end)' * sum(Gv, 2)/100; % Cost of grid import 
z = FiT(Nt_begin:Nt_end)' * sum(Gfeedv, 2)/100; % Revenue of grid export 
z1 = value(sum((Xp2p_price(Nt_begin:Nt_end)' * sum(Total_Xp2p, 2))))/100; % P2P selling profit to prosumers

title({
    'Scenario-1 P2P + 100% PV + 100% BESS',
    ['Community cost (Euro) = ' num2str(w)], 
    ['Cost/house (Euro) = ' num2str(x)], 
    ['Peak demand (kW) = ' num2str(y)], 
    ['Grid Import Cost(Euro) = ' num2str(y1)], 
    ['Export Revenue (Euro) = ' num2str(z)],
    ['P2P Selling profit(Euro) = ' num2str(z1)]
    }, 'FontSize', 8);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 32 20]);

saveas(gcf, 'fig9.tif')

% Plot hourly prices
plot(Pg, 'LineWidth', 1.15, 'Color', 'r')
hold on
plot(FiT, 'LineWidth', 1.15, 'Color', 'g')
hold on
plot(Ip2p_price, 'LineWidth', 1.15, 'Color', 'b')
hold on
plot(Xp2p_price, 'LineWidth', 1.15, 'Color', 'c')
xlabel('Time(hours)')
ylabel('Price (cent)')
legend({'Retail Price','FiT','P2P Buy Price','P2P Sell Price'},'Location','northwest','NumColumns',3)
