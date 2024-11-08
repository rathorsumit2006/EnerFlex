clc;
clear all;
close all;
tic;

% filename_in='Data Input_Irish pilot_Jan_2020_55Prosumers_AllPVnBESS Version.xlsx';% Put the input file name 
filename_in='Data_electricity.xlsx'; % Put the input file name 
filename_load='Load_.xlsx'; % Put the input file name 
filename_PV='PV.xlsx'; % Put the input file name 
filename_tariff='tariffs.xlsx'; % Put the input file name
filename_out='Case1.xlsx'; % Put the output filename 

dem1 = xlsread(filename_load,'Sheet1','A1:BC24'); % Selects Demand sheet of the input filename and corresponding cell area
dem = dem1; % Convert into kWhr

Pg = xlsread(filename_tariff,'tariff_minutes','A1:A24'); % Selects electricity price timeseries in relevant sheet of the input filename and corresponding cell area
FiT = xlsread(filename_tariff,'tariff_minutes','B1:B24'); % Selects feed-in-tariff timeseries in relevant sheet of the input filename and corresponding cell area

Solar1 = xlsread(filename_PV,'PV_2','A1:BC24'); % Selects PV generation sheet of the input filename and corresponding cell area
Solar = Solar1; % Convert into kWhr
res = Solar;

alpha = 3.3; % Maximum charge rate of battery 3.3 kW or 5kW
beta = 3.3; % Maximum discharge rate of battery 3.3 kW or 5kW

% Battery capacity
Su = 5; % Upper bounds of storage level in battery 5kWh or 13.4 kWh battery size
Sl = 1; % Lower bounds of storage level in battery 
eta_ch = 0.95; % Battery charging efficiency 
eta_dis = 0.95; % Battery discharging efficiency

Nh = length(dem(1,:)); % No of houses
psip2p = 1-0.0001; % Distribution network losses and conversion of DG for P2P sale % Original paper value: 0.076

% Battery placement
BatteryPlace = xlsread(filename_in,'Battery Placement','B3:BD3');  % Battery location B1-0%, B2-25%, B3-50%, B4-75%, B5-100%
Sp0 = 2*ones(1,Nh); % Initial battery state
Spp0 = BatteryPlace .* Sp0; % initial storage level in battery
S0 = Spp0(1,BatteryPlace==1);

Nt_end = 24;
Nt_begin = 1;
Nt = 24;
% dt = 1 / 60; % Time step in hours (1 minute)

% Creating decision variable sets
yalmip('clear')
G = sdpvar(Nt,Nh,'full'); % Grid import decision variable
Gfeed = sdpvar(Nt,Nh,'full'); % Grid feed-in decision variable
D = sdpvar(Nt,sum(BatteryPlace),'full'); % Discharging decision variable
C = sdpvar(Nt,sum(BatteryPlace),'full'); % Charging decision variable
S = sdpvar(Nt,sum(BatteryPlace),'full'); % SoC decision variable
O_ch = binvar(Nt,sum(BatteryPlace)); % BESS operating mode decision variable

% Objective function
Objective = (Pg(Nt_begin:Nt_end)' * (sum(G,2)) - FiT(Nt_begin:Nt_end)' * (sum(Gfeed,2)));

% Constraints YALMIP
Constraints = S(1,:) == S0 + eta_ch * C(1,:) - (1 / eta_dis) * D(1,:);
for i = 2:Nt
    Constraints = [Constraints; S(i,:) == S(i-1,:) + eta_ch * C(i,:)  - (1 / eta_dis) * D(i,:)];
end

k = 0;
for i = 1:Nh
    if BatteryPlace(i) == 1
        k = k + 1;
        Constraints = [Constraints, res(Nt_begin:Nt_end,i) + G(:,i) + D(:,k) >= dem(Nt_begin:Nt_end,i) + C(:,k) + Gfeed(:,i)];
    else
        Constraints = [Constraints, res(Nt_begin:Nt_end,i) + G(:,i) >= dem(Nt_begin:Nt_end,i) + Gfeed(:,i)];
    end
end

for t1 = 1:Nt
    Constraints = [Constraints, 0 <= G(t1,:) <= 1000];
    Constraints = [Constraints, 0 <= Gfeed(t1,:) <= 1000];
    Constraints = [Constraints, 0 <= D(t1,:) <= beta * (1 - O_ch(t1,:))];
    Constraints = [Constraints, 0 <= C(t1,:) <= alpha * O_ch(t1,:)];
end

for t1 = 1:Nt
    if (t1 == 1)
        Constraints = [Constraints, S(t1,:) == S0];     
    elseif (t1 < Nt)
        Constraints = [Constraints, Sl <= S(t1,:) <= Su];  
    else
        Constraints = [Constraints, S(t1,:) == S0];  
    end
end

ops = sdpsettings('verbose',1,'solver','gurobi','debug',1);
sol = optimize(Constraints,Objective,ops)

% Storing the sdp variable as values
Gv = value(G);
Gfeedv = value(Gfeed);
Cv = value(C);
Dv = value(D);
Sv = value(S);
Ochv = value(O_ch);
Pinj = Gv - Gfeedv; % Net injection in a node; Net injection = Grid import - Grid Export

% Summation of variable timeseries to extract KPI
Tot_Grid_import = sum(Gv,'all'); % total grid import 
Tot_Grid_export = sum(Gfeedv,'all'); % total grid feed-in
Peak_demand = max(sum(Gv,2));
Avg_demand = mean(sum(Gv,2));
Cost_Grid_import = Pg(Nt_begin:Nt_end)' * sum(Gv,2); % Cost of grid import
Revenue_Grid_export = FiT(Nt_begin:Nt_end)' * sum(Gfeedv,2); % Cost of grid export
KPI = [value(Objective(1,1)) Tot_Grid_import Tot_Grid_export 0 Peak_demand Cost_Grid_import Revenue_Grid_export Cost_Grid_import - Revenue_Grid_export];

% Writing results in excel files
xlswrite('Case2_without_KPI.xlsx',KPI,'KPIs_X_All_Jan20','D25') % writing KPIs in excel
xlswrite(filename_out,Gv,'G','C2')
xlswrite(filename_out,Gfeedv,'Gfeed','C2')
xlswrite(filename_out,Sv,'SoC','C2')
xlswrite(filename_out,Cv,'Battery Charge','C2')
xlswrite(filename_out,Dv,'Battery Discharge','C2')
xlswrite(filename_out,Pinj,'Net Pinj','C2')

toc;
fprintf('Electricity cost for community per day: %f\n', value(Objective(1,1)) / 100);
fprintf('Average Electricity cost per House per day: %f\n', value(Objective(1,1)) / (100 * 55));
fprintf('Peak Demand: %f\n', max(sum(Gv,2)));
fprintf('Average Demand: %f\n', mean(sum(Gv,2)));
fprintf('Cost of Grid Import per day: %f\n', Pg(Nt_begin:Nt_end)' * sum(Gv,2) / 100);
fprintf('Revenue from Grid Export day: %f\n', FiT(Nt_begin:Nt_end)' * sum(Gfeedv,2) / 100);

plot(sum(res'), 'LineWidth', 1.15, 'Color', 'r')
hold on
plot(sum(dem'), 'LineWidth', 1.15, 'Color', 'g')
hold on
plot(sum(Cv'), 'LineWidth', 1.15, 'Color', 'b')
hold on
plot(sum(Dv'), 'LineWidth', 1.15, 'Color', 'c')
hold on
plot(sum(value(Gfeed')), 'LineWidth', 1.15, 'Color', 'm')
hold on
plot(sum(value(G')), 'LineWidth', 1.15, 'Color', 'k')
legend({'Solar','Demand','Charging','Discharging','Export','Import'},'Location','northwest','NumColumns',3)
w = value(Objective(1,1)) / 100;
x = value(Objective(1,1)) / (100 * 55);
y = max(sum(Gv,2));
z = FiT(Nt_begin:Nt_end)' * sum(Gfeedv,2) / 100;
title({
    'Scenario NO P2P + 100% PV + 100% BESS',
    ['Community cost (Euro) = ' num2str(w)], 
    ['Cost/house (Euro) = ' num2str(x)], 
    ['Peak demand (kW) = ' num2str(y)], 
    ['Export Revenue (Euro) = ' num2str(z)]
    }, 'FontSize', 8); % Adjust FontSize as needed
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 32 20]); % x_width=10cm y_width=15cm
saveas(gcf,'fig9.tif')

C_Gv = Pg .* Gv;
C_Gfeed = FiT .* Gfeedv;
Indi_Cost = (sum(C_Gv) - sum(C_Gfeed));
Total_community_cost = sum(Indi_Cost);

balance = sum(res') + sum(value(G')) + sum(Dv') - sum(value(Gfeed')) - sum(Cv') - sum(dem'); % check power balance
