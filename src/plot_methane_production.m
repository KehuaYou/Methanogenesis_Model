% Plot of methanogenesis model (Introductory Carbon Balance Model)
% Find the ICBM model in Andren & Katterer (1997)
C_stable_old = 0.99;
C_labile_old = 0.01;
r_temperature = 1;
k1= 0.3e-7;     % per second
k2 = k1/0.5e4;  % per second
h=0.7;  % fraction of labile organic carbon humidification into stable pool

delta_t = 86400*0.1;
dt = (0:delta_t:86400*365*100);

C_labile = C_labile_old.*exp(-k1*r_temperature*dt);
C_stable= (C_stable_old-h*k1*r_temperature*C_labile_old/(r_temperature*(k2-k1))).*exp(-k2*r_temperature.*dt) + ...
        h*k1*r_temperature*C_labile_old/(r_temperature*(k2-k1)).*exp(-k1*r_temperature.*dt);
qg_biogenic = 0.5*((1-h)*k1*r_temperature.*C_labile + k2*r_temperature.*C_stable);

cumulative_biogenic = zeros(length(dt),1);
for i = 2:length(dt)
    cumulative_biogenic(i) = cumulative_biogenic(i-1) + (qg_biogenic(i-1)+qg_biogenic(i))/2*delta_t;
end

figure
plot(dt./(86400*365), C_stable./C_stable_old.*100,'r','LineWidth', 1)
hold on
plot(dt./(86400*365), C_labile./C_labile_old.*100,'g','LineWidth', 1)
set(gca,'FontSize',18, 'LineWidth', 1)
xlabel('Time, (years)', 'FontSize', 18)
ylabel('Percentage of remaining organic carbon, (%)','FontSize', 18)
legend('stable','labile')


figure
plot(dt./(86400*365), qg_biogenic.*1e6.*86400,'b','LineWidth', 1)
set(gca,'FontSize',18, 'LineWidth', 1)
xlabel('Time, (years)', 'FontSize', 18)
ylabel('Methane generation rate, (mg CH4-C / kg C /day)','FontSize', 18)


figure
plot(dt./(86400*365), cumulative_biogenic*1e3,'g','LineWidth', 1)
set(gca,'FontSize',18, 'LineWidth', 1)
xlabel('Time, (years)', 'FontSize', 18)
ylabel('Cummulative methane generation, (g CH4-C / kg C) ','FontSize', 18)

% data from Knoblauch(2018): Figure 2b
% time in days, cumulative methane production in g CH4-C/kg C 
data = [4.18613	0.0196521
8.40237	0.082611
15.3948	0.137674
20.9963	0.192748
32.1612	0.247778
47.493	0.294902
60.0568	0.361732
96.27	0.43625
138.068	0.54222
172.893	0.620686
206.328	0.699163
285.714	0.856018
334.428	0.906816
365.045	0.934134
664.206	1.12863
875.706	1.26476
1080.62	1.9285
];

hold on
scatter(data(:,1)./365, data(:,2))





