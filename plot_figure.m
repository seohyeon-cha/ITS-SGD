figure
line=2;
close all

% MISO
load('MISO.mat','tx_power','rate_w','N_IRS');
plot(tx_power, rate_w,'ro-','LineWidth',line);
grid on
hold on

% Perfect CSI
load('perfect_CSI.mat','tx_power','rate_w','N_IRS');
plot(tx_power, rate_w,'m^-','LineWidth',line);

% SGD - adam
load('mini_adam.mat','tx_power','rate_w','N_IRS');
plot(tx_power, rate_w,'bv-','LineWidth',line);

% Random phase
load('random.mat','tx_power','rate_w','N_IRS');
plot(tx_power, rate_w,'ks-','LineWidth',line);


%%
xlim([-10,30]);
ylabel('Achievable rate [bps/Hz]'); xlabel('Transmit power [dBm]');
legend('Conventional MISO','Perfect','SGD','Random','Location','best')