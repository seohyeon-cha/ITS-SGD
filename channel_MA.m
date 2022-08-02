rng(123);

% Components 
N_IRS = 16; % number of elements in IRS
N_BS = 1; % number of BS antenna
L = 4; % number of channel paths
N_active = 4;
% Path loss 
l = 100; %[m] distance between transmitter and receiver
eta = 2; % pathloss exponent ------------------------------------------------------------------2.2/2.8/3.5 도 해보기

% noise power = bandwidth W x noise PSD No x noise figure Nf
% W = 100 MHz, No = - 174 dBm/Hz, Nf = 6 dB
noise = (100*1e6)*(1e-3*10^(-174/10))*(10^(6/10)); % [W] 계산해보면 -118 dB임 
lambda = 3e8/(28e9); % f = 28GHz
pathloss = (lambda/(4*pi))^2/l^2;
path_gain = sqrt(pathloss/noise);
d = lambda/2;

% Channel gain
iteration = 200;
% AoD
phi_L = zeros(L,1,iteration);
theta_L = zeros(L,1,iteration);

% IRS phase shift matrix 
IRS_phase_vector = exp(1j.*rand(N_IRS,1,iteration).*2.*pi);
W = zeros(N_active, 1, iteration);

% Make batch size of realizations 
for j0 = 1:iteration
    phi_L(:,:,j0) = 2/3*pi*(-1+2*rand(L,1));
    theta_L(:,:,j0) = 1/3*pi*(-1+2*rand(L,1));
    W_j = rand(N_active, 1);
    W(:,:,j0) = W_j./norm(W_j);
end

save('channel_MA.mat','phi_L','theta_L','path_gain','IRS_phase_vector','iteration','N_IRS','N_BS','N_active','W','L','lambda','d');


