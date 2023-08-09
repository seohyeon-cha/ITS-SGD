rng(123);

% Components 
N_y = 2;
N_z = 2;
N_UPA = N_y * N_z; % number of elements in IRS
N_BS = 1; % number of BS antenna
L = 4; % number of channel paths

% Path loss 
l = 100; %[m] distance between transmitter and receiver
eta = 2; % pathloss exponent 

% noise power = bandwidth W x noise PSD No x noise figure Nf
% W = 100 MHz, No = - 174 dBm/Hz, Nf = 6 dB
noise = (100*1e6)*(1e-3*10^(-174/10))*(10^(6/10)); 
lambda = 3e8/(28e9); % f = 28GHz
pathloss = (lambda/(4*pi))^2/l^2;
path_gain = sqrt(pathloss/noise);
d = lambda/2;
% d = 1; 

% Channel gain
iteration = 200;
% AoD
phi_L = zeros(L,1,iteration);
theta_L = zeros(L,1,iteration);


% IRS phase shift matrix 
IRS_phase_vector = exp(1j.*rand(N_UPA,1,iteration).*2.*pi);

% Make batch size of realizations 
for j0 = 1:iteration
    phi_L(:,:,j0) = 2/3*pi*(-1+2*rand(L,1));
    theta_L(:,:,j0) = 1/2*pi*(-1+2*rand(L,1));
end



save('channel_MISO.mat','phi_L','theta_L','path_gain','IRS_phase_vector','iteration','N_y','N_z','N_BS','L','lambda','d');



