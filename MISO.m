clear;
load('channel_MISO.mat','phi_L','theta_L','path_gain','IRS_phase_vector','iteration','N_IRS','N_BS','L','lambda','d');

%% SNR 따라서, rate 구하는 거 (iteration 100번) 
tx_power = linspace(-10, 30, 10); %dBm
rate_w = zeros(1,length(tx_power));
% active antenna 일때는 d 바꿔야 하는 거 아닌가? (해보니까 결과는 크게 안 다르긴 함) 
% 일반 MISO 설정을 어떻게 할지 생각해보기 

for p0 = 1:length(tx_power)
    fprintf('tx_power=%d\n',tx_power(p0));
    pow = 1e-3*10.^(tx_power(p0)/10); %dBm -> W
    rate = zeros(1,iteration);
    G_BS = 2; % Antenna gain 
    for ite = 1:iteration
        theta_l_epoch = theta_L(:,:,ite); %(L,1)
        phi_l_epoch = phi_L(:,:,ite);  %(L,1)
        batch_size = 50;
        ratec = zeros(1,batch_size);
        for b0 = 1:batch_size
            % Channel vector between IRS and user
            h = zeros(N_IRS,1);
            for l=1:L
                hl = sqrt(G_BS)*path_gain*sqrt(0.5).*(randn(1,1)+1j*randn(1,1));
                h = h + hl.* channel_ht(theta_l_epoch(l,1),phi_l_epoch(l,1),N_IRS,d,lambda);
            end
            h = h./sqrt(L);
            ratec(b0) = log(1+pow/N_IRS.*real(h'*h));
        end
        rate(ite) = mean(ratec);
    end
    rate_w(p0)=mean(rate);
end
figure
plot(tx_power,rate_w,'ro-');
save('MISO.mat','tx_power','rate_w','N_IRS');
