clc;
load('channel.mat','phi_L','theta_L','path_gain','IRS_phase_vector','iteration','N_IRS','N_BS','L');

%% SNR 따라서, rate 구하는 거 (iteration 100번) 
tx_power = linspace(-10, 30, 10); %dBm
rate_w = zeros(1,length(tx_power));

for p0 = 1:length(tx_power)
    fprintf('tx_power=%d\n',tx_power(p0));
    pow = 1e-3*10.^(tx_power(p0)/10); %dBm -> W
    rate = zeros(1,iteration);
    for ite = 1:iteration
        ratec = zeros(1,batch_size);
        for b0 = 1:batch_size
            % Channel vector between IRS and user
            h = zeros(N_IRS,1);
            for l=1:L
                hl = path_gain*sqrt(0.5).*(randn(1,1)+1j*randn(1,1));
                h = h + hl.* channel_ht(theta_l_epoch(l,1),phi_l_epoch(l,1),N_IRS,d,lambda);
            end
            h = h./sqrt(L);
            ratec(b0) = log(1+pow.*real(h'*h));
        end
        rate(ite) = mean(ratec);
    end
    rate_w(p0)=mean(rate);
end
figure
plot(tx_power,rate_w,'ro-');
save('MISO.mat','tx_power','rate_w','N_IRS');