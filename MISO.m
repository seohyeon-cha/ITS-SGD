clear;
load('channel_MISO.mat','phi_L','theta_L','path_gain','IRS_phase_vector','iteration','N_y','N_z','N_BS','L','lambda','d');


tx_power = linspace(-10, 30, 10); %dBm
rate_w = zeros(1,length(tx_power));

N_UPA = N_y * N_z;
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
            h = zeros(N_UPA,1);
            for l=1:L
                hl = sqrt(G_BS)*path_gain*sqrt(0.5).*(randn(1,1)+1j*randn(1,1));
                h = h + hl.* channel_ht_any(theta_l_epoch(l,1),phi_l_epoch(l,1),N_y,N_z,d,lambda);
%                 h = h + hl.* channel_ht(theta_l_epoch(l,1),phi_l_epoch(l,1),N_UPA,d,lambda);
            end
            h = h./sqrt(L);
            ratec(b0) = log2(1+pow.*real(h'*h));
        end
        rate(ite) = mean(ratec);
    end
    rate_w(p0)=mean(rate);
end
figure
plot(tx_power,rate_w,'ro-');
save('MISO.mat','tx_power','rate_w','N_UPA');
