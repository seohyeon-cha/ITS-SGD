clc;
load('channel.mat','phi_L','theta_L','path_gain','IRS_phase_vector','iteration','N_IRS','N_BS','L');

%% T matrix generation
lambda = 3e8/(28e9); 
d = lambda/2;
Rd = (d*sqrt(N_IRS))/sqrt(pi); 
theta_0 = atan((sqrt(N_IRS)*d/2)/(Rd));

rho_srf = 1;
G_BS = 2/(1-cos(theta_0));
G_IRS = 2;
[r_avg, r_matrix] = average_distance(Rd, d, N_IRS);
r_m = reshape(r_matrix,[N_IRS,1]);
k0 = lambda * sqrt(rho_srf*G_BS*G_IRS)/(4*pi*r_avg);

T = k0*exp(-1j*2*pi*r_m/lambda);
% T = 10*T;

tx_power = linspace(-10, 30, 9); %dBm
rate_w = zeros(1,length(tx_power));

for p0 = 1:length(tx_power)
    fprintf('tx_power=%d\n',tx_power(p0));
    pow = 1e-3*10.^(tx_power(p0)/10); %dBm -> W
    rate = zeros(1,iteration);
    batch_size = 50;
    beta1 = 0.9;
    beta2 = 0.999;
    mt = zeros(N_IRS, 1);
    vt = zeros(N_IRS, 1);

    for ite = 1:iteration
        Theta = IRS_phase_vector(:,:,ite);
        phi = angle(Theta);
        theta_l_epoch = theta_L(:,:,ite); %(L,1)
        phi_l_epoch = phi_L(:,:,ite);  %(L,1)
        t = 0;

        while 1
            t=t+1;
            grad_in = zeros(N_IRS,1,batch_size);
            
            for b0 = 1:batch_size
                
                % Channel vector between IRS and user
                h = zeros(N_IRS,1);
                for l=1:L
                    hl = path_gain*sqrt(0.5).*(randn(1,1)+1j*randn(1,1));
                    h = h + hl.* channel_ht(theta_l_epoch(l,1),phi_l_epoch(l,1),N_IRS,d,lambda);
                end
                h = h./sqrt(L);

                % Channel realization
                H = diag(h)*T;
                grad_in(:, :, b0) = real((2.*pow.*(H*H')*Theta)...
                    ./(1+pow.*Theta'*(H*H')*Theta).*(-1j*conj(Theta))); 

            end
            grad = mean(grad_in, 3);
            mt = beta1 .* mt + (1-beta1)*grad;
            vt = beta2 .* vt + (1-beta2)*grad.^2;
            lr = 0.1;
            eps = 1e-8;
            update = lr * mt ./ (1-beta1^t)./(sqrt(vt./(1-beta2^t))+eps);
            phi = phi + update;

            % New theta
            Theta = exp(1j*phi);

            if t>iteration
                break
            end
        end
        rate(ite) = test_rate(Theta, T, pow, path_gain, theta_l_epoch, phi_l_epoch, d, lambda, N_IRS, L);
    end
    rate_w(p0)=mean(rate);
end
figure
plot(tx_power,rate_w,'ro-');
save('mini_adam.mat','tx_power','rate_w','N_IRS');
