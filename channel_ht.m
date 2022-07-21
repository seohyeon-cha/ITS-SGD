function [h] = channel_ht(theta, phi, N_IRS, d, lambda)
    n = round(sqrt(N_IRS));
    h_theta = exp(1j*2*pi*d*sin(theta).*(0:n-1)'/lambda); %(n,1)
    h_phi_given_theta = exp(1j*2*pi*d*cos(theta).*sin(phi).*(0:n-1)'/lambda); %(n,1)
    h = kron(h_theta, h_phi_given_theta);
end