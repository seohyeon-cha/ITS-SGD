function [h] = channel_ht_any(theta, phi, N_y, N_z, d, lambda)
    h_theta = exp(1j*2*pi*d*sin(theta).*(0:N_z-1)'/lambda); %(n,1)
    h_phi_given_theta = exp(1j*2*pi*d*cos(theta).*sin(phi).*(0:N_y-1)'/lambda); %(n,1)
    h = kron(h_theta, h_phi_given_theta);
end