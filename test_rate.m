function f = test_rate(IRS_phase_up, T, pow, path_gain, theta_l_epoch, phi_l_epoch, d, lambda, N_IRS, L)
    rng(123);
    iteration = 1000;
    ratec = zeros(1,iteration);
    for ite = 1:iteration
        h = zeros(N_IRS, 1);
        for l=1:L
            hl = path_gain*sqrt(0.5).*(randn(1,1)+1j*randn(1,1));
            h = h + hl.* channel_ht(theta_l_epoch(l,1),phi_l_epoch(l,1),N_IRS,d,lambda);
        end
        h = h./sqrt(L);

        % Channel realization
        H = diag(h)*T;
        f0 = real(IRS_phase_up'*(H*H')*IRS_phase_up);
        ratec(ite) = log2(1+pow*f0);
    end
    f= mean(ratec);
end
