function r_avg = average_distance_MA(Rd, Rr, d, N_IRS)
    n = round(sqrt(N_IRS));
    AA_pos = [-sqrt(2)*Rr/2,-sqrt(2)*Rr/2,0];
    pos = zeros(n/2,n/2,3);

    for i = 1:n/2
        for j=1:n/2
            pos(i,j,:) = [-((n-1)/2)*d+(i-1)*d, -((n-1)/2)*d+(j-1)*d, Rd]-AA_pos;
        end
    end
    r_matrix = sqrt(pos(:,:,1).^2 + pos(:,:,2).^2 + pos(:,:,3).^2);

    r_sum = sum(r_matrix,'all');
    r_avg = r_sum/(N_IRS/4);
end
