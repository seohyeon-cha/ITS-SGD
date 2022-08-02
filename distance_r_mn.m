function r_matrix = distance_r_mn(Rd, Rr, d, N_IRS)
    n = round(sqrt(N_IRS));
    AA_pos = [-sqrt(2)*Rr/2,-sqrt(2)*Rr/2,0;-sqrt(2)*Rr/2,sqrt(2)*Rr/2,0 ;sqrt(2)*Rr/2,-sqrt(2)*Rr/2,0 ;sqrt(2)*Rr/2,sqrt(2)*Rr/2,0];
    r_matrix = zeros(n,n,4);
    for id_act=1:4
        pos = zeros(n,n,3);
        for i = 1:n
            for j=1:n
                pos(i,j,:) = [-((n-1)/2)*d+(i-1)*d, -((n-1)/2)*d+(j-1)*d, Rd]-AA_pos(id_act, :);
            end
        end
        r_matrix(:,:,id_act) = sqrt(pos(:,:,1).^2 + pos(:,:,2).^2 + pos(:,:,3).^2);
    end
end
