function [r_avg, r_matrix] = average_distance(Rd, d, N_IRS)
    n = round(sqrt(N_IRS));
    pos = zeros(n,n,3);
    for i = 1:n
        for j=1:n
            pos(i,j,:) = [-((n-1)/2)*d+(i-1)*d, -((n-1)/2)*d+(j-1)*d, Rd];
        end
    end
    
    r_matrix = sqrt(pos(:,:,1).^2 + pos(:,:,2).^2 + pos(:,:,3).^2);
    r_sum = sum(r_matrix,'all');
    r_avg = r_sum/N_IRS;
end
