function [G,G_image] = Greens_Function(position,position_prime,params,free_boundary)
    % Greens function for forces in a half elastic space, see 
    % "Elastic interactions of active cells with soft materials" for
    % details of implementation... THIS IS NOT CORRECT FOR A POISSON RATIO
    % THAT IS NOT 1/2

    % Positions of the image forces
    position_im_prime = position_prime;
    position_im_prime(3) = -position_prime(3);
    
    poisson_ratio = params.poisson_ratio; % I am pretty sure this is always 1/2 for our cases...
    delta = @(i,j) (i==j) + 0; % hacky kronecker delta function
    
    % With free/clamped boundary conditions
    if free_boundary
        M = 3 - 4*poisson_ratio;
        J = -1;
        C = 2*(1-poisson_ratio);
        B = 2*(1-2*poisson_ratio);
    else
        M = -1;
        J = 1/(3-4*poisson_ratio);
        C = 0;
        B = 0;
    end
    
    % A sort of distance vector, given in paper. Distance between where we
    % measure and the image force
    s_vector = -position_im_prime + position; % column vector
    s = norm(s_vector); % scalar
    
    % distance vector between where we are measuring and the true force
    r_vector = position_prime - position;
    r = norm(r_vector);
    
    G_infinity = 1/(8*pi*r)*(eye(3) + r_vector*r_vector'./r^2); % Greens Function for infinite space, same as always
    
    G_im_infinity = 1/(8*pi*s)*(eye(3) + s_vector*s_vector'./s^2); % Greens Function for infinite space, with s
    
    %DPHI = NaN*ones(3,1); % phi_i from paper, first derivative (doesnt need to be computed for incompressible...)
    %D2PSI = NaN*ones(3); % psi_ij from paper, second derivative (doesnt need to be computed for incompressible...)
    sij3 = NaN*ones(3); % third derivative of s = |r' - r|
    Ds_inv = -s_vector/s^3; % derivative of 1/s = 1/|r'-r|
    D2s_inv = -eye(3)/s^3 + 3*(s_vector*s_vector')/s^5; % second derivative of 1/s
    G_image = NaN*ones(3);
    
%     for i = 1:3
%         DPHI(i) = (delta(i,3)+fix(i)*s_vector(i)/s)/(position(3) + position_prime(3) +s);
%     end
    
    for i=1:3
        for j=1:3
            sij3(i,j) = -(delta(i,j)*s_vector(3)+delta(i,3)*s_vector(j)+delta(3,j)*s_vector(i))/s^3 + 3*s_vector(i)*s_vector(j)*s_vector(3)/s^5;
        end
    end
    
    % this part doesnt contribute for poisson_ratio = 1/2
%     for i = 1:3
%         for j = 1:3
%              D2PSI(i,j) = delta(j,3)*(delta(i,3)+s_vector(i)/s)/(position(3)+position_prime(3)+s)...
%                  - (position(3)+position_prime(3))*(delta(i,3)+s_vector(i)/s)*(delta(j,3)+s_vector(j)/s)/(position(3) + position_prime(3) + s)^2 ...
%                  +(position(3) + position_prime(3))/(position(3) + position_prime(3) + s)*(delta(i,j)/s - s_vector(i)*s_vector(j)/s^3) ...
%                  + delta(i,3)*DPHI(j) - delta(i,j)/s + s_vector(i)*s_vector(j)/s^3;
%         end
%     end
    
    % Combining this all, we get (I ommited a few terms proportional to (1-2*poisson_ratio)
    for i = 1:3
        for j = 1:3
             G_image(i,j) = M*G_im_infinity(i,j) ...
                 +J*position_prime(3)/(4*pi)*(sij3(i,j) - 2*delta(j,3)*sij3(i,3)-4*(1-poisson_ratio)*delta(i,3)*(Ds_inv(j) - 2*delta(j,3)*Ds_inv(3)))...
                 -J*position_prime(3)^2/(4*pi)*(D2s_inv(i,j)-2*delta(j,3)*D2s_inv(i,3));   
        end
    end
   
    % 
    G = G_infinity + G_image;
    
end

