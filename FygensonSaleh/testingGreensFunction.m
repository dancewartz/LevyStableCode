xprime = 0;
yprime = 0;
zprime = 1;
position_prime = [xprime yprime zprime]';
params.poisson_ratio = 1/2;
for theta = [0 pi/4 pi/2]
    bx = sin(theta);
    by = cos(theta);
    b = [bx 0 by]';
    d = 0.01;
    free_boundary = false;

    xvals = linspace(-1,1,11);
    zvals = linspace(0,2,11);
    [X,Z]= meshgrid(xvals,zvals);
    U = zeros(size(X));
    V = zeros(size(Z));
    for i = 1:length(xvals)
        for j = 1:length(zvals)
            x = X(i,j);
            z = Z(i,j);
            position = [x 0 z]';

            [G_minus,G_image_minus] = Greens_Function(position,position_prime-b*d/2,params,free_boundary);

            [G_plus,G_image_plus] = Greens_Function(position,position_prime+b*d/2,params,free_boundary);
            
            response = (G_image_minus - G_image_plus)*b; % think about this a lot...
            U(i,j) = response(1);
            V(i,j) = response(3);
        end
    end


    figure
    hold on
    yline(0,'LineWidth',3);
    title(['\theta = ' num2str(theta) ', Free boundary = '  num2str(free_boundary)])
    R = sqrt(U.^2 + V.^2);
    Un = U./R;
    Vn = V./R;
    quiver(X,Z,Un,Vn);
    quiver(xprime-bx*d/2,zprime-by*d/2,bx/5,by/5,'linewidth',5);
    quiver(xprime+bx*d/2,zprime+by*d/2,-bx/5,-by/5,'Linewidth',5);
    
end