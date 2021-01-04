function [Fxs,ts] = tensor_sum(params)
%% Initializing parameters and setting up simulation
    tmax = params.tmax; % maximum time of simulation in seconds
    dt = params.dt; % time step in seconds
    Npoints = params.Npoints; % number of force dipoles
    L = params.L; % Box size, box spans -L/2 to L/2
    m_monopole = params.m_monopole; % fall off exponent for monopole, this is 1
    d = params.d; % Dipole size
    tau = params.tau; % sets our timescale -- the random switching time
    ts = 0:dt:tmax-dt;%linspace(0,tmax,round(tmax/dt));
    F = params.F;
    cutoff = d/2;
    % x-projection of trajectory
    Fxs = zeros(size(ts));
    
    % Time at which dipole will turn on, and how long that turn on event
    % will last
    nexton = -(tau)*log(rand(1,Npoints));
    durations = -(tau)*log(rand(1,Npoints));

    xs = L*(rand(1,Npoints)-0.5);
    ys = L*(rand(1,Npoints)-0.5);
    zs = L*(rand(1,Npoints)-0.5);
    rs = sqrt(xs.^2+ys.^2+zs.^2);
    
    % This doesnt allow dipoles to get too close to the origin
    while any(rs < cutoff)
        num = sum(rs<cutoff);
        xs(rs<cutoff) = L*(rand(1,num)-0.5);
        ys(rs<cutoff) = L*(rand(1,num)-0.5);
        zs(rs<cutoff) = L*(rand(1,num)-0.5);
        rs = sqrt(xs.^2+ys.^2+zs.^2);
    end
    
    phis = 2*pi*rand(1,Npoints);
    thetas = acos((2*rand(1,Npoints)-1));
    bx = sin(thetas).*cos(phis);
    by = sin(thetas).*sin(phis);
    bz = cos(thetas);
    uxc = zeros(1,Npoints);

    tic;

    %% Simulation Lööp
    for s = 1:length(ts)
        % find those points which have turned on this step
        turnon = (ts(s) - nexton <= dt) & (ts(s) - nexton > 0); % forces that are turning on this time step
        turnoff = (ts(s) >= nexton+durations); % forces that are turning off on this time step
        status = turnon - turnoff; % status change of force 1 for turned on, -1 for turned off 
        ind_change = find(status); % find indicies of all the forces that changed

        % Positions of the individual forces in the dipole
        xs1 = xs(ind_change) - bx(ind_change)*d/2;
        ys1 = ys(ind_change) - by(ind_change)*d/2;
        zs1 = zs(ind_change) - bz(ind_change)*d/2;
        rs1 = sqrt(xs1.^2+ys1.^2+zs1.^2);
        
        xs2 = xs(ind_change) + bx(ind_change)*d/2;
        ys2 = ys(ind_change) + by(ind_change)*d/2;
        zs2 = zs(ind_change) + bz(ind_change)*d/2;
        rs2 = sqrt(xs2.^2+ys2.^2+zs2.^2);

        % Calculate the change in the response due to the turn on or turn
        % off event
        uxc(ind_change) = uxc(ind_change) + status(ind_change).*((rs1.^(-m_monopole)).*(F*bx(ind_change) + (xs1.*xs1./rs1.^2).*F.*bx(ind_change) + (xs1.*ys1./rs1.^2).*F.*by(ind_change) + (xs1.*zs1./rs1.^2).*F.*bz(ind_change)) ...
            -(rs2.^(-m_monopole)).*(F*bx(ind_change) + (xs2.*xs2./rs2.^2).*(F.*bx(ind_change)) + (xs2.*ys2./rs2.^2).*(F.*by(ind_change)) + (xs2.*zs2./rs2.^2).*(F.*bz(ind_change))));

        % sum the contribution due to each force
        Fxs(s) = sum(uxc);


        % update all the state of those that switch
        % need to change x,y,z,r, Fmag, duration,nexton
        usi = size(find(turnoff)); % update size
        nexton(turnoff) = ts(s) + (-tau)*log(rand(usi)); % Finding the next time dipole will turn on
        durations(turnoff) = (-tau)*log(rand(usi)); % finding the duration of those events

        % Redistributing and reorienting dipoles that turned off
        xs(turnoff) = L*(rand(usi)-0.5);
        ys(turnoff) = L*(rand(usi)-0.5);
        zs(turnoff) = L*(rand(usi)-0.5);
        rs(turnoff) = sqrt(xs(turnoff).^2+ys(turnoff).^2+zs(turnoff).^2);
        
        % I believe we can copy and paste this here 
        while any(rs < cutoff)
            num = sum(rs<cutoff);
            xs(rs<cutoff) = L*(rand(1,num)-0.5);
            ys(rs<cutoff) = L*(rand(1,num)-0.5);
            zs(rs<cutoff) = L*(rand(1,num)-0.5);
            rs = sqrt(xs.^2+ys.^2+zs.^2);
        end
        
        phis = 2*pi*rand(usi);
        thetas = acos((2*rand(usi)-1));
        bx(turnoff) = sin(thetas).*cos(phis);
        by(turnoff) = sin(thetas).*sin(phis);
        bz(turnoff) = cos(thetas);
    
        % Track progress...
        if(rem(s,1e5)==0)
            fprintf('s = %d, %3.3g percent done \n',s,100*s/length(ts));
            toc;
        end


    end