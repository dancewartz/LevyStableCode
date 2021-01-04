function [Fxs,ts] =  tensor_sum_2d(params)
    tmax = params.tmax;
    dt = params.dt;
    Npoints = params.Npoints; 
    L = params.L;
    m_monopole = params.m_monopole;
    d = params.d;
    tau = params.tau; % sets our timescale -- the random switching time
    ts = linspace(0,tmax,round(tmax/dt));
    F = params.F;
    saff_length = params.saff_length;
    
    Fxs = zeros(size(ts));
    
    nexton = -tau*log(rand(1,Npoints));
    durations = -tau*log(rand(1,Npoints));
    
    xs = L*(rand(1,Npoints)-0.5);
    ys = L*(rand(1,Npoints)-0.5);
    rs = sqrt(xs.^2 + ys.^2);
    
    phis = 2*pi*rand(1,Npoints);
    
    bx = cos(phis);
    by = sin(phis);
    
    uxc = zeros(1,Npoints);
    
    tic;
    for s = 1:length(ts)
        turnon = (ts(s) - nexton <= dt) & (ts(s) - nexton > 0); % forces that are turning on this time step
        turnoff = (ts(s) >= nexton+durations); % forces that are turning off on this time step
        status = turnon - turnoff; % status change of force 1 for turned on, -1 for turned off 
        ind_change = find(status); % find indicies of all the forces that changed

        
        xs1 = xs(ind_change) - bx(ind_change)*d/2;
        ys1 = ys(ind_change) - by(ind_change)*d/2;
        rs1 = sqrt(xs1.^2+ys1.^2);
        
        xs2 = xs(ind_change) + bx(ind_change)*d/2;
        ys2 = ys(ind_change) + by(ind_change)*d/2;
        rs2 = sqrt(xs2.^2+ys2.^2);

        kr1 = rs1./saff_length;
        kr2 = rs2./saff_length;
        
        f1 = StruveH0(kr1)-StruveH1(kr1)./kr1-1/2*(bessely(0,kr1)-bessely(2,kr1))+2./(pi*kr1.^2);
        f2 = StruveH0(kr2)-StruveH1(kr2)./kr2-1/2*(bessely(0,kr2)-bessely(2,kr2))+2./(pi*kr2.^2);
        
        g1 = StruveH0(kr1) - 2*StruveH1(kr1)./kr1 + bessely(2,kr1) + 4./(pi*kr1.^2);
        g2 = StruveH0(kr2) - 2*StruveH1(kr2)./kr2 + bessely(2,kr2) + 4./(pi*kr2.^2);
        
        uxc(ind_change) = uxc(ind_change) + F*status(ind_change).*(f1.*bx(ind_change) - g1.*(xs1.*xs1./rs1.^2.*bx(ind_change) + xs1.*ys1./rs1.^2.*by(ind_change))) ... 
            -F*status(ind_change).*(f2.*bx(ind_change) - g2.*(xs2.*xs2./rs2.^2.*bx(ind_change) + xs2.*ys2./rs2.^2.*by(ind_change)));

        Fxs(s) = sum(uxc);
        
        if isnan(Fxs(s))
            s
            break
        end
        
        % update all the state of those that switch
        % need to change x,y,z,r, Fmag, duration,nexton
        usi = size(nexton(turnoff)); % update size
        nexton(turnoff) = ts(s) + (-tau)*log(rand(usi));
        durations(turnoff) = (-tau)*log(rand(usi));
        
        xs(turnoff) = L*(rand(usi)-0.5);
        ys(turnoff) = L*(rand(usi)-0.5);
    
        rs(turnoff) = sqrt(xs(turnoff).^2+ys(turnoff).^2);
        phis = 2*pi*rand(usi);
        bx(turnoff) = cos(phis);
        by(turnoff) = sin(phis);

        if(rem(s,1e5)==0)
            fprintf('s = %d, %3.3g percent done \n',s,100*s/length(ts));
            toc;
        end
    
        
    end
end