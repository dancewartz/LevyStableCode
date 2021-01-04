function [Ps,deltas,tlag,xilag,msd,fds,alphas,betas,gams,deltasd] = van_hoves(ts,us,maxlagtime,numpoints,dodist)
dt = ts(2)-ts(1);
%tlag = (0:maxn)*dt;
tlag_est = logspace(log10(dt),log10(maxlagtime),numpoints);
nlag = [0 unique(round(tlag_est/dt))];
tlag = nlag*dt;
Ps = cell(size(tlag));
deltas = cell(size(tlag));
msd = NaN*ones(size(tlag));
xilag = NaN*ones(size(tlag)); 
fds = cell(size(tlag));
alphas = NaN*ones(size(tlag));
for j = 1:length(nlag)
    du = us - circshift(us,nlag(j));
    du(1:nlag(j)) = [];
    msd(j) = mean(du.^2);
    if(dodist)
        if(j>1)
            fds{j} = fitdist(du(1:100:end).','stable');
            alphas(j) = fds{j}.alpha;
            betas(j) = fds{j}.beta;
            gams(j) = fds{j}.gam;
            deltasd(j) = fds{j}.delta;
            fprintf('Time lag %3.3g alpha = %3.3g \n',tlag(j),alphas(j));
        end
    end
    
    clf
    H1 = histogram(abs(du),'normalization','pdf','linestyle','none');
    xilag(j) = exp(mean(log(abs(du))));
    set(gca,'xscale','log','yscale','log');
    xl = xlim;
    xf = logspace(log10(xl(1)),log10(xl(2)),100);
    H = histogram(abs(du),'normalization','pdf','binedges',xf,'linestyle','none');
    deltas{j} = (H.BinEdges(1:end-1)+H.BinEdges(2:end))/2;
    Ps{j} = H.Values;
    
end