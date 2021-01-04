%% Generate "force" trajectory
tmax = 1e4; %1e4 usually enough
dt = 0.01;  % default 0.01
rho = 10^(-1);    % number density of force dipoles (fillaments) per cubic micron; default 1
L = 10;      % system size in microns
G0 = 1; 
Npoints = round(rho*L^3);
m_monopole = 1;  % scaling for the monopole force response - should always be 1
d = 0.4;     % physical size of dipoles in micronsf
beta = 0.17; % rheological exponent

tau = 1/2; % time scale in seconds
F = 0.25; % Force strength in picoNewtons

A = 0.66; % Fit parameters from Hoffman Crocker
B = 0.009; % ^ 

% Putting all these parameters in the params struct to store 
params = struct;
params.tmax = tmax;
params.dt = dt;
params.rho = rho;
params.L = L;
params.Npoints = Npoints;
params.d = d;
params.beta = beta;
params.tau = tau;
params.F = F;
params.m_monopole = m_monopole;
params.A = A;
params.B = B;
params.G0 = G0;

% % Parameterrs from old analysis, not necessary for generating trajectories
% thin_by = 1;
% burn = round(10/dt);  % burn-in time

[Fxs,ts] = tensor_sum(params);

save('simulation_data.mat','Fxs','ts','params')
%% Convolve with kernel to find displacement comment to do elastic case
% u = convolve_fft(ts(1:thin_by:end),Fxs(1:thin_by:end),params);
% u = real(u); % UHHHHHH Seems necessary because sometimes u has an imaginary part?? Idk 

% %% Fit
% thin_again = 100;
% stabledist = fitdist(u(burn/thin_by:thin_again:end).','stable')
% histogram(u,'normalization','pdf','linestyle','none');
% set(gca,'xscale','log','yscale','log');
% xl = xlim;
% %hold on
% xf = logspace(log10(xl(1)),log10(xl(2)),100);
% clf
% histogram(u(burn/thin_by:thin_again:end),'normalization','pdf','linestyle','none','binedges',real(xf));
% hold on
% %xf = logspace(log10(xl(1)),log10(xl(2)),100);
% stable_theory = makedist('stable');
% stable_theory.alpha = 3/(m_monopole+1);
% stable_theory.beta = 0;
% stable_theory.gam = stabledist.gam;
% stable_theory.delta = 0;
% 
% plot(xf,pdf(stabledist,xf),'LineWidth',4);
% set(gca,'xscale','log','yscale','log');
% plot(xf,pdf(stable_theory,xf),'--','LineWidth',3);
% 
% %% Find MSD + exponent
% maxlagtime = 100; % going much beyond 3-5 is going to lead to turnover
% %[msd,tlag] = mean_squared_displacement(ts(burn:thin_by:end),u(burn/thin_by:end),maxlagtime,100);
% figure
% numpoints = 100;
% dodist = false;
% [Ps,deltas,tlag,xilag,msd,fds,alphas] = van_hoves(ts,u(burn/thin_by:end),maxlagtime,numpoints,dodist);
% g = (tlag<0.1)&(tlag>0);
% p = polyfit(log(tlag(g)),log(msd(g)),1);
% 
% plot(tlag,msd,'ko-','LineWidth',4,'MarkerSize',10);
% hold on
% plot(tlag,exp(polyval(p,log(tlag))),'LineWidth',4);
% set(gca,'xscale','log','yscale','log');
% xlabel('Lag time'); ylabel('MSD');
% set(gca,'FontSize',24)
% title(sprintf('MSD; fit exponent = %3.3g',p(1)));
% axis tight;
% yl = ylim;
% plot(ones(1,100),linspace(yl(1),yl(2),100),'--','color',[0.8 0.8 0.8],'LineWidth',2)
% drawnow
% figure
% clf;
% hold on;
% c = 1;
% hs = {};
% legs = {};
% for j = 7:10:length(Ps) %skip lag time zero
%     hs{c} = plot(deltas{j}/xilag(j),Ps{j}*xilag(j),'LineWidth',1);
%     legs{c} = sprintf('Lagtime %3.3g',tlag(j));
%     c = c + 1;
% end
% set(gca,'xscale','log','yscale','log');
% xlabel('\Delta x / \xi'); ylabel('P(\Delta x / \xi)');
% legend(legs)
% axis tight
% set(gca,'FontSize',24)
% 
% drawnow
% 
% %% Find alphas:
% figure
% numpoints2 = 10;
% maxlagtime2 = 100;
% dodist = true;
% [Ps2,deltas2,tlag2,xilag2,msd2,fds2,alphas2,betas2,gams2,deltasd] = van_hoves(ts,u(burn:end),maxlagtime2,numpoints2,dodist);
% plot(tlag2,alphas2);
% set(gca,'xscale','log')
% xlabel('Lag time'); ylabel('Levy \alpha exponent')
% 
% %%
% figure
% Ndist = length(Ps2); %Ndr = ceil(sqrt(length(Ps2)));
% Ndr = 3;
% for k = 2:Ndist-1
%     subplot(Ndr,Ndr,k-1)
%     levydist = makedist('stable','beta',0,'alpha',fds2{k}.alpha,'delta',0,'gam',fds2{k}.gam);
%     plot(deltas2{k},Ps2{k});
%     hold on
%     plot(deltas2{k},2*pdf(levydist,deltas2{k})); %factor of 2 here is from symmetry
%     set(gca,'xscale','log','yscale','log')
%     title(sprintf('Lag time %3.3g',tlag2(k)));
%     %axis tight
%     xlim([1e-2 1e3])
%     ylim([1e-8 100])
% end