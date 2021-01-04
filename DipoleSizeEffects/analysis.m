%% Making Figure
clear
fontsize = 20;
figure(1)
hold on
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Stability Parameter \alpha')
legend

dvals = 4*logspace(-3,-1,5);

for d = dvals([1 4 5])
    %% Getting parameters from simulation
    load(['simulation_d' num2str(d) '.mat'])
    u = Fxs;
    dt = params.dt;
    tau = params.tau;
    tlag = logspace(-2,2,25); 
    nlag = unique(round(tlag/dt)); % converting to indicies
    tlag = dt*nlag;

    alphas = zeros(size(nlag));
    
    %% Making Trajectory
    i = 0;
    for n = nlag
        i = i+1;
    
        du = u - circshift(u,n); % Difference in u at two different lag times
        du(1:n) = []; % First n values are useless... 

        stable_fit = fitdist(du(1:100:end)','stable');
        alphas(i) = stable_fit.alpha;
    end
    
    figure(1)
    plot(tlag,alphas,'o-','LineWidth',4,'MarkerSize',10,'DisplayName',['d = ' num2str(round(d,3))])
    
end

saveas(gca,'Figures/dipolesize.png')