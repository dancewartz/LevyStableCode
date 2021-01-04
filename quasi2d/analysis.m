%% Load in data and set up parameters
clear
load('simulation_data.mat')
dt = params.dt;
tau = params.tau;
tlag = logspace(-2,2); 
nlag = unique(round(tlag/dt)); % converting to indicies
tlag = dt*nlag; % lag times we will measure at
lag_time_van_hove = tau; % lag time at which we compute the Van Hove Correlation Plot

msd = zeros(size(nlag)); % Mean squared displacement at each lag time
xis = zeros(size(nlag)); % Rescaling factor as a function of t ~ gamma
gams = zeros(size(nlag)); % Value of gamma at each lag time
alphas = zeros(size(nlag)); % Value of levu alpha at each lag time 

fontsize = 20;
%% Make all the Figures

% Plot for a sample trajectory
figure(1)
hold on
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
title('Sample Trajectory Projected Along x-Axis')
xlabel('Time (sec)')
ylabel('Displacement u_x (\mum)')
plot(ts(1:10:10000),u(1:10:10000),'LineWidth',2)

% Plot of MSD
figure(2)
hold on
set(gca,'FontSize',fontsize,'xscale','log','yscale','log')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Mean Squared Displacement (\mum^2)')
xline(tau,'--b','\tau','LineWidth',3);


% Van Hove plot
figure(3)
hold on
set(gca,'xscale','log','yscale','log','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
title(['Van Hove Correlation Plot: Lag Time = ' num2str(lag_time_van_hove) 's'])
xlabel('Displacement u_x (\mum)')
ylabel('Probability P(u_x)')
legend

% Plot of xi vs time
figure(4)
hold on
set(gca,'xscale','log','yscale','log','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Rescaling Parameter \xi')

% Plot of rescaled probabilities
figure(5)
hold on
set(gca,'xscale','log','yscale','log','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
title('Rescaled Van Hove Correlation Plots')
xlabel('Displacement u_x/\xi ')
ylabel('Probability P(u_x/\xi)')
legend

% Plot of alpha vs time
figure(6)
hold on
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Stability Parameter \alpha')

%% Actaully Calculating the curves to plot
i = 0;
for n = nlag
    i = i+1;
    
    du = u - circshift(u,n); % Difference in u at two different lag times
    du(1:n) = []; % First n values are useless... 
    
    % Calculating MSD
    msd(i) = mean(du.^2);
    
    xis(i) = exp(mean(log(abs(du))));
    
%     stable_fit = fitdist(du(1:100:end)','stable');
%     alphas(i) = stable_fit.alpha;
    
    if mod(i,10) == 0
        figure(5)
        histogram(abs(du)/xis(i),'Normalization','pdf','Displaystyle','stairs','lineWidth',3,'DisplayName',['tlag = ' num2str(dt*n)])
    end
    
end

figure(2)
g = tlag < 0.1;
msdfit = polyfit(log(tlag(g)),log(msd(g)),1);
pfit = plot(tlag,exp(polyval(msdfit,log(tlag))),'Linewidth',4,'DisplayName',['Power Law With Exponent ' num2str(msdfit(1))]);
pmsd = plot(tlag,msd,'ko-','LineWidth',4,'MarkerSize',10,'DisplayName','MSD Extracted From Simulation');
legend([pfit,pmsd],['Power Law With Exponent ' num2str(msdfit(1))],'MSD Extracted From Simulation');

figure(3)
hold on
n_van_hove = round(lag_time_van_hove/dt);
du = u - circshift(u,n_van_hove);
du(1:n_van_hove) = [];
stabledist = fitdist(du(1:100:end)','stable');
[N,edges] = histcounts(du,logspace(min(log10(du(du>0))),max(log10(du(du>0))),60),'Normalization','pdf');
%plot_range = (plot_range(1:end-1) + plot_range(2:end))/2;
plot(edges(1:end-1),N,'.','MarkerSize',12);
histogram(du,'Displaystyle','stairs','lineWidth',4,'Normalization','pdf','DisplayName','Trajectory Data')
plot(edges,pdf(stabledist,edges),'--r','LineWidth',2,'DisplayName',['Fit with \alpha = ' num2str(stabledist.alpha)]);

figure(4)
plot(tlag,xis,'LineWidth',3)

figure(6)
plot(tlag,alphas,'o-','LineWidth',4,'MarkerSize',10)

figure(7);
hold on
i = 0;
saffs = logspace(-3,5,9);
for saff_length = saffs
    i = i+1;
    load(['saff_length' num2str(saff_length) '.mat']);
    du = u - circshift(u,n_van_hove);
    du(1:n_van_hove) = [];
    stable_dist = fitdist(du(1:100:end)','stable');
    
    [N,edges] = histcounts(du,logspace(min(log10(du(du>0))),max(log10(du(du>0))),40),'Normalization','pdf');
    edges = (edges(1:end-1) + edges(2:end))/2;
    %edges = exp(edges(1:end-1));
    subplot(3,3,i)
    hold on
    title(['\kappa^{-1} = ' num2str(round(saffs(i),4)) ' (\mum); \alpha = ' num2str(stable_dist.alpha)])
    xlabel('Displacement u_x (\mum)')
    ylabel('P(u_x)')
    plot(edges,N,'.','MarkerSize',6)
    %histogram(du,'Displaystyle','stairs','lineWidth',4,'Normalization','pdf')
    set(gca,'xscale','log','yscale','log')
    
    %stablefit = fitdist(du(1:100:end)','stable');
    plot(edges,pdf(stable_dist,edges),'--','LineWidth',2)
    
end
saveas(gca,'Figures/saffvanhove.png')
