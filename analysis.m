%% Load in data and set up parameters
clear
load('simulation_data.mat')
u = Fxs;
%% CHANGE THIS TO MAKE ANALYSIS FOR VISCOELASTIC CASE
convolve = true;
%%
if convolve==true 
    u = convolve_fft(ts,u,params); % Takes data from elastic -> viscoelastic case
end
u = real(u);
dt = params.dt;
tau = params.tau;
tlag = logspace(-2,2); 
nlag = unique(round(tlag/dt)); % converting to indicies
tlag = dt*nlag; % lag times we will measure at
lag_time_van_hove = 2*tau; % lag time at which we compute the Van Hove Correlation Plot

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
ylabel('Probability Density P(u_x)')
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
ylabel('Probability Density P(u_x/\xi)')
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
    
    %du = (abs(du) - mean(abs(du)))/sqrt(var(abs(du)));
    
    % Calculating MSD
    msd(i) = mean(du.^2);
    
    xis(i) = exp(mean(log(abs(du))));
    
    stable_fit = fitdist(du(1:100:end)','stable');
    alphas(i) = stable_fit.alpha;
    
    if mod(i,10) == 0
        figure(5)
        histogram(du/xis(i),'Normalization','pdf','Displaystyle','stairs','lineWidth',3,'DisplayName',['lag time = ' num2str(dt*n) 's'])
    end
    
end

figure(1)
if convolve == false
    saveas(gca,'Figures/sample_trajectory.png')
else
    saveas(gca,'Figures/visco_sample_trajectory.png')
end

figure(2)
g = tlag < 0.1;
msdfit = polyfit(log(tlag(g)),log(msd(g)),1);
pmsd = plot(tlag,msd,'ko-','LineWidth',4,'MarkerSize',10,'DisplayName','MSD Extracted From Simulation');
pfit = plot(tlag,exp(polyval(msdfit,log(tlag))),'Linewidth',4,'DisplayName',['Power Law With Exponent ' num2str(msdfit(1))]);
legend([pfit,pmsd],['Best fit to t^\Delta with \Delta = ' num2str(msdfit(1))],'MSD Extracted From Simulation');
if convolve == false
    saveas(gca,'Figures/MSD.png')
else
    saveas(gca,'Figures/viscoMSD.png')
end

figure(3)
hold on
n_van_hove = round(lag_time_van_hove/dt);
du = u - circshift(u,n_van_hove);
stabledist = fitdist(du(1:1000:end)','stable');
[~,edges] = histcounts(log10(du(du>0)));
plot_range = logspace(-4,1,300);
xlim([plot_range(1) plot_range(end)])
histogram(du,plot_range,'Displaystyle','stairs','lineWidth',3,'Normalization','pdf','DisplayName','Trajectory Data')
plot(plot_range,pdf(stabledist,plot_range),'--r','LineWidth',4,'DisplayName',['Fit with \alpha = ' num2str(stabledist.alpha)]);
if convolve == false
    saveas(gca,'Figures/vanhovelag.png')
else
    saveas(gca,'Figures/viscovanhovelag.png')
end

figure(4)
plot(tlag,xis,'LineWidth',3)
if convolve == false
    saveas(gca,'Figures/xivstime.png')
else
    saveas(gca,'Figures/viscoxivstime.png')
end

figure(5)
if convolve == false
    saveas(gca,'Figures/rescaledprobs.png')
else
    saveas(gca,'Figures/viscorescaledprobs.png')
end

figure(6)
plot(tlag,alphas,'o-','LineWidth',4,'MarkerSize',10)
if convolve == false
    saveas(gca,'Figures/alphavstime.png')
else
    saveas(gca,'Figures/viscoalphavstime.png')
end
