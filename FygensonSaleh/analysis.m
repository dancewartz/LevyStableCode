%% Defining Parameters
% Importing simulation data and parameters from the simulation
load('sim.mat')
dt = params.dt;
tau = params.tau;
tlag = logspace(-2,2,100);
nlag = unique(round(tlag/dt));
tlag = nlag*dt; 

msd = zeros(size(tlag));
gammas = zeros(size(tlag));
alphas = zeros(size(tlag));
xis = zeros(size(tlag));

%% Setting up figures

% Plot of the running value of alpha vs lag time
figure(1)
hold on
set(gca,'FontSize',24)
set(gca,'xscale','log')
xline(tau,'--r',{'\tau'},'DisplayName','t = \tau');
xlabel('lag time (sec)')
ylabel('\alpha')
legend

% Plog of the running value of gamma vs lag time
figure(2) 
hold on
set(gca,'FontSize',24)
legend;
xlabel('lag time')
ylabel('\gamma')
set(gca,'FontSize',24)
set(gca,'Xscale','log','yscale','log')
xline(tau,'--r',{'\tau'},'DisplayName','t = \tau','linewidth',3);

% Plot of the mean squared displacement vs lag time
figure(3)
hold on;
set(gca,'xscale','log','yscale','log');
xlabel('Lag Time (sec)'); ylabel('MSD (\mum^2)');
set(gca,'FontSize',24)
axis tight;
yl = ylim;
xline(tau,'--r',{'\tau'},'DisplayName','t = \tau','linewidth',3);
legend


% Van Hove Correlation Plot
figure(4)
hold on
set(gca,'FontSize',24)
title('Van Hove Correlation Plot')
xlabel('Displacement u (\mum)')
ylabel('Probability Density P(u)')
set(gca,'yscale','log','xscale','linear')
legend

% Plot of rescaled Van Hoves,
figure(5)
hold on
set(gca,'FontSize',24)
xlabel('u')
ylabel('P(u)')
title('Rescaled Van Hoves')
set(gca,'xscale','log','yscale','log')
legend

% My choice for the fit distribution is ARBITRARY
u = Fxs;
n = nlag(10); % This could be any nlag, I chose the 10th because it is the first to show up in figure 5
du = u - circshift(u,n);
du(1:n) = [];
xi = exp(mean(log(abs(du))));
stabledist = fitdist(du(1:1000:end)'/xi,'Stable');

stable_theory = makedist('Stable');
stable_theory.alpha = stabledist.alpha;
stable_theory.beta = 0;
stable_theory.gam = stabledist.gam;
stable_theory.delta = 0;

% plot(logspace(-1.5,1.5),2*pdf(stable_theory,logspace(-1.5,1.5)),'LineWidth',2,'DisplayName',['Theory with \alpha = ' num2str(stabledist.alpha)]);

% Plot of the running xi vs time
figure(6)
hold on
set(gca,'FontSize',24)
xlabel('lag time');
ylabel('\xi');
xline(tau,'--r',{'\tau'},'DisplayName','t = \tau','LineWidth',3);
set(gca,'xscale','log','yscale','log')
legend

% Plot of gamma vs xi, looking for a linear relationship
figure(7)
hold on
set(gca,'FontSize',24)
xlabel('\gamma');
ylabel('\xi');
title('\xi vs gamma')
%xline(tau,'--r',{'\tau'},'DisplayName','t = \tau');
set(gca,'xscale','linear','yscale','linear')
legend
ind_tau = find(tlag == tau);


% Loading in data
results = [Fxs;Fys;Fzs];

%% Analysis of Data

for q = [1 2 3]
    % Loading in data for one direction
    u = results(q,1:end);
    if q == 1
        label = 'x';
    elseif q == 2
        label = 'y';
    elseif q == 3
        label = 'z';
    end

    % convolving with viscoelastic kernel 
%     u = convolve_fft(ts,u,params); % commenting this line out will give the elastic results
    u = real(u);
    i = 0;
    
    for n = nlag
        % calculating differences in displacement
        i = i+1;
        du = u - circshift(u,n);
        du(1:n) = [];
        msd(i) = mean(du.^2);
        
        % This is a rescaling factor from Yu Shi's work. Seems to
        % accomplish the goal
        xi = exp(mean(log(abs(du))));
        xis(i) = xi;

        if mod(i,20) == 0 && q==3
            figure(5)
            histogram(abs(du)/xi,'Normalization','pdf','DisplayStyle','Stairs','Linewidth',3,'DisplayName',[label ', tlag = ' num2str(tlag(i))])
        end
    % This is commented because it is slow uncomment to calculate alpha nad
    % gamma vs lag time 
    
%         stabledist = fitdist(du(1:1000:end)','stable');
%         gammas(i) = stabledist.gam;
%         alphas(i) = stabledist.alpha;
    end

    % Plotting the running alpha vs lag time
    figure(1)
    plot(tlag,alphas,'Displayname',label)
    
    % Plotting the gammas and fitting it to a power law
    g = (tlag<0.1)&(tlag>0);
    figure(2) 
    plot(tlag,gammas,'o-','LineWidth',4,'MarkerSize',10,'Displayname',label)
    p = polyfit(log10(tlag(g)),log10(gammas(g)),1);
    plot(tlag,10.^(polyval(p,log10(tlag))),'Displayname',['slope = '  num2str(p(1))]);


    % Plotting the MSD and fitting it to a power law 
    p = polyfit(log(tlag(g)),log(msd(g)),1);
    figure(3)
    if q == 1
        plot(tlag,msd,'o-','LineWidth',4,'MarkerSize',10,'DisplayName',label);
    elseif q == 2
        plot(tlag,msd,'+-','LineWidth',4,'MarkerSize',10,'DisplayName',label);
    elseif q == 3
        plot(tlag,msd,'*-','LineWidth',4,'MarkerSize',10,'DisplayName',label);
    end
    plot(tlag,exp(polyval(p,log(tlag))),'LineWidth',4,'DisplayName',['slope = ' num2str(round(p(1),2))]);


    % Making a characteristic Van Hove, just like in Fygenson Saleh
    tlag2 = 2*tau;

    du = u - circshift(u,round(tlag2/dt));
    du(1:n) = [];
    dist = fitdist(du(1:100:end)','Stable');
    figure(4)
    title(['Van Hove; Lag Time is ' num2str(tlag2) 's'])
    histogram(du,'Normalization','pdf','Displaystyle','stairs','DisplayName',[label ': \alpha = ' num2str(dist.alpha)],'LineWidth',4)
    [N,edges] = histcounts(du,'Normalization','pdf');
    plot(edges,pdf(dist,edges),'--','LineWidth',3,'DisplayName',[label ' : Best Fit'])
    
    % Plotting xi vs lag time, looks like gamma....
    figure(6)
    plot(tlag,xis,'o-','LineWidth',4,'MarkerSize',10,'DisplayName',['\xi, ' label]);
    xilinearfit = polyfit(log10(tlag(g)),log10(xis(g)),1);
    plot(tlag,10.^(polyval(xilinearfit,log10(tlag))),'Displayname',['fit with power law exponent ' num2str(xilinearfit(1))],'Linewidth',3);
    
    % Plots gamma vs xi for only the x direction, change q == 2 for y and q
    % == 3 for z
    if q == 1
        figure(7)
        scatter(gammas,xis)
    end
    
end

figure(1)
saveas(gca,'figures/fygensonalphavslagtime.png')
figure(2)
saveas(gca,'figures/fygensongammavslagtime.png')
figure(3)
saveas(gca,'figures/fygensonMSD.png')
figure(4)
saveas(gca,'figures/fygensonvanhovelag.png')
figure(5)
saveas(gca,'figures/fygensonrescaledprobs.png')
figure(6)
saveas(gca,'figures/fygensonxivslagtime.png')
figure(7)
saveas(gca,'figures/fygensongammasvsxis.png')

