figure;
hold on;
title('Experimental Data');
xlabel('u (nm)');
ylabel('Probability Density P(u)');
set(gca,'xscale','log','yscale','log','FontSize',15);
legend

for q = 1
    if q == 1
        T = readtable('FygensonSalehx.csv');
        label = 'x';
    elseif q == 2
        T = readtable('FygensonSalehy.csv');
        label = 'y';
    else
        T = readtable('FygensonSalehz.csv');
        label = 'z';
    end
    x = table2array(T(1:end,1));
    px = table2array(T(1:end,2));

    g = (x > 10^1.1) & (x<10^2.2);

    fit = polyfit(log10(x(g)),log10(px(g)),1);

    plot_range = logspace(1.1,2.2);

    scatter(x(x>0),px(x>0),'filled','DisplayName',['Experimental Distribution of u_' label])
    plot(plot_range,1.5*plot_range.^fit(1)*10^fit(2),'--','LineWidth',4, 'DisplayName',['u^{-(1+\alpha)}, \alpha = ' num2str(-fit(1)-1)])
end

saveas(gca,'figures/ExperimentalVanHov.png')