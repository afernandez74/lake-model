figure

plot(dates_clim,daily_LDv,'k-.','LineWidth',1)
ylabel('Lake Depth (cm)', 'FontSize', 20)
grid on
title('20th Century Modeled and Observed Lake Depth', 'FontSize', 25)
ax=gca;
ax.FontSize=20;
% ax.XLim=[1 12];
% ax.XTickLabel={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ax.XMinorTick = 'off';
hold
plot(dates_obs,LLobs(:,4),'k-','LineWidth',1.3)
plot(dates_sat,LLsat(:,4),'ks','Markersize',10)
plot(dates_sat_uncert,LLsat_uncert(:,4),'.','Markersize',10)
legend({'Modeled Daily Water Depth','Measured Lake Depth (levelogger)','Estimated Lake Level (satelite images)'},...
    'FontSize', 15,'Location','southeast')
% xlim([datetime(2005,7,1) dates_obs(end)])

%%
figure

% plot(dates_clim,daily_d18O_aragonite,'k--','LineWidth',0.1)
hold
plot(dates_core,movmean(d18O_core(:,2),[10 0]),'b.-','MarkerSize',20,'LineWidth',0.3)
plot(dates_clim_years,year_summer_arag,'r--','LineWidth',0.1)
plot(dates_clim_years,movmean(year_summer_arag,[10 0]),'r--','LineWidth',3)
legend({'Modeled Daily Aragonite \delta ^{18}O','Measured Sediment \delta ^{18}O','Averaged modeled \delta ^{18}O for summer','20-year moving average modeled summer \delta ^{18} O'},...
    'FontSize', 15,'Location','southeast')
ax=gca;
ax.FontSize=20;
ax.XMinorTick = 'off';


%%
figure 
plot(dates_clim,daily_dl,'k-','LineWidth',0.2)
hold
plot(dates_samples,d18O_samples(:,4),'ks','Markersize',10)
ax=gca;
ax.FontSize=20;
ax.XMinorTick = 'off';
title('20th Century Modeled and Observed Surface Lake Isotopic Composition', 'FontSize', 25)
ax.YLabel.String='Water VSMOW \delta ^{18} O ';
% xlim([datetime(2005,7,1) dates_obs(end)])

%%
figure
plot(dates_clim,movmean(daily_LDv,[365/2 365/2]),'k-.','LineWidth',1)
ylabel('Lake Depth (cm)', 'FontSize', 20)
grid on
title('20th Century Modeled and Observed Lake Depth', 'FontSize', 25)
ax=gca;
ax.FontSize=20;
% ax.XLim=[1 12];
% ax.XTickLabel={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ax.XMinorTick = 'off';
hold
plot(dates_obs,movmean(LLobs(:,4),[365/2 365/2]),'k-','LineWidth',1.3)
plot(dates_sat,LLsat(:,4),'ks','Markersize',10)
plot(dates_sat,LLsat(:,4),'.','Markersize',10)

legend({'Modeled Daily Water Depth','Measured Lake Depth (levelogger)','Estimated Lake Level (satelite images)'},...
    'FontSize', 15,'Location','southeast')
xlim([datetime(2005,7,1) dates_obs(end)])