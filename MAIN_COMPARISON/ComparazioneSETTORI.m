clc; clearvars; close all
%% addapath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_OT')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
%% load results
load RisultatiOT_SETTORI_ITC_VS_ENERGY_LAG_PESATI_REVISION.mat
P_forecast_containerR = P_forecast_container;
load RisultatiOT_SETTORI_ITC_VS_ENERGY_CONT_PESATI_REVISION.mat
P_forecast_containerC = P_forecast_container;
load RisultatiOT_PESATI_REVISION.mat
load RisultatiOT_SETTORI_ITC_VS_ENERGY_LAG_PESATI_REVISION_BECH.mat
load('DatiPerOT_SETTORI_PESATI_REVISION.mat','reti','NodiTutti')
anni = 2003:2017;
reti = reti(:,:,:,8);
%% caption
leg = num2str(anni(size(TestLag_container,2)+1:size(TestLag_container,1))');
leg1 = ["Value-";"Optimum-"];
leg1 = repmat(leg1,length(leg),1);
legx = string(repelem(leg,2,1));
leg2 = strcat(leg1,legx);
leg = string(leg);

%% show optimal lag
posmin=zeros(length(size(TestLag_container,2)+1:size(TestLag_container,1)),1);
x = 1;
for k = size(TestLag_container,2)+1:size(TestLag_container,1)
    [a,b]=min(TestLag_container(k,:));
    posmin(x)=b;
    x = x+1;
end

%% plot forecasted network
Forecast = P_forecast_container(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
ForecastR = P_forecast_containerR(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
ForecastC = P_forecast_containerC(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);


ForecastOT = zeros(size(Forecast,1),size(Forecast,2));
ForecastOT_R = zeros(size(Forecast,1),size(Forecast,2));
ForecastOT_C = zeros(size(Forecast,1),size(Forecast,2));
Reali = zeros(size(Forecast,1),size(Forecast,2));

for k = 1:size(Forecast,2)
    ForecastOT(:,k)=Forecast(:,k,posmin(k));
    ForecastOT_R(:,k)=ForecastR(:,k,posmin(k)); 
    ForecastOT_C(:,k)=ForecastC(:,k,posmin(k)); 
end

%% scores
rmseModel = zeros(size(Forecast,3)-1,3);
rmseModelR = zeros(size(Forecast,3)-1,3);
rmseModelC = zeros(size(Forecast,3)-1,3);

for z = 1:size(Forecast,2)-1
    Reali(:,z)=reshape(reti(:,:,5+z),size(reti,1)*size(reti,2),1);
    Reali(:,z) = Reali(:,z)/sum(Reali(:,z));
    
    rmseModel(z,1) = sqrt(mean((Reali(:,z+1)-ForecastOT(:,z)).^2));
    rmseModel(z,2) = (mean((Reali(:,z+1)-ForecastOT(:,z)).^2));
    rmseModel(z,3) = (mean(abs((ForecastOT(:,z)-Reali(:,z+1)))));

    rmseModelR(z,1) = sqrt(mean((Reali(:,z+1)-ForecastOT_R(:,z)).^2));
    rmseModelR(z,2) = (mean((Reali(:,z+1)-ForecastOT_R(:,z)).^2));
    rmseModelR(z,3) = (mean(abs((ForecastOT_R(:,z)-Reali(:,z+1)))));

    rmseModelC(z,1) = sqrt(mean((Reali(:,z+1)-ForecastOT_C(:,z)).^2));
    rmseModelC(z,2) = (mean((Reali(:,z+1)-ForecastOT_C(:,z)).^2));
    rmseModelC(z,3) = (mean(abs((ForecastOT_C(:,z)-Reali(:,z+1)))));

end

figure
tiledlayout(2,1, 'TileSpacing', 'compact')

nexttile
bar([rmseModel(:,1),rmseModelR(:,1),rmseModelC(:,1)])
title('RMSE')
xticks(1:length(rmseModel))
xticklabels(leg(2:end))
set(gca,'FontSize',10,'fontWeight','bold')
axis tight
grid on


nexttile
bar([rmseModel(:,2),rmseModelR(:,2),rmseModelC(:,2)])
title('MSE')
xticks(1:length(rmseModel))
xticklabels(leg(2:end))
set(gca,'FontSize',10,'fontWeight','bold')
axis tight
grid on

lgd = legend('Wasserstein','Wasserstein Reg.','Wasserstein Cont.');
title(lgd,'Models')
lgd.Layout.Tile = 'east';