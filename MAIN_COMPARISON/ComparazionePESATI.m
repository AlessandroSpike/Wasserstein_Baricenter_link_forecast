clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_OT')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_ALTERNATIVES')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
%% load results
load RisultatiOT_PESATI_REVISION_GDP_TRADE.mat
P_forecast_containerGDP = P_forecast_container;
load RisultatiOT_PESATI_REVISION_TRADE.mat
load RisultatiNMF_PESATI_REVISION_TRADE.mat
load('DatiPerOT_PESATI_REVISION_TRADE.mat','reti','NodiTutti')
anni = 2003:2021;
%% colori
rgb = vals2colormap(size(TestLag_container,2)+1:size(TestLag_container,1), 'turbo');
%% show optimal lag
leg = num2str(anni(size(TestLag_container,2)+1:size(TestLag_container,1))');
leg1 = ["Value-";"Optimum-"];
leg1 = repmat(leg1,length(leg),1);
legx = string(repelem(leg,2,1));
leg2 = strcat(leg1,legx);
leg = string(leg);
posmin=zeros(length(size(TestLag_container,2)+1:size(TestLag_container,1)),1);
figure
x = 1;
for k = size(TestLag_container,2)+1:size(TestLag_container,1)
    plot(TestLag_container(k,:),'color',rgb(x,:),'linewidth',2)
    hold on
    [a,b]=min(TestLag_container(k,:));
    posmin(x)=b;
    plot(b, a, 'o','color',rgb(x,:), 'MarkerSize', 10);
    x = x+1;
end
axis tight
grid on
ylabel('Wasserstein Distance')
xlabel('Lag')
xticks(1:size(TestLag_container,2))
xticklabels(1:size(TestLag_container,2))
title('Model Selection')
set(gca,'fontsize',10,'Fontweight','bold')
legend(leg2,'location','eastoutside');

save('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_ALTERNATIVES\ModelSelectionPESATI_REVISION.mat','posmin')
%% plot forecasted network
Forecast = P_forecast_container(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
ForecastGDP = P_forecast_containerGDP(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastk = P_forecast_containerK_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastk2 = P_forecast_containerK2_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastnmf = P_forecast_containerNMF_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastnmf2 = P_forecast_containerNMF2_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastrwr = P_forecast_containerRWR_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastrwr2 = P_forecast_containerRWR2_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastaa = P_forecast_containerAA_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastaa2 = P_forecast_containerAA2_S(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);


ForecastOT = zeros(size(Forecast,1),size(Forecast,2));
ForecastOT_GDP = zeros(size(Forecast,1),size(Forecast,2));
ForecastK = zeros(size(Forecast,1),size(Forecast,2));
ForecastK2 = zeros(size(Forecast,1),size(Forecast,2));
ForecastNMF = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastNMF2 = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastRWR = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastRWR2 = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastAA = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastAA2 = zeros(size(Forecastnmf,1),size(Forecast,2));

Reali = zeros(size(Forecastnmf,1),size(Forecast,2));

for k = 1:size(Forecast,2)
    ForecastOT(:,k)=Forecast(:,k,posmin(k));
    ForecastOT_GDP(:,k)=ForecastGDP(:,k,posmin(k));
    ForecastNMF(:,k)=Forecastnmf(:,k,posmin(k));
    ForecastNMF2(:,k)=Forecastnmf2(:,k,posmin(k));
    ForecastK(:,k)=Forecastk(:,k,posmin(k));
    ForecastK2(:,k)=Forecastk2(:,k,posmin(k));
    ForecastRWR(:,k)=Forecastrwr(:,k,posmin(k));
    ForecastRWR2(:,k)=Forecastrwr2(:,k,posmin(k));
    ForecastAA(:,k)=Forecastaa(:,k,posmin(k));
    ForecastAA2(:,k)=Forecastaa2(:,k,posmin(k));
end

%% scores
rmseModel = zeros(size(Forecast,3)-1,3);
rmseModelGDP = zeros(size(Forecast,3)-1,3);
rmseNMF = zeros(size(Forecast,3)-1,3);
rmseNFM2 = zeros(size(Forecast,3)-1,3);
rmseK = zeros(size(Forecast,3)-1,3);
rmseK2 = zeros(size(Forecast,3)-1,3);
rmseRWR = zeros(size(Forecast,3)-1,3);
rmseRWR2 = zeros(size(Forecast,3)-1,3);
rmseAA = zeros(size(Forecast,3)-1,3);
rmseAA2 = zeros(size(Forecast,3)-1,3);

for z = 1:size(Forecast,2)-1
    Reali(:,z)=reshape(reti(:,:,5+z),size(reti,1)*size(reti,2),1);
    Reali(:,z) = Reali(:,z)/sum(Reali(:,z));
    
    rmseModel(z,1) = sqrt(mean((Reali(:,z+1)-ForecastOT(:,z)).^2));
    rmseModel(z,2) = (mean((Reali(:,z+1)-ForecastOT(:,z)).^2));
    rmseModel(z,3) = (mean(abs((ForecastOT(:,z)-Reali(:,z+1)))));

    rmseModelGDP(z,1) = sqrt(mean((Reali(:,z+1)-ForecastOT_GDP(:,z)).^2));
    rmseModelGDP(z,2) = (mean((Reali(:,z+1)-ForecastOT_GDP(:,z)).^2));
    rmseModelGDP(z,3) = (mean(abs((ForecastOT_GDP(:,z)-Reali(:,z+1)))));

    rmseNMF(z,1) = sqrt(mean((Reali(:,z+1)-ForecastNMF(:,z)).^2));
    rmseNMF(z,2) = (mean((Reali(:,z+1)-ForecastNMF(:,z)).^2));
    rmseNMF(z,3) = (mean(abs((ForecastNMF(:,z)-Reali(:,z+1)))));

    rmseNFM2(z,1) = sqrt(mean((Reali(:,z+1)-ForecastNMF2(:,z+1)).^2));
    rmseNFM2(z,2) = (mean((Reali(:,z+1)-ForecastNMF2(:,z)).^2));
    rmseNFM2(z,3) = (mean(abs((ForecastNMF2(:,z+1)-Reali(:,z+1)))));

    rmseK(z,1) = sqrt(mean((Reali(:,z+1)-ForecastK(:,z)).^2));
    rmseK(z,2) = (mean((Reali(:,z+1)-ForecastK(:,z)).^2));
    rmseK(z,3) = (mean(abs((ForecastK(:,z)-Reali(:,z+1)))));

    rmseK2(z,1) = sqrt(mean((Reali(:,z+1)-ForecastK2(:,z+1)).^2));
    rmseK2(z,2) = (mean((Reali(:,z+1)-ForecastK2(:,z)).^2));
    rmseK2(z,3) = (mean(abs((ForecastK2(:,z+1)-Reali(:,z+1)))));


    rmseRWR(z,1) = sqrt(mean((Reali(:,z+1)-ForecastRWR(:,z)).^2));
    rmseRWR(z,2) = (mean((Reali(:,z+1)-ForecastRWR(:,z)).^2));
    rmseRWR(z,3) = (mean(abs((ForecastRWR(:,z)-Reali(:,z+1)))));

    rmseRWR2(z,1) = sqrt(mean((Reali(:,z+1)-ForecastRWR2(:,z+1)).^2));
    rmseRWR2(z,2) = (mean((Reali(:,z+1)-ForecastRWR2(:,z)).^2));
    rmseRWR2(z,3) = (mean(abs((ForecastRWR2(:,z+1)-Reali(:,z+1)))));

    rmseAA(z,1) = sqrt(mean((Reali(:,z+1)-ForecastAA(:,z)).^2));
    rmseAA(z,2) = (mean((Reali(:,z+1)-ForecastAA(:,z)).^2));
    rmseAA(z,3) = (mean(abs((ForecastAA(:,z)-Reali(:,z+1)))));

    rmseAA2(z,1) = sqrt(mean((Reali(:,z+1)-ForecastAA2(:,z+1)).^2));
    rmseAA2(z,2) = (mean((Reali(:,z+1)-ForecastAA2(:,z)).^2));
    rmseAA2(z,3) = (mean(abs((ForecastAA2(:,z+1)-Reali(:,z+1)))));

end

figure
tiledlayout(2,1, 'TileSpacing', 'compact')

nexttile
b=bar([rmseModel(:,1),rmseModelGDP(:,1),rmseNMF(:,1),rmseNFM2(:,1),rmseK(:,1),rmseK2(:,1),rmseRWR(:,1),rmseRWR2(:,1),rmseAA(:,1),rmseAA2(:,1)]);
b(8).FaceColor = [.2 .6 .5];
b(9).FaceColor = [238/256	121/256 159/256];
b(10).FaceColor = [100/256	149/256	237/256];
title('RMSE')
xticks(1:length(rmseNFM2))
xticklabels(leg(2:end))
set(gca,'FontSize',10,'fontWeight','bold')
axis tight
grid on


nexttile
b=bar([rmseModel(:,2),rmseModelGDP(:,2),rmseNMF(:,2),rmseNFM2(:,2),rmseK(:,2),rmseK2(:,2),rmseRWR(:,2),rmseRWR2(:,2),rmseAA(:,2),rmseAA2(:,2)]);
b(8).FaceColor = [.2 .6 .5];
b(9).FaceColor = [238/256	121/256 159/256];
b(10).FaceColor = [100/256	149/256	237/256];
title('MSE')
xticks(1:length(rmseNFM2))
xticklabels(leg(2:end))
set(gca,'FontSize',10,'fontWeight','bold')
axis tight
grid on

% nexttile
% bar([rmseModel(:,3),rmseModelGDP(:,3),rmseNMF(:,3),rmseNFM2(:,3),rmseK(:,3),rmseK2(:,3)])
% title('MAE')
% xticks(1:length(rmseNFM2))
% xticklabels(leg(2:end))
% set(gca,'FontSize',10,'fontWeight','bold')
% axis tight
% grid on

lgd = legend('Wasserstein','Wasserstein GDP','Aver. NMF','Exp. NMF','Aver. Katz','Exp. Katz','Aver. RWR','Exp. RWR','Aver. AA','Exp. AA');
title(lgd,'Models')
lgd.Layout.Tile = 'east';

%% figura network
periodo = 10;
net1a = reshape(ForecastOT(:,periodo),sqrt(size(ForecastOT,1)),sqrt(size(ForecastOT,1)));
deg = sum(net1a);
deg2 = sum(net1a,2);
net1a(net1a<prctile(net1a,98))=0;
net1a = net1a-diag(diag(net1a));
net1 = digraph(net1a);
LWidths = 3*net1.Edges.Weight/max(net1.Edges.Weight);

net1b = reshape(ForecastOT_GDP(:,periodo),sqrt(size(ForecastOT,1)),sqrt(size(ForecastOT,1)));
degb = sum(net1b);
deg2b = sum(net1b,2);
net1b(net1b<prctile(net1b,98))=0;
net1b = net1b-diag(diag(net1b));
net2 = digraph(net1b);
LWidths2 = 3*net2.Edges.Weight/max(net2.Edges.Weight);

figure
subplot(1,2,1)
p = plot(net1,'NodeLabel',NodiTutti,'LineWidth',LWidths,'Layout','force','UseGravity',true,...
    'EdgeColor',[190/256 190/256 190/256]);
p.NodeCData = (deg);
colormap('hot')
c = colorbar('southoutside');
c.Label.String = 'Out-strength';
axis off
p.MarkerSize=exp(4*deg2+2);
axis off
title('Distance Cost')
bubblelim([0 round(max(deg2),2)])
bubblelegend('In-strength')
set(gca,'FontSize',10,'fontWeight','bold')

subplot(1,2,2)
p = plot(net2,'NodeLabel',NodiTutti,'LineWidth',LWidths2,'Layout','force',...
    'UseGravity',true,'EdgeColor',[190/256 190/256 190/256]);
p.NodeCData = degb;
c = colorbar('southoutside');
c.Label.String = 'Out-strength';
axis off
p.MarkerSize=exp(4*deg2b+2);
axis off
title('Distance GDP')
bubblelim([0 round(max(deg2b),2)])
bubblelegend('In-strength')
set(gca,'FontSize',10,'fontWeight','bold')
