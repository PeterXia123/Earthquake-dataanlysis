%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Data analysis of tsunami hazard and loss information for modelling CATBond triggering   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   Earthquake source simulations   %%%
% Source simulations are conducted using 'slipsimu_Tohoku_v3.m'.
% 5001-5500: Mw7.6 scenarios - 'SourceModel_PTHATohokuM76'
% 5501-6000: Mw7.8 scenarios - 'SourceModel_PTHATohokuM78'
% 6001-6500: Mw8.0 scenarios - 'SourceModel_PTHATohokuM80'
% 6501-7000: Mw8.2 scenarios - 'SourceModel_PTHATohokuM82'
% 7001-7500: Mw8.4 scenarios - 'SourceModel_PTHATohokuM84'
% 7501-8000: Mw8.6 scenarios - 'SourceModel_PTHATohokuM86'
% 8001-8500: Mw8.8 scenarios - 'SourceModel_PTHATohokuM88'
% 8501-9000: Mw9.0 scenarios - 'SourceModel_PTHATohokuM90'

%%%   Tsunami hazard simulations   %%%
% Tsunami hazard simulations are conducted using 'run_TOHOKUMCv3.m' for Region1-50m, Region3-50m, and Region4-50m.
% 5001-5500: Mw7.6 scenarios - e.g. 'PTHARegion50m1-5001to5050'
% 5501-6000: Mw7.8 scenarios - e.g. 'PTHARegion50m1-5501to5550'
% 6001-6500: Mw8.0 scenarios - e.g. 'PTHARegion50m1-6001to6050'
% 6501-7000: Mw8.2 scenarios - e.g. 'PTHARegion50m1-6501to6550'
% 7001-7500: Mw8.4 scenarios - e.g. 'PTHARegion50m1-7001to7050'
% 7501-8000: Mw8.6 scenarios - e.g. 'PTHARegion50m1-7501to7550'
% 8001-8500: Mw8.8 scenarios - e.g. 'PTHARegion50m1-8001to8050'
% 8501-9000: Mw9.0 scenarios - e.g. 'PTHARegion50m1-8501to8550'

%%%   Earthquake-tsunami loss estimation   %%%
% Earthquake-tsunami loss estimation is performed using 'PETLE_TOHOKUMCv3.m'.
% Region1: PETLE_Iwanuma_base.mat
% Region3: PETLE_Onagawa_base.mat, PETLE_Ishinomaki_base.mat
% Region4: PETLE_Shizugawa_base.mat, PETLE_Kesennuma_base.mat, 

clear
close all

load japan_coastline_gshhg

% Output from 'DataExtract_MultihazardCATBond_TOHOKUMCv3.m'
load CATbond_MultihazardLoss_Iwanuma

%% Preliminary figures

figure (1); 
contourf(Lon4,Lat4,-DEP4,100,'LineStyle','none'); 
hold on; shading flat; view(0,90); caxis([-50 50]);
plot(japan_gshhg(:,1),japan_gshhg(:,2),'k-'); hold on;
plot(FOCUSAREA([3 4 4 3 3]),FOCUSAREA([1 1 2 2 1]),'w-','LineWidth',2); hold on;
axis equal; axis(RANGE); colorbar;

figure (2); 
plot(japan_gshhg(:,1),japan_gshhg(:,2),'k-'); hold on;
plot(MLITdata(Oid,22),MLITdata(Oid,21),'kv','MarkerSize',3); hold on;
plot(MLITdata(Wid,22),MLITdata(Wid,21),'g^','MarkerSize',3); hold on;
plot(MLITdata(Sid,22),MLITdata(Sid,21),'rs','MarkerSize',3); hold on;
plot(MLITdata(Cid,22),MLITdata(Cid,21),'bo','MarkerSize',3); hold on;
axis equal; axis([FOCUSAREA(3:4) FOCUSAREA(1:2)]);

TLC_ProbLevel(:,1) = [0.005 0.002 0.001]';

for ii = 1:3
    [~,tmp1] = min(abs(TLC(:,1)-TLC_ProbLevel(ii,1)));
    [~,tmp2] = min(abs(TLC(:,2)-TLC_ProbLevel(ii,1)));
    [~,tmp3] = min(abs(TLC(:,3)-TLC_ProbLevel(ii,1)));
    TLC_ProbLevel(ii,2) = TLrange(tmp1);
    TLC_ProbLevel(ii,3) = TLrange(tmp2);
    TLC_ProbLevel(ii,4) = TLrange(tmp3);
end

figure (3);
loglog(TLrange,TLC(:,1),'b-',TLrange,TLC(:,2),'r-',TLrange,TLC(:,3),'g-','LineWidth',1); hold on; 
loglog(TLC_ProbLevel(:,2),TLC_ProbLevel(:,1),'bo',TLC_ProbLevel(:,3),TLC_ProbLevel(:,1),'ro',TLC_ProbLevel(:,4),TLC_ProbLevel(:,1),'go','LineWidth',1); hold on; 
title('Loss'); axis square; axis([1 TC_MUSD(2) 0.00001 0.1]); 
grid on; xlabel('Loss (million USD)'); ylabel('Annual exceedance prob');

figure (4); 
contourf(Lon2,Lat2,-DEP2,100,'LineStyle','none'); 
hold on; shading flat; view(0,90); caxis([-5000 5000]);
plot(japan_gshhg(:,1),japan_gshhg(:,2),'k-'); hold on;
plot(PntInfo(IDgpsb,8),PntInfo(IDgpsb,7),'ko','MarkerFaceColor','y'); hold on;
plot(PntInfo(IDcity,8),PntInfo(IDcity,7),'ko','MarkerFaceColor','w'); hold on;
plot(PntInfo(IDsnet,8),PntInfo(IDsnet,7),'ko','MarkerFaceColor','c'); hold on;
axis equal; axis([140.5 144.5 35 41]); colorbar;

%%

%%%   DATA   %%%
% SCENARIO: Source information (magnitude and 8 kinds of source location) - 500 by 17 by 8
% PGV     : Peak ground velocity at SoilInfoJSHIS (columns 1 and 2)
% WHLabs  : Maximum wave height/depth information (1 min interval) at 119 recording locations - 120 by 119 by 500 by 8
% TL1     : Tsunami loss based on depth information - 500 by 8 

BldgLoc = [mean(MLITdata(:,13)) mean(MLITdata(:,14))]

PntInfo(:,10) = deg2km(distance(BldgLoc(1),BldgLoc(2),PntInfo(:,7),PntInfo(:,8)));
% Find the closest PGV location to BldgLoc 
[~,PGVLocID] = min(abs(BldgLoc(1) - SoilInfoJSHIS(:,2)) + abs(BldgLoc(2) - SoilInfoJSHIS(:,1)));

% Reduce PV to the observation station
PGVobs(:,:) = PGV(PGVLocID,:,:);

pgv_loc_pick = [1000 3000];
pgv_loc_pick(3) = deg2km(distance(SoilInfoJSHIS(pgv_loc_pick(1),2),SoilInfoJSHIS(pgv_loc_pick(1),1),SoilInfoJSHIS(pgv_loc_pick(2),2),SoilInfoJSHIS(pgv_loc_pick(2),1)));
pgv_loc_pick(4) = SoilInfoJSHIS(pgv_loc_pick(1),3);%?
pgv_loc_pick(5) = SoilInfoJSHIS(pgv_loc_pick(2),3);%?
PGVtmp1 = PGV(pgv_loc_pick(1),:,:); PGVtmp1 = PGVtmp1(:);
PGVtmp2 = PGV(pgv_loc_pick(2),:,:); PGVtmp2 = PGVtmp2(:);

figure (50)
plot(PGVtmp1,PGVtmp2,'b.'); axis square; xlabel(['Site ',num2str(pgv_loc_pick(1))]); ylabel(['Site ',num2str(pgv_loc_pick(2))]);
title([num2str(pgv_loc_pick(3)),' km : ',num2str(pgv_loc_pick(4)),' m/s vs ',num2str(pgv_loc_pick(5)),' m/s']);

% Calculate the distance and azimuth between the city and the source
for ii = 1:500
    for jj = 1:8
        for kk = 1:8
            DIST_CITYtoSOURCE(ii,jj,kk) = deg2km(distance(BldgLoc(1),BldgLoc(2),SCENARIO(ii,1+2*(kk-1)+1,jj),SCENARIO(ii,1+2*(kk-1)+2,jj)));
            AZI_CITYtoSOURCE(ii,jj,kk)  = azimuth(BldgLoc(1),BldgLoc(2),SCENARIO(ii,1+2*(kk-1)+1,jj),SCENARIO(ii,1+2*(kk-1)+2,jj));
        end
    end
end

% Reduce WHLabs to the maximum tsunami height/depth over 2 hours duration
for ii = 1:119
    for jj = 1:500
        for kk = 1:8
            WHLabs2(ii,jj,kk) = WHLabs(end,ii,jj,kk);
        end
    end
end

offshore = find(PntInfo(:,9) ~= 2);
onshore  = find(PntInfo(:,9) == 2);

WH_off = WHLabs2(offshore,:,:);
WH_on  = WHLabs2(onshore,:,:);
%?? what

PntInfo_off = PntInfo(offshore,:);
PntInfo_on  = PntInfo(onshore,:);

[tmp1,tmp2] = min(PntInfo_off(:,10));
disp(['Closest station: ',num2str(tmp2),' at ',num2str(tmp1),' km at depth of ',num2str(PntInfo_off(tmp2,6)),' m']);

pick_station = [1 50 71]
disp(['Picked station 1: ',num2str(pick_station(1)),' at ',num2str(PntInfo_off(pick_station(1),10)),' km at depth of ',num2str(PntInfo_off(pick_station(1),6)),' m']);
disp(['Picked station 2: ',num2str(pick_station(2)),' at ',num2str(PntInfo_off(pick_station(2),10)),' km at depth of ',num2str(PntInfo_off(pick_station(2),6)),' m']);
disp(['Picked station 3: ',num2str(pick_station(3)),' at ',num2str(PntInfo_off(pick_station(3),10)),' km at depth of ',num2str(PntInfo_off(pick_station(3),6)),' m']);
%why it is 0
[tmp3,tmp4] = min(PntInfo_on(:,10));
disp(['Closest city: ',num2str(tmp4),' at ',num2str(tmp3),' km']);

figure (4);
plot(SoilInfoJSHIS(PGVLocID,1),SoilInfoJSHIS(PGVLocID,2),'bs','MarkerFaceColor','b'); hold on;
plot(PntInfo_off(tmp2,8),PntInfo_off(tmp2,7),'ko','MarkerFaceColor','r'); hold on;
plot(PntInfo_off(pick_station(1),8),PntInfo_off(pick_station(1),7),'ko','MarkerFaceColor','k'); hold on;
plot(PntInfo_off(pick_station(2),8),PntInfo_off(pick_station(2),7),'ko','MarkerFaceColor','k'); hold on;
plot(PntInfo_off(pick_station(3),8),PntInfo_off(pick_station(3),7),'ko','MarkerFaceColor','k'); hold on;
plot(PntInfo_on(tmp4,8),PntInfo_on(tmp4,7),'ko','MarkerFaceColor','r'); hold on;

figure (100+6); plot(japan_gshhg(:,1),japan_gshhg(:,2),'b-'); hold on; axis equal;
figure (100+7); plot(japan_gshhg(:,1),japan_gshhg(:,2),'b-'); hold on; axis equal;

for iii = length(MW):-1:1
% for iii = 1:length(MW)   

%     if iii == 1;
%         symbol = 'bo';
%     elseif iii == 2;
%         symbol = 'ro';
%     elseif iii == 3;
%         symbol = 'go';
%     elseif iii == 4;
%         symbol = 'ko';
%     elseif iii == 5;
%         symbol = 'bs';
%     elseif iii == 6;
%         symbol = 'rs';
%     elseif iii == 7;
%         symbol = 'gs';
%     elseif iii == 8;
%         symbol = 'ks';
%     end
    
    if iii == 1
        symbol = 'bo';
    elseif iii == 2
        symbol = 'bo';
    elseif iii == 3
        symbol = 'rs';
    elseif iii == 4
        symbol = 'rs';
    elseif iii == 5
        symbol = 'g^';
    elseif iii == 6
        symbol = 'g^';
    elseif iii == 7
        symbol = 'kv';
    elseif iii == 8
        symbol = 'kv';
    end
    
    figure (100+1);
    subplot(121); plot(SCENARIO(:,1,iii),TL_E(:,iii),symbol); hold on; axis square;
    xlabel('Magnitude'); ylabel('EQ loss');

    subplot(122); plot(SCENARIO(:,1,iii),TL_T(:,iii),symbol); hold on; axis square;
    xlabel('Magnitude'); ylabel('Tsunami loss');    
    
    figure (100+2);
    subplot(121); plot(DIST_CITYtoSOURCE(:,iii,1),TL_E(:,iii),symbol); hold on; axis square;
    xlabel('Distance (geometrical centroid)'); ylabel('EQ loss');
    
    subplot(122); plot(DIST_CITYtoSOURCE(:,iii,1),TL_T(:,iii),symbol); hold on; axis square;
    xlabel('Distance (geometrical centroid)'); ylabel('Tsunami loss');
    
    figure (100+3);
    subplot(121); plot(DIST_CITYtoSOURCE(:,iii,8),TL_E(:,iii),symbol); hold on; axis square;
    xlabel('Distance (slip centroid)'); ylabel('EQ loss');

    subplot(122); plot(DIST_CITYtoSOURCE(:,iii,8),TL_T(:,iii),symbol); hold on; axis square;
    xlabel('Distance (slip centroid)'); ylabel('Tsunami loss');
    
    figure (100+4);
    subplot(121); plot(AZI_CITYtoSOURCE(:,iii,1),TL_E(:,iii),symbol); hold on; axis square;
    xlabel('Azimuth (geometrical centroid)'); ylabel('EQ loss');
    
    subplot(122); plot(AZI_CITYtoSOURCE(:,iii,1),TL_T(:,iii),symbol); hold on; axis square;
    xlabel('Azimuth (geometrical centroid)'); ylabel('Tsunami loss');
    
    figure (100+5);
    subplot(121); plot(AZI_CITYtoSOURCE(:,iii,8),TL_E(:,iii),symbol); hold on; axis square;
    xlabel('Azimuth (slip centroid)'); ylabel('EQ loss');
    
    subplot(122); plot(AZI_CITYtoSOURCE(:,iii,8),TL_T(:,iii),symbol); hold on; axis square;
    xlabel('Azimuth (slip centroid)'); ylabel('Tsunami loss');

%     figure (100+6);
    for ii = 1:500
        polygontmp(:,1) = SCENARIO(ii,2,iii) + 0.1*unitcircle(:,1);
        polygontmp(:,2) = SCENARIO(ii,3,iii) + 0.1*unitcircle(:,2);
        PGtmp0(ii+(iii-1)*500,1)    = TL_E(ii,iii);
        PGtmp0(ii+(iii-1)*500,2)    = TL_T(ii,iii);
        PGtmp1(ii+(iii-1)*500,1:73) = polygontmp(:,1);
        PGtmp2(ii+(iii-1)*500,1:73) = polygontmp(:,2);
%         subplot(121); fill(polygontmp(:,2),polygontmp(:,1),TL_E(ii,iii),'EdgeColor','none'); hold on;
%         subplot(122); fill(polygontmp(:,2),polygontmp(:,1),TL_T(ii,iii),'EdgeColor','none'); hold on;
        clear polygontmp
    end
    
%     figure (100+7);
    for ii = 1:500
        polygontmp(:,1) = SCENARIO(ii,16,iii) + 0.1*unitcircle(:,1);
        polygontmp(:,2) = SCENARIO(ii,17,iii) + 0.1*unitcircle(:,2);
        PGtmp3(ii+(iii-1)*500,1:73) = polygontmp(:,1);
        PGtmp4(ii+(iii-1)*500,1:73) = polygontmp(:,2);
%         subplot(121); fill(polygontmp(:,2),polygontmp(:,1),TL_E(ii,iii),'EdgeColor','none'); hold on;
%         subplot(122); fill(polygontmp(:,2),polygontmp(:,1),TL_T(ii,iii),'EdgeColor','none'); hold on;
        clear polygontmp
    end
    
    figure (200+0);
    plot(PGVobs(:,iii),TL_E(:,iii),symbol); hold on; axis square;
    xlabel(['PGV at ',num2str(PGVLocID)]); ylabel('EQ loss');

    figure (200+1);
    plot(WH_off(tmp2,:,iii),TL_T(:,iii),symbol); hold on; axis square;
    xlabel(['Wave height at ',num2str(tmp2)]); ylabel('Tsunami loss          ./÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷//÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷');
    %this one is same as this one
    figure (200+2);
    plot(WH_off(pick_station(1),:,iii),TL_T(:,iii),symbol); hold on; axis square;
    xlabel(['Wave height at ',num2str(pick_station(1))]); ylabel('Tsunami loss');

    figure (200+3);
    plot(WH_off(pick_station(2),:,iii),TL_T(:,iii),symbol); hold on; axis square;
    xlabel(['Wave height at ',num2str(pick_station(2))]); ylabel('Tsunami loss');
    
    figure (200+4);
    plot(WH_off(pick_station(3),:,iii),TL_T(:,iii),symbol); hold on; axis square;
    xlabel(['Wave height at ',num2str(pick_station(3))]); ylabel('Tsunami loss');
    
end

TMP = [PGtmp0 PGtmp1 PGtmp2 PGtmp3 PGtmp4];
clear PGtmp0 PGtmp1 PGtmp2 PGtmp3 PGtmp4

TMPE = sortrows(TMP,1);
PGtmpE0 = TMPE(:,1);
PGtmpE1 = TMPE(:,3:75);
PGtmpE2 = TMPE(:,76:148);
PGtmpE3 = TMPE(:,149:221);
PGtmpE4 = TMPE(:,222:294);

TMPT = sortrows(TMP,2);
PGtmpT0 = TMPT(:,2);
PGtmpT1 = TMPT(:,3:75);
PGtmpT2 = TMPT(:,76:148);
PGtmpT3 = TMPT(:,149:221);
PGtmpT4 = TMPT(:,222:294);

clear TMP TMPE TMPT

figure (100+6); 
subplot(121);
for ii = 1:4000
    fill(PGtmpE2(ii,:),PGtmpE1(ii,:),PGtmpE0(ii,1),'EdgeColor','none'); hold on;
end
axis([140 145 34.5 41.5]); colormap(jet); grid on; colorbar; caxis([0 TC_MUSD(1)/4]); title('EQ loss - Geometrical centroid');
subplot(122);
for ii = 1:4000
    fill(PGtmpT2(ii,:),PGtmpT1(ii,:),PGtmpT0(ii,1),'EdgeColor','none'); hold on;
end
axis([140 145 34.5 41.5]); colormap(jet); grid on; colorbar; caxis([0 TC_MUSD(1)]); title('Tsunami loss - Geometrical centroid');

figure (100+7); 
subplot(121);
for ii = 1:4000
    fill(PGtmpE4(ii,:),PGtmpE3(ii,:),PGtmpE0(ii,1),'EdgeColor','none'); hold on;
end
axis([140 145 34.5 41.5]); colormap(jet); grid on; colorbar; caxis([0 TC_MUSD(1)/4]); title('EQ loss - Slip centroid');
subplot(122);
for ii = 1:4000
    fill(PGtmpT4(ii,:),PGtmpT3(ii,:),PGtmpT0(ii,1),'EdgeColor','none'); hold on;
end
axis([140 145 34.5 41.5]); colormap(jet); grid on; colorbar; caxis([0 TC_MUSD(1)]); title('Tsunami loss - Slip centroid');

%% Preliminary investigations of tsunami loss prediction ability
TL = TL_T(:);
non_zeroTL = find(TL > 0);

Sigma_max = [100 0];

for ii = 1:length(offshore)
    
    data_wh(:,:) = WH_off(ii,:,:);
    data_wh = data_wh(:);
    
    [b,bint,r] = regress(log10(TL(non_zeroTL)),[ones(length(non_zeroTL),1) log10(data_wh(non_zeroTL))]);
    
    tl_pred = 10.^(b(1) + b(2)*log10([0.01 10]));
    
    Sigma(ii,1) = std(r);
    
    if Sigma(ii,1) < Sigma_max(1)
        
        figure (1000);
        loglog(data_wh(non_zeroTL),TL(non_zeroTL),'bo',[0.01 10],tl_pred,'r-');
        title(['StationID: ',num2str(ii),' - sigma = ',num2str(Sigma(ii,1))]);
        
        Sigma_max(1) = Sigma(ii,1);
        Sigma_max(2) = ii;
        
    end
    
    clear data_wh b b_int r
    
end

figure (4);
plot(PntInfo_off(Sigma_max(2),8),PntInfo_off(Sigma_max(2),7),'ko','MarkerFaceColor','m'); hold on;

figure (200+5); 
plot(japan_gshhg(:,1),japan_gshhg(:,2),'b-'); hold on; axis equal;
axis([140 145 34.5 41.5]); colormap(jet); grid on; colorbar; caxis([min(Sigma) max(Sigma)]); title('Sigma');

for ii = 1:99
    
    polygontmp(:,1) = PntInfo_off(ii,7) + 0.1*unitcircle(:,1);
    polygontmp(:,2) = PntInfo_off(ii,8) + 0.1*unitcircle(:,2);
    fill(polygontmp(:,2),polygontmp(:,1),Sigma(ii,1),'EdgeColor','none'); hold on;
    clear polygontmp
    
end

