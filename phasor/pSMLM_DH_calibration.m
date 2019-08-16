function [curveAngle curveDistAngle magnratiorange wobbleMatrix] = pSMLM_DH_calibration(input,zposcali,imagegen)
%% Phasor double-helix calibration
%pSMLM_DH_calibration(locdata, zposcalib, imagegen)
%Input: 
%     - frame-x-y data of 2D phasor localizations 
%     - array of z-positions at which calibration images are taken. Should
%       be same length as max number of frames of input.
%     - Boolean operation for whether or not image should be generated.
%Output: 
%     - Curve for DH calib angle over Z
%     - Curve for DH calib distance over angle
%     - Range for ratio of magn ratios (magn_ratio_top/magn_ratio_bot) over
%     Z
%     - wobble file to correct for wobble in xy

%% Required variables
%Very rough min and max distance between points, in pixels. Preferably this
%should be well outside real value
minDistApprox = 5;
maxDistApprox = 20;

maxWobbleDistance = 1; %distance in pixels the centre of a DHPSF can wobble between 2 frames
maxWobbleMemory = 3; %memory of tracking centres of PSFs
movmeansize = 15; %nr of frames to movmean over (odd) for wobble correction
wobbleIntervalLen = 5; %nr of frames to average the wobble over

%% Pairing of localizations
[locationsDH, ~] = DH_phasor_linking(input,[],[minDistApprox maxDistApprox],[],[]);
%Calculate mean and std of distance and angle, based on frame
for i = 1:max(locationsDH(:,1))
    distancemean(i) = mean(locationsDH(locationsDH(:,1)==i,6));
    distancestd(i) = std(locationsDH(locationsDH(:,1)==i,6),0,1);
    anglemean(i) = mean(locationsDH(locationsDH(:,1)==i,5));
    anglestd(i) = std(locationsDH(locationsDH(:,1)==i,5),0,1);
end

%Give error warning if min/max approximate distance were too close
if min(abs(minDistApprox-min(distancemean)),abs(maxDistApprox-max(distancemean))) < 2
   warning('Approximated min and max distance is close to found min or max distance!')
   warning('Consider changing min/max distance approximation in pSMLM_DH_calibration')
end

%% Fit calibration curve to angle
%Fitcurve is an arbitrary x^3-curve
ft = @(a,x) a(1).*x.^3+a(2).*x.^2+a(3).*x+a(4);
weightarr = 1./anglestd;
curveAngle = fitnlm(zposcali',anglemean,ft,[1 1 1 1],'Weight',weightarr);

%Make Y-cali data
% ycalidata = zeros(size(zposcali,2),1);
afitdata = table2array(curveAngle.Coefficients(:,1));
ycalidata = afitdata(1).*zposcali.^3+afitdata(2).*zposcali.^2+afitdata(3).*zposcali+afitdata(4);

curveDistAngle = fitnlm(anglemean,distancemean,ft,[1 1 1 1],'Weight',weightarr);
afitdata2 = table2array(curveDistAngle.Coefficients(:,1));
ycalidata2 = afitdata2(1).*anglemean.^3+afitdata2(2).*anglemean.^2+afitdata2(3).*anglemean+afitdata2(4);

%% Wobble correction
%Easiest way: just perform tracking, then look at tracks to determine
%wobble
%Initiate tracking parameters
param.maxDisp = maxWobbleDistance;
param.mem = maxWobbleMemory;
param.dim = 2; %2D
param.good = 1; %reject short tracks?
param.quiet = 1; %no verbose
wobbleTracks = trackWithDummy([locationsDH(:,2:3) locationsDH(:,1)],param); %perform tracking
%wobbleTracks now is a x,y,frame,trackID matrix

%wobbleMatrix needs to be minz-maxz-offsetx-offsety matrix, where offsetx/y
%are the average offset between minz and maxz values.
%First, normalize the tracks to so that middle of tracks is 0,0

%loop over tracks
movmeanarrtot = [];
for trackid = 1:max(wobbleTracks(:,4))
clear movmeanarr
    tempTracks = wobbleTracks(wobbleTracks(:,4)==trackid,:);
    if (size(tempTracks,1)>1)
        xmovmean = movmean(tempTracks(:,1),movmeansize);
        ymovmean = movmean(tempTracks(:,2),movmeansize);
        midpos = round((max(tempTracks(:,3))-min(tempTracks(:,3)))/2);
%         offsetx = xmovmean(tempTracks(:,3)==midpos);
%         offsety = ymovmean(tempTracks(:,3)==midpos);
        offsetx = xmovmean(midpos);
        offsety = ymovmean(midpos);
        movmeanarr(:,1) = xmovmean-offsetx;
        movmeanarr(:,2) = ymovmean-offsety;
        movmeanarr(:,3) = tempTracks(:,3);
        movmeanarrtot = [movmeanarrtot; movmeanarr];
    end
end

%take avg tracks after normalization
for frameid = min(wobbleTracks(:,3)):max(wobbleTracks(:,3))
    avgWobbleTracksMovMean(frameid,1:2) = mean(movmeanarrtot(movmeanarrtot(:,3)==frameid,1:2));
    stdWobbleTracksMovMean(frameid,1:2) = std(movmeanarrtot(movmeanarrtot(:,3)==frameid,1:2),0,1);
end

%Make more spaced points, i.e. 5-10 frames - so that every zpos interval
%has a x,y value pair associated with it
wobbleIntervalBins = round(max(size(zposcali)-1)/wobbleIntervalLen);
wobbleMatrix = [];
for wobbleInterval = 1:wobbleIntervalBins
    valstart = round((wobbleInterval-1)*(max(size(zposcali))-1)/wobbleIntervalBins)+1;
    valend = round((wobbleInterval)*(max(size(zposcali))-1)/wobbleIntervalBins)+1;
    meanXwobble = mean(avgWobbleTracksMovMean(valstart:valend,1));
    meanYwobble = mean(avgWobbleTracksMovMean(valstart:valend,2));
    wobbleMatrix(wobbleInterval,:) = [zposcali(valstart),zposcali(valend),meanXwobble,meanYwobble];
end


%% Magnitude min/max ratios, based on 4x std
%Get mean and stdev of every zpos
%same bins as wobbleMatrix
%Make more spaced points, i.e. 5-10 frames - so that every zpos interval
%has a x,y value pair associated with it
magnratioIntervalBins = round(max(size(zposcali)-1)/wobbleIntervalLen);
magnratiorange = zeros(magnratioIntervalBins,4);
meanmagnrat = [];
stdmagnrat = [];
for maggratioInterval = 1:magnratioIntervalBins
    valstart = round((maggratioInterval-1)*(max(size(zposcali))-1)/magnratioIntervalBins)+1;
    valend = round((maggratioInterval)*(max(size(zposcali))-1)/magnratioIntervalBins)+1;
    %Get original locationsDH array with only values in these frames
    tt = locationsDH(locationsDH(:,1)>=valstart,:);
    tt = tt(tt(:,1)<valend,:);
    meanmagnrat = mean(tt(:,9));
    stdmagnrat = std(tt(:,9),0,1);
    magnratiorange(maggratioInterval,:) = [zposcali(valstart),zposcali(valend),meanmagnrat-3*stdmagnrat,meanmagnrat+3*stdmagnrat];
end
% keyboard
%% Show data and calibration
if imagegen
    figure(1);clf(1);    
    subplot('Position',[0.1 0.1 0.35 0.8])
    hold on
    yyaxis left
    fill([zposcali fliplr(zposcali)],0.1*[(distancemean+distancestd*2) fliplr(distancemean-distancestd*2)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(zposcali,0.1*distancemean,'k-');
    %     errorbar(zposcali,distancemean,distancestd);
    ylabel('Distance between lobes (µm)')
    axis([-1 1 1.1-.25 1.1+.25])
    ax = gca;
    ax.YColor = 'k';
    yyaxis right
    fill([zposcali fliplr(zposcali)],[(anglemean+anglestd*2) fliplr(anglemean-anglestd*2)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
%     plot(zposcali,anglemean,'k:');
%     errorbar(zposcali,anglemean,anglestd);
    plot(zposcali,ycalidata,'k-')
    % plot(zposcali,curveCI,'k:')
    ylabel('Angle (rad)')
    ax = gca;
    ax.YColor = 'k';
    xlabel('Axial position (µm)')
    title({'Calibration phasor-DoubleHelix',''})
    set(gca,'TickDir','out')
    grid on
    box on
    
%     subplot('Position',[0.6 0.7 0.3 0.2])
%     hold on
%     plot(anglemean,distancemean)
%     plot(anglemean,ycalidata2,'k-')
%     ylabel('Distance between points (px)')
%     xlabel('Angle (rad)')
%     title('Fit angle-distance')
    
    subplot('Position',[0.6 0.4 0.3 0.2])
    hold on
%     errorbar(zposcali(1:100),avgWobbleTracksMovMean(1:100,1),stdWobbleTracksMovMean(1:100,1),'b-')
%     errorbar(zposcali(1:100),avgWobbleTracksMovMean(1:100,2),stdWobbleTracksMovMean(1:100,2),'r-')
    fill([zposcali fliplr(zposcali)],...
        100*[(avgWobbleTracksMovMean(:,1)+stdWobbleTracksMovMean(:,1))' flipud(avgWobbleTracksMovMean(:,1)-stdWobbleTracksMovMean(:,1))'],...
        [0 0 1],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
    a=plot(zposcali,100*avgWobbleTracksMovMean(:,1),'b')
    fill([zposcali fliplr(zposcali)],...
        100*[(avgWobbleTracksMovMean(:,2)+stdWobbleTracksMovMean(:,2))' flipud(avgWobbleTracksMovMean(:,2)-stdWobbleTracksMovMean(:,2))'],...
        [1 0 0],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
    b=plot(zposcali,100*avgWobbleTracksMovMean(:,2),'r')
%     for i = 1:size(wobbleMatrix,1)
%         scatter((wobbleMatrix(i,1))*.5+(wobbleMatrix(i,2))*.5,wobbleMatrix(i,3),'k*')
%         scatter((wobbleMatrix(i,1))*.5+(wobbleMatrix(i,2))*.5,wobbleMatrix(i,4),'k*')
%     end
    leg =legend('x','y','Location','best')
    leg.ItemTokenSize = [10,18];
    xlabel('Axial position (µm)')
    ylabel('Wobble (nm)')
    title({'Wobble correction based on z position',''})
    set(gca,'TickDir','out')
    axis([-.8 .8 -25 25])
    box on
    grid on
    
    
%     subplot('Position',[0.6 0.1 0.3 0.2])
%     hold on
%     for i = 1:size(magnratiorange,1)
%         scatter((magnratiorange(i,1))*.5+(magnratiorange(i,2))*.5,magnratiorange(i,3),7,'ko','filled')
%         scatter((magnratiorange(i,1))*.5+(magnratiorange(i,2))*.5,magnratiorange(i,4),7,'ko','filled')
%     end
%     xlabel('Zpos (µm)')
%     ylabel('Magnitude ratio (3x std)')
%     title('Magnitude ratio based on z pos')
end
end