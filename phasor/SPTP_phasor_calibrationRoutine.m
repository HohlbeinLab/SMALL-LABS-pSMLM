function [curveXnlm,curveYnlm] = SPTP_phasor_calibrationRoutine(calibstack,zpositions)
%% Calibration procedure of SP or TP via phasor
%required input: calibstack: a N-by-N-by-frames-by-nrlocs stack
%               zpositions: a frames-by-1 list of z-positions of all frames
%               in µm
%output: nlm-like coefficients to fitted curves, based on Z-µm
% ROIsize = size(calibstack,1);
% %Get middle pixel for every PSF
% midpx = (size(pixelatedPSF,1)+2)/2;
% minpos = midpx-(ROIsize-1)/2;
% maxpos = midpx+(ROIsize-1)/2;
% iterations = size(calibstack,4);
distX = zeros(size(calibstack,3),size(calibstack,4));
distY = zeros(size(calibstack,3),size(calibstack,4));
posX1 = zeros(size(calibstack,3),size(calibstack,4));
posX2 = zeros(size(calibstack,3),size(calibstack,4));
posY1 = zeros(size(calibstack,3),size(calibstack,4));
posY2 = zeros(size(calibstack,3),size(calibstack,4));
posXtot = zeros(size(calibstack,3),size(calibstack,4));
posYtot = zeros(size(calibstack,3),size(calibstack,4));

for itpos = 1:size(calibstack,4)
    if sum(sum(sum(calibstack(:,:,:,itpos)))) > 0
[distX(:,itpos) distY(:,itpos) posX1(:,itpos) posX2(:,itpos) posY1(:,itpos) posY2(:,itpos) posXtot(:,itpos) posYtot(:,itpos)] = ctPhasor(calibstack(:,:,:,itpos));
    end
end
for i = 1:size(distX,1)
    meandistX(i) = mean(distX(i,~isnan(distX(i,:))));
    stddistX(i) = std(distX(i,~isnan(distX(i,:))),0,2);
    meandistY(i) = mean(distY(i,~isnan(distY(i,:))));
    stddistY(i) = std(distY(i,~isnan(distY(i,:))),0,2);
end
meandistX = meandistX';
stddistX = stddistX';
meandistY = meandistY';
stddistY = stddistY';
meandistX(isnan(meandistX))=0;
stddistX(isnan(stddistX))=0;
meandistY(isnan(meandistY))=0;
stddistY(isnan(stddistY))=0;

% distX(isnan(distX)) = 0;
% distY(isnan(distY)) = 0;
% distX(:,sum(distX,1)==0) = [];
% distY(:,sum(distY,1)==0) = [];
% 
% 
% 
% meandistX = mean(distX,2);
% meandistY = mean(distY,2);
% stddistX = std(distX,0,2);
% stddistY = std(distY,0,2);

%Fitparameters
fo = fitoptions('Method','NonlinearLeastSquares',...
    'StartPoint',[1,1,1,1]);
ft = fittype('real(a*x^3+b*x^2+c*x+d)','options',fo);
ftnlm = @(a,x) real(a(1).*x.^3+a(2).*x.^2+a(3).*x+a(4));

%Fit x distance curve
fitarraydistX = meandistX(meandistX>(0.1))';%50/100
fitarraydistXstd = stddistX(meandistX>(0.1))';%50/100
[fitarraydistX_xpos,~] = find(meandistX == fitarraydistX);
fitarraydistX_xpos=zpositions(fitarraydistX_xpos)';
if sum(fitarraydistXstd)>0
[curveXnlm] = fitnlm((fitarraydistX_xpos),fitarraydistX',ftnlm,[1 1 1 1],'Weight',1./fitarraydistXstd);
else
[curveXnlm] = fitnlm((fitarraydistX_xpos),fitarraydistX',ftnlm,[1 1 1 1]);
end
[curveXnlmpred,curveXnlmconf] = predict(curveXnlm,fitarraydistX_xpos,'Simultaneous',true);

%Fit Y distance curve
fitarraydistY = meandistY(meandistY>(0.1))';%50/100
fitarraydistYstd = stddistY(meandistY>(0.1))';%50/100
[fitarraydistY_xpos,~] = find(meandistY == fitarraydistY);
fitarraydistY_xpos=zpositions(fitarraydistY_xpos)';
if sum(fitarraydistYstd)>0
[curveYnlm] = fitnlm((fitarraydistY_xpos),fitarraydistY',ftnlm,[1 1 1 1],'Weight',1./fitarraydistYstd);
else
[curveYnlm] = fitnlm((fitarraydistY_xpos),fitarraydistY',ftnlm,[1 1 1 1]);
end
[curveYnlmpred,curveYnlmconf] = predict(curveYnlm,fitarraydistY_xpos,'Simultaneous',true);

figure(6);clf(6);
hold on
a1=errorbar(zpositions,meandistX,stddistX,'r--o','MarkerSize',2,'Color',[1 0.4 0.2])
a2=errorbar(zpositions,meandistY,stddistY,'b--o','MarkerSize',2,'Color',[0.2 0.4 1])
k = plot(fitarraydistX_xpos,curveXnlmpred,'r-')
set(k,'LineWidth',2)
k = plot(fitarraydistY_xpos,curveYnlmpred,'b-')
set(k,'LineWidth',2)
k=plot(fitarraydistX_xpos,curveXnlmconf,'r:')
set(k,'LineWidth',1.5)
k=plot(fitarraydistY_xpos,curveYnlmconf,'b:')
set(k,'LineWidth',1.5)
axis([-inf inf 0 12])
grid on
legend({'X-real distance','Y-real distance','X-fit','Y-fit'},'location','north')
xlabel('z position (µm)')
ylabel('Distance between lobes (px)')
end