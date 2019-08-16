function [fitcurveXY,fitcurveYX] = astig_phasor_calibrationRoutine(calibstack,zpositions)
%% Calibration procedure of astigmatism via phasor
%Perform phasor on the stack. Input stack should have correct dimensions
%already
%calibstack is x-y-frame-repeatOnFrame 4D matrix
for frame = 1:size(zpositions,2)
    try
    [~,magnitudes(:,:,frame),~,~,~] = Phasor_localization_SMALLLABS(permute(calibstack(:,:,frame,:),[1 2 4 3]),0,'');
    end
end
avgMagnRatioXY = zeros(size(zpositions,2),1);
avgMagnRatioYX = zeros(size(zpositions,2),1);
stdMagnRatioXY = zeros(size(zpositions,2),1);
stdMagnRatioYX = zeros(size(zpositions,2),1);
%Get avg and std of magnitude ratio
for frame = 1:size(zpositions,2)
    %if there are any calculated magnitudes whatsoever
    if (frame <= size(magnitudes,3))
        if (~isempty(magnitudes(magnitudes(:,1,frame)>0,:,frame)))
            t = magnitudes(magnitudes(:,1,frame)>0,:,frame);
            avgMagnRatioXY(frame,1) = mean(t(:,1)./t(:,2));
            stdMagnRatioXY(frame,1) = std((t(:,1)./t(:,2))',0,2);
            avgMagnRatioYX(frame,1) = mean(t(:,2)./t(:,1));
            stdMagnRatioYX(frame,1) = std((t(:,2)./t(:,1))',0,2);
        end
    end
end
%cleanup
zpos = [1:size(zpositions,2)];
zposXY = (zpos(((avgMagnRatioXY > 0 )&(avgMagnRatioXY>avgMagnRatioYX))))';
stdXY = stdMagnRatioXY(((avgMagnRatioXY > 0 )&(avgMagnRatioXY>avgMagnRatioYX)));
magnXY = avgMagnRatioXY(((avgMagnRatioXY > 0 )&(avgMagnRatioXY>avgMagnRatioYX)));
zposYX = (zpos(((avgMagnRatioYX > 0 )&(avgMagnRatioYX>avgMagnRatioXY))))';
stdYX = stdMagnRatioYX(((avgMagnRatioYX > 0 )&(avgMagnRatioYX>avgMagnRatioXY)));
magnYX = avgMagnRatioYX(((avgMagnRatioYX > 0 )&(avgMagnRatioYX>avgMagnRatioXY)));

%Find cross-point
if magnXY(1) > 1
    crosspoint = max(zposXY(magnXY>1))
else
    crosspoint = max(zposYX(magnYX>1))
end
% crosspoint = 40;

zposdistance = (zpositions(2)-zpositions(1));

zposXY = (zposXY - crosspoint)*zposdistance;
zposYX = (zposYX - crosspoint)*zposdistance;

fo = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[1 1 1 1]);
ft = fittype('a*x^3+b*x^2+c*x+d','options',fo);

%Fit equations (a1*(z-c1)^2+b1 and a2*(z-c2)^2+b2)
weightsarr = (1./stdXY).*(magnXY);
weightsarr(isinf(weightsarr)) = 0;
[fitcurveXY,gofXY] = fit(zposXY,magnXY,ft,'weight',weightsarr)
weightsarr = (1./stdYX).*(magnYX);
weightsarr(isinf(weightsarr)) = 0;
[fitcurveYX,gofYX] = fit(zposYX,magnYX,ft,'weight',weightsarr)


figure(6);clf(6);
hold on
errorbar(zposXY,magnXY,stdXY,'b.')
plot(fitcurveXY,'b-')
errorbar(zposYX,magnYX,stdYX,'r.')
plot(fitcurveYX,'r-')
grid on
axis([min(zpositions) max(zpositions) 0 inf])
ylabel('Phasor magnitude ratio')
xlabel('z position (µm)')
legend('off')

end