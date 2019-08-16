%Use the circle-tangent method with phasor to get distance of foci in X and
%Y of astigmatic PSFs
%input: 3D PSF stack of 2D PSF
%output: distance in X and Y, positions in X and Y (redundant atm), and
%central position in X and Y coordinates (8 vars in total)

function [distX distY PositionXPSF1 PositionXPSF2 PositionYPSF1 PositionYPSF2 PositionX PositionY] = ctPhasor(PSFarr)

    pixelsizephasor = size(PSFarr,1);
    
    %Do fourier transformation of images
    phasorfft = fft2(PSFarr);
    %Calculate the angle of the X-phasor by taking the first harmonic in X
    angX = permute(angle(phasorfft(1,2,:)),[3 1 2]);
    %Correct the angle
    angX(angX>0)=angX(angX>0)-2*pi;
    %Calculate X position
    PositionX = (abs(angX)/(2*pi/pixelsizephasor));
    FWHMY = abs(phasorfft(2,1,:));
    %Calculate the angle of the Y-phasor by taking the first harmonic in Y
    angY = permute(angle(phasorfft(2,1,:)),[3 1 2]);
    %Correct the angle
    angY(angY>0)=angY(angY>0)-2*pi;
    %Calculate Y position
    PositionY = (abs(angY)/(2*pi/pixelsizephasor));
    FWHMX = abs(phasorfft(1,2,:));
    maxMagn = max(FWHMX,FWHMY);
%     totmaxMagn = max(max(max(maxMagn)));
%     for i = 1:size(maxMagn,1)
%     for j = 1:size(maxMagn,1)
%         maxMagn(i,1,j) = totmaxMagn;
%     end
%     end
%     keyboard

    %% ---------------Y pos calculation ----------------------------
    %Make equation for line perpendicular to phasor point
    %First for Y pos
    
    %Orig point is phasorfft(2,1,:). The slope of the line perpendicular to this is
    %one over the slope of the line for phasorfft(2,1,:)
    %The orig slope is real/imag
    origSlopeY = imag(phasorfft(2,1,:))./real(phasorfft(2,1,:));
    perpendicularSlopeY = -1./origSlopeY;
    %perp.slopeX is the 'a' of y=ax+b
    %here we calculate the b since we know a crossing-point
    perpendicularIntersectionpointY = imag(phasorfft(2,1,:))-(real(phasorfft(2,1,:)).*perpendicularSlopeY);
    
    %now we have to check where this lines crosses with the circle of the
    %original (big) magnitude of a single PSF. This should have 2 intersection
    %points!
    %Mathematical description of circle: y^2 + x^2 = r^2, assuming the circle
    %passes through the origin (0,0). r should be the magnitude of a single
    %PSF.
%     a = (perpendicularSlopeY.^2+1);
%     b = 2*perpendicularSlopeY.*perpendicularIntersectionpointY;
%     c = -1*maxMagn.^2+perpendicularIntersectionpointY.^2;
%     D = max(0,b.^2-4.*a.*c);%abs(b.^2-4.*a.*c);%
%     
%     x1 = (-b+sqrt(D))./(2.*a);
%     y1 = perpendicularSlopeY.*x1+perpendicularIntersectionpointY;
%     x2 = (-b-sqrt(D))./(2.*a);
%     y2 = perpendicularSlopeY.*x2+perpendicularIntersectionpointY;

    %new method (17-10-2018):
    %solving y^2+x^2=r^2 with y=ax+b result in:
    a = perpendicularSlopeY;
    b = perpendicularIntersectionpointY;
    r = maxMagn;
    x1 = -1*((sqrt((a.^2+1).*r.^2-b.^2)+a.*b)./(a.^2+1));
    x2 = ((sqrt((a.^2+1).*r.^2-b.^2)-a.*b)./(a.^2+1));
    y1 = x1.*a+b;
    y2 = x2.*a+b;
    
%     counterXDbelowzero = 0;
%     for i = 1:size(PSFarr,3)
%         if D(i) <= 0
%             x1(i) = 0.1;
%             y1(i) = 0.1;
%             x2(i) = 0.1;
%             y2(i) = 0.1;
%             counterXDbelowzero = counterXDbelowzero+1;
%         end
%     end
    x1 = real(x1);
    x2 = real(x2);
    y1 = real(y1);
    y2 = real(y2);
    %Set the intersection points back to complex numbers
    origYPSF1 = permute(complex(x1,y1),[1 3 2]);
    origYPSF2 = permute(complex(x2,y2),[1 3 2]);
    
    %For both,
    %Calculate the angle of the Y-phasor by taking the first harmonic in Y
    angY = angle(origYPSF1);
    % Correct the angle
    angY(angY>0)=angY(angY>0)-2*pi;
    % Calculate Y position
    PositionYPSF1 = (abs(angY)/(2*pi/pixelsizephasor));
    %Calculate the angle of the Y-phasor by taking the first harmonic in Y
    angY = angle(origYPSF2);
    % Correct the angle
    angY(angY>0)=angY(angY>0)-2*pi;
    % Calculate Y position
    PositionYPSF2 = (abs(angY)/(2*pi/pixelsizephasor));
       
    %% ---------------X pos calculation ----------------------------
    %Make equation for line perpendicular to phasor point
    %Now for X pos
    
    %Orig point is phasorfft(1,2,:). The slope of the line perpendicular to this is
    %one over the slope of the line for phasorfft(1,2,:)
    %The orig slope is real/imag
    origSlopeX = imag(phasorfft(1,2,:))./real(phasorfft(1,2,:));
    perpendicularSlopeX = -1./origSlopeX;
    %perp.slopeX is the 'a' of y=ax+b
    %here we calculate the b since we know a crossing-point
    perpendicularIntersectionpointX = imag(phasorfft(1,2,:))-(real(phasorfft(1,2,:)).*perpendicularSlopeX);
    
    %now we have to check where this lines crosses with the circle of the
    %original (big) magnitude of a single PSF. This should have 2 intersection
    %points!
    %Mathematical description of circle: y^2 + x^2 = r^2, assuming the circle
    %passees through the origin (0,0). r should be the magnitude of a single
    %PSF.
%     a = (perpendicularSlopeX.^2+1);
%     b = 2.*perpendicularSlopeX.*perpendicularIntersectionpointX;
%     c = -1.*maxMagn.^2+perpendicularIntersectionpointX.^2;
%     D = max(0,b.^2-4.*a.*c);%b.^2-4.*a.*c;%
%     x1 = (-b+sqrt(D))./(2.*a);
%     y1 = perpendicularSlopeX.*x1+perpendicularIntersectionpointX;
%     x2 = (-b-sqrt(D))./(2.*a);
%     y2 = perpendicularSlopeX.*x2+perpendicularIntersectionpointX;
    %new method (17-10-2018):
    %solving y^2+x^2=r^2 with y=ax+b result in:
    a = perpendicularSlopeX;
    b = perpendicularIntersectionpointX;
    r = maxMagn;
    x1 = -1*((sqrt((a.^2+1).*r.^2-b.^2)+a.*b)./(a.^2+1));
    x2 = ((sqrt((a.^2+1).*r.^2-b.^2)-a.*b)./(a.^2+1));
    y1 = x1.*a+b;
    y2 = x2.*a+b;
    x1 = real(x1);
    x2 = real(x2);
    y1 = real(y1);
    y2 = real(y2);
    
%     for i = 1:size(PSFarr,3)
%         if D(i) <= 0
%             x1(i) = 0.01;
%             y1(i) = 0.01;
%             x2(i) = 0.01;
%             y2(i) = 0.01;
%         end
%     end
%     x1 = real(x1);
%     x2 = real(x2);
%     y1 = real(y1);
%     y2 = real(y2);
    
    %Set the intersection points back to complex numbers
    origXPSF1 = permute(complex(x1,y1),[1 3 2]);
    origXPSF2 = permute(complex(x2,y2),[1 3 2]);
    
    %For both,
    %Calculate the angle of the Y-phasor by taking the first harmonic in Y
    angX = angle(origXPSF1);
    % Correct the angle
    angX(angX>0)=angX(angX>0)-2*pi;
    % Calculate Y position
    PositionXPSF1 = (abs(angX)/(2*pi/pixelsizephasor));
    %Calculate the angle of the Y-phasor by taking the first harmonic in Y
    angX = angle(origXPSF2);
    % Correct the angle
    angX(angX>0)=angX(angX>0)-2*pi;
    % Calculate Y position
    PositionXPSF2 = (abs(angX)/(2*pi/pixelsizephasor));
    
    %% Calculate distX and distY for all and plot
    % FOR NOW, PERFORM X, Y DISTANCE CALCULATION ONLY IF CENTRAL 5x5 PIXELS
    % CONTAIN >= 90% OF PHOTONS - ELSE IT'S HIGH CHANCE THERE IS A SECOND
    % EMITTER PRESENT
    startpos = (size(PSFarr,1)-3)/2+1;
    endpos = startpos+3;
    distX = zeros(size(PSFarr,3),1);
    distY = zeros(size(PSFarr,3),1);
    for i = 1:size(PSFarr,3)
%         medval = min(min(PSFarr(:,:,i)));
%         medvalfull = medval*(size(PSFarr,1)^2);
%         medvalpart = medval*25;
%         tempPSF = PSFarr(:,:,i)-medval;
%         if sum(sum(tempPSF(startpos:endpos,startpos:endpos))) >= 0.000000001*(sum(sum(tempPSF(:,:))))
            distX(i) = abs((PositionXPSF1(1,i)-PositionXPSF2(1,i)));
            distY(i) = abs((PositionYPSF1(1,i)-PositionYPSF2(1,i)));
%         else
%             distX(i) = 9;
%             distY(i) = 9;
%         end
    end
    
    %%
%     keyboard
%     dataarr = [FWHMX FWHMY angX angY perpendicularSlopeX perpendicularSlopeY perpendicularIntersectionpointX perpendicularIntersectionpointY angle(origXPSF1) angle(origXPSF2) angle(origYPSF1) angle(origYPSF2)];
end