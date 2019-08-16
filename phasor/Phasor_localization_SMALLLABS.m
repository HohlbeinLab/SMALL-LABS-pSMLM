function [loc,magnitudes,badfit,amplitude,mediannoise] = Phasor_localization_SMALLLABS(PSFdata,enable3D,calfile)
global verbose
if verbose
disp([char(datetime),'   Running phasor localization'])
end

loc = zeros(size(PSFdata,3),3);
magnitudes = zeros(size(PSFdata,3),2);
badfit = zeros(size(PSFdata,3),1);
tic;
    pixelsizephasor = size(PSFdata,1);
    %Do fourier transformation of images
    phasorfft = fft2(PSFdata);
    %Calculate the angle of the X-phasor by taking the first harmonic in X
    angX = permute(angle(phasorfft(1,2,:)),[3 1 2]);
    %Correct the angle
    angX(angX>0)=angX(angX>0)-2*pi;
    %Calculate X position
    PositionX = (abs(angX)/(2*pi/pixelsizephasor));
    FWHMX = abs(phasorfft(1,2,:));
    %Calculate the angle of the Y-phasor by taking the first harmonic in Y
    angY = permute(angle(phasorfft(2,1,:)),[3 1 2]);
    %Correct the angle
    angY(angY>0)=angY(angY>0)-2*pi;
    %Calculate Y position
    PositionY = (abs(angY)/(2*pi/pixelsizephasor));
    FWHMY = abs(phasorfft(2,1,:));
    %save data
    loc(:,1)=PositionX;
    loc(:,2)=PositionY;
    magnitudes(:,1) = FWHMX;
    magnitudes(:,2) = FWHMY;
    
    %Find the noise or signal based on an aperture method
    %First, copy-translate from java
    centersize = pixelsizephasor-4;
    ringsize = 2;
    noisearray = ones(pixelsizephasor,pixelsizephasor);
    signalarray = zeros(pixelsizephasor,pixelsizephasor);
    for xx = 1:pixelsizephasor
        for yy = 1:pixelsizephasor
            %get distance to center
            dist = sqrt(abs(xx-(pixelsizephasor+1)/2)^2+abs(yy-(pixelsizephasor+1)/2)^2);
            distarr(xx,yy) = dist;
            if dist > (centersize+1)/2
                if dist <= ringsize+(centersize+1)/2
                    noisearray(xx,yy) = 1;
                else
                    noisearray(xx,yy) = 0;
                end
            else
                noisearray(xx,yy) = 0;
            end
            if dist <= (centersize+1)/2
                signalarray(xx,yy) = 1;
            else
                signalarray(xx,yy) = 0;
            end
        end
    end
    
    %Create 1D arrays of noise, signal (both Boolean), and data
    noisearr1d = reshape(noisearray,pixelsizephasor^2,1);
    signalarr1d = reshape(signalarray,pixelsizephasor^2,1);
    PSFdata1d = reshape(PSFdata,pixelsizephasor^2,size(PSFdata,3));
    
    %Get median, std of noise, and sum of signal, only using the pixels of
    %the 1d arrays of noise and signal.
    mediannoise = median(PSFdata1d(noisearr1d==1,:))';
    stdnoise = std(PSFdata1d(noisearr1d==1,:),0,1)';
    sumsignal = sum(PSFdata1d(signalarr1d==1,:),1)';
    amplitude = max(PSFdata1d(signalarr1d==1,:))'-mediannoise;
    
    %Check goodness of fit:
    %if sumsignal is likely to be the same as the median noise level, its a
    %bad psf
    %real ratio: sumsignal should be larger than
    %(mediannoise+stdnoise*.5)*nrpixelssignal
    intensityratio = sumsignal./((mediannoise+stdnoise*.5)*sum(sum(signalarray)));
    
    magnratio = permute(FWHMX./FWHMY,[3 1 2]);
    
    badfit = magnratio>2 | magnratio<0.5 | ...;  %badfit if ratio of magnitudes is way off;
            intensityratio < 1; %or if intensity isn't actually significantly(ish) higher than median of noise
    if enable3D == 1
        if string(calfile) ~=''
            loc(:,3) = phasorcalcZpos(permute([FWHMX FWHMY],[3 2 1]),calfile);
        else
            fprintf('No phasor calibration file found!\n');
        end
    end
time = toc;
end