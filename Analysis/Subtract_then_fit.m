function  fits=Subtract_then_fit(mov_fname,mov,movsz,...
    off_frames,moloffwin,guesses,roinum,dfrlmsz,MLE_fit,stdtol,...
    maxerr,do_avgsub,which_gaussian,fit_ang,usegpu,usephasor,phasorrad,allparams)
%% Subtract_mol_off_frames
% subtracts the average (or median) intensity of off frames for each guess
% stored in Mol_off_frames_fname.
%
% If you just want to do fitting, and not do background subtraction, set
% off_frames = 'nobgsub'. The program will take care of everything else.
%
%%%% Inputs %%%%
% mov_fname the filename of the movie
%
% mov is the movie data as a 3D array where the third dimension is the
% frame number.
%
% movsz is the output of size(mov)
%
% off_frames is the ouput from Mol_off_frames.mat which contains the off
% frames list for all guesses.
%
% moloffwin is the number of frames around the current frame to use for the
% BGSUB
%
% guesses is guesses array from Guessing.mat
%
% dfrlmsz is the  size of a diffraction limited spot in pixels. It's the
% nominal diameter, NOT the FWHM or something similar. Integer please
%
% MLE_fit  a Boolean determining whether or not MLE fitting is used. Set to
% 1 to use MLE and to 0 to use least squares. Note that MLE is quite slow,
% and so its not recommended for a large number of guesses
%
% stdtol is tolerance on fit Gaussian STD.
%
% maxerr is the maximum error of the fit for MLE fit, using variance
% default 0.1 (can't be above this) for LSQR fit, using the 95% confidence
% interval on the position
%
% do_avgsub is a Boolean determining whether or not to subtract the mean of
% the off frames. Set to 1 to subtract the mean and to 0 to subtract the
% median.
%
% which_gaussian determines what functional form of Gaussian function the
% molecules will be fit to if using least-squares fitting (MLE fitting only
% fits symmetric Gaussian). Set to 1 to use a symmetric Gaussian. Set to 2
% to use an asymmetric Gaussian (with angle determined by fit_ang). Set to
% 3 to use a freely rotating asymmetric Gaussian.
%
% fit_ang is the angle in degrees for an asymmetric Gaussian fit, see above
%
% usegpu is Boolean determining whether or not to use a CUDA enabled GPU
% for fitting if available.
%
%%%% Output %%%%
% a .mat file, importantly containing the fits structure that has fields
%frame number of the fit:
% fits.frame
%row coordinate of the fit:
% fits.row
%column coordinate of the fit:
% fits.col
%standard deviation in the row dimension of the Gaussian fit (if using a
%symmetric Gaussian this will be the same as the other width):
% fits.widthr
%standard deviation in the column dimension of the Gaussian fit (if using a
%symmetric Gaussian this will be the same as the other width):
% fits.widthc
%angle of asymmetric Gaussian fit:
% fits.ang
%offset of Gaussian fit:
% fits.offset
%amplitude of Gaussian fit:
% fits.amp
%error on fit (for MLE fitting, this is the variance, for least squares
%fitting, this is the mean 95% confidence interval on the position):
% fits.err
%sum of pixels in ROI around guess:
% fits.sum
%goodfit boolean:
% fits.goodfit
%
%%%% Dependencies %%%%
% TIFFStack 
% MLEwG (for MLE fitting) 
% gaussfit (for least squares fitting)
% gpufit
%
%     Copyright (C) 2018  Benjamin P Isaacoff
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
subnfit=tic;%for measuring the time to run the entire program

[pathstr,fname] = fileparts([mov_fname '.mat']);

disp([char(datetime),'   Fitting ',fname])

% plot_on is for debugging
plot_on=0;
% check if a GPU is available
if usegpu
    usegpu=parallel.gpu.GPUDevice.isAvailable;
end
if usegpu
    mov=single(mov);
end

if usephasor
    MLE_fit = 1; %for some archaic code - offset related
end

%check that the bgsub is actually happening
bgsub=1;
if strcmp(off_frames,'nobgsub');bgsub=0;end
if bgsub
    %check number of fits vs length of off frames
    if size(guesses,1)~=numel(off_frames);error('Unequal number of fits and number of off frames lists');end
end
%% The Averaging and Subtraction
%the conversion between dfrlmsz and the STD of the Gaussian, reccomended
%using the full width at 20% max given by (2*sqrt(2*log(5)))
dfD2std=(2*sqrt(2*log(5)));
%the guessed std
gesss=dfrlmsz/dfD2std;

% Use a CUDA enabled GPU to perform the fitting with GPUfit for significant
% improvement in fitting speed. Otherwise use the CPU to fit.

dataset=single(NaN((dfrlmsz*2+1)^2,size(guesses,1)));
initial_parameters=single(NaN(5,size(guesses,1)));
molr=guesses(:,2);
molc=guesses(:,3);
fits.frame=guesses(:,1);
fits.molid=(1:size(guesses,1))';
framelist=guesses(:,1);
offset=NaN(1,size(guesses,1));
if MLE_fit && usegpu
offset=NaN(1,size(guesses,1));
end
if usephasor
    shiftonedge = zeros(size(guesses,1),2); %list of shifts in pixels, if a localization is found on the edge. Used later by phasor to zoom in correctly.
end
%looping through all the guesses
if bgsub
    minxpos = NaN(1,size(guesses,1));
    maxxpos = NaN(1,size(guesses,1));
    minypos = NaN(1,size(guesses,1));
    maxypos = NaN(1,size(guesses,1));
    for ii=1:size(guesses,1)
        %current frame and molecule position
        curfrmnum=framelist(ii);
        curmolr=molr(ii);
        curmolc=molc(ii);
        frmlst=off_frames{ii};
        
        %the average (or median) frame
        %get x,y pos where the molecule is
        minxpos(ii) = curmolr-dfrlmsz;
        maxxpos(ii) = curmolr+dfrlmsz;
        minypos(ii) = curmolc-dfrlmsz;
        maxypos(ii) = curmolc+dfrlmsz;
        %if this is below zero or above max val, correct for it
        if minxpos(ii) <= 0
            if usephasor
                shiftonedge(ii,:) = [max(minxpos(ii)-1,phasorrad-((dfrlmsz*2+1)+1)/2+1), 0];
            end
			maxxpos(ii) = maxxpos(ii)+minxpos(ii)*-1+1;
			minxpos(ii) = 1;
        end
		if minypos(ii) <= 0
            if usephasor
                shiftonedge(ii,:) = [0, max(minypos(ii)-1,phasorrad-((dfrlmsz*2+1)+1)/2+1)];
            end
			maxypos(ii) = maxypos(ii)+minypos(ii)*-1+1;
			minypos(ii) = 1;
		end
		if maxxpos(ii) > movsz(1)
            if usephasor
                shiftonedge(ii,:) = [min((maxxpos(ii)-movsz(1)),((dfrlmsz*2+1)+1)/2-phasorrad-1), 0];
            end
			minxpos(ii) = minxpos(ii)-(maxxpos(ii)-movsz(1));
			maxxpos(ii) = movsz(1);
        end
		if maxypos(ii) > movsz(2)
            if usephasor
                %Shift on y-axis, with max =
                %(ROIsize+1)/2)-phasorrad-1, i.e. so that
                %midpoint+shiftonedge+phasorrad cannot be larger than
                %movsz(2). Same for others
                shiftonedge(ii,:) = [0, min((maxypos(ii)-movsz(2)),((dfrlmsz*2+1)+1)/2-phasorrad-1)];
            end
			minypos(ii) = minypos(ii)-(maxypos(ii)-movsz(2));
			maxypos(ii) = movsz(2);
        end
        if do_avgsub
			mean_mov=mean(single(mov(minxpos(ii):maxxpos(ii),minypos(ii):maxypos(ii),frmlst)),3);
            %mean_mov=mean(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3); Legacy
        else
            try
			mean_mov=median(single(mov(minxpos(ii):maxxpos(ii),minypos(ii):maxypos(ii),frmlst)),3);
            catch
                keyboard
            end
            %mean_mov=median(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3); Legacy
        end
        %the molecule image
        molim=single(mov(minxpos(ii):maxxpos(ii),minypos(ii):maxypos(ii),curfrmnum));
%         molim=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum)); Legacy
        %the subtracted image
        data=molim-mean_mov;
        data=reshape(data,[],1);
        gessb=min(data(:));
        gessN=range(data(:));
        if MLE_fit && ~usegpu
            %the guessed amplitude, using the formula in MLEwG
            gessN=range(data(:))*(4*pi*gesss^2);
        end
        if usegpu
            params0=[gessN;dfrlmsz;dfrlmsz;gesss;gessb];
        else
            params0=[dfrlmsz,dfrlmsz,gesss,gessb,gessN];
        end
        if MLE_fit && usegpu
            dataset(:,ii)=data+2*abs(min(data));
            offset(1,ii)=2*abs(min(data));
            params0(5)=abs(min(data));
            initial_parameters(:,ii)=params0;
        else
            initial_parameters(:,ii)=params0;
            dataset(:,ii)=data;
        end
    end
else
    for ii=1:size(guesses,1)
        %current frame and molecule position
        curfrmnum=framelist(ii);
        curmolr=molr(ii);
        curmolc=molc(ii);
        
        %the average (or median) frame
        %get x,y pos where the molecule is
        minxpos(ii) = curmolr-dfrlmsz;
        maxxpos(ii) = curmolr+dfrlmsz;
        minypos(ii) = curmolc-dfrlmsz;
        maxypos(ii) = curmolc+dfrlmsz;
        %if this is below zero or above max val, correct for it
        if minxpos(ii) <= 0
            if usephasor
                shiftonedge(ii,:) = [max(minxpos(ii)-1,phasorrad-((dfrlmsz*2+1)+1)/2+1), 0];
            end
			maxxpos(ii) = maxxpos(ii)+minxpos(ii)*-1+1;
			minxpos(ii) = 1;
        end
		if minypos(ii) <= 0
            if usephasor
                shiftonedge(ii,:) = [0, max(minypos(ii)-1,phasorrad-((dfrlmsz*2+1)+1)/2+1)];
            end
			maxypos(ii) = maxypos(ii)+minypos(ii)*-1+1;
			minypos(ii) = 1;
		end
		if maxxpos(ii) > movsz(1)
            if usephasor
                shiftonedge(ii,:) = [min((maxxpos(ii)-movsz(1)),((dfrlmsz*2+1)+1)/2-phasorrad-1), 0];
            end
			minxpos(ii) = minxpos(ii)-(maxxpos(ii)-movsz(1));
			maxxpos(ii) = movsz(1);
        end
		if maxypos(ii) > movsz(2)
            if usephasor
                %Shift on y-axis, with max =
                %(ROIsize+1)/2)-phasorrad-1, i.e. so that
                %midpoint+shiftonedge+phasorrad cannot be larger than
                %movsz(2). Same for others
                shiftonedge(ii,:) = [0, min((maxypos(ii)-movsz(2)),((dfrlmsz*2+1)+1)/2-phasorrad-1)];
            end
			minypos(ii) = minypos(ii)-(maxypos(ii)-movsz(2));
			maxypos(ii) = movsz(2);
        end
        
        try
            data = single(mov(minxpos(ii):maxxpos(ii),minypos(ii):maxypos(ii),curfrmnum));
%         data=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
        catch
            keyboard
        end
        data=reshape(data,[],1);
        gessb=min(data(:));
        gessN=range(data(:));
        if MLE_fit && ~usegpu
        %the guessed amplitude, using the formula in MLEwG
        gessN=range(data(:))*(4*pi*gesss^2);
        end
        if usegpu
            params0=[gessN;dfrlmsz;dfrlmsz;gesss;gessb];
        else
            params0=[dfrlmsz;dfrlmsz;gesss;gessb;gessN];
        end
        initial_parameters(:,ii)=params0;
        dataset(:,ii)=data;
    end
end
%% %% Fitting %%%%
%% If Gaussian (default)
if usephasor == false
    if usegpu
        %fitting with gpufit LSE fitting
        tolerance = 1e-4;
        % maximum number of iterations
        max_n_iterations = 1e4;
        % estimator id
        if MLE_fit
            estimator_id = EstimatorID.MLE;
        else
            estimator_id = EstimatorID.LSE;
        end
        % model ID
        if which_gaussian==1
            model_id = ModelID.GAUSS_2D;
            params_to_fit=[];
        elseif which_gaussian==2
            model_id=ModelID.GAUSS_2D_ROTATED;
            initial_parameters(6,:)=fit_ang;
            initial_parameters(5:7,:)=initial_parameters(4:6,:);
            params_to_fit=[1,1,1,1,1,1,0]';
        elseif which_gaussian==3
            model_id=ModelID.GAUSS_2D_ROTATED;
            params_to_fit=[];
            initial_parameters(6,:)=fit_ang;
            initial_parameters(5:7,:)=initial_parameters(4:6,:);
        end
        [parameters, states, chi_squares,~,~] = gpufit(dataset, [], ...
            model_id, initial_parameters, tolerance, max_n_iterations, params_to_fit, estimator_id, []);
        fits.amp=parameters(1,:)';
        if MLE_fit
            fits.offset=(parameters(5,:)-offset)';
        else
            fits.offset=parameters(5,:)';
        end
        fits.row=parameters(2,:)'-dfrlmsz+molr;
        fits.col=parameters(3,:)'-dfrlmsz+molc;
        fits.widthr=parameters(4,:)';
        fits.widthc=parameters(4,:)';
        if which_gaussian==1
            fits.ang=zeros(size(guesses,1),1);
        else
            fits.ang=parameters(6,:);
        end
        fits.err=(1-(chi_squares)./(sum((dataset-mean(dataset,1)).^2)))';
        fits.chi_squares=chi_squares';    
        if MLE_fit
            fits.err=(1-chi_squares./sum(2.*((mean(dataset,1)-dataset)-dataset.*log(mean(dataset,1)./dataset))))';
            errbad=fits.err<maxerr | states~=0;
        else
            errbad=fits.err<maxerr;
        end
        if MLE_fit && bgsub
            fits.sum=sum(dataset-offset,1)';
        else
            fits.sum=sum(dataset,1)';
        end
        fits.rowCI=sqrt(((fits.widthr.^2+1/12)./fits.sum)+(4*sqrt(pi()).*fits.widthr.^3.*fits.chi_squares)./(fits.sum.^2)); %Localization error based on Thompson, Larson, and Webb Biophys J. 2002 82 27752783. Equation 14
        fits.colCI=sqrt(((fits.widthc.^2+1/12)./fits.sum)+(4*sqrt(pi()).*fits.widthc.^3.*fits.chi_squares)./(fits.sum.^2)); %Where s is the gaussian width, a is the pixel size, N is the integrated intensity, and b^2 is the fit error (chi-squares), all spatial units are in pixels

        %determining if it's a goodfit or not (remember this field was
        %initialized to false)
        fits.goodfit=false(size(guesses,1),1);
        for ii=1:size(guesses,1)
            if (mean([fits.widthr(ii),fits.widthc(ii)])<=(stdtol*gesss) && mean([fits.widthr(ii),fits.widthc(ii)])>=(gesss/stdtol)) && ... %Compare width with diffraction limit
                    ~errbad(ii) && ... %too much error on fit?
                    fits.amp(ii)<fits.sum(ii) && ... %the amplitude of the fit shouldn't be bigger than the integral
                    ~any([fits.row(ii),fits.col(ii),fits.amp(ii),fits.sum(ii)]<0) && ... %none of the fitted parameters should be negative, except the offset!
                    fits.rowCI(ii)<=dfrlmsz && fits.colCI(ii)<=dfrlmsz %none of the localization errors are larger than the gaussian widths
                fits.goodfit(ii)=true;%goodfit boolean
            end
        end
        fits.states=states';
        %plotting for debugging/tests
        if plot_on
            for kk=1:size(guesses,1)
                curfrmnum=framelist(kk);
                curmolr=molr(kk);
                curmolc=molc(kk);
                frmlst=off_frames{kk};

                if bgsub
                    %the average (or median) frame
                    if do_avgsub
                        mean_mov=mean(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3);
                    else
                        mean_mov=median(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3);
                    end
                    %the molecule image
                    molim=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
                    %the subtracted image
                    data=molim-mean_mov;

                    h12=figure(12);
                    subplot(1,4,1)
                    imshow(mean_mov,[])
                    title('Mean BG')
                    subplot(1,4,2)
                    imshow(molim,[])
                    title('Raw Molecule')
                    subplot(1,4,3)
                    imshow(data,[])
                    title('BGSUB')

                    [x, y] = ndgrid(0:size(data,1)-1,0:size(data,2)-1);
                    fitim=gaussian_2d(x, y, parameters(:,kk));
                    subplot(1,4,4)
                    imshow(fitim,[])
                    title('Fit Profile')
                    annotation('textbox', [0 0.9 1 0.1], ...
                        'String', ['Frame # ',num2str(curfrmnum),'   Guess # ',num2str(kk),' R^2=',num2str(fits.err(kk))], ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center')

                    keyboard
                    try
                        close(h12)
                    catch
                    end
                else
                    data=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
                    h12=figure(12);
                    subplot(1,2,1)
                    imshow(data,[])
                    title('Molecule Image')
                    [x, y] = ndgrid(size(data));
                    fitim=gaussian_2d(x, y, parameters(:,kk));
                    subplot(1,2,2)
                    imshow(fitim,[])
                    title('Fit Profile')
                    annotation('textbox', [0 0.9 1 0.1], ...
                        'String', ['Frame # ',num2str(curfrmnum),'   Guess # ',num2str(kk),' R^2=',num2str(fits.err(kk))], ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center')

                    keyboard
                    try
                        close(h12)
                    catch
                    end
                end
            end
        end


    else
        %initializing the fits structure
        fits.row=NaN(size(guesses,1),1);%row coordinate of the fit
        fits.col=NaN(size(guesses,1),1);%column coordinate of the fit
        fits.widthr=NaN(size(guesses,1),1);%standard deviation in the row dimension of the Gaussian fit
        fits.widthc=NaN(size(guesses,1),1);%standard deviation in the column dimension of the Gaussian fit
        fits.ang=NaN(size(guesses,1),1);%angle of asymmetric Gaussian fit
        fits.offset=NaN(size(guesses,1),1);%offset
        fits.amp=NaN(size(guesses,1),1);%amplitude of Gaussian fit
        fits.err=NaN(size(guesses,1),1);%error on fit
        fits.sum=sum(dataset,1)';%sum of pixels in ROI around guess
        fits.goodfit=false(size(guesses,1),1);%goodfit boolean
        goodfit=false(size(guesses,1),1);%goodfit boolean
        sumsum=fits.sum;
        fit_sd_r=NaN(size(guesses,1),1);
        fit_sd_rCI=NaN(size(guesses,1),1);
        fit_sd_c=NaN(size(guesses,1),1);
        fit_sd_cCI=NaN(size(guesses,1),1);
        fit_off=NaN(size(guesses,1),1);
        fit_offCI=NaN(size(guesses,1),1);
        fit_amp=NaN(size(guesses,1),1);
        fit_ampCI=NaN(size(guesses,1),1);
        fit_err=NaN(size(guesses,1),1);
        act_r=NaN(size(guesses,1),1);
        act_rCI=NaN(size(guesses,1),1);
        fit_ang=NaN(size(guesses,1),1);
        fit_angCI=NaN(size(guesses,1),1);
        act_c=NaN(size(guesses,1),1);
        act_cCI=NaN(size(guesses,1),1);
%         keyboard
        for ii=1:size(guesses,1) %
            if MLE_fit
                %fitting with MLE
                [paramsF,varianceF] = MLEwG (reshape(dataset(:,ii),[2*dfrlmsz+1,2*dfrlmsz+1]),double(initial_parameters(:,ii)'),1,plot_on,1);
                %shifting
                paramsF([1,2])=paramsF([1,2])+0.5;
                fit_r=paramsF(1);fit_c=paramsF(2);
                fit_sd_r(ii)=paramsF(3);fit_sd_c(ii)=paramsF(3);
                %recalculating the values based on their equations to match
                paramsF(5)=paramsF(5)/(2*pi*paramsF(3)^2);
                if paramsF(4)>=0
                    paramsF(4)=sqrt(paramsF(4));
                else
                    paramsF(4)=-sqrt(-paramsF(4));
                end
                fit_off(ii)=paramsF(4);
                fit_amp(ii)=paramsF(5);
                fit_ang(ii)=0;
                fit_err(ii)=varianceF;
                errbad=varianceF>maxerr;%too much error on fit?
            else
                %fitting with least squares
                [fitPars,conf95,~,~,resid]=gaussFit(double(reshape(dataset(:,ii),[2*dfrlmsz+1,2*dfrlmsz+1])),'searchBool',0,'nPixels',2*dfrlmsz+1,...
                    'checkVals',0,'ffSwitch',which_gaussian);
                %converting the variables to match the output of MLEwG, and
                %arranging for each particular Gaussian fit
                fit_r=fitPars(1);fit_c=fitPars(2);
                if which_gaussian==1
                    fit_sd_r(ii)=fitPars(3);fit_sd_c(ii)=fitPars(3);
                    fit_off(ii)=fitPars(5);
                    fit_amp(ii)=fitPars(4);
                    fit_ang(ii)=0;
                    fit_sd_rCI(ii)=conf95(3);fit_sd_cCI(ii)=conf95(3);
                    fit_offCI(ii)=conf95(5);
                    fit_ampCI(ii)=conf95(4);
                    fit_angCI(ii)=0;
                elseif which_gaussian==2
                    fit_sd_r(ii)=fitPars(3);fit_sd_c(ii)=fitPars(4);
                    fit_sd_rCI(ii)=conf95(3);fit_sd_cCI(ii)=conf95(4);
                    fit_off(ii)=fitPars(6);
                    fit_offCI(ii)=conf95(6);
                    fit_amp(ii)=fitPars(5);
                    fit_ampCI(ii)=conf95(5);
                    fit_ang(ii)=0;
                    fit_angCI(ii)=0;
                elseif   which_gaussian==3
                    fit_sd_r(ii)=fitPars(4);fit_sd_c(ii)=fitPars(5);
                    fit_sd_rCI(ii)=conf95(4);fit_sd_cCI(ii)=conf95(5);
                    fit_off(ii)=fitPars(7);
                    fit_offCI(ii)=conf95(7);
                    fit_amp(ii)=fitPars(6);
                    fit_ampCI(ii)=conf95(6);
                    fit_ang(ii)=fitPars(3);
                    fit_angCI(ii)=conf95(3);
                end
                fit_err(ii)=1-(sum(resid.^2)/sum((dataset(:,ii)-mean(dataset(:,ii))).^2));
                errbad=fit_err(ii)<maxerr;%too much error on fit?
            end
            %Convert back into full frame coordinates, NOTE the -1!
            act_r(ii)=fit_r-dfrlmsz-1+molr(ii);
            act_c(ii)=fit_c-dfrlmsz-1+molc(ii);
            try
                act_rCI(ii)=conf95(1);
                act_cCI(ii)=conf95(2);
            catch
                act_rCI(ii)=0;
                act_cCI(ii)=0;
            end
            if (mean([fit_sd_r(ii),fit_sd_c(ii)])<=(stdtol*gesss) && mean([fit_sd_r(ii),fit_sd_c(ii)])>=(gesss/stdtol)) && ... %Compare width with diffraction limit
                    ~errbad && ... %too much error on fit?
                    fit_amp(ii)<sumsum(ii) && ... %the amplitude of the fit shouldn't be bigger than the integral
                    ~any([fit_r,fit_c,fit_amp(ii),sumsum(ii)]<0) %none of the fitted parameters should be negative, except the offset!

                goodfit(ii)=true;%goodfit boolean
            else
                goodfit(ii)=false;
            end
            %The sum(:) of the the data
        end

        %putting the fit results into the fits structure
        fits.row=act_r;%row coordinate of the fit
        fits.rowCI=act_rCI;%row coordinate confidence interval of the fit
        fits.col=act_c;%column coordinate of the fit
        fits.colCI=act_cCI;%column coordinate of the fit
        fits.widthr=fit_sd_r;%standard deviation in the row dimension of the Gaussian fit
        fits.widthrCI=fit_sd_rCI;%standard deviation in the row dimension of the Gaussian fit
        fits.widthc=fit_sd_c;%standard deviation in the column dimension of the Gaussian fit
        fits.widthcCI=fit_sd_cCI;%standard deviation in the column dimension of the Gaussian fit
        fits.ang=fit_ang;%angle of asymmetric Gaussian fit
        fits.angCI=fit_angCI;%Confidence interval of angle of asymmetric Gaussian fit
        fits.offset=fit_off;%offset
        fits.offsetCI=fit_offCI;%Confidence interval of offset
        fits.amp=fit_amp;%amplitude of Gaussian fit
        fits.ampCI=fit_ampCI;%Confidence interval of amplitude of Gaussian fit
        fits.err=fit_err;%error on fit
        fits.goodfit=goodfit;%determining if it's a goodfit or not (remember this field was
        %initialized to false)


        %plotting for debugging/tests
        if plot_on
            for ii=1:size(guesses,1)
                curfrmnum=framelist(ii);
                curmolr=molr(ii);
                curmolc=molc(ii);
                frmlst=off_frames{ii};

                if bgsub
                    %the average (or median) frame
                    if do_avgsub
                        mean_mov=mean(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3);
                    else
                        mean_mov=median(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3);
                    end
                    %the molecule image
                    molim=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
                    %the subtracted image
                    data=molim-mean_mov;

                    h12=figure(12);
                    subplot(1,4,1)
                    imshow(mean_mov,[])
                    title('Mean BG')
                    subplot(1,4,2)
                    imshow(molim,[])
                    title('Raw Molecule')
                    subplot(1,4,3)
                    imshow(data,[])
                    title('BGSUB')
                    [x, y] = ndgrid(1:size(data,1),1:size(data,2));
                    fitim=gaussian_2d(x, y, [fits.amp(ii),fits.row(ii),fits.col(ii),fits.widthc(ii),fits.offset(ii)]');
                    subplot(1,4,4)
                    imshow(fitim,[])
                    title('Fit Profile')
                    annotation('textbox', [0 0.9 1 0.1], ...
                        'String', ['Frame # ',num2str(curfrmnum),'   Guess # ',num2str(ii),' R^2=',num2str(fits.err)], ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center')

                    keyboard
                    try
                        close(h12)
                    catch
                    end
                else
                    data=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
                    h12=figure(12);
                    subplot(1,2,1)
                    imshow(data,[])
                    title('Molecule Image')
                    [x, y] = ndgrid(size(data));
                    fitim=gaussian_2d(x, y, parameters(:,ii));
                    subplot(1,2,2)
                    imshow(fitim,[])
                    title('Fit Profile')
                    annotation('textbox', [0 0.9 1 0.1], ...
                        'String', ['Frame # ',num2str(curfrmnum),'   Guess # ',num2str(ii),' R^2=',num2str(fits.err(ii))], ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center')

                    keyboard
                    try
                        close(h12)
                    catch
                    end
                end
            end
        end
    end
%% Phasor fitting
else
    tic
    %Parameters is a 5xN array with values:
    %amplitude
    %posX
    %posY
    %widthXY (identical appearantly?)
    %offset
    %States is a 1xN array, 0 is good fit, 1 is bad fit??? - appears to
    %be 1 if multiple psfs are detected
    %Chi_squares is a 1xN array with X2-values? Used to determine bad
    %fits as well, if it's too high
    %Make dataset smaller - zoom in on center based on phasorrad and
    %shiftonedge
    phasordataset = reshape(dataset,sqrt(size(dataset(:,1),1)),sqrt(size(dataset(:,1),1)),size(dataset,2));
    midval = ones(size(phasordataset,3),2).*[(sqrt(size(dataset(:,1),1))+1)/2 (sqrt(size(dataset(:,1),1))+1)/2] ...
        + shiftonedge;
    phasorradsp = allparams.phasor_SPTPradius;
    phasordatasetsm = zeros(phasorrad*2+1,phasorrad*2+1,size(phasordataset,3));
    phasordatasetsp = zeros(phasorradsp*2+1,phasorradsp*2+1,size(phasordataset,3));
    %Next 3 lines might be simplified to 1 line? Can't get it to work
    %atm
    for ii = 1:size(phasordataset,3)
        try
            phasordatasetsm(:,:,ii) = phasordataset([midval(ii,1)-phasorrad:midval(ii,1)+phasorrad],[midval(ii,2)-phasorrad:midval(ii,2)+phasorrad],ii);phasordatasetsm(:,:,ii) = phasordataset([midval(ii,1)-phasorrad:midval(ii,1)+phasorrad],[midval(ii,2)-phasorrad:midval(ii,2)+phasorrad],ii);
        catch
            keyboard
        end
    end
    %For now, give posX and posY, leave rest at approximations
    %Reshape-function makes 1D array to 2D array
    if allparams.phasor_2d || allparams.phasor_DH
        [phasorlocs,phasormagnitudes,badfits,phasoramp,phasoroffset] = Phasor_localization_SMALLLABS(phasordatasetsm,0,''); %2D
    elseif allparams.phasor_astig
        [phasorlocs,phasormagnitudes,badfits,phasoramp,phasoroffset] = Phasor_localization_SMALLLABS(phasordatasetsm,1,allparams.phasor_astigcal); %astig
    elseif allparams.phasor_SP || allparams.phasor_TP
        [phasorlocs,phasormagnitudes,badfits,phasoramp,phasoroffset] = Phasor_localization_SMALLLABS(phasordatasetsm,0,'');
%         SP_TP_phasor_linking([fits.frame phasorlocs],11.5,1)
                % cleanup below here
                parameters(2:3,:) = phasorlocs(:,[2 1])'; %removing 3rd column from phasorxy - this is z pos
                %fits.amp will be max of the intensity via phasor, corrected for
                %median of background
                fits.amp=phasoramp;
                %fits.offset is based on MLE fitting. There, the offset array is
                %subtracted from the offsetMLE. offsetMLE is linearly related to
                %offsetPhasor, and on average 1.28x higher. Thus, I correct for
                %that (no idea why)
                fits.offset = phasoroffset/1.2836-offset';
                %x, y pos
%                 fits.row=parameters(2,:)'+minxpos'+(dfrlmsz-phasorrad)+shiftonedge(:,1);%-phasorrad+molr-shiftonedge(:,1)*0;
%                 fits.col=parameters(3,:)'+minypos'+(dfrlmsz-phasorrad)+shiftonedge(:,2);%-phasorrad+molc-shiftonedge(:,2)*0;
                fits.row=parameters(2,:)'+minxpos'+(dfrlmsz-phasorrad*2)+shiftonedge(:,1);%-phasorrad+molr-shiftonedge(:,1)*0;
                fits.col=parameters(3,:)'+minypos'+(dfrlmsz-phasorrad*2)+shiftonedge(:,2);%-phasorrad+molc-shiftonedge(:,2)*0;
%                 figure(1);clf(1);
%                 scatter(fits.frame,phasormagnitudes(:,1)./phasormagnitudes(:,2))
%                 axis([0 inf 0 3])
                phasordatapointstp = SP_TP_phasor_linking([fits.frame fits.row fits.col],11.5);
                fits.row = [];
                fits.col = [];
                fits.frame = [];
                fits.frame = phasordatapointstp(:,1);
                fits.row = phasordatapointstp(:,3);%switching x and y?
                fits.col = phasordatapointstp(:,2);
                %set some values to 0 because they're unused and array size
                %changes
                shiftonedge = zeros(size(fits.frame,1),2);
                minxpos = zeros(1,size(fits.frame,1));
                minypos = zeros(1,size(fits.frame,1));
                fits.molid = [1:size(fits.frame,1)]';
                fits.amp = zeros(1,size(fits.frame,1));
                fits.offset = zeros(1,size(fits.frame,1));
                clear phasordatasetsp
                %Get bigger loc sizes
                %loop over localizations
%                 keyboard
%                 tic
%                 goodlist = [];
%                 for locs = 1:size(phasordatapointstp,1)
%                     try
%                     phasordatasetsp(:,:,locs) = mov([round(phasordatapointstp(locs,2))-phasorradsp:round(phasordatapointstp(locs,2))+phasorradsp],...
%                     [round(phasordatapointstp(locs,3))-phasorradsp:round(phasordatapointstp(locs,3))+phasorradsp],fits.frame(locs));
%                     goodlist = [goodlist; locs];
%                     catch
% %                         fits.frame(locs) = 1
%                     end
%                 end
%                 toc
                
                goodlist = [1:size(phasordatapointstp,1)]';
                for locs = 1:size(phasordatapointstp,1)
                    try
                    phasordatasetsp(:,:,locs) = mov([round(phasordatapointstp(locs,2))-phasorradsp:round(phasordatapointstp(locs,2))+phasorradsp],...
                    [round(phasordatapointstp(locs,3))-phasorradsp:round(phasordatapointstp(locs,3))+phasorradsp],fits.frame(locs));
%                     goodlist = [goodlist; locs];
                    catch
                        % it's too much on an edge, so remove it
                        goodlist(locs) = -1;
%                         fits.frame(locs) = 1
                    end
                end
                goodlist(goodlist==-1) = [];
                fits.frame = fits.frame(goodlist);
                fits.molid = [1:size(goodlist,1)]';
                fits.amp = fits.amp(goodlist);
                fits.offset = fits.offset(goodlist);
                fits.row = fits.row(goodlist);
                fits.col = fits.col(goodlist);
                clear parameters
                parameters = zeros(3,size(fits.frame,1));
                parameters(1,:) = fits.frame;
                phasorlocs = [fits.row fits.col];
                minxpos = minxpos(goodlist);
                minypos = minypos(goodlist);
                shiftonedge = shiftonedge(goodlist,:);
                
                phasordatasetsp(:,:,size(phasordatapointstp,1)+1:end) = [];
                phasordatapointstp(phasordatasetsp(1,1,:)==0,:) = [];
                phasordatasetsp(:,:,phasordatasetsp(1,1,:)==0) = [];
                
                [distX distY posX1 posX2 posY1 posY2 posXtot posYtot] = ctPhasor(phasordatasetsp);
                if allparams.phasor_SP
                    load(allparams.phasor_SPcal);
                elseif allparams.phasor_TP
                    load(allparams.phasor_TPcal);
                end
                [fits.zpos] = SPTP_phasor_zposcalculationRoutine(distX, distY, SPTP_phasor_curveX, SPTP_phasor_curveY);
    end
    parameters(2:3,:) = phasorlocs(:,[2 1])'; %removing 3rd column from phasorxy - this is z pos
    %x, y pos
    fits.row=parameters(2,:)'+minxpos'+(dfrlmsz-phasorrad)+shiftonedge(:,1);%-phasorrad+molr-shiftonedge(:,1)*0;
    fits.col=parameters(3,:)'+minypos'+(dfrlmsz-phasorrad)+shiftonedge(:,2);%-phasorrad+molc-shiftonedge(:,2)*0;
    if allparams.phasor_astig
        fits.zpos = phasorlocs(:,3)/1000; %nm to um
    end
     if allparams.phasor_SP || allparams.phasor_TP
        fits.zpos = fits.zpos/1000; %nm to um
    end
    fits.MagnX = phasormagnitudes(:,1);
    fits.MagnY = phasormagnitudes(:,2);
    %badfits = 1 if bad fit
    %fits.amp will be max of the intensity via phasor, corrected for
    %median of background
    fits.amp=phasoramp;
    %fits.offset is based on MLE fitting. There, the offset array is
    %subtracted from the offsetMLE. offsetMLE is linearly related to
    %offsetPhasor, and on average 1.28x higher. Thus, I correct for
    %that (no idea why)
    fits.offset = phasoroffset/1.2836-offset';
%         keyboard
    
    %set fwhm x and y to 0 - no info from phasor
    fits.widthr=zeros(size(guesses,1),1);
    fits.widthc=zeros(size(guesses,1),1);
    %set angle to 0 - no info from phasor
    fits.ang=zeros(size(guesses,1),1);
    %fits.err and fits.chi_squares also to 0 via phasor
    fits.err=zeros(size(guesses,1),1);
    fits.chi_squares=zeros(size(guesses,1),1);
    %errbad is fully in badfits from phasor
    errbad = badfits;
    if bgsub
        fits.sum=sum(dataset-offset,1)';
    else
        fits.sum=sum(dataset,1)';
    end
    %CI =0 for phasor
    fits.rowCI=zeros(size(guesses,1),1);
    fits.colCI=zeros(size(guesses,1),1);
    
    %determining if it's a goodfit or not (remember this field was
    %initialized to false)
    fits.goodfit=false(size(guesses,1),1);
    for ii=1:size(guesses,1)
        if ~errbad(ii)
            fits.goodfit(ii)=true;%goodfit boolean
        end
    end
    fits.states=errbad;
    
    
    %% SP phasor - after 2D phasor
    if allparams.phasor_SP || allparams.phasor_TP
        %Re-make all of fits structure I guess?
        fits.molid = [1:size(parameters,2)]';
        particles = fits.molid.*0;
        fits.amp = zeros(size(particles,1),1);
        fits.offset = zeros(size(particles,1),1);
        fits.widthr = zeros(size(particles,1),1);
        fits.widthc = zeros(size(particles,1),1);
        fits.ang = zeros(size(particles,1),1);
        fits.err = zeros(size(particles,1),1);
        fits.chi_squares = zeros(size(particles,1),1);
        fits.sum = zeros(size(particles,1),1);
        fits.rowCI = zeros(size(particles,1),1);
        fits.colCI = zeros(size(particles,1),1);
        fits.goodfit = boolean(ones(size(particles,1),1));
        fits.states = zeros(size(particles,1),1);
    end
    %% DH phasor - after 2D phasor
    if allparams.phasor_DH %Double Helix fitting
%         keyboard
%         pSMLM_DH_calibration(input,zposcali,imagegen)
% keyboard
%         pSMLM_DH_calibration_SL([fits.frame fits.col fits.row fits.MagnY fits.MagnX],[-0.75:0.01:0.75],1,'H:\MATLAB_BIPNAS\SMALLLABS\phasor\DHcalib_ph5rad_fullim.mat')
        %Load calibration data
        load(allparams.phasor_DHcal);
        loc_list = [fits.frame fits.col fits.row fits.MagnY fits.MagnX]; %Transposed x,y due to matlab shenanigans
%         keyboard
%         tic
        [particles, possibilities] = DH_phasor_linking(loc_list,curveAngle,curveDistAngle,[],wobbleMatrix);
%         toc
%         [t] = DH_phasor_linking(loc_list,curveAngle,curveDistAngle,[],wobbleMatrix)
        %% Testing some things
        %loop over frames
        testingmode = 0;
        if testingmode
            clear fwhms
            pxsize = 23;
            for frame = 1:movsz(3)
                %get ROI to encompass both locs + half diff-lim-spot
                tt = loc_list(loc_list(:,1)==frame,:);
                tp = possibilities{frame};
                tl = particles(particles(:,1)==frame,:);
                c = 1;
                for p = 1:size(tp,1)
                    if ~isempty(cell2mat(tp(p,2)))
                        if cell2mat(tp(p,1))<cell2mat(tp(p,2))
                            pos1 = [tt(cell2mat(tp(p,1)),2), tt(cell2mat(tp(p,1)),3)];
                            pos2 = [tt(cell2mat(tp(p,2)),2), tt(cell2mat(tp(p,2)),3)];
                            minx = (min(pos1(1),pos2(1))-dfrlmsz*1.5);
                            maxx = (max(pos1(1),pos2(1))+dfrlmsz*1.5);
                            miny = (min(pos1(2),pos2(2))-dfrlmsz*1.5);
                            maxy = (max(pos1(2),pos2(2))+dfrlmsz*1.5);
                            meanx = round(minx+(maxx-minx)/2);
                            meany = round(miny+(maxy-miny)/2);
                            linkedparticle = [];
                            if size(particles(particles(:,1)==frame,:),1) > 0
                                if size(particles(particles(:,1)==frame,:),1) > 1
                                    [a,b] = min(abs(particles(particles(:,1)==frame,2:3)-[meanx meany]));
                                    b = min(b);
                                    t = particles(particles(:,1)==frame,:);
                                    linkedparticle = t(b,:);
                                else
                                    t = particles(particles(:,1)==frame,:);
                                    linkedparticle = t(1,:);
                                end
                            end
                            fullROI = mov(:,:,frame);
                            try
                                DHROI = mov((meany-(pxsize-1)/2):(meany+(pxsize-1)/2),(meanx-(pxsize-1)/2):(meanx+(pxsize-1)/2),frame);
                                pxsize(frame,c) = size(DHROI,1);
                                [loc,fwhms(frame,c,:)] = Phasor_localization_SMALLLABS(DHROI,0,'');
                                zposarr(frame,c) = linkedparticle(1,4);
                                c = c+1;
                                %                             figure(2);clf(2);
                                %                             imagesc(DHROI);
                            catch
                            end
                        end
                    end
                end
            end
            %         figure(4);clf(4);
            %         hold on
            %         for i = 1:6
            %             scatter([1:151],fwhms(:,i,1)./fwhms(:,i,2))
            %         end
            bins = [-1:0.01:1];
            zz = reshape(zposarr,size(zposarr,1)*size(zposarr,2),1);
            zf = reshape(fwhms(:,:,1)./fwhms(:,:,2),size(zposarr,1)*size(zposarr,2),1);
            zz(isnan(zf)) = [];
            zf(isnan(zf)) = [];
            arrcount = zeros(size(bins,2),1);
            zfarr = zeros(size(bins,2),1);
            zfarrcomp = zeros(size(bins,2),size(zz,1));
            for iii = 1:size(zz,1)
                [~,curbin] = find(zz(iii)>bins);
                curbin = max(curbin);
                if (curbin >= 0) & (curbin <= size(bins,2))
                    zfarr(curbin) = zfarr(curbin)+zf(iii);
                    zfarrcomp(curbin,iii) = zf(iii);
                    arrcount(curbin) = arrcount(curbin)+1;
                end
            end
            zfarrmean = zfarr./arrcount;
            for iii = 1:size(bins,2)
                ttt = zfarrcomp(iii,:);
                ttt(ttt==0) = [];
                zfarrstd(iii) = std(ttt,0,2);
                zfarrsize(iii) = size(ttt,2);
            end
            bins(isnan(zfarrmean)) = [];
            zfarrstd(isnan(zfarrmean)) = [];
            zfarrmean(isnan(zfarrmean)) = [];
            
            figure(5);clf(5);
            hold on
            errorbar(bins,zfarrmean,zfarrstd)
            axis([-0.8 0.8 0.0 12])
        end
%%
        %Re-make all of fits structure I guess?
        fits.frame = particles(:,1);
        fits.molid = [1:size(particles,1)]';
        fits.amp = zeros(size(particles,1),1);
        fits.offset = zeros(size(particles,1),1);
        fits.row = particles(:,3); %back-transpose x and y due to matlab shenanigans
        fits.col = particles(:,2); %back-transpose x and y due to matlab shenanigans
        fits.zpos = particles(:,4);
        fits.widthr = zeros(size(particles,1),1);
        fits.widthc = zeros(size(particles,1),1);
        fits.ang = zeros(size(particles,1),1);
        fits.err = zeros(size(particles,1),1);
        fits.chi_squares = zeros(size(particles,1),1);
        fits.sum = zeros(size(particles,1),1);
        fits.rowCI = zeros(size(particles,1),1);
        fits.colCI = zeros(size(particles,1),1);
        fits.goodfit = boolean(ones(size(particles,1),1));
        fits.states = zeros(size(particles,1),1);
    end
    toc;
end
%% Resume fitting code
tictoc=toc(subnfit);%the time to run the entire program
%save the data
fits.roinum=roinum;
if bgsub
    fname=[pathstr,filesep,fname,'_AccBGSUB_fits.mat'];
else
    fname=[pathstr,filesep,fname,'_fits.mat'];
end
save(fname,'fits','MLE_fit','stdtol','maxerr','dfrlmsz','movsz','moloffwin',...
    'tictoc','do_avgsub','which_gaussian','-v7.3')
end

function g = gaussian_2d(x, y, p)
% Generates a 2D Gaussian peak.
% http://gpufit.readthedocs.io/en/latest/api.html#gauss-2d
%
% x,y - x and y grid position values p - parameters (amplitude, x,y center
% position, width, offset)

g = p(1) * exp(-((x - p(2)).^2 + (y - p(3)).^2) / (2 * p(4)^2)) + p(5);

end
