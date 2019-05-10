
%     %% Perform cross-correlation to account for drift effects
%         if verbose
%             disp([char(datetime),'   Starting cross-correlation drift correction'])
%         end
%         %load fits
%         if params.bgsub
%             %try loading in the fits & tracking results
%             load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits');
%         else
%             load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits');
%         end
%         %change params so it appears to be 2D
%         paramscc = params;
%         paramscc.avgshifthist_axialcolnrs = 1;
%         totframes = max((fits.frame));
%         %Number of bins used for cross-correlation - should be user-entered
%         nrccbins = 10;
%         clear ccim
%         % Create an avg shift image consisting of all positions within the
%         % specified frames (specified by amount of cc-bins)
%     for cc_bins = 1:nrccbins %for all the bins (frametime)
%         %Calculate min and max frame
%         minframe = floor(totframes/nrccbins)*(cc_bins-1)+1;
%         maxframe = floor(totframes/nrccbins)*(cc_bins);
%         %Find number of localizations
%         nrlocsincc = size(fits.frame(fits.frame>=minframe & fits.frame<=maxframe),1);
%         %Make structures with x,y,goodfit params
%         fitscc = {};
%         fitscc.row = fits.row(fits.frame>=minframe & fits.frame<=maxframe);
%         fitscc.col = fits.col(fits.frame>=minframe & fits.frame<=maxframe);
%         fitscc.goodfit = ones(nrlocsincc,1);
%         %perform avgshifthist without image gen to obtain partial images
%         ccim(:,:,:,cc_bins) = avgshifthist(fitscc,paramscc,movsz,[dlocs{ii},filesep,dnames{ii}],0);
%     end
%     
%     cc_points = [];
%     subpixloc = zeros(nrccbins-1,3);
%     totpix = zeros(nrccbins-1,2);
%     clear ypeak xpeak
%     clear stats cccim maxpix
%     %Perform cross-correlation for all bins compared to bin 1
%     for cc_bins = 1:nrccbins %for all the bins (frametime)
%         clear cc
%         cc = xcorr2(ccimwhite(:,:,1),ccimwhite(:,:,cc_bins));
%         cccim(:,:,cc_bins) = cc;
%         %find maximum pixel
%         [C,mp] = max(cc(:));
%         [maxpix(cc_bins,1),maxpix(cc_bins,2)] = ind2sub(size(cc),mp);
%         %Perform subpixel localization on 5x5-area around max pixel
%         [subpixloc(cc_bins,:)] = Phasor_localization_SMALLLABS( ...
%              cc(maxpix(cc_bins,1)-2:maxpix(cc_bins,1)+2, ...
%              maxpix(cc_bins,2)-2:maxpix(cc_bins,2)+2), ...
%              0,'');
%          %Store subpixel location
%          totpix(cc_bins,:) = maxpix(cc_bins,:)+subpixloc(cc_bins,1:2);
%     end
%     %Normalize subpixel location of cc-dift
%     totpix = totpix-totpix(1,:);
%     totpix = totpix./paramscc.avgshifthist_lateralsubpix;
%     %Interpolate with spline with required sampling interval
%     vqx = interp1([1:10],totpix(:,1),[1:0.1:10],'spline');
%     vqy = interp1([1:10],totpix(:,2),[1:0.1:10],'spline');
%     
%     %for now, plot drift over time
%     figure(2);clf(2);
%     hold on
%     plot([1:10],totpix(:,1),'kx')
%     plot([1:0.1:10],vqx,'k-')
%     plot([1:10],totpix(:,2),'r*')
%     plot([1:0.1:10],vqy,'r-')
%     grid on
%     