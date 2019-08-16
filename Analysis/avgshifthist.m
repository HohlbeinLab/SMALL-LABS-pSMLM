function [totimarr] = avgshifthist(fits,params,movsz,prefilename,imagegen)
    %Make 0 zpos if it doesn't exist
    if ~isfield(fits,'zpos')
        fits.zpos = zeros(size(fits.col));
    end
    histsubpix = params.avgshifthist_lateralsubpix; %number of subpixels corresponding to every pixel
    if sum(fits.zpos) == 0
        histzstacks = 1; %number of z-stacks for z-information
    else
        histzstacks = params.avgshifthist_axialcolnrs;
    end
    %Get min and max z pos, based on a x% pctle
    zstackprctlefilter = 1; % Percentage (in actual %s) to filter out of zpos for image generation, on each side
    prctlieoutput = prctile(fits.zpos,[zstackprctlefilter 100-zstackprctlefilter]);
    zlb = prctlieoutput(1);
    zub = prctlieoutput(2);
    zposedges = linspace(zlb,zub,(histzstacks+1));

    latshifts = params.avgshifthist_lateralshifts*2; %number of lateral shifts - should be even!
    axshifts = 0; %number of axial shifts - to be added later?
    edgeval = 0;%params.egdesz;
    avgshhist = zeros((movsz(1)-edgeval*2)*histsubpix, (movsz(2)-edgeval*2)*histsubpix, histzstacks);
    %loop over all localizations
    try
        for locs = 1:size(fits.row,1)
            %         fprintf('loc %.0f of %.0f \r\n',locs,size(fits.row,1));
            if fits.goodfit(locs)
                %find correct xbin, ybin, zbin
                xbin = max(min(round((fits.row(locs,1)-edgeval)*histsubpix,0),(movsz(1)-edgeval*2)*histsubpix),1);
                ybin = max(min(round((fits.col(locs,1)-edgeval)*histsubpix,0),(movsz(2)-edgeval*2)*histsubpix),1);
                zbin = min(max(min(find(zposedges > fits.zpos(locs,1)))-1,1),histzstacks);
                if isempty(zbin)
                    zbin = 1;
                end
                %add value in bin itself
%                 avgshhist(xbin,ybin,zbin) = avgshhist(xbin,ybin,zbin)+1;
                %add in bins around
                for xxbin = max(xbin-latshifts/2,1):min(xbin+latshifts/2,(movsz(1)-edgeval*2)*histsubpix)
                    for yybin = max(ybin-latshifts/2,1):min(ybin+latshifts/2,(movsz(2)-edgeval*2)*histsubpix)
                        
                        avgshhist(xxbin,yybin,zbin) = avgshhist(xxbin,yybin,zbin)+latshifts-(sqrt((xxbin-xbin)^2)+sqrt((yybin-ybin)^2));
%                         tic
%                         for i = 1:10000
%                             x=latshifts-(abs(xxbin-xbin)+abs(yybin-ybin));
%                         end
%                         toc
%                         tic
%                         for i = 1:10000
%                             x=latshifts-(sqrt((xxbin-xbin)^2)+sqrt((yybin-ybin)^2));
%                         end
%                         toc
%                         tic
%                         for i = 1:10000
%                         v=latshifts-(pdist2(xxbin,xbin,'cityblock')+pdist2(yybin,ybin,'cityblock'));
%                         end
%                         toc
%                         fprintf('%.0f - %.0f\n',x,v)
                    end
                end
            end
        end
        %Clean up edges
        avgshhist(1:latshifts,:,:) = 0;
        avgshhist(:,1:latshifts,:) = 0;
        avgshhist(end-latshifts:end,:,:) = 0;
        avgshhist(:,end-latshifts:end,:) = 0;

        %If wanted, make single z-slices
        %     figure(2);clf(2);
        %     for zst = 1:histzstacks
        %         subplot(ceil(sqrt(histzstacks)),ceil(sqrt(histzstacks)),zst)
        %         imagesc(avgshhist(:,:,zst));
        %         axis off
        %     end

        clear coloverz imarr ax gca
        figure(1);clf(1);
        for zst = 1:histzstacks
            im = imagesc(avgshhist(:,:,zst));
            imarr(:,:,zst) = im.CData;
            if sum(sum(im.CData)) >0
                imarr(:,:,zst) = imarr(:,:,zst)./max(max(imarr(:,:,zst)));
            end
        end
        coloverz = colormap(jet(histzstacks));

        totimarr = zeros(size(imarr,1),size(imarr,2),3);
        %Loop over all z bins
        for zz = 1:size(imarr,3)
            t = zeros(64,64,3);
            %make temp col array
            for xx = 1:64
                for yy = 1:64
                    if size(imarr,3) > 1
                        t(xx,yy,:) = coloverz(zz,:);
                    else
                        t(xx,yy,:) = [1 1 1];
                    end
                end
            end
            %convolute col array with image array and sum
            totimarr(:,:,:) = totimarr(:,:,:)+t(xx,yy,:).*imarr(:,:,zz);
        end
        %Certain % of pixels may be saturated in at least 1 color
        saturationfraction = 0.025; %2.5%
        %reshape image
        temparr = reshape(totimarr,3*size(totimarr,1)*size(totimarr,2),1);
        %remove 0-values
        temparr(temparr==0) = [];
        %get value at satfraction
        satvalue = prctile(temparr,100-saturationfraction*100);
        totimarr = totimarr./satvalue;

        %% Make image
        clf(1);
        % plot image
        hold on
        imagesc('CData',totimarr);
        % axis stuff
        axis tight;
        axis square;
        axis off;
        ax = gca;
        ax.YDir = 'reverse';
        % pixel size, create some kind of scalebar
        fullimsize = movsz(2)*params.avgshifthist_pixelsize/1000;
        optimalscalebarsize = fullimsize/10;
        %round scalebar to nearest 1 µm
        scalebarsize = round(optimalscalebarsize);
        %draw scalebar in bottom right (high x, high y), with some spacing
        spacingx = 0.17*histsubpix; %in % of image
        spacingy = spacingx/1.2; %in % of image
        sbsizex = (scalebarsize/fullimsize)*movsz(2);%size scalebar in x
        sbsizey = 0.10*sbsizex;
        fill([movsz(2)-spacingx movsz(2)-spacingx-sbsizex movsz(2)-spacingx-sbsizex movsz(2)-spacingx]*histsubpix,...
             [movsz(1)-spacingy movsz(1)-spacingy movsz(1)-spacingy-sbsizey movsz(1)-spacingy-sbsizey]*histsubpix,'w')
        figpos = get(gcf,'Position');
        axpos = get(gca,'Position');
        figwidth = figpos(3)*axpos(3);
        text(((movsz(2)-spacingx-sbsizex*.5)*histsubpix),(movsz(1)-spacingy-sbsizey-0.2)*histsubpix,[num2str(scalebarsize) ' µm'], ...
        'Color','white','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figwidth/50)
        %Zpos colorbar
        if histzstacks > 1
            cbar = colorbar();
            cbar.Ticks = [0 1];
            cbar.TickLabels = {[num2str(round(zlb*1000,0)) ' nm'],[num2str(round(zub*1000,0)) ' nm']};
        end
        %% Export image
        if imagegen
        resstr = ['-r' num2str(params.avgshifthist_res)];
        export_fig([prefilename,'AvgShiftedHist.tif'], resstr,'-transparent');
        end
    catch
        fprintf('Some error occured in avgshifthist making\n')
    end
end