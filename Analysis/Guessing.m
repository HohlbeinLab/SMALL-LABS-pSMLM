function [guesses,roinum]=Guessing(mov_fname,mov,movsz,goodframe,dfrlmsz,...
    bpthrsh,egdesz,pctile_frame,debugmode,mask_fname,make_guessmovie,guessing_filter,allparams)
%% Guessing
% make a list of guesses for sinlge molecules. Using a bandpass filter to
% filter pixel noise first, then uses bwpropfilt to find blobs of the
% correct size
%
%%%% Inputs %%%%
% mov_fname is the full filename of the movie to be analyzed for output
% file naming purposes
%
% mov is the movie data as a 3D array where the third dimension is the
% frame number.
%
% movsz is the output of size(mov)
%
% goodframe is an optional logical vector indicating which frames are to be
% ignored. The length of the goodframe vector should be the number of
% frames in mov. To ignore a frame set the corresponding element in
% goodframe to false.
%
% dfrlmsz is the  size of a diffraction limited spot in pixels. It's the
% nominal diameter, NOT the FWHM or something similar. Integer please!
% KM: No longer integer necessary at input, it'll use the integer for ROI
% size selection, but the float value for checking filtered image for
% circle-equivalent sizes.
%
% bpthrsh is the the percentile of brightnesses of the bandpassed image
% below which those pixels will be ignored.
%
% edgesz is the number of pixels on the edge of the image that will be
% ignored.
%
% pctile_frame is a boolean determining whether bpthrsh will be applied
% frame by frame, or to the entire movie. Using the entire movie (setting
% to 0) is more sensitive to low frequency noise and background changes,
% but is a more robust guessing method. Using each frame tends to produce a
% constant number of guesses per frame, regardless of their absolute
% brightness.
%
% debugmode is a boolean to determine if you want to go through and look at
% the guesses.
%
% mask_fname is the filename of a mask to use for guessing. If no mask is
% being used just leave it empty. If mask_fname is set 1, then the program
% will look for a file in the same directory as the movie with '_PhaseMask'
% appened to the name of the movie. The mask is a .mat file which has a
% logical array (or at least where nonzero entries will be converted to 1s)
% called PhaseMask that is the same size as a frame in the current movie.
%
% make_guessmovie is a Boolean determining whether or not to make a .avi
% movie of the guesses. This can be helpful to determine how successful the
% guessing was.
%
%%%% Outputs %%%%
% guesses is an array with columns 1. frame #, 2. row #, 3. column # of the
% guesses
%
% The program currently also writes a .mat file with guesses and all of the
% user parameters saved.
%
%
%%%% Dependencies %%%%
% bpass
%
%     Copyright (C) 2017  Benjamin P Isaacoff
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
%
%
%did you not set dfrlmsz to an integer?
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end

%pad size for the bandpass function
pdsz=50;
global verbose
tic;%for measuring the time to run the entire program
% last updated 3/10/18 BPI
%% Setup

[pathstr,fname] = fileparts(mov_fname);
if verbose
    disp([char(datetime),'   Making guesses for ',fname])
end
%intializing the guess indices cell array
guesses=zeros(1,3);
roinum=0;
%% Guessing
%making the phasemask logical map
if ~isempty(mask_fname)
    if ischar(mask_fname)
        %the strrep is to get rid of the avgsub, note that this shouldn't
        %do anything if bgsub=0
        try
            load([pathstr,filesep,mask_fname,'.mat'],'PhaseMask')
        catch
            [datalist,dataloc,~]=uigetfile([pathstr,filesep,'*.*']);
            if ~iscell(datalist); datalist={datalist}; end
            datalist=[dataloc datalist];
            [dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);
            load([dlocs{1,1},filesep,dnames{1,2}],'PhaseMask')
        end
    else
        load([pathstr,filesep,strrep(fname,'_avgsub',[]),'_PhaseMask.mat'],'PhaseMask')
        
    end
    PhaseMasklg=PhaseMask;
    PhaseMasklg(PhaseMasklg~=0)=1;
    PhaseMasklg=logical(PhaseMasklg);
else
    PhaseMasklg=true(movsz([1,2]));
    PhaseMask=true(movsz([1,2]));
end
if debugmode
    figure();
end
if strcmp(guessing_filter,'bpass')
    %using the percentiles on the entire movie
    if ~pctile_frame
        %initializing the bandpassed movie
        bimgmov=zeros(movsz);
        goodfrmmov=false(movsz);
        %looping through and making the bandpassed movie
        if isempty(gcp('nocreate'))
            numWorkers = 0;
        else
            parpoolinfo = gcp;
            numWorkers = parpoolinfo.NumWorkers;
        end
        parfor(ll=1:movsz(3),numWorkers)
            if goodframe(ll)
                goodfrmmov(:,:,ll)=true;
                %padding the current frame to avoid the Fourier ringing
                %associated with the edges of the image
                curfrm=padarray(mov(:,:,ll),[pdsz,pdsz],'symmetric');

                %bandpass parameters
                LP=2;%lnoise, should always be 1
                HP=round(dfrlmsz*1.5);%lobject, set by diffraction limit
                T=0;%threshold, now always zero
                lzero=egdesz;%how many pixels around the edge should be ignored, optional
                %bandpass it
                bimg=bpass(curfrm,LP,HP,T,lzero+pdsz);
                %removed the padded pixels around the edge
                bimgmov(:,:,ll)=bimg((pdsz+1):(movsz(1)+pdsz),(pdsz+1):(movsz(2)+pdsz));
            end
        end

        %convert it to a logical movie by thresholding with the bpthrsh
        %percentile of the brightnesses for nonzero pixels
        bimgmov=logical(bimgmov.*(bimgmov>prctile(bimgmov(bimgmov>0 & ...
            goodfrmmov & repmat(PhaseMasklg,[1,1,movsz(3)])),bpthrsh)).*repmat(PhaseMasklg,[1,1,movsz(3)]));
    end

    if make_guessmovie
        figure();
        v = VideoWriter([pathstr,filesep,fname,'_Guesses.avi'],'Uncompressed AVI');
        open(v);

        disp(['Making guesses movie for ',fname]);
    end

    %If not parrallel pool
    if isempty(gcp('nocreate'))
        [roinum, guesses] = nonparrGuessingbpass(movsz,goodframe,pctile_frame,mov,pdsz,dfrlmsz,bimgmov,allparams,debugmode,make_guessmovie,guesses,roinum,PhaseMask,1,movsz(3)); 
    %If there is a parallel pool
    else
        parpoolinfo = gcp;
        numWorkers = parpoolinfo.NumWorkers;
        framestarts = round(linspace(1,movsz(3),numWorkers+1));
        for k = 1:numWorkers
            frameends(k) = framestarts(k+1)-1;
        end
        framestarts(end) = [];
        frameends(end) = movsz(3);

        parfor k = 1:numWorkers
            [roinump{k}, guessesp{k}] = nonparrGuessingbpass(movsz,goodframe,pctile_frame,mov,pdsz,dfrlmsz,bimgmov,allparams,debugmode,make_guessmovie,guesses,roinum,PhaseMask,framestarts(k),frameends(k)); 
        end

        totlistlength = 0;
        for k = 1:numWorkers
            totlistlength = totlistlength+size((guessesp{k}),1)-1;
        end
        guesses = zeros(totlistlength+1,3);
        roinum = zeros(totlistlength+1,1);
        counter = 2;
        for k = 1:numWorkers
            guesses(counter:counter+size(cell2mat(guessesp(k)),1)-1-1,:) = guessesp{k}(2:end,:);
            roinum(counter:counter+size(cell2mat(guessesp(k)),1)-1-1,:) = roinump{k}(2:end,:);
            counter = counter+size(cell2mat(guessesp(k)),1)-1;
        end
     %end ifstatement parallel pool
    end
elseif strcmp(guessing_filter,'wavelet')
    %Do wavelet filter on all frames - roughly 2x faster as per-frame
    [allfrmwvlt,V1,V2]=mywaveletfilteratrousNdim(single(mov)-mean(mean(single(mov))),1);
    
    if make_guessmovie
        figure();
        v = VideoWriter([pathstr,filesep,fname,'_Guesses.avi'],'Uncompressed AVI');
        open(v);

        disp(['Making guesses movie for ',fname]);
    end
    
    %If not parrallel pool
    if isempty(gcp('nocreate'))
        [roinum, guesses] = nonparrGuessingwvlt(allfrmwvlt,movsz,goodframe,pctile_frame,mov,pdsz,dfrlmsz,allparams,debugmode,make_guessmovie,guesses,roinum,PhaseMask,1,movsz(3)); 
    %If there is a parallel pool
    else
        parpoolinfo = gcp;
        numWorkers = parpoolinfo.NumWorkers;
        framestarts = round(linspace(1,movsz(3),numWorkers+1));
        for k = 1:numWorkers
            frameends(k) = framestarts(k+1)-1;
        end
        framestarts(end) = [];
        frameends(end) = movsz(3);
        parfor k = 1:numWorkers
            [roinump{k}, guessesp{k}] =nonparrGuessingwvlt(allfrmwvlt,movsz,goodframe,pctile_frame,mov,pdsz,dfrlmsz,allparams,debugmode,make_guessmovie,guesses,roinum,PhaseMask,framestarts(k),frameends(k)); 
        end

        totlistlength = 0;
        for k = 1:numWorkers
            totlistlength = totlistlength+size((guessesp{k}),1)-1;
        end
        guesses = zeros(totlistlength+1,3);
        roinum = zeros(totlistlength+1,1);
        counter = 2;
        for k = 1:numWorkers
            guesses(counter:counter+size(cell2mat(guessesp(k)),1)-1-1,:) = guessesp{k}(2:end,:);
            roinum(counter:counter+size(cell2mat(guessesp(k)),1)-1-1,:) = roinump{k}(2:end,:);
            counter = counter+size(cell2mat(guessesp(k)),1)-1;
        end
    %end ifstatement parallel pool
    end
end
    
guesses=guesses(2:end,:);%get rid of first row of zeros
roinum(1)=[];
if make_guessmovie
    close(v)
end

tictoc=toc;%the time to run the entire program

if allparams.saveloadmat
[pathstr,name,~] = fileparts([mov_fname '.mat']);
save([pathstr,filesep,name,'_guesses.mat'],'guesses','goodframe','dfrlmsz','egdesz','pctile_frame','bpthrsh',...
    'movsz','tictoc','mask_fname','roinum','-v7.3');
end

end

function [roinum, guesses] = nonparrGuessingbpass(movsz,goodframe,pctile_frame,mov,pdsz,dfrlmsz,bimgmov,allparams,debugmode,make_guessmovie,guesses,roinum,PhaseMask,beginframe,endframe)
    for ll=beginframe:endframe
        if goodframe(ll)
            %using the percentile on each frame
            if pctile_frame
                %padding the current frame to avoid the Fourier ringing
                %associated with the edges of the image
                curfrm=mov(:,:,ll);

                curfrmbp=padarray(curfrm,[pdsz,pdsz],'symmetric');

                %bandpass parameters
                LP=1;%lnoise, should always be 1
                HP=round(dfrlmsz*1.5);%lobject, set by diffraction limit
                T=0;%threshold, now always zero
                lzero=egdesz;%how many pixels around the edge should be ignored, optional
                %bandpass it
                bimg=bpass(curfrmbp,LP,HP,T,lzero+pdsz);
                %pull out the actual data
                bimg=bimg((pdsz+1):(movsz(1)+pdsz),(pdsz+1):(movsz(2)+pdsz));

                %threshold with the bpthrsh percentile of the brightnesses for
                %nonzero pixels, then turn it into a logical array
                logim=logical(bimg.*(bimg>prctile(bimg(bimg>0 & PhaseMasklg),bpthrsh)).*PhaseMasklg);
            else
                logim=bimgmov(:,:,ll);
            end

            %search for shapes with an EquivDiameter of floor(dfrlmsz/2) to
            %2*dfrlmsz
            bw2=bwpropfilt(logim,'EquivDiameter',[floor(allparams.dfrlmsz_float/2),2*allparams.dfrlmsz_float]);
            rgps=regionprops(bw2,'centroid');% find the centroids of those shapes
            centroids = cat(1, rgps.Centroid);%just rearraging the array

            %filling the array for this frame
            if ~isempty(centroids)
                guesses=cat(1,guesses,[repmat(ll,size(centroids(:,2))),round(centroids(:,2)),round(centroids(:,1))]);
                roinum=[roinum;diag(PhaseMask(round(centroids(:,2)),round(centroids(:,1))))];
            end

            if debugmode || make_guessmovie %plot the guesses, for checking parameters
                if ~pctile_frame
                    curfrm=mov(:,:,ll);
                end
                imshow(double(curfrm),prctile(double(curfrm(curfrm>0)),[.1,99.8]))
                if ~isempty(centroids)
                    %viscircles is reversed
                    vcs=viscircles([centroids(:,1),centroids(:,2)],repmat(dfrlmsz,[length(centroids(:,2)),1]));
                    set(vcs.Children,'LineWidth',1)
                end
                if debugmode
                    title([fname,'   frame ',num2str(ll)],'Interpreter','none')
                    drawnow
                elseif make_guessmovie
                    frame = getframe;
                    writeVideo(v,frame);
                end
            end
        end
    end
end






function [roinum, guesses] = nonparrGuessingwvlt(allfrmwvlt,movsz,goodframe,pctile_frame,mov,pdsz,dfrlmsz,allparams,debugmode,make_guessmovie,guesses,roinum,PhaseMask,beginframe,endframe)
    for ll=beginframe:endframe
        if goodframe(ll)
            %Get wavelet of this frame
            curfrmwvlt = allfrmwvlt(:,:,ll);
            %get std
            stdwvlt = std(reshape(curfrmwvlt,movsz(1)*movsz(2),1),0,1);
            %make boolean im where val is larger than Xx std?
            
            %Replaced original code with ~ 16% speed increase (KM)
            logimwvlt = logical(curfrmwvlt>=stdwvlt*allparams.wvltfilter_std);
            sizerangeeqvdiam = 1.5; %2 standard, 1.5 seems to work better?
            bw2=bwpropfilt(logimwvlt,'EquivDiameter',[floor(allparams.dfrlmsz_float/sizerangeeqvdiam),sizerangeeqvdiam*allparams.dfrlmsz_float]);
            %findpeaks method
            cent = FastPeakFind(single(bw2.*curfrmwvlt));
            centroids = [cent(1:2:end),cent(2:2:end)];
            
            %Original code here for keepsakes------------
%             logimwvlt = logical(curfrmwvlt>=stdwvlt*allparams.wvltfilter_std);
%             %search for shapes with an EquivDiameter of floor(dfrlmsz/X) to
%             %X*dfrlmsz
%             sizerangeeqvdiam = 1.5; %2 standard, 1.5 seems to work better?
%             bw2=bwpropfilt(logimwvlt,'EquivDiameter',[floor(allparams.dfrlmsz_float/sizerangeeqvdiam),sizerangeeqvdiam*allparams.dfrlmsz_float]);
%             rgps=regionprops(bw2,'centroid');% find the centroids of those shapes
%             centroids = cat(1, rgps.Centroid);%just rearraging the array
            %End original code-------------------
             
            %filling the array for this frame
            if ~isempty(centroids)
                guesses=cat(1,guesses,[repmat(ll,size(centroids(:,2))),round(centroids(:,2)),round(centroids(:,1))]);
                roinum=[roinum;diag(PhaseMask(round(centroids(:,2)),round(centroids(:,1))))];
            end
            
            
            if debugmode || make_guessmovie %plot the guesses, for checking parameters
                curfrm=mov(:,:,ll);
                imshow(double(curfrm),prctile(double(curfrm(curfrm>0)),[.1,99.8]))
                if ~isempty(centroids)
                    %viscircles is reversed
                    vcs=viscircles([centroids(:,1),centroids(:,2)],repmat(dfrlmsz,[length(centroids(:,2)),1]));
                    set(vcs.Children,'LineWidth',1)
                end
                if debugmode
                    title([fname,'   frame ',num2str(ll)],'Interpreter','none')
                    drawnow
                elseif make_guessmovie
                    frame = getframe;
                    writeVideo(v,frame);
                end
            end
        end
    end
end






