function SMALLLABS_main(file_or_directory_name,dfrlmsz,avgwin,moloffwin,varargin)
%% SMALLLABS_main
%
%%%% If you obtained this code from anywhere other than the Biteen Lab
%%%% GitHub please visit https://github.com/BiteenMatlab/SMALL-LABS to
%%%% obtain the most up-to-date version of this code.
%
% SMALLLABS_main is the wrapper function for the SMALL-LABS algorithm which
% does accurate background substraction of single molecule imaging movies
% so that the molecules can be fit and their intensity accurately measured.
%
% Note: that by setting bgsub=false SMALLLABS_main will do fitting without doing
% background subtraction. Doing this obviates a lot of parameters, it
% doesn't matter what they're set to. See the User Guide for more details.
%
%
%%%% Inputs %%%%
%%% required
%   file_or_directory_name   is the name of the directory where the movies
%   will be selected OR the filename of a single movie OR a list of movies
%   to be fit. See the User Guide for more details and/or Movie2mat for
%   details about the file formats currently supported.
%
%   dfrlmsz   is the size of a diffraction limited spot in pixels. It's the
%   nominal diameter, NOT the FWHM or something similar. Must be an
%   integer! For an expected diffraction limited standard deviation, std,
%   using the full width at 20% max, dfrlmsz = std*(2*sqrt(2*log(5)))
%
%   avgwin   is the length of the temporal window (in frames) to be used
%   for the average subtraction.
%
%   moloffwin   is the length of the temporal window (in frames) to be
%   checked to determine in which frames that molecule was off and are thus
%   safe to subtract. Needs to be an even integer.
%
%%% optional (varargin)
% See the list of optional parameters in parameters section below, details
% can be found in the User Guide. They are called with a name value pair as
% is standard for Matlab functions. e.g., 'bpthrsh',83 would set the
% parameter bpthrsh=83
%
%
%%%% Outputs %%%%
% This function will output a series of movies & .mat files depending on
% which functions are called. See each function or the User Guide for
% details
%
%%%% Dependencies %%%%
% AVGSUB_moves
% TIFFStack 
% Guessing
% bpass
% Mol_off_frames
% MLEwG
% gaussFit
% Subtract_then_fit
% Track_3D2
% Tracking
% Track_filter
% hungarian
% ViewFits
% ViewFitsTracking
% gpufit
% Movie2mat
%
%%%%
% Written by Benjamin P Isaacoff at the University of Michigan last update
% 3/10/18 BPI & SAL
%
% Updated by Koen Martens at the University of Wageningen in 2019 (adding
% phasor-based localization, GUI, wavelet filter, avg shifted hist,
% drift-corection)
%
%     Copyright (C) 2018  Benjamin P Isaacoff
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

global verbose
verbose=true; %Set verbose to true to have matlab print the current step and date/time to the command window. Otherwise set to false
%% Parameter Defaults
% You are of course welcome to change the default values, but I would
% strongly urge you to instead set them as inputs to the function using a
% name value pair (e.g., 'bpthrsh',83). Please see the User Guide for
% details about the meaning and function of these parameters. The default
% parameter values are:

%Save and load intermediate steps as mat
params.saveloadmat = true;

%Check if GUI is used. Normally 0, calling SMALLLABS from GUI sets it to 1
params.useGUI = 0;

%%% Actions %%%
% Do background subtraction
params.bgsub = true;
% Make the average subtracted movie
params.makeAvgsub = true;
% Do guessing
params.guessing = true;
% Check the guesses
params.checkGuesses = true;
% Make the off frames list
params.makeOffFrames = true;
% Do fitting
params.fitting = true;
% Do tracking
params.tracking = true;
% Make the ViewFits movie
params.makeViewFits = true;

%%% AVGSUB parameters %%%
% Subtract the temporal average? otherwise use median
params.do_avg = false;
% AVGSUB offset
params.offset = 5000; %default 1000

%%% Guessing parameters %%%
% the threshold percentile of the banpassed movie
params.bpthrsh = 97; %standard 95
% how many pixels to ignore around the edge of the frame
params.egdesz = dfrlmsz;
% compare brightnesses in each frame? otherwise use entire movie
params.pctile_frame = true;
% use a mask for guessing? Enter it's filename or set to true to look
% for .mat file w/ '_PhaseMask' appended to movie name
params.mask_fname=[];
% make a movie (.avi) of the guesses to check parameters
params.make_guessmovie = false;

%%% Fitting parameters %%%
% do MLE fitting? If not least squares fitting will be used
params.MLE_fit = true;
% Goodfit parameters. See Subtract_mol_off_frames for the details
params.stdtol = 1.5;
params.maxerr = 0.10; %if you change this please also change the if statement after the next loop
% subtract the mean of off frames? If not use median
params.do_avgsub = false;
% Which Gaussian function to fit to when using LSQ fitting? 1. symmetric,
% 2. fixed angle asymmetric, 3. free angle asymmetric
params.which_gaussian = 1;
% the angle to fit to for a fixed angle Gaussian
params.fit_ang= 0;
% use a GPU if one is a available?
params.usegpu= true;


%Make ThunderSTORM csv output file?
params.makeThSTORMoutput = true;

%% Phasor fitting addition
% Koen Martens, 08-04-2019, based on extension of pSMLM-3D
params.phasor_fit = true; %true for pSMLM fitting. Overrides Guassian fitting
params.phasor_rad = 3;
params.phasor_2d = true;
params.phasor_DH = false; %true if DH phasor is wanted. Requires correct calibration file
params.phasor_DHcal = ''; %Full path to DH calibration file here. Can be made via 
params.phasor_astig = false; %true if astig phasor is wanted. Requires correct calibration file
params.phasor_astigcal = ''; %Full path to astig calibration file here. Can be made via 
params.phasor_SP = false; %true if SP phasor is wanted. Requires correct calibration file
params.phasor_SPcal = ''; %Full path to SP calibration file here. Can be made via 
params.phasor_TP = false; %true if TP phasor is wanted. Requires correct calibration file
params.phasor_TPcal = ''; %Full path to TP calibration file here. Can be made via 
params.phasor_SPTPradius = 10; %'Large' radius for SP/TP phasor.

%For now, here is some quick and dirty calibration code, this should be
%done seperately later.
%     zposcali = ([1:151]*0.01-0.76); %in µm
%     %Input (now from csv from ThStorm, should be from SMALLLABS later)
%     input = csvread('H:\Data\DoubleHelix\sequence-as-stack-Beads-DH-Exp-as-stack\AllBeads_phasor2D_phasor3.csv',1,0);
%     %Required input is frame-x-y
%     fxyinput = input(:,[2:4 9 10]);
%     %Set x-y to pixels rather than nm
%     fxyinput(:,2:3) = fxyinput(:,2:3)./100;
%     %Perform calibration
%     pSMLM_DH_calibration_SL(fxyinput,zposcali,1,'H:\Data\DoubleHelix\sequence-as-stack-Beads-DH-Exp-as-stack\SLpSMLM_DH_calib_phasor3.mat');

%%% Filter parameters %%%
params.filtermethod = 'wavelet'; %Choose between 'wavelet' and 'bpass'
params.wvltfilter_std = 2; %multiplicator of wavlet filter std.    

%%% Partial analysis %%%
params.fullmovie = 1;
params.nrframes = 500;
params.specificframeanalysis = 0;
params.specificframe = 100;
    
%% Resume original code
%%% Tracking parameters %%%
% default is [0.01,200,0.5,3,3,1,0]
% save the separate tracks .mat file?
params.savetracks = true;
% minimum merit
params.trackparams(1)=0.01;
% Integration time (ms)
params.trackparams(2)=200;
% gamma
params.trackparams(3)=0.5;
% maximum step size
params.trackparams(4)=3;
% minimum track length
params.trackparams(5)=3;
% speed estimation window halfsize
params.trackparams(6)=1;
% time delay between consecutive frames (ms)
params.trackparams(7)=0;
%%% ViewFits parameters %%%
% use the original movie? if not, use the avgsub movie
params.orig_movie = true;
% diameter of the circles showing the fits
params.circ_D = dfrlmsz;
% linewidth of the circles
params.linewidth = 1;%1
% write a .avi movie showing the fits. If not, goes to debug mode
params.write_mov=true;
% autoscale frame by frame?
params.autoscale_on = false;
% use the tracking viewfits instead
params.trackingVF = false;
%% Average shifted histogram

%Make average shifted histogram?
params.makeavgshifthist = true;
params.avgshifthist_res = 300; %DPI of avgshifthist
params.avgshifthist_lateralsubpix = 15; %axial subpixels per pixel
params.avgshifthist_lateralshifts = 3; %axial subpixels per pixel
params.avgshifthist_axialcolnrs = 9; %lateral colors per pixel
params.avgshifthist_pixelsize = 100; %pixel size in nm

%% Drift correction via cross-correlation
params.driftcorr_cc = true; %performing driftcorr_cc?
params.driftcorrcc_lateralsubpix = 10; %nr of subpixels for xy
params.dirftcorrcc_temporalbins = 20; %nr of temporal bins for cc

%% Evaluating the inputs

paramsnames=fieldnames(params);
% if any parameters are included as inputs, change the parameter mentioned
if nargin>4
    for ii=1:2:nargin-5
        whichField = strcmp(paramsnames,varargin{ii});    
        if strcmp('handles',varargin{ii})
            GUIhandles = varargin{ii+1};
            axis(GUIhandles.axes1);
            plot([1,2],[1,2]);
            axis off
        else
            try
                eval(['params.' paramsnames{whichField} ' = varargin{ii+1};'])
            catch
                error([varargin{ii}, '  is not an input parameter. Check the spelling.'])
            end
        end
    end
end
%changing the default error if doing MLE fitting
if params.MLE_fit && params.maxerr==3
    params.maxerr=0.1;
end

% verify that a GPU is available to use if usegpu==true
if params.usegpu
    try
        params.usegpu=parallel.gpu.GPUDevice.isAvailable;
    catch
         params.usegpu=false;
    end
end
%% Checking the required inputs
% Erroring out if dfrlmsz isn't an integer, because it's a "strong"
% parameter. You should really input what you mean.
params.dfrlmsz_float = dfrlmsz;
dfrlmsz = round(dfrlmsz);
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end

%Rounding moloffwin and reseting it to proper values. Not
%erroring because it's a "weak" parameter. Only needed if doing bgsub
if params.bgsub && params.makeOffFrames
    if moloffwin~=(ceil(moloffwin/2)*2)
        moloffwin=(ceil(moloffwin/2)*2);
        warning(['moloffwin must be an even integer. It has been reset to avgwin = ',num2str(moloffwin)])
    end
end

%% Select the movies
% check if it's a directory or a file
if exist(file_or_directory_name)==7
    %Select movies with uigetfile. If you make an error in specifying the
    %directory, it opens in the current directory.
    disp('Select the movie(s)')
    try
        [datalist,dataloc,findex]=uigetfile([file_or_directory_name,filesep,'*.*'],'multiselect','on');
    catch
        curdir=pwd;
        [datalist,dataloc,findex]=uigetfile([curdir,filesep,'*.*'],'multiselect','on');
    end
    if findex==0
        error('No movies selected')
    end
    %convert to a list of directories and filenames
    if ~iscell(datalist); datalist={datalist}; end
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
    [dlocs,dnames,exts]=cellfun(@fileparts,datalist,'uniformoutput',false);
elseif exist(file_or_directory_name)==2
    %get the directory and filename, and format into the cell list as above
    [dname,fname,ext] = fileparts(file_or_directory_name);
    %if it's a .txt file then assume it's a list of filenames
    if strcmp(ext,'.txt')
        %open the file for reading
        fid=fopen(file_or_directory_name,'r');
        %initializing loop variables
        linetxt='foobar';
        ii=0;
        while all(linetxt~=-1)
            ii=ii+1;
            linetxt=fgetl(fid);
            if linetxt~=-1
                [dname,fname,ext]=fileparts(linetxt);
                dlocs{ii}=dname;
                dnames{ii}=fname;
                exts{ii}=ext;
            end
        end
        fclose(fid);%close the file
        %remove any leading or trailing blank spaces
        dlocs=strtrim(dlocs);
        exts=strtrim(exts);
    else
        dlocs{1}=dname;
        dnames{1}=fname;
        exts{1}=ext;
    end
    clear dname fname ext
else
    error('Please input either a directory name or a filename.')
end

%% Loop through all the movies

% turn off the warning if a variable isn't found in a .mat file
warning('off','MATLAB:load:variableNotFound');

% try to add gpufit to path
try
    addpath(genpath('gpufit'))
end

if params.useGUI
    drawGUIloadbar(0.01,'Starting',GUIhandles)
else
h2=waitbar(0);
set(findall(h2,'type','text'),'Interpreter','none');
waitbar(0,h2,'Starting SMALLLABS');
end
wholeshabang=tic;
for ii=1:numel(dlocs)
    clear goodframe
    if params.useGUI
        drawGUIloadbar((7*ii-6)/numel(dlocs)/7,{['Initializing']},GUIhandles) 
    else
        try; waitbar((7*ii-6)/numel(dlocs)/7,h2); end
    end
    %% Convert the movies to .mat
    % Movies must be version 7.3 mat files with the movie data saved in the
    % 'mov' variable. This function converts several standard scientific movie
    % data types into this format.
    if params.fullmovie
        mov = Movie2mat([dlocs{ii},filesep,dnames{ii},exts{ii}],params);
    elseif params.specificframeanalysis %only specific frame - keep in mind the frames required for temporal subtraction, should also be loaded
        if params.bgsub
            startframe = params.specificframe-max(avgwin,moloffwin);
            endframe = params.specificframe+max(avgwin,moloffwin);
            if startframe < 0
                error('Impossible! Background subtraction is bigger than frame, use a different frame for preview!')
            end
        else
            startframe = params.specificframe;
            endframe = params.specificframe+1;
        end
        mov = Movie2mat([dlocs{ii},filesep,dnames{ii},exts{ii}],params,[startframe,endframe]);
    else %only first X frames
        mov = Movie2mat([dlocs{ii},filesep,dnames{ii},exts{ii}],params,[1,params.nrframes]);
    end
%     if params.fullmovie == 0
%         mov = mov(:,:,1:params.nrframes);
%     end
    if params.saveloadmat
        load([dlocs{ii},filesep,dnames{ii},'.mat'],'mov');
    end
    mov=single(mov);
    movsz=size(mov);
    try
        warning('off','MATLAB:load:variableNotFound');
        load([dlocs{ii},filesep,dnames{ii},'.mat'],'goodframe');
        warning('on','MATLAB:load:variableNotFound');
    catch
    end
    if ~exist('goodframe','var')
        goodframe=true(movsz(3),1);
    end
    %% The Average Subtraction
    % Only if doing bgsub.     
    % Because this step can be somewhat slow, first there is a check to
    % determine if the avgsub movie has already been created. If so, then
    % there's a further check to determine whether or not any of the
    % parameters have changed. If they're identical, this step is skipped.
    % To bypass this process, simply delete the avgsub movie from the
    % directory.
    
    % AVGSUB_movs will save an average subtracted movie .mat file, called
    % moviename_avgsub.mat 
    if params.makeAvgsub && params.bgsub
        if exist([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'file')~=0
            avgsubparams=load([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'subwidth','offset','do_avg','goodframe');
            % KM: Prevented skipping avgsub due to partial analysis
            % shenanigans
            %             if avgsubparams.subwidth==avgwin && avgsubparams.offset==params.offset && avgsubparams.do_avg==params.do_avg && sum(avgsubparams.goodframe)==sum(goodframe)
%                 disp(sprintf(['Avgsub parameters unchanged from previous calculation.\r Skipping avgsub on ',dnames{ii}],1))
%             else
                if params.useGUI
                    drawGUIloadbar((7*ii-5)/numel(dlocs)/7,{['Running AVGSUB on ',dnames{ii}]},GUIhandles) 
                else
                    try; waitbar((7*ii-5)/numel(dlocs)/7,h2,{['Running AVGSUB on ',dnames{ii}],'Overall Progress'}); end
                end
                bgsub_mov=AVGSUB_movs([dlocs{ii},filesep,dnames{ii}],mov,goodframe,...
                    params.do_avg,avgwin,params.offset,params.saveloadmat);
%                 AVGSUB_movs([dlocs{ii},filesep,dnames{ii}],mov,goodframe,...
%                     params.do_avg,avgwin,params.offset,params.saveloadmat);
%             end
        else
            %do the avgsub
            if params.useGUI
                drawGUIloadbar((7*ii-5)/numel(dlocs)/7,{['Running AVGSUB on ',dnames{ii}]},GUIhandles) 
            else
                try; waitbar((7*ii-5)/numel(dlocs)/7,h2,{['Running AVGSUB on ',dnames{ii}],'Overall Progress'}); end
            end
            bgsub_mov=AVGSUB_movs([dlocs{ii},filesep,dnames{ii}],mov,goodframe,...
                params.do_avg,avgwin,params.offset,params.saveloadmat);
        end        
    end
    %% Make Guesses
    % loop through each movie and make the guesses, which will be saved in a
    % .mat file called moviename_guesses.mat
    if params.guessing
        if params.useGUI
            drawGUIloadbar((7*ii-4)/numel(dlocs)/7,{['Making guesses for ',dnames{ii}]},GUIhandles) 
        else
            try; waitbar((7*ii-4)/numel(dlocs)/7,h2,{['Making guesses for ',dnames{ii}],'Overall Progress'}); end
        end
        if params.bgsub
            % try loading in the bgsub movie
            if params.saveloadmat
                try
                    bgsub_mov=load([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'mov');
                    bgsub_mov=bgsub_mov.mov;
                catch
                    warning('No avgsub file was found, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                    keyboard
                    params.bgsub=0;
                end
            end
            if params.bgsub
                [guesses,roinum] = Guessing([dlocs{ii},filesep,dnames{ii},'_avgsub'],bgsub_mov,movsz,goodframe,dfrlmsz,...
                    params.bpthrsh,params.egdesz,params.pctile_frame,params.checkGuesses,...
                    params.mask_fname,params.make_guessmovie,params.filtermethod,params);
            else
                [guesses,roinum] = Guessing([dlocs{ii},filesep,dnames{ii}],mov,movsz,goodframe,dfrlmsz,...
                    params.bpthrsh,params.egdesz,params.pctile_frame,params.checkGuesses,...
                    params.mask_fname,params.make_guessmovie,params.filtermethod,params);
            end
        else
            [guesses,roinum] = Guessing([dlocs{ii},filesep,dnames{ii}],mov,movsz,goodframe,dfrlmsz,...
                params.bpthrsh,params.egdesz,params.pctile_frame,params.checkGuesses,...
                params.mask_fname,params.make_guessmovie,params.filtermethod,params);
        end
        if params.saveloadmat
            clear bgsub_mov guesses
        end
    end   
    %% Make the off frames list
    % Only if doing bgsub. 
    % Loop through all of the movies and using the guesses .mat file, will write
    % the off frames list to a .mat file, called guessesname_Mol_off_frames.mat
    if params.makeOffFrames && params.bgsub
        if params.useGUI
            drawGUIloadbar((7*ii-3)/numel(dlocs)/7,{['Making off-frames lists for ',dnames{ii}]},GUIhandles) 
        else
            try; waitbar((7*ii-3)/numel(dlocs)/7,h2,{['Making off-frames lists for ',dnames{ii}],'Overall Progress'}); end
        end
            %try loading in the guesses
        if params.saveloadmat
            try
                load([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],'guesses','dfrlmsz')
            catch
                warning('No avgsub_guesses file was found, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                keyboard
                params.bgsub=0;
            end
        end
        if params.bgsub
            off_frames=Mol_off_frames([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],...
                guesses,goodframe,movsz,dfrlmsz,moloffwin,params.saveloadmat);
        end
        if params.saveloadmat
            clear guesses
        end
    end
    %% Subtract and fit
    % If not doing bgsub then a string ('nobgsub') is sent to
    % Subtract_then_fit which will proceed accordingly.     
    % Loop through the guesses and fit the subtracted images using the off
    % frames list. If doing bgsub outputs a .mat file, called
    % moviename_AccBGSUB_fits.mat, otherwise if not doing bgsub, it's
    % called moviename_fits.mat.
    if params.fitting
        if params.useGUI
            drawGUIloadbar((7*ii-2)/numel(dlocs)/7,{['Fitting ',dnames{ii}]},GUIhandles) 
        else
            try; waitbar((7*ii-2)/numel(dlocs)/7,h2,{['Fitting ',dnames{ii}],'Overall Progress'}); end
        end
        if params.bgsub
            %try loading in the mol_off_frames
            if params.saveloadmat
                try
                    load([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses_Mol_off_frames.mat'],'off_frames','dfrlmsz','moloffwin')
                catch
                    warning('No Mol_off_frames file was found, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                    keyboard
                    params.bgsub=0;
                end            
            end
            if params.bgsub
                % load in the guesses
                if params.saveloadmat
                    load([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],'guesses','roinum');
                end
                % fit it
                fits = Subtract_then_fit([dlocs{ii},filesep,dnames{ii}],mov,movsz,...
                    off_frames,moloffwin,guesses,roinum,dfrlmsz,params.MLE_fit,params.stdtol,...
                    params.maxerr,params.do_avgsub,params.which_gaussian,params.fit_ang,params.usegpu,...
                    params.phasor_fit,params.phasor_rad,params);
            else
                % load in the guesses
                if params.saveloadmat
                    load([dlocs{ii},filesep,dnames{ii},'_guesses.mat'],'guesses','rounum','dfrlmsz')
                end
                % fit it
                fits = Subtract_then_fit([dlocs{ii},filesep,dnames{ii}],mov,movsz,...
                    'nobgsub','nobgsub',guesses,roinum,dfrlmsz,params.MLE_fit,params.stdtol,...
                    params.maxerr,params.do_avgsub,params.which_gaussian,params.fit_ang,params.usegpu,...
                    params.phasor_fit,params.phasor_rad,params);
            end
        else
            % load in the guesses
            if params.saveloadmat
                load([dlocs{ii},filesep,dnames{ii},'_guesses.mat'],'guesses','roinum','dfrlmsz')
            end
            % fit it
            fits = Subtract_then_fit([dlocs{ii},filesep,dnames{ii}],mov,movsz,...
                'nobgsub','nobgsub',guesses,roinum,dfrlmsz,params.MLE_fit,params.stdtol,...
                params.maxerr,params.do_avgsub,params.which_gaussian,params.fit_ang,params.usegpu,...
                    params.phasor_fit,params.phasor_rad,params);
        end
        if params.saveloadmat
            clear guesses off_frames
        end
    end
    %% Tracking
    % Does tracking and appends the tracks to the fits .mat file. Also will
    % append a logical vector called trk_filt which indicates if the fit
    % passed was successfully tracked and wasn't the first or last frame in
    % a track.
    if params.tracking
        if params.useGUI
            drawGUIloadbar((7*ii-1)/numel(dlocs)/7,{['Tracking ',dnames{ii}]},GUIhandles) 
        else
            try; waitbar((7*ii-1)/numel(dlocs)/7,h2,{['Tracking ',dnames{ii}],'Overall Progress'}); end
        end
        if params.bgsub
            % try loading in the fits
            if params.saveloadmat
                try
                    load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits')
                catch
                    warning('Fitting file was without AccBGSUB, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                    keyboard
                    params.bgsub=0;
                    load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits')
                end
            end
            % track it
            if params.bgsub
                Track_filter([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],fits,...
                    params.trackparams,params.savetracks);
            else
                Track_filter([dlocs{ii},filesep,dnames{ii},'_fits.mat'],fits,...
                    params.trackparams,params.savetracks);
            end
        else            
            if params.saveloadmat
                load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits')
            end
            Track_filter([dlocs{ii},filesep,dnames{ii},'_fits.mat'],fits,...
                params.trackparams,params.savetracks);
        end
        if params.saveloadmat
            clear fits
        end
    end
    %% Preview option
    if params.specificframeanalysis
        pause(0.5);
        %load fits
        if params.saveloadmat
            if params.bgsub
                %try loading in the fits & tracking results
                load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits');
            else
                load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits');
            end
        end
        if params.bgsub
            if params.specificframe > (max(avgwin,moloffwin)+1)
                chosenframe = (max(avgwin,moloffwin)+1);
            else
                chosenframe = params.specificframe;
            end
        else
            chosenframe = 1;
        end
%         fitsonframe = [fits.row(fits.frame==(max(avgwin,moloffwin)+1)) fits.col(fits.frame==(max(avgwin,moloffwin)+1))];
        fitsonframe = [fits.row(fits.frame==chosenframe) fits.col(fits.frame==chosenframe)];
        figure(10);clf(10);
%         imagesc(mov(:,:,max(avgwin,moloffwin)+1))
        imagesc(mov(:,:,chosenframe));
        hold on
        if ~isempty(fitsonframe)
            scatter(fitsonframe(:,2),fitsonframe(:,1),130,'ro','filled')
        end
        set(gca,'YDir','reverse')
%         axis off
        colorbar
        colormap parula
%         caxis([prctile(reshape(mov(:,:,max(avgwin,moloffwin)+1),[movsz(1)*movsz(2),1]),10),...
%             max(reshape(mov(:,:,max(avgwin,moloffwin)+1),[movsz(1)*movsz(2),1]))*.8]);
    end
    %% Drift corr
    if (params.driftcorr_cc && params.specificframeanalysis == 0)
        if verbose
            disp([char(datetime),'   Starting cross-correlation drift correction'])
        end
        %load fits
        if params.saveloadmat
            if params.bgsub
                %try loading in the fits & tracking results
                load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits');
            else
                load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits');
            end
        end
        totframes = max((fits.frame));
        nrccbins = params.dirftcorrcc_temporalbins;
        % XY drift
        coords = [fits.row, fits.col, fits.frame];
        [coordscorr, finaldrift, A,b] = RCC(coords, max(fits.frame)/nrccbins, movsz(1), params.avgshifthist_pixelsize, params.avgshifthist_pixelsize/params.driftcorrcc_lateralsubpix, 0.2);
        fits.row = coordscorr(:,1);
        fits.col = coordscorr(:,2);
        if isfield(fits,'zpos')
            % Zdrift
            paramscc = params;
            %Number of bins used for cross-correlation - should be user-entered
            clear ccim meanframe
            % Create an avg shift image consisting of all positions within the
            % specified frames (specified by amount of cc-bins)
            for cc_bins = 1:nrccbins %for all the bins (frametime)
                %Calculate min and max frame
                minframe = floor(totframes/nrccbins)*(cc_bins-1)+1;
                maxframe = floor(totframes/nrccbins)*(cc_bins);
                meanframe(cc_bins) = minframe*.5+maxframe*.5;
                %Find number of localizations
                nrlocsincc = size(fits.frame(fits.frame>=minframe & fits.frame<=maxframe),1);
                %Make structures with x,y,goodfit params
                fitscc = {};
                fitscc.row = fits.row(fits.frame>=minframe & fits.frame<=maxframe);
                fitscc.col = (fits.zpos(fits.frame>=minframe & fits.frame<=maxframe)-min(fits.zpos))*(movsz(2)/(max(fits.zpos)-min(fits.zpos)));
                fitscc.goodfit = ones(nrlocsincc,1);
%                 paramscc.driftcorrcc_lateralsubpix = paramscc.avgshifthist_lateralsubpix;
                paramscc.avgshifthist_lateralsubpix = paramscc.driftcorrcc_lateralsubpix;
                %perform avgshifthist without image gen to obtain partial images
                ccim(:,:,:,cc_bins) = avgshifthist(fitscc,paramscc,movsz,[dlocs{ii},filesep,dnames{ii}],0);
            end
            
            ccimwhite = zeros(size(ccim,1),size(ccim,2),size(ccim,4));
            for i = 1:size(ccim,4)
                for j = 1:size(ccim,3)
                    ccimwhite(:,:,i) = ccimwhite(:,:,i) + ccim(:,:,j,i);
                end
            end
            %
            %     figure(6);clf(6);
            %     for i = 1:9
            %         subplot(3,3,i)
            %         imagesc(ccimwhite(:,:,i))
            %         axis square
            %         axis off
            %     end
            %
            subpixloc = zeros(nrccbins-1,3);
            totpix = zeros(nrccbins-1,2);
            clear ypeak xpeak
            clear stats cccim maxpix
            %Perform cross-correlation for all bins compared to bin 1
            for cc_bins = 1:nrccbins %for all the bins (frametime)
                clear cc
                cc = xcorr2(ccimwhite(:,:,1),ccimwhite(:,:,cc_bins));
                cccim(:,:,cc_bins) = cc;
                %find maximum pixel
                [C,mp] = max(cc(:));
                [maxpix(cc_bins,1),maxpix(cc_bins,2)] = ind2sub(size(cc),mp);
                %Perform subpixel localization on 5x5-area around max pixel
                [subpixloc(cc_bins,:)] = Phasor_localization_SMALLLABS( ...
                    cc(maxpix(cc_bins,1)-2:maxpix(cc_bins,1)+2, ...
                    maxpix(cc_bins,2)-2:maxpix(cc_bins,2)+2), ...
                    0,'');
                %Store subpixel location
                totpix(cc_bins,:) = maxpix(cc_bins,:)+subpixloc(cc_bins,1:2);
            end
            %Normalize subpixel location of cc-dift
            totpix = (totpix-totpix(1,:))*-1;
            totpix = totpix./paramscc.driftcorrcc_lateralsubpix;
            %Interpolate with spline with required sampling interval
            driftZinterpolated = interp1(meanframe,totpix(:,2),[1:totframes],'spline');
            
            minzposbeforecc = min(fits.zpos);
            maxzposbeforecc = max(fits.zpos);
            fits.zpos = fits.zpos-(driftZinterpolated(fits.frame)'-driftZinterpolated(1))/(movsz(2)/(maxzposbeforecc-minzposbeforecc));
        end
        %for now, plot drift over time
        figure(2);clf(2);
        hold on
        if isfield(fits,'zpos')
            yyaxis left
        end
        plot([1:totframes],finaldrift(:,1),'k-')
        plot([1:totframes],finaldrift(:,2),'c-')
        ylabel('XY drift (pixels)')
        grid on
        if isfield(fits,'zpos')
            yyaxis right
            plot(meanframe,1000*totpix(:,2)/(movsz(2)/(maxzposbeforecc-minzposbeforecc)),'r*')
            plot([1:totframes],1000*driftZinterpolated/(movsz(2)/(maxzposbeforecc-minzposbeforecc)),'r-')
            ylabel('Z drift (nm)')
        end
        
        avgshifthist(fits,params,movsz,[dlocs{ii},filesep,dnames{ii}],1);
        % %%
        %     figure(2);clf(2);
        %     hold on
        %     scatter([1:max(fits.frame)], finaldrift(:,1),1,'filled')
        %     scatter([1:max(fits.frame)], finaldrift(:,2),1,'filled')
        % %     scatter([1:max(fits.frame)], finaldriftz(:,1),1,'filled')
        % %     scatter([1:max(fits.frame)], finaldriftz(:,2),1,'filled')
        %     grid on
        %     fits.zpos = coordscorrz(:,2)/10;
        %     avgshifthist(fits,params,movsz,[dlocs{ii},filesep,dnames{ii}],1);
        
        if params.bgsub
            save([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits_driftcorrcc.mat'],'fits','movsz');
        else
            save([dlocs{ii},filesep,dnames{ii},'_fits_driftcorrcc.mat'],'fits','movsz');
        end
    end
    %% Make a average-shifted histograms (adapted from ThunderSTORM)
    %added by KM, 2019-04-10
    if (params.makeavgshifthist && params.specificframeanalysis == 0)
        if verbose
            disp([char(datetime),'   Making average-shifted histogram'])
        end
        if params.saveloadmat
            if params.bgsub
                %try loading in the fits & tracking results
                if params.driftcorr_cc
                    load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits_driftcorrcc.mat'],'fits');
                else
                    load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits');
                end
            else
                if params.driftcorr_cc
                    load([dlocs{ii},filesep,dnames{ii},'_fits_driftcorrcc.mat'],'fits');
                else
                    load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits');
                end
            end
        end
        avgshifthist(fits,params,movsz,[dlocs{ii},filesep,dnames{ii}],1);
    end
    %% Make ThunderSTORM-like csv output
    
    if( params.makeThSTORMoutput && params.specificframeanalysis == 0)
        if verbose
        disp([char(datetime),'   Making CSV output'])
        end
        try
            if params.saveloadmat
                if params.bgsub
                    %try loading in the fits & tracking results
                    if params.driftcorr_cc
                        load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits_driftcorrcc.mat'],'fits');
                    else
                        load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits');
                    end
                else
                    if params.driftcorr_cc
                        load([dlocs{ii},filesep,dnames{ii},'_fits_driftcorrcc.mat'],'fits');
                    else
                        load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits');
                    end
                end
            end
            array = zeros(size(fits.row,1),5);
            array(:,1) = [1:size(fits.row,1)]'; % id
            array(:,2) = fits.frame; %frame

    %         array(:,3:4) = [fits.col fits.row]*100; %x,y pos (nm) - transposed due to matlab shenanigans

    %         array(:,3:4) = [fits.col-.1440477+0.0184402 fits.row+.1256075+0.0184402]*params.avgshifthist_pixelsize; %x,y pos (nm) - transposed due to matlab shenanigans - corrected for fitting on beads smlmch
            array(:,3:4) = [fits.col fits.row]*params.avgshifthist_pixelsize; %x,y pos (nm) - transposed due to matlab shenanigans - corrected for fitting on beads smlmch
            try
            array(:,5) = fits.zpos*1000; % zpos (nm)
            catch
                array(:,5) = zeros(size(fits.frame,1),1);
            end
            %export to OUT file in same folder, same name
                temp = num2cell('"id","frame","x [nm]","y [nm]","z [nm]"');
                %,"intensity [photon]","offset [photon]","bkgstd [photon]","sigma1 [nm]","sigma2 [nm]"
                dlmwrite([dlocs{ii},filesep,dnames{ii},'_ThSTORMinput.csv'],temp, '')
                dlmwrite([dlocs{ii},filesep,dnames{ii},'_ThSTORMinput.csv'],array,'precision',10,'-append');
        catch
            fprintf('Error writing CSV file!\n')
        end
    end
    %% Make the ViewFits movie
    % Make a ViewFits movie, or just go into debug mode, to look at the
    % results. Outpits an avi file called moviename_ViewFits.avi
    if( params.makeViewFits && params.specificframeanalysis == 0)
        if params.useGUI
            drawGUIloadbar((7*ii)/numel(dlocs)/7,{['Making Viewfits ',dnames{ii}]},GUIhandles) 
        else
            try; waitbar((7*ii)/numel(dlocs)/7,h2,{['Making Viewfits ',dnames{ii}],'Overall Progress'}); end
        end
        if params.bgsub
            %try loading in the fits & tracking results
            if params.saveloadmat
                load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits','trk_filt','tracks');
            end
            
            %set trk_filt to empty if it doesn't exist            
            if ~exist('trk_filt','var')
                trk_filt=[];
            end          
            if params.orig_movie
                if ~params.trackingVF
                    ViewFits([dlocs{ii},filesep,dnames{ii},'.mat'],...
                        mov,trk_filt,movsz,goodframe,fits,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                else
                    if ~exist('bgsub_mov','var')
                        try; bgsubmov=load([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'mov'); end
                        bgsub_mov=single(bgsubmov.mov);
                        clear bgsubmov
                    end
                    ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'.mat'],...
                        mov,movsz,tracks,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                end
            else
                if ~params.trackingVF
                    ViewFits([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],...
                        bgsub_mov,trk_filt,movsz,goodframe,fits,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                else
                    ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],...
                        bgsub_mov,movsz,tracks,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                end
            end
        else
            %try loading in the fits & tracking results
            if params.saveloadmat
                load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits','trk_filt','tracks');
            end
            
            %set trk_filt to empty if it doesn't exist            
            if ~exist('trk_filt','var')
                trk_filt=[];
            end
            
            if ~params.trackingVF
%                 ViewFits([dlocs{ii},filesep,dnames{ii},'.mat'],...
%                     mov,trk_filt,movsz,goodframe,fits,...
%                     params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                fits.row = fits.row*4;
                fits.col = fits.col*4;
                ViewFits([dlocs{ii},filesep,dnames{ii},'.mat'],...
                    imresize(mov,4,'nearest'),trk_filt,movsz,goodframe,fits,...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            else
                ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'.mat'],...
                    mov,movsz,tracks,...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            end
        end
        if params.saveloadmat
            clear fits trk_filt tracks
        end
    end
end

if params.useGUI
    drawGUIloadbar(1,{['Finished']},GUIhandles) 
else
    try
        delete(h2)
    end
end
tictoc=toc(wholeshabang);
% disp(num2str(tictoc))

%turn warning back on
warning('on','MATLAB:load:variableNotFound')

% try to remove gpufit path
try
    rmpath(genpath('gpufit'))
end
if verbose
    disp([char(datetime),'   Finished in ', num2str(round(tictoc,5)) ' seconds!'])
end
end