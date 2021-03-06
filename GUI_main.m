function varargout = GUI_main(varargin)
% GUI_MAIN MATLAB code for SMALLLABS
% Made by Koen J.A. Martens, Hohlbein Lab
% Wageningen University, The Netherlands
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_main

% Last Modified by GUIDE v2.5 17-Sep-2019 18:16:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_main_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_main is made visible.
function GUI_main_OpeningFcn(hObject, eventdata, handles, varargin)

%First get the path of GUI_main
fullpath = mfilename('fullpath');
fullpath = extractBefore(fullpath,'GUI_main');
addpath(genpath(fullpath));
handles.GUIsettingssavelocation = [fullpath 'GUIStoredSettings.mat'];
%hide loadbar axis
% set(handles.axes1,'visible', 'off');
axis(handles.axes1);
axis off
try
    %Following line throws errors, might be fixed via setappdata, getappdata
    load(handles.GUIsettingssavelocation,'fvals_to_save');
    %loop over handles
    structhandles = struct2cell(handles);
    for ii = 1:size(structhandles,1)
        for jj = 1:size(fvals_to_save)
            if strcmp(structhandles{ii,1}.Tag,fvals_to_save{jj,1}.Tag)
                    eval(strcat('handles.',structhandles{ii,1}.Tag,'.Value = fvals_to_save{',...
                    num2str(jj),',1}.Value;'));
                    eval(strcat('handles.',structhandles{ii,1}.Tag,'.String = fvals_to_save{',...
                    num2str(jj),',1}.String;'));
            end
        end
    end
catch
end
% Choose default command line output for GUI_main
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
%Update active/inactiveness of parts of the code
enable_disable_phasor_gauss(handles);
enable_disable_filtering(handles);
phasor_dropdown_visibility(handles);
enable_disable_bgsub(handles);
enable_disable_driftcorr(handles);
enable_disable_ASH(handles);
enable_disable_tracking(handles);
% UIWAIT makes GUI_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function varargout = GUI_main_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
disp(' ');
disp(' ');
disp([char(datetime),'   Starting SMALLLABS-pSMLM from GUI'])
fnames = fieldnames(handles);
fvals = struct2cell(handles);
fvals_to_save = {};
cntr = 1;
for ii = 1:size(fvals,1)
    if isa(fvals{ii,1},'matlab.ui.control.UIControl')
        fvals_to_save{cntr,1} = fvals{ii,1};
        cntr=cntr+1;
    end
end
save(handles.GUIsettingssavelocation,'fvals_to_save')

filterstr = 'wavelet';
if handles.filter_bpass_radio.Value
    filterstr = 'bpass';
end
fullmovie = 1;
if (handles.firstframes_rdio.Value || handles.spcfc_frame_rdio.Value)
    fullmovie = 0;
end

gaussiantype = 1;
if handles.which_gauss_2_radio.Value
    gaussiantype = 2;
elseif handles.which_gauss_3_radio.Value
    gaussiantype = 3;
end

phasorType = get(handles.popupmenu_phasor,'Value'); %Order 2D-astig-DH-SP-TP
switch phasorType
    case 1
        phasor_2d = 1;
        phasor_ast = 0;
        phasor_DH = 0;
        phasor_SP = 0;
        phasor_TP = 0;
    case 2         
        phasor_2d = 0;
        phasor_ast = 1;
        phasor_DH = 0;
        phasor_SP = 0;
        phasor_TP = 0;   
    case 3         
        phasor_2d = 0;
        phasor_ast = 0;
        phasor_DH = 1;
        phasor_SP = 0;
        phasor_TP = 0;   
    case 4         
        phasor_2d = 0;
        phasor_ast = 0;
        phasor_DH = 0;
        phasor_SP = 1;
        phasor_TP = 0;   
    case 5         
        phasor_2d = 0;
        phasor_ast = 0;
        phasor_DH = 0;
        phasor_SP = 0;
        phasor_TP = 1;   
end
% run program - batch if necessary
%Check if batch
if ~isempty(handles.uipanel1.UserData) %if this is true, then it's batch
    batchsize = size(handles.uipanel1.UserData{1,1},2);
    for k = 1:batchsize
        fullfilepath{1,k} =  [char(handles.uipanel1.UserData{1,2}) char(handles.uipanel1.UserData{1,1}(k))];
    end
    
else
    batchsize = 1;
    fullfilepath{1,1} = handles.filepath.String;
end
for batchiteration = 1:batchsize
SMALLLABS_main(fullfilepath{1,batchiteration},str2num(handles.dfrlmsz_edit.String),...
    str2num(handles.avgwin_edit.String), str2num(handles.moloffwin_edit.String),...
    'bgsub',handles.bgsub_chk.Value,...
    'makeAvgsub',handles.bgsub_chk.Value,...%handles.makeavgsubmovie_chk.Value,...
    'guessing',handles.performguess_chk.Value,...
    'checkGuesses',handles.checkguess_chk.Value,...
    'makeOffFrames',handles.bgsub_chk.Value,...%handles.makeoffframelist_chk.Value,...
    'fitting',handles.performguess_chk.Value,...
    'tracking',handles.performtracking_chk.Value,...
    'make_guessmovie',handles.makeguessmov_chk.Value,...
    'makeViewFits',handles.makeviewfitsmov_chk.Value,...
    'do_avg',handles.bgsub_type_mean_rdiobtn.Value,...
    'offset',str2num(handles.avgsub_offset_edit.String),...
    'filtermethod',filterstr,...
    'wvltfilter_std',str2num(handles.wvltfilter_std_edit.String),...
    'phasor_fit',handles.fitting_phasor_radio.Value,...
    'phasor_rad',str2num(handles.phasor_rad_edit.String),...
    'phasor_2d',phasor_2d,...%handles.phasor2d_radiobtn.Value,...
    'phasor_astig',phasor_ast,...%handles.phasorastig_radiobtn.Value,...
    'phasor_astigcal',handles.phasor_ast_calib.String,...
    'phasor_DH',phasor_DH,...%handles.phasordh_radiobtn.Value,...
    'phasor_DHcal',handles.phasor_DH_calib.String,...
    'phasor_SPTPradius',str2num(handles.phasor_SP_rad_edit.String)*phasor_SP+str2num(handles.phasor_TP_rad_edit.String)*phasor_TP,...
    'phasor_SP',phasor_SP,...%handles.phasorSP_radiobtn.Value,...
    'phasor_SPcal',handles.phasor_SP_calib.String,...
    'phasor_TP',phasor_TP,...%,handles.phasorTP_radiobtn.Value,...
    'phasor_TPcal',handles.phasor_TP_calib.String,...
    'makeavgshifthist',handles.makeavgshifthist_chk.Value,...
    'avgshifthist_res',str2num(handles.avgshifthist_res_edit.String),...
    'avgshifthist_lateralsubpix',str2num(handles.avgshfthist_latsubpx_edit.String),...
    'avgshifthist_lateralshifts',str2num(handles.avgshifthist_latshift_edit.String),...
    'avgshifthist_axialcolnrs',20,...%str2num(handles.avgshfthist_axcolnr_edit.String),...
    'avgshifthist_pixelsize',str2num(handles.pxsize_edit.String),...
    'fullmovie',fullmovie,...
    'nrframes',str2num(handles.edit_nrframes.String),...
    'specificframeanalysis',handles.spcfc_frame_rdio.Value,...
    'specificframe',str2num(handles.edit_specificframe.String),...
    'bpthrsh',str2num(handles.bpthrsh_edit.String),...
    'egdesz',str2num(handles.edgesz_edit.String),...
    'pctile_frame',handles.pctile_frame.Value,...
    'which_gaussian',gaussiantype,...
    'usegpu',handles.gaussGPU_checkbox.Value,...
    'MLE_fit',handles.gauss_mle_radio.Value,...
    'fit_ang',str2num(handles.fit_ang_edit.String),...
    'driftcorr_cc',handles.do_driftcorr_cc.Value,...
    'driftcorrcc_lateralsubpix',str2num(handles.edit_driftcorrcc_latsubpix.String),...
    'dirftcorrcc_temporalbins',str2num(handles.edit_driftcorrcc_tempbins.String),...
    'trackparams',[str2num(handles.edit_tracking1.String) str2num(handles.edit_tracking2.String) str2num(handles.edit_tracking3.String) str2num(handles.edit_tracking4.String) str2num(handles.edit_tracking5.String) str2num(handles.edit_tracking6.String) str2num(handles.edit_tracking7.String)],...
    'makeThSTORMoutput',handles.checkbox_ThStorm.Value,...
    'saveloadmat',handles.checkbox_savemat.Value,...
    'useGUI',1,...
    'handles',handles);
end

% --- Executes on button press in bgsub_chk.
function bgsub_chk_Callback(hObject, eventdata, handles)
enable_disable_bgsub(handles)


% --- Executes on button press in performguess_chk.
function performguess_chk_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);


% --- Executes on button press in performfitting_chk.
function performfitting_chk_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
%Store everything
try
    fnames = fieldnames(handles);
    fvals = struct2cell(handles);
    fvals_to_save = {};
    cntr = 1;
    for ii = 1:size(fvals,1)
        if isa(fvals{ii,1},'matlab.ui.control.UIControl')
            fvals_to_save{cntr,1} = fvals{ii,1};
            cntr=cntr+1;
        end
    end
    save(handles.GUIsettingssavelocation,'fvals_to_save')
catch
    fprintf('\n Could not store current vars!\n')
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in phasordh_radiobtn.
function phasordh_radiobtn_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function phasorastig_radiobtn_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function phasorSP_radiobtn_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function phasorTP_radiobtn_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function phasor2d_radiobtn_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function which_gauss_3_radio_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function which_gauss_1_radio_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function which_gauss_2_radio_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);

function filter_wvlt_radio_Callback(hObject, eventdata, handles)
enable_disable_filtering(handles);

function filter_bpass_radio_Callback(hObject, eventdata, handles)
enable_disable_filtering(handles);


% --- Executes on button press in browse_btn.
function browse_btn_Callback(hObject, eventdata, handles)
folderselect(handles,'filepath',1);



% --- Executes on button press in make_ASH.
function make_ASH_Callback(hObject, eventdata, handles)
%Get fits 
fldrstrt = pwd;
if isfolder(handles.filepath.String)
    fldrstrt = handles.filepath.String;
else
    lastslash = max(strfind(handles.filepath.String,'\'));
    if ~isempty(lastslash)
    prevfolder = extractBefore(handles.filepath.String,lastslash+1);
    if isfolder(prevfolder)
        fldrstrt = prevfolder;
    end
    end
end
curfldr = pwd;
cd(fldrstrt)
[fileASH,pathASH] = uigetfile({'*_fits.mat;*_fits_driftcorrcc.mat'});
cd(curfldr)

%Load fits
try
load([pathASH fileASH],'fits','movsz');
catch
    fprintf('Please choose a *_fits.mat file!\n')
end
% Get necessary parameters
p.avgshifthist_res = str2num(handles.avgshifthist_res_edit.String);
p.avgshifthist_lateralsubpix = str2num(handles.avgshfthist_latsubpx_edit.String);
p.avgshifthist_lateralshifts = str2num(handles.avgshifthist_latshift_edit.String);
p.avgshifthist_axialcolnrs = 50;%str2num(handles.avgshfthist_axcolnr_edit.String);
p.avgshifthist_pixelsize = str2num(handles.pxsize_edit.String);
%Run code
disp('Making average shifted histogram...')
avgshifthist(fits,p,movsz,[pathASH fileASH],1);
disp('Completed average shifted histogram...')


function pxsize_edit_Callback(hObject, eventdata, handles)
AiryDiscDiameter_pxunits = round((1.15*2*(1.22*str2num(handles.emm_wavelength_edit.String))/...
    (2*str2num(handles.objNA_edit.String)))/str2num(handles.pxsize_edit.String),2);
handles.dfrlmsz_edit.String = num2str(AiryDiscDiameter_pxunits);


function emm_wavelength_edit_Callback(hObject, eventdata, handles)
AiryDiscDiameter_pxunits = round((1.15*2*(1.22*str2num(handles.emm_wavelength_edit.String))/...
    (2*str2num(handles.objNA_edit.String)))/str2num(handles.pxsize_edit.String),2);
handles.dfrlmsz_edit.String = num2str(AiryDiscDiameter_pxunits);


function objNA_edit_Callback(hObject, eventdata, handles)
AiryDiscDiameter_pxunits = round((1.15*2*(1.22*str2num(handles.emm_wavelength_edit.String))/...
    (2*str2num(handles.objNA_edit.String)))/str2num(handles.pxsize_edit.String),2);
handles.dfrlmsz_edit.String = num2str(AiryDiscDiameter_pxunits);


function makeavgshifthist_chk_Callback(hObject, eventdata, handles)
enable_disable_ASH(handles)

function performtracking_chk_Callback(hObject, eventdata, handles)
enable_disable_tracking(handles)

% --- Executes on button press in fitting_gauss_radio.
function fitting_gauss_radio_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);


% --- Executes on button press in fitting_gauss_radio.
function fitting_phasor_radio_Callback(hObject, eventdata, handles)
enable_disable_phasor_gauss(handles);


function popupmenu_phasor_Callback(hObject, eventdata, handles)
phasor_dropdown_visibility(handles)


function do_driftcorr_cc_Callback(hObject, eventdata, handles)
enable_disable_driftcorr(handles)


% --- Executes on button press in pushbutton_phasor_ast_calib.
function pushbutton_phasor_ast_calib_Callback(hObject, eventdata, handles)
folderselect(handles,'phasor_ast_calib',2);


% --- Executes on button press in pushbutton_phasor_dh_calib.
function pushbutton_phasor_dh_calib_Callback(hObject, eventdata, handles)
folderselect(handles,'phasor_DH_calib',2);


% --- Executes on button press in pushbutton_phasor_sp_calib.
function pushbutton_phasor_sp_calib_Callback(hObject, eventdata, handles)
folderselect(handles,'phasor_SP_calib',2);


% --- Executes on button press in pushbutton_phasor_tp_calib.
function pushbutton_phasor_tp_calib_Callback(hObject, eventdata, handles)
folderselect(handles,'phasor_TP_calib',2);


% --- Executes on button press in pushbutton_perform_DH_calib.
function pushbutton_perform_DH_calib_Callback(hObject, eventdata, handles)
%Ask for calibration movie
[file_calibmovie,path_calibmovie] = fileselecttoopen(handles,'phasor_DH_calib',{'*.tiff;*.tif';'*.*'},'Select calibration movie');
%Ask for location to write
if (file_calibmovie~=0)
    [file,path] = folderselecttowrite(handles,'phasor_DH_calib','*.mat','Select storage location');
    if file ~= 0
        prompt = {'Start z-position (�m):','End z-position (�m):','Z-position step size (�m):'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'-0.75','0.75','0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        %Perform 2D phasor via SMALL-LABS and the current settings (if
        %applicable)
        filterstr = 'wavelet';
        if handles.filter_bpass_radio.Value
            filterstr = 'bpass';
        end
        SMALLLABS_main([path_calibmovie file_calibmovie],str2num(handles.dfrlmsz_edit.String),...
            str2num(handles.avgwin_edit.String), str2num(handles.moloffwin_edit.String),...
            'bgsub',0,...
            'makeAvgsub',0,...%handles.makeavgsubmovie_chk.Value,...
            'guessing',1,...
            'checkGuesses',0,...
            'makeOffFrames',0,...%handles.makeoffframelist_chk.Value,...
            'fitting',1,...
            'tracking',0,...
            'makeViewFits',0,...
            'do_avg',0,...
            'offset',str2num(handles.avgsub_offset_edit.String),...
            'filtermethod',filterstr,...
            'wvltfilter_std',str2num(handles.wvltfilter_std_edit.String),...
            'phasor_fit',1,...
            'phasor_rad',str2num(handles.phasor_rad_edit.String),...
            'phasor_2d',1,...%handles.phasor2d_radiobtn.Value,...
            'phasor_astig',0,...%handles.phasorastig_radiobtn.Value,...
            'phasor_astigcal','',...
            'phasor_DH',0,...%handles.phasordh_radiobtn.Value,...
            'phasor_DHcal','',...
            'phasor_SP',0,...%handles.phasorSP_radiobtn.Value,...
            'phasor_SPcal','',...
            'phasor_TP',0,...%,handles.phasorTP_radiobtn.Value,...
            'phasor_TPcal','',...
            'makeavgshifthist',0,...
            'avgshifthist_res',0,...
            'avgshifthist_lateralsubpix',0,...
            'avgshifthist_lateralshifts',0,...
            'avgshifthist_axialcolnrs',0,...
            'avgshifthist_pixelsize',str2num(handles.pxsize_edit.String),...
            'fullmovie',1,...
            'nrframes',0,...
            'specificframeanalysis',0,...
            'specificframe',0,...
            'bpthrsh',str2num(handles.bpthrsh_edit.String),...
            'egdesz',str2num(handles.edgesz_edit.String),...
            'pctile_frame',handles.pctile_frame.Value,...
            'which_gaussian',1,...
            'usegpu',0,...
            'MLE_fit',1,...
            'fit_ang',0,...
            'driftcorr_cc',0,...
            'driftcorrcc_lateralsubpix',0,...
            'dirftcorrcc_temporalbins',0,...
            'useGUI',1,...
            'handles',handles);
        %Read resulting .mat data
        load([path_calibmovie file_calibmovie(1:end-4) '_fits.mat']);
        input = [fits.frame,fits.col,fits.row zeros(size(fits.frame,1),2)];
        pSMLM_DH_calibration_SL(input,[str2num(answer{1,1}):str2num(answer{3,1}):str2num(answer{2,1})],1,[path file]);
        delete_SMALLLABS_files(path_calibmovie,1,1,1,1,1);
    end
end

function pushbutton_perform_SP_calib_Callback(hObject, eventdata, handles)
% Ask for calibration movie
[file_calibmovie,path_calibmovie] = fileselecttoopen(handles,'phasor_SP_calib',{'*.tiff;*.tif';'*.*'},'Select calibration movie');
%Ask for location to write
if (file_calibmovie~=0)
    [file,path] = folderselecttowrite(handles,'phasor_SP_calib','*.mat','Select storage location');
    if file ~= 0
        prompt = {'Start z-position (�m):','End z-position (�m):','Z-position step size (�m):'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'-0.75','0.75','0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        %Perform 2D phasor via SMALL-LABS and the current settings (if
        %applicable)
        filterstr = 'wavelet';
        if handles.filter_bpass_radio.Value
            filterstr = 'bpass';
        end
        SMALLLABS_main([path_calibmovie file_calibmovie],str2num(handles.dfrlmsz_edit.String),...
            str2num(handles.avgwin_edit.String), str2num(handles.moloffwin_edit.String),...
            'bgsub',0,...
            'makeAvgsub',0,...%handles.makeavgsubmovie_chk.Value,...
            'guessing',1,...
            'checkGuesses',0,...
            'makeOffFrames',0,...%handles.makeoffframelist_chk.Value,...
            'fitting',1,...
            'tracking',0,...
            'makeViewFits',0,...
            'do_avg',0,...
            'offset',str2num(handles.avgsub_offset_edit.String),...
            'filtermethod',filterstr,...
            'wvltfilter_std',str2num(handles.wvltfilter_std_edit.String),...
            'phasor_fit',1,...
            'phasor_rad',str2num(handles.phasor_rad_edit.String),...
            'phasor_2d',1,...%handles.phasor2d_radiobtn.Value,...
            'phasor_astig',0,...%handles.phasorastig_radiobtn.Value,...
            'phasor_astigcal','',...
            'phasor_DH',0,...%handles.phasordh_radiobtn.Value,...
            'phasor_DHcal','',...
            'phasor_SP',0,...%handles.phasorSP_radiobtn.Value,...
            'phasor_SPcal','',...
            'phasor_TP',0,...%,handles.phasorTP_radiobtn.Value,...
            'phasor_TPcal','',...
            'makeavgshifthist',0,...
            'avgshifthist_res',0,...
            'avgshifthist_lateralsubpix',0,...
            'avgshifthist_lateralshifts',0,...
            'avgshifthist_axialcolnrs',0,...
            'avgshifthist_pixelsize',str2num(handles.pxsize_edit.String),...
            'fullmovie',1,...
            'nrframes',0,...
            'specificframeanalysis',0,...
            'specificframe',0,...
            'bpthrsh',str2num(handles.bpthrsh_edit.String),...
            'egdesz',str2num(handles.edgesz_edit.String),...
            'pctile_frame',handles.pctile_frame.Value,...
            'which_gaussian',1,...
            'usegpu',0,...
            'MLE_fit',1,...
            'fit_ang',0,...
            'driftcorr_cc',0,...
            'driftcorrcc_lateralsubpix',0,...
            'dirftcorrcc_temporalbins',0,...
            'useGUI',1,...
    'handles',handles);
        %Read resulting .mat data
        fileextstartpos = (strfind(file_calibmovie,'.tif'));
        load([path_calibmovie file_calibmovie(1:fileextstartpos-1) '_fits.mat']);
        
        %open movie file again
        mov = TiffLoader_SL([path_calibmovie,file_calibmovie]);
        
        %some variables
        ROIsize = (str2num(handles.phasor_SP_rad_edit.String)*2+1);
        halfROI = (ROIsize-1)/2;
        
        %Link 2D points to get midpoints
        SPcalibmidpoints = SP_TP_phasor_linking([fits.frame fits.row fits.col],15);
        
        %Make imagestack
        SPcalibimagestack = zeros(ROIsize,ROIsize,max(SPcalibmidpoints(:,1)),1);
        
        %Fill imagestack
        for i = 1:size(SPcalibmidpoints,1)
            if i > 2
                if SPcalibmidpoints(i,1) > SPcalibmidpoints(i-1,1)
                    counter = 1;
                end
            end
            try
                SPcalibimagestack(:,:,SPcalibmidpoints(i,1),counter) = mov(round(SPcalibmidpoints(i,2))-halfROI:round(SPcalibmidpoints(i,2))+halfROI,...
                    round(SPcalibmidpoints(i,3))-halfROI:round(SPcalibmidpoints(i,3))+halfROI,...
                    round(SPcalibmidpoints(i,1)));
                counter = counter+1;
            end
        end
        SPcalibimagestack(isnan(SPcalibimagestack)) = 0;
        
        SPTP_phasor_calibrationRoutine_SL(SPcalibimagestack,[str2num(answer{1,1}):str2num(answer{3,1}):str2num(answer{2,1})],[path file])
        delete_SMALLLABS_files(path_calibmovie,1,1,1,1,1);
    end
end


function pushbutton_perform_TP_calib_Callback(hObject, eventdata, handles)
% Ask for calibration movie
[file_calibmovie,path_calibmovie] = fileselecttoopen(handles,'phasor_TP_calib',{'*.tiff;*.tif';'*.*'},'Select calibration movie');
%Ask for location to write
if (file_calibmovie~=0)
    [file,path] = folderselecttowrite(handles,'phasor_TP_calib','*.mat','Select storage location');
    if file ~= 0
        prompt = {'Start z-position (�m):','End z-position (�m):','Z-position step size (�m):'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'-0.75','0.75','0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        %Perform 2D phasor via SMALL-LABS and the current settings (if
        %applicable)
        filterstr = 'wavelet';
        if handles.filter_bpass_radio.Value
            filterstr = 'bpass';
        end
        SMALLLABS_main([path_calibmovie file_calibmovie],str2num(handles.dfrlmsz_edit.String),...
            str2num(handles.avgwin_edit.String), str2num(handles.moloffwin_edit.String),...
            'bgsub',0,...
            'makeAvgsub',0,...%handles.makeavgsubmovie_chk.Value,...
            'guessing',1,...
            'checkGuesses',0,...
            'makeOffFrames',0,...%handles.makeoffframelist_chk.Value,...
            'fitting',1,...
            'tracking',0,...
            'makeViewFits',0,...
            'do_avg',0,...
            'offset',str2num(handles.avgsub_offset_edit.String),...
            'filtermethod',filterstr,...
            'wvltfilter_std',str2num(handles.wvltfilter_std_edit.String),...
            'phasor_fit',1,...
            'phasor_rad',str2num(handles.phasor_rad_edit.String),...
            'phasor_2d',1,...%handles.phasor2d_radiobtn.Value,...
            'phasor_astig',0,...%handles.phasorastig_radiobtn.Value,...
            'phasor_astigcal','',...
            'phasor_DH',0,...%handles.phasordh_radiobtn.Value,...
            'phasor_DHcal','',...
            'phasor_SP',0,...%handles.phasorSP_radiobtn.Value,...
            'phasor_SPcal','',...
            'phasor_TP',0,...%,handles.phasorTP_radiobtn.Value,...
            'phasor_TPcal','',...
            'makeavgshifthist',0,...
            'avgshifthist_res',0,...
            'avgshifthist_lateralsubpix',0,...
            'avgshifthist_lateralshifts',0,...
            'avgshifthist_axialcolnrs',0,...
            'avgshifthist_pixelsize',str2num(handles.pxsize_edit.String),...
            'fullmovie',1,...
            'nrframes',0,...
            'specificframeanalysis',0,...
            'specificframe',0,...
            'bpthrsh',str2num(handles.bpthrsh_edit.String),...
            'egdesz',str2num(handles.edgesz_edit.String),...
            'pctile_frame',handles.pctile_frame.Value,...
            'which_gaussian',1,...
            'usegpu',0,...
            'MLE_fit',1,...
            'fit_ang',0,...
            'driftcorr_cc',0,...
            'driftcorrcc_lateralsubpix',0,...
            'dirftcorrcc_temporalbins',0,...
            'useGUI',1,...
            'handles',handles);
        %Read resulting .mat data
        load([path_calibmovie file_calibmovie(1:end-5) '_fits.mat']);
        
        %open movie file again
        mov = TiffLoader_SL([path_calibmovie,file_calibmovie]);
        
        %some variables
        ROIsize = (str2num(handles.phasor_TP_rad_edit.String)*2+1);
        halfROI = (ROIsize-1)/2;
        
        %Link 2D points to get midpoints
        TPcalibmidpoints = SP_TP_phasor_linking([fits.frame fits.row fits.col],15);
        
        %Make imagestack
        TPcalibimagestack = zeros(ROIsize,ROIsize,max(TPcalibmidpoints(:,1)),1);
        
        %Fill imagestack
        for i = 1:size(TPcalibmidpoints,1)
            if i > 2
                if TPcalibmidpoints(i,1) > TPcalibmidpoints(i-1,1)
                    counter = 1;
                end
            end
            try
                TPcalibimagestack(:,:,TPcalibmidpoints(i,1),counter) = mov(round(TPcalibmidpoints(i,2))-halfROI:round(TPcalibmidpoints(i,2))+halfROI,...
                    round(TPcalibmidpoints(i,3))-halfROI:round(TPcalibmidpoints(i,3))+halfROI,...
                    round(TPcalibmidpoints(i,1)));
                counter = counter+1;
            end
        end
        TPcalibimagestack(isnan(TPcalibimagestack)) = 0;
        
        SPTP_phasor_calibrationRoutine_SL(TPcalibimagestack,[str2num(answer{1,1}):str2num(answer{3,1}):str2num(answer{2,1})],[path file])
        delete_SMALLLABS_files(path_calibmovie,1,1,1,1,1);
    end
end
%% EMPTY CALLBACKS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)

function edit_specificframe_Callback(hObject, eventdata, handles)
function edit_driftcorrcc_tempbins_Callback(hObject, eventdata, handles)
function edit_driftcorrcc_latsubpix_Callback(hObject, eventdata, handles)
function fit_ang_edit_Callback(hObject, eventdata, handles)
function gaussGPU_checkbox_Callback(hObject, eventdata, handles)
function pctile_frame_Callback(hObject, eventdata, handles)
function edgesz_edit_Callback(hObject, eventdata, handles)
function bpthrsh_edit_Callback(hObject, eventdata, handles)
function edit_nrframes_Callback(hObject, eventdata, handles)
function phasor_ast_calib_Callback(hObject, eventdata, handles)
function phasor_SP_calib_Callback(hObject, eventdata, handles)
function phasor_TP_calib_Callback(hObject, eventdata, handles)
function avgshifthist_res_edit_Callback(hObject, eventdata, handles)
function avgshfthist_latsubpx_edit_Callback(hObject, eventdata, handles)
function avgshfthist_axcolnr_edit_Callback(hObject, eventdata, handles)
function avgshifthist_latshift_edit_Callback(hObject, eventdata, handles)
function wvltfilter_std_edit_Callback(hObject, eventdata, handles)
function phasor_rad_edit_Callback(hObject, eventdata, handles)
function phasor_DH_calib_Callback(hObject, eventdata, handles)
function avgsub_offset_edit_Callback(hObject, eventdata, handles)
function dfrlmsz_edit_Callback(hObject, eventdata, handles)
function avgwin_edit_Callback(hObject, eventdata, handles)
function moloffwin_edit_Callback(hObject, eventdata, handles)
function radiobutton4_Callback(hObject, eventdata, handles)
function radiobutton5_Callback(hObject, eventdata, handles)
function makeviewfitsmov_chk_Callback(hObject, eventdata, handles)
function filepath_Callback(hObject, eventdata, handles)
function checkguess_chk_Callback(hObject, eventdata, handles)
function makeoffframelist_chk_Callback(hObject, eventdata, handles)
function makeavgsubmovie_chk_Callback(hObject, eventdata, handles)

function pushbutton_perform_astig_calib_Callback(hObject, eventdata, handles)
% Ask for calibration movie
[file_calibmovie,path_calibmovie] = fileselecttoopen(handles,'phasor_ast_calib',{'*.tiff;*.tif';'*.*'},'Select calibration movie');
%Ask for location to write
if (file_calibmovie~=0)
    [file,path] = folderselecttowrite(handles,'phasor_ast_calib','*.mat','Select storage location');
    if file ~= 0
        prompt = {'Start z-position (�m):','End z-position (�m):','Z-position step size (�m):'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'-0.75','0.75','0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        %Perform 2D phasor via SMALL-LABS and the current settings (if
        %applicable)
        filterstr = 'wavelet';
        if handles.filter_bpass_radio.Value
            filterstr = 'bpass';
        end
        SMALLLABS_main([path_calibmovie file_calibmovie],str2num(handles.dfrlmsz_edit.String),...
            str2num(handles.avgwin_edit.String), str2num(handles.moloffwin_edit.String),...
            'bgsub',0,...
            'makeAvgsub',0,...%handles.makeavgsubmovie_chk.Value,...
            'guessing',1,...
            'checkGuesses',0,...
            'makeOffFrames',0,...%handles.makeoffframelist_chk.Value,...
            'fitting',1,...
            'tracking',0,...
            'makeViewFits',0,...
            'do_avg',0,...
            'offset',str2num(handles.avgsub_offset_edit.String),...
            'filtermethod',filterstr,...
            'wvltfilter_std',str2num(handles.wvltfilter_std_edit.String),...
            'phasor_fit',1,...
            'phasor_rad',str2num(handles.phasor_rad_edit.String),...
            'phasor_2d',1,...%handles.phasor2d_radiobtn.Value,...
            'phasor_astig',0,...%handles.phasorastig_radiobtn.Value,...
            'phasor_astigcal','',...
            'phasor_DH',0,...%handles.phasordh_radiobtn.Value,...
            'phasor_DHcal','',...
            'phasor_SP',0,...%handles.phasorSP_radiobtn.Value,...
            'phasor_SPcal','',...
            'phasor_TP',0,...%,handles.phasorTP_radiobtn.Value,...
            'phasor_TPcal','',...
            'makeavgshifthist',0,...
            'avgshifthist_res',0,...
            'avgshifthist_lateralsubpix',0,...
            'avgshifthist_lateralshifts',0,...
            'avgshifthist_axialcolnrs',0,...
            'avgshifthist_pixelsize',str2num(handles.pxsize_edit.String),...
            'fullmovie',1,...
            'nrframes',0,...
            'specificframeanalysis',0,...
            'specificframe',0,...
            'bpthrsh',str2num(handles.bpthrsh_edit.String),...
            'egdesz',str2num(handles.edgesz_edit.String),...
            'pctile_frame',handles.pctile_frame.Value,...
            'which_gaussian',1,...
            'usegpu',0,...
            'MLE_fit',1,...
            'fit_ang',0,...
            'driftcorr_cc',0,...
            'driftcorrcc_lateralsubpix',0,...
            'dirftcorrcc_temporalbins',0,...
            'useGUI',1,...
            'handles',handles);
        %Read resulting .mat data
        fileextstartpos = (strfind(file_calibmovie,'.tif'));
        load([path_calibmovie file_calibmovie(1:fileextstartpos-1) '_fits.mat']);
        
        %open movie file again
        mov = TiffLoader_SL([path_calibmovie,file_calibmovie]);
        
        %some variables
        ROIsize = (str2num(handles.phasor_rad_edit.String)*2+1);
        halfROI = (ROIsize-1)/2;
        
        %Link 2D points to get midpoints
        astigcalibmidpoints = [fits.frame fits.row fits.col];
        
        %Make imagestack
        astigcalibimagestack = zeros(ROIsize,ROIsize,max(astigcalibmidpoints(:,1)),1);
%         keyboard
        %Fill imagestack
        for i = 1:size(astigcalibmidpoints,1)
            if i > 2
                if astigcalibmidpoints(i,1) > astigcalibmidpoints(i-1,1)
                    counter = 1;
                end
            end
            try
                astigcalibimagestack(:,:,astigcalibmidpoints(i,1),counter) = mov(round(astigcalibmidpoints(i,2))-halfROI:round(astigcalibmidpoints(i,2))+halfROI,...
                    round(astigcalibmidpoints(i,3))-halfROI:round(astigcalibmidpoints(i,3))+halfROI,...
                    round(astigcalibmidpoints(i,1)));
                counter = counter+1;
            end
        end
        astigcalibimagestack(isnan(astigcalibimagestack)) = 0;
        
        astig_phasor_calibrationRoutine_SL(astigcalibimagestack,[str2num(answer{1,1}):str2num(answer{3,1}):str2num(answer{2,1})],[path file])
        delete_SMALLLABS_files(path_calibmovie,1,1,1,1,1);
    end
end

%% Additional functions -----------------------------------------------------------------------------------------------------------------------------------------------
%  Additional functions -----------------------------------------------------------------------------------------------------------------------------------------------
function enable_disable_phasor_gauss(handles)
genarronoff{1} = 'off';
genarronoff{2} = 'on';
if handles.performguess_chk.Value
    if handles.fitting_phasor_radio.Value
        colGauss = [0.651 0.651 0.651];
        colPhasor = [0 0 0];
        gaussen = 'off';
        phasoren = 'on';
    else
        colPhasor = [0.651 0.651 0.651];
        colGauss = [0 0 0];
        gaussen = 'on';
        phasoren = 'off';
    end
else
    gaussen = 'off';
    phasoren = 'off';
    colPhasor = [0.651 0.651 0.651];
    colGauss = [0.651 0.651 0.651];
end

switch handles.performguess_chk.Value
    case 0
        handles.checkguess_chk.Enable = 'off';
    case 1
        handles.checkguess_chk.Enable = 'on';
end
        
handles.gauss_mle_radio.Enable = gaussen;
handles.gauss_ls_radio.Enable = gaussen;
handles.gaussGPU_checkbox.Enable = gaussen;
handles.uibuttongroup7.ForegroundColor = colGauss;
handles.uibuttongroup8.ForegroundColor = colGauss;
handles.which_gauss_1_radio.Enable = gaussen;
handles.which_gauss_2_radio.Enable = gaussen;
handles.which_gauss_3_radio.Enable = gaussen;
handles.fit_ang_edit.Enable = genarronoff{1+(handles.performguess_chk.Value*handles.fitting_gauss_radio.Value*handles.which_gauss_3_radio.Value)};
handles.text33.Enable = genarronoff{1+(handles.performguess_chk.Value*handles.fitting_gauss_radio.Value*handles.which_gauss_3_radio.Value)};
handles.uibuttongroup5.ForegroundColor = colPhasor;
handles.text9.Enable = phasoren;
handles.phasor_rad_edit.Enable = phasoren;
handles.phasor_ast_calib.Enable = phasoren;
handles.phasor_DH_calib.Enable = phasoren;
handles.phasor_SP_calib.Enable = phasoren;
handles.phasor_TP_calib.Enable = phasoren;
handles.pushbutton_phasor_ast_calib.Enable = phasoren;
handles.pushbutton_phasor_dh_calib.Enable = phasoren;
handles.pushbutton_phasor_sp_calib.Enable = phasoren;
handles.pushbutton_phasor_tp_calib.Enable = phasoren;
handles.popupmenu_phasor.Enable = phasoren;
handles.text22.Enable = phasoren;
handles.text38.Enable = phasoren;
handles.fitting_phasor_radio.Enable = genarronoff{1+handles.performguess_chk.Value};
handles.fitting_gauss_radio.Enable = genarronoff{1+handles.performguess_chk.Value};


function enable_disable_filtering(handles)
genarronoff{1} = 'off';
genarronoff{2} = 'on';
handles.text21.Enable = genarronoff{1+handles.filter_wvlt_radio.Value};
handles.wvltfilter_std_edit.Enable = genarronoff{1+handles.filter_wvlt_radio.Value};
handles.pctile_frame.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.text30.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.text31.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.bpthrsh_edit.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.edgesz_edit.Enable = genarronoff{1+handles.filter_bpass_radio.Value};

function enable_disable_ASH(handles)
genarronoff{1} = 'off';
genarronoff{2} = 'on';
handles.text15.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.text16.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.text17.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.text18.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.text19.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.text20.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.avgshifthist_res_edit.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.avgshfthist_latsubpx_edit.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.avgshifthist_latshift_edit.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};
handles.avgshfthist_axcolnr_edit.Enable = genarronoff{1+handles.makeavgshifthist_chk.Value};

function enable_disable_tracking(handles)
genarronoff{1} = 'off';
genarronoff{2} = 'on';
handles.text40.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.text41.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.text42.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.text43.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.text44.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.text45.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.text46.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking1.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking2.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking3.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking4.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking5.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking6.Enable = genarronoff{1+handles.performtracking_chk.Value};
handles.edit_tracking7.Enable = genarronoff{1+handles.performtracking_chk.Value};

function enable_disable_driftcorr(handles)
genarronoff{1} = 'off';
genarronoff{2} = 'on';
handles.text34.Enable = genarronoff{1+handles.do_driftcorr_cc.Value};
handles.text35.Enable = genarronoff{1+handles.do_driftcorr_cc.Value};
handles.edit_driftcorrcc_latsubpix.Enable = genarronoff{1+handles.do_driftcorr_cc.Value};
handles.edit_driftcorrcc_tempbins.Enable = genarronoff{1+handles.do_driftcorr_cc.Value};

function enable_disable_bgsub(handles)
genarronoff{1} = 'off';
genarronoff{2} = 'on';
switch handles.bgsub_chk.Value
    case 0
        handles.uibuttongroup1.ForegroundColor = [0.651 0.651 0.651];
    case 1
        handles.uibuttongroup1.ForegroundColor = [0 0 0];
end        
handles.bgsub_type_median_rdiobtn.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.bgsub_type_mean_rdiobtn.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.text3.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.text5.Enable = genarronoff{1+handles.bgsub_chk.Value};
% handles.text4.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.text6.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.avgsub_offset_edit.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.avgwin_edit.Enable = genarronoff{1+handles.bgsub_chk.Value};
% handles.dfrlmsz_edit.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.text28.Enable = genarronoff{1+handles.bgsub_chk.Value};
handles.moloffwin_edit.Enable = genarronoff{1+handles.bgsub_chk.Value};

function phasor_dropdown_visibility(handles)
%set location of all calib text boxes to that of astigmatism (so no manual
%alignemnt is necessary)
handles.phasor_DH_calib.Position = handles.phasor_ast_calib.Position;
handles.phasor_SP_calib.Position = handles.phasor_ast_calib.Position;
handles.phasor_TP_calib.Position = handles.phasor_ast_calib.Position;
handles.pushbutton_phasor_dh_calib.Position = handles.pushbutton_phasor_ast_calib.Position;
handles.pushbutton_phasor_sp_calib.Position = handles.pushbutton_phasor_ast_calib.Position;
handles.pushbutton_phasor_tp_calib.Position = handles.pushbutton_phasor_ast_calib.Position;
handles.pushbutton_perform_DH_calib.Position = handles.pushbutton_perform_astig_calib.Position;
handles.pushbutton_perform_SP_calib.Position = handles.pushbutton_perform_astig_calib.Position;
handles.pushbutton_perform_TP_calib.Position = handles.pushbutton_perform_astig_calib.Position;
handles.phasor_TP_rad_edit.Position = handles.phasor_SP_rad_edit.Position;
%Display the correct calibration file text, hide the others
switch handles.popupmenu_phasor.Value
    case 1 %2D
        handles.text22.Visible = 'off';
        handles.phasor_ast_calib.Visible = 'off';
        handles.pushbutton_phasor_ast_calib.Visible = 'off';
        handles.pushbutton_perform_astig_calib.Visible = 'off';
        handles.phasor_DH_calib.Visible = 'off';
        handles.pushbutton_phasor_dh_calib.Visible = 'off';
        handles.pushbutton_perform_DH_calib.Visible = 'off';
        handles.phasor_SP_calib.Visible = 'off';
        handles.pushbutton_phasor_sp_calib.Visible = 'off';
        handles.pushbutton_perform_SP_calib.Visible = 'off';
        handles.phasor_TP_calib.Visible = 'off';
        handles.pushbutton_phasor_tp_calib.Visible = 'off';
        handles.pushbutton_perform_TP_calib.Visible = 'off';
        handles.phasor_SP_rad_edit.Visible = 'off';
        handles.phasor_TP_rad_edit.Visible = 'off';
    case 2 %ast
        handles.text22.Visible = 'on';
        handles.phasor_ast_calib.Visible = 'on';
        handles.pushbutton_perform_astig_calib.Visible = 'on';
        handles.pushbutton_phasor_ast_calib.Visible = 'on';
        handles.phasor_DH_calib.Visible = 'off';
        handles.pushbutton_phasor_dh_calib.Visible = 'off';
        handles.pushbutton_perform_DH_calib.Visible = 'off';
        handles.phasor_SP_calib.Visible = 'off';
        handles.pushbutton_phasor_sp_calib.Visible = 'off';
        handles.pushbutton_perform_SP_calib.Visible = 'off';
        handles.phasor_TP_calib.Visible = 'off';
        handles.pushbutton_phasor_tp_calib.Visible = 'off';
        handles.pushbutton_perform_TP_calib.Visible = 'off';
        handles.phasor_SP_rad_edit.Visible = 'off';
        handles.phasor_TP_rad_edit.Visible = 'off';
    case 3 %DH
        handles.text22.Visible = 'on';
        handles.phasor_ast_calib.Visible = 'off';
        handles.pushbutton_phasor_ast_calib.Visible = 'off';
        handles.pushbutton_perform_astig_calib.Visible = 'off';
        handles.phasor_DH_calib.Visible = 'on';
        handles.pushbutton_phasor_dh_calib.Visible = 'on';
        handles.pushbutton_perform_DH_calib.Visible = 'on';
        handles.phasor_SP_calib.Visible = 'off';
        handles.pushbutton_phasor_sp_calib.Visible = 'off';
        handles.pushbutton_perform_SP_calib.Visible = 'off';
        handles.phasor_TP_calib.Visible = 'off';
        handles.pushbutton_phasor_tp_calib.Visible = 'off';
        handles.pushbutton_perform_TP_calib.Visible = 'off';
        handles.phasor_SP_rad_edit.Visible = 'off';
        handles.phasor_TP_rad_edit.Visible = 'off';
    case 4 %SP
        handles.text22.Visible = 'on';
        handles.phasor_ast_calib.Visible = 'off';
        handles.pushbutton_phasor_ast_calib.Visible = 'off';
        handles.pushbutton_perform_astig_calib.Visible = 'off';
        handles.phasor_DH_calib.Visible = 'off';
        handles.pushbutton_phasor_dh_calib.Visible = 'off';
        handles.pushbutton_perform_DH_calib.Visible = 'off';
        handles.phasor_SP_calib.Visible = 'on';
        handles.pushbutton_phasor_sp_calib.Visible = 'on';
        handles.pushbutton_perform_SP_calib.Visible = 'on';
        handles.phasor_TP_calib.Visible = 'off';
        handles.pushbutton_phasor_tp_calib.Visible = 'off';
        handles.pushbutton_perform_TP_calib.Visible = 'off';
        handles.phasor_SP_rad_edit.Visible = 'on';
        handles.phasor_TP_rad_edit.Visible = 'off';
    case 5 %TP
        handles.text22.Visible = 'on';
        handles.phasor_ast_calib.Visible = 'off';
        handles.pushbutton_phasor_ast_calib.Visible = 'off';
        handles.pushbutton_perform_astig_calib.Visible = 'off';
        handles.phasor_DH_calib.Visible = 'off';
        handles.pushbutton_phasor_dh_calib.Visible = 'off';
        handles.pushbutton_perform_DH_calib.Visible = 'off';
        handles.phasor_SP_calib.Visible = 'off';
        handles.pushbutton_phasor_sp_calib.Visible = 'off';
        handles.pushbutton_perform_SP_calib.Visible = 'off';
        handles.phasor_TP_calib.Visible = 'on';
        handles.pushbutton_phasor_tp_calib.Visible = 'on';
        handles.pushbutton_perform_TP_calib.Visible = 'on';
        handles.phasor_SP_rad_edit.Visible = 'off';
        handles.phasor_TP_rad_edit.Visible = 'on';
end

function folderselect(handles,handlename,getorput)
fullhandlename = strcat('handles.',string(handlename),'.String');
fldrstrt = pwd;
if isfolder(fullhandlename)
    fldrstrt = fullhandlename;
else
    lastslash = max(strfind(eval(fullhandlename),'\'));
    prevfolder = '';
    if ~isempty(lastslash)
        prevfolder = extractBefore(eval(fullhandlename),lastslash+1);
    end
    if isfolder(prevfolder)
        fldrstrt = prevfolder;
    end
end
curfldr = pwd;
cd(fldrstrt)
[file,path] = uigetfile('*.*','MultiSelect', 'on');
cd(curfldr) 
if ~isequal(file,0)
    if (size(file,2) > 1 && iscell(file)) %if batching
        eval(['handles.',handlename,'.String = "Multiple Files";']) ;
        handles.uipanel1.UserData = {file path};
    else
        handles.uipanel1.UserData = {};
        eval(['handles.',handlename,'.String = ''',path,file,''';']) ;
        %     strcat("handles.",string(handlename),".String{1,1}") = [path file];%strcat('handles.',string(handlename),'.String{1,1}') = [path file];
    end
else
    handles.uipanel1.UserData = {};
end



function [file,path]=folderselecttowrite(handles,handlename,ext,title)
fullhandlename = strcat('handles.',string(handlename),'.String');
fldrstrt = pwd;
if isfolder(fullhandlename)
    fldrstrt = fullhandlename;
else
    lastslash = max(strfind(eval(fullhandlename),'\'));
    prevfolder = '';
    if ~isempty(lastslash)
        prevfolder = extractBefore(eval(fullhandlename),lastslash+1);
    end
    if isfolder(prevfolder)
        fldrstrt = prevfolder;
    end
end
curfldr = pwd;
cd(fldrstrt)
[file,path] = uiputfile(ext,title);
cd(curfldr) 
   
if ~isequal(file,0)
    eval(['handles.',handlename,'.String = ''',path,file,''';']) ; 
end


function [file,path]=fileselecttoopen(handles,handlename,ext,title)
fullhandlename = strcat('handles.',string(handlename),'.String');
fldrstrt = pwd;
if isfolder(fullhandlename)
    fldrstrt = fullhandlename;
else
    lastslash = max(strfind(eval(fullhandlename),'\'));
    prevfolder = '';
    if ~isempty(lastslash)
        prevfolder = extractBefore(eval(fullhandlename),lastslash+1);
    end
    if isfolder(prevfolder)
        fldrstrt = prevfolder;
    end
end
curfldr = pwd;
cd(fldrstrt)
[file,path] = uigetfile(ext,title);
cd(curfldr) 
   
if ~isequal(file,0)
    eval(['handles.',handlename,'.String = ''',path,file,''';']) ; 
end


% --- Executes on button press in checkbox_ThStorm.
function checkbox_ThStorm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ThStorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ThStorm


% --- Executes on button press in checkbox_savemat.
function checkbox_savemat_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_savemat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_savemat
