function varargout = GUI_main(varargin)
% GUI_MAIN MATLAB code for GUI_main.fig
%      GUI_MAIN, by itself, creates a new GUI_MAIN or raises the existing
%      singleton*.
%
%      H = GUI_MAIN returns the handle to a new GUI_MAIN or the handle to
%      the existing singleton*.
%
%      GUI_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MAIN.M with the given input arguments.
%
%      GUI_MAIN('Property','Value',...) creates a new GUI_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_main

% Last Modified by GUIDE v2.5 09-May-2019 15:44:12

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
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_main (see VARARGIN)
try
    %Following line throws errors, might be fixed via setappdata, getappdata
    load('GUI\StoredSettings.mat','fvals_to_save');
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
enable_disable_filtering(handles)
% UIWAIT makes GUI_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(' ');
disp(' ');
disp([char(datetime),'   Starting SMALLLABS_JH from GUI'])

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
save('GUI\StoredSettings.mat','fvals_to_save')

filterstr = 'wavelet';
if handles.filter_bpass_radio.Value
    filterstr = 'bpass';
end
fullmovie = 1;
if handles.firstframes_rdio.Value
    fullmovie = 0;
end

gaussiantype = 1;
if handles.which_gauss_2_radio.Value
    gaussiantype = 2;
elseif handles.which_gauss_3_radio.Value
    gaussiantype = 3;
end
% run program
SMALLLABS_main(handles.filepath.String{1,1},str2num(handles.dfrlmsz_edit.String),...
    str2num(handles.avgwin_edit.String), str2num(handles.moloffwin_edit.String),...
    'bgsub',handles.bgsub_chk.Value,...
    'makeAvgsub',handles.makeavgsubmovie_chk.Value,...
    'guessing',handles.performguess_chk.Value,...
    'checkGuesses',handles.checkguess_chk.Value,...
    'makeOffFrames',handles.makeoffframelist_chk.Value,...
    'fitting',handles.performfitting_chk.Value,...
    'tracking',handles.performtracking_chk.Value,...
    'makeViewFits',handles.makeviewfitsmov_chk.Value,...
    'do_avg',handles.bgsub_type_mean_rdiobtn.Value,...
    'offset',str2num(handles.avgsub_offset_edit.String),...
    'filtermethod',filterstr,...
    'wvltfilter_std',str2num(handles.wvltfilter_std_edit.String),...
    'phasor_fit',handles.fitting_phasor_radio.Value,...
    'phasor_rad',str2num(handles.phasor_rad_edit.String),...
    'phasor_2d',handles.phasor2d_radiobtn.Value,...
    'phasor_astig',handles.phasorastig_radiobtn.Value,...
    'phasor_astigcal',handles.phasor_ast_calib.String{1,1},...
    'makeavgshifthist',handles.makeavgshifthist_chk.Value,...
    'avgshifthist_res',str2num(handles.avgshifthist_res_edit.String),...
    'avgshifthist_lateralsubpix',str2num(handles.avgshfthist_latsubpx_edit.String),...
    'avgshifthist_lateralshifts',str2num(handles.avgshifthist_latshift_edit.String),...
    'avgshifthist_axialcolnrs',str2num(handles.avgshfthist_axcolnr_edit.String),...
    'avgshifthist_pixelsize',str2num(handles.pxsize_edit.String),...
    'fullmovie',fullmovie,...
    'nrframes',str2num(handles.edit_nrframes.String),...
    'bpthrsh',str2num(handles.bpthrsh_edit.String),...
    'egdesz',str2num(handles.edgesz_edit.String),...
    'pctile_frame',handles.pctile_frame.Value,...
    'which_gaussian',gaussiantype,...
    'usegpu',handles.gaussGPU_checkbox.Value,...
    'MLE_fit',handles.gauss_mle_radio.Value,...
    'fit_ang',str2num(handles.fit_ang_edit.String),...
    'driftcorr_cc',handles.do_driftcorr_cc.Value,...
    'driftcorrcc_lateralsubpix',str2num(handles.edit_driftcorrcc_latsubpix.String),...
    'dirftcorrcc_temporalbins',str2num(handles.edit_driftcorrcc_tempbins.String));

% --- Executes on button press in bgsub_chk.
function bgsub_chk_Callback(hObject, eventdata, handles)
% hObject    handle to bgsub_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bgsub_chk


% --- Executes on button press in makeavgsubmovie_chk.
function makeavgsubmovie_chk_Callback(hObject, eventdata, handles)
% hObject    handle to makeavgsubmovie_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeavgsubmovie_chk


% --- Executes on button press in performguess_chk.
function performguess_chk_Callback(hObject, eventdata, handles)
% hObject    handle to performguess_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of performguess_chk


% --- Executes on button press in checkguess_chk.
function checkguess_chk_Callback(hObject, eventdata, handles)
% hObject    handle to checkguess_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkguess_chk


% --- Executes on button press in makeoffframelist_chk.
function makeoffframelist_chk_Callback(hObject, eventdata, handles)
% hObject    handle to makeoffframelist_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeoffframelist_chk


% --- Executes on button press in performfitting_chk.
function performfitting_chk_Callback(hObject, eventdata, handles)
% hObject    handle to performfitting_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
enable_disable_phasor_gauss(handles);

% Hint: get(hObject,'Value') returns toggle state of performfitting_chk


% --- Executes on button press in performtracking_chk.
function performtracking_chk_Callback(hObject, eventdata, handles)
% hObject    handle to performtracking_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of performtracking_chk


% --- Executes on button press in makeviewfitsmov_chk.
function makeviewfitsmov_chk_Callback(hObject, eventdata, handles)
% hObject    handle to makeviewfitsmov_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeviewfitsmov_chk



function filepath_Callback(hObject, eventdata, handles)
% hObject    handle to filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filepath as text
%        str2double(get(hObject,'String')) returns contents of filepath as a double


% --- Executes during object creation, after setting all properties.
function filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
    save('GUI\StoredSettings.mat','fvals_to_save')
catch
    fprintf('\n Could not store current vars!\n')
end

% Hint: delete(hObject) closes the figure
delete(hObject);



function avgsub_offset_edit_Callback(hObject, eventdata, handles)
% hObject    handle to avgsub_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgsub_offset_edit as text
%        str2double(get(hObject,'String')) returns contents of avgsub_offset_edit as a double


% --- Executes during object creation, after setting all properties.
function avgsub_offset_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgsub_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dfrlmsz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dfrlmsz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dfrlmsz_edit as text
%        str2double(get(hObject,'String')) returns contents of dfrlmsz_edit as a double


% --- Executes during object creation, after setting all properties.
function dfrlmsz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfrlmsz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgwin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to avgwin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgwin_edit as text
%        str2double(get(hObject,'String')) returns contents of avgwin_edit as a double


% --- Executes during object creation, after setting all properties.
function avgwin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgwin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function moloffwin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to moloffwin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of moloffwin_edit as text
%        str2double(get(hObject,'String')) returns contents of moloffwin_edit as a double


% --- Executes during object creation, after setting all properties.
function moloffwin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to moloffwin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


function phasorastig_radiobtn_Callback(hObject, eventdata, handles)
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


function phasor_rad_edit_Callback(hObject, eventdata, handles)
% hObject    handle to phasor_rad_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phasor_rad_edit as text
%        str2double(get(hObject,'String')) returns contents of phasor_rad_edit as a double


% --- Executes during object creation, after setting all properties.
function phasor_rad_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phasor_rad_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makeavgshifthist_chk.
function makeavgshifthist_chk_Callback(hObject, eventdata, handles)
% hObject    handle to makeavgshifthist_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeavgshifthist_chk


% --- Executes on button press in browse_btn.
function browse_btn_Callback(hObject, eventdata, handles)
% hObject    handle to browse_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fldrstrt = pwd;
if isfolder(handles.filepath.String{1,1})
    fldrstrt = handles.filepath.String{1,1};
else
    lastslash = max(strfind(handles.filepath.String{1,1},'\'));
    prevfolder = extractBefore(handles.filepath.String{1,1},lastslash+1);
    if isfolder(prevfolder)
        fldrstrt = prevfolder;
    end
end
curfldr = pwd;
cd(fldrstrt)
[file,path] = uigetfile('*.*');
cd(curfldr)
if ~isequal(file,0)
    handles.filepath.String{1,1} = [path file];
end



function avgshifthist_res_edit_Callback(hObject, eventdata, handles)
% hObject    handle to avgshifthist_res_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgshifthist_res_edit as text
%        str2double(get(hObject,'String')) returns contents of avgshifthist_res_edit as a double


% --- Executes during object creation, after setting all properties.
function avgshifthist_res_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgshifthist_res_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgshfthist_latsubpx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to avgshfthist_latsubpx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgshfthist_latsubpx_edit as text
%        str2double(get(hObject,'String')) returns contents of avgshfthist_latsubpx_edit as a double


% --- Executes during object creation, after setting all properties.
function avgshfthist_latsubpx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgshfthist_latsubpx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgshfthist_axcolnr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to avgshfthist_axcolnr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgshfthist_axcolnr_edit as text
%        str2double(get(hObject,'String')) returns contents of avgshfthist_axcolnr_edit as a double


% --- Executes during object creation, after setting all properties.
function avgshfthist_axcolnr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgshfthist_axcolnr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgshifthist_latshift_edit_Callback(hObject, eventdata, handles)
% hObject    handle to avgshifthist_latshift_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgshifthist_latshift_edit as text
%        str2double(get(hObject,'String')) returns contents of avgshifthist_latshift_edit as a double


% --- Executes during object creation, after setting all properties.
function avgshifthist_latshift_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgshifthist_latshift_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wvltfilter_std_edit_Callback(hObject, eventdata, handles)
% hObject    handle to wvltfilter_std_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wvltfilter_std_edit as text
%        str2double(get(hObject,'String')) returns contents of wvltfilter_std_edit as a double


% --- Executes during object creation, after setting all properties.
function wvltfilter_std_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wvltfilter_std_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in make_ASH.
function make_ASH_Callback(hObject, eventdata, handles)
% hObject    handle to make_ASH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get fits 
fldrstrt = pwd;
if isfolder(handles.filepath.String{1,1})
    fldrstrt = handles.filepath.String{1,1};
else
    lastslash = max(strfind(handles.filepath.String{1,1},'\'));
    prevfolder = extractBefore(handles.filepath.String{1,1},lastslash+1);
    if isfolder(prevfolder)
        fldrstrt = prevfolder;
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
p.avgshifthist_axialcolnrs = str2num(handles.avgshfthist_axcolnr_edit.String);
p.avgshifthist_pixelsize = str2num(handles.pxsize_edit.String);
%Run code
disp('Making average shifted histogram...')
avgshifthist(fits,p,movsz,[pathASH fileASH],1)
disp('Completed average shifted histogram...')



function phasor_ast_calib_Callback(hObject, eventdata, handles)
% hObject    handle to phasor_ast_calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phasor_ast_calib as text
%        str2double(get(hObject,'String')) returns contents of phasor_ast_calib as a double


% --- Executes during object creation, after setting all properties.
function phasor_ast_calib_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phasor_ast_calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pxsize_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pxsize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pxsize_edit as text
%        str2double(get(hObject,'String')) returns contents of pxsize_edit as a double

% Update diffractionlimited spot size here upon updating
% Same code as in emm_wavelength_edit_Callback, objNA_edit, and
% pxsize_edit_Callback
AiryDiscDiameter_pxunits = round((1.15*2*(1.22*str2num(handles.emm_wavelength_edit.String))/...
    (2*str2num(handles.objNA_edit.String)))/str2num(handles.pxsize_edit.String),2);
handles.dfrlmsz_edit.String = num2str(AiryDiscDiameter_pxunits);

% --- Executes during object creation, after setting all properties.
function pxsize_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pxsize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function emm_wavelength_edit_Callback(hObject, eventdata, handles)
% hObject    handle to emm_wavelength_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of emm_wavelength_edit as text
%        str2double(get(hObject,'String')) returns contents of emm_wavelength_edit as a double
% Update diffractionlimited spot size here upon updating
% Same code as in emm_wavelength_edit_Callback, objNA_edit, and
% pxsize_edit_Callback
AiryDiscDiameter_pxunits = round((1.15*2*(1.22*str2num(handles.emm_wavelength_edit.String))/...
    (2*str2num(handles.objNA_edit.String)))/str2num(handles.pxsize_edit.String),2);
handles.dfrlmsz_edit.String = num2str(AiryDiscDiameter_pxunits);


% --- Executes during object creation, after setting all properties.
function emm_wavelength_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to emm_wavelength_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function objNA_edit_Callback(hObject, eventdata, handles)
% hObject    handle to objNA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of objNA_edit as text
%        str2double(get(hObject,'String')) returns contents of objNA_edit as a double
% Update diffractionlimited spot size here upon updating
% Same code as in emm_wavelength_edit_Callback, objNA_edit, and
% pxsize_edit_Callback
AiryDiscDiameter_pxunits = round((1.15*2*(1.22*str2num(handles.emm_wavelength_edit.String))/...
    (2*str2num(handles.objNA_edit.String)))/str2num(handles.pxsize_edit.String),2);
handles.dfrlmsz_edit.String = num2str(AiryDiscDiameter_pxunits);


% --- Executes during object creation, after setting all properties.
function objNA_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to objNA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nrframes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nrframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nrframes as text
%        str2double(get(hObject,'String')) returns contents of edit_nrframes as a double


% --- Executes during object creation, after setting all properties.
function edit_nrframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nrframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bpthrsh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bpthrsh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bpthrsh_edit as text
%        str2double(get(hObject,'String')) returns contents of bpthrsh_edit as a double


% --- Executes during object creation, after setting all properties.
function bpthrsh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bpthrsh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edgesz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to edgesz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edgesz_edit as text
%        str2double(get(hObject,'String')) returns contents of edgesz_edit as a double


% --- Executes during object creation, after setting all properties.
function edgesz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgesz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pctile_frame.
function pctile_frame_Callback(hObject, eventdata, handles)
% hObject    handle to pctile_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pctile_frame


% --- Executes on button press in gaussGPU_checkbox.
function gaussGPU_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to gaussGPU_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gaussGPU_checkbox



function fit_ang_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fit_ang_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_ang_edit as text
%        str2double(get(hObject,'String')) returns contents of fit_ang_edit as a double


% --- Executes during object creation, after setting all properties.
function fit_ang_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_ang_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fitting_gauss_radio.
function fitting_gauss_radio_Callback(hObject, eventdata, handles)
% hObject    handle to fitting_gauss_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fitting_gauss_radio
enable_disable_phasor_gauss(handles);


% --- Executes on button press in fitting_gauss_radio.
function fitting_phasor_radio_Callback(hObject, eventdata, handles)
% hObject    handle to fitting_gauss_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fitting_gauss_radio
enable_disable_phasor_gauss(handles);

function enable_disable_phasor_gauss(handles)

genarronoff{1} = 'off';
genarronoff{2} = 'on';
if handles.performfitting_chk.Value
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


handles.gauss_mle_radio.Enable = gaussen;
handles.gauss_ls_radio.Enable = gaussen;
handles.gaussGPU_checkbox.Enable = gaussen;
handles.uibuttongroup7.ForegroundColor = colGauss;
handles.uibuttongroup8.ForegroundColor = colGauss;
handles.which_gauss_1_radio.Enable = gaussen;
handles.which_gauss_2_radio.Enable = gaussen;
handles.which_gauss_3_radio.Enable = gaussen;
handles.fit_ang_edit.Enable = genarronoff{1+(handles.performfitting_chk.Value*handles.fitting_gauss_radio.Value*handles.which_gauss_3_radio.Value)};
handles.text33.Enable = genarronoff{1+(handles.performfitting_chk.Value*handles.fitting_gauss_radio.Value*handles.which_gauss_3_radio.Value)};

handles.uibuttongroup5.ForegroundColor = colPhasor;
handles.text9.Enable = phasoren;
handles.phasor_rad_edit.Enable = phasoren;
handles.phasor2d_radiobtn.Enable = phasoren;
handles.phasorastig_radiobtn.Enable = phasoren;
handles.phasor_ast_calib.Enable = genarronoff{1+(handles.performfitting_chk.Value*handles.fitting_phasor_radio.Value*handles.phasorastig_radiobtn.Value)};


handles.fitting_phasor_radio.Enable = genarronoff{1+handles.performfitting_chk.Value};
handles.fitting_gauss_radio.Enable = genarronoff{1+handles.performfitting_chk.Value};


function enable_disable_filtering(handles)

genarronoff{1} = 'off';
genarronoff{2} = 'on';

handles.text21.Enable = genarronoff{1+handles.filter_wvlt_radio.Value};
handles.wvltfilter_std_edit.Enable = genarronoff{1+handles.filter_wvlt_radio.Value};


handles.pctile_frame.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.text30.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.text31.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.text32.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.bpthrsh_edit.Enable = genarronoff{1+handles.filter_bpass_radio.Value};
handles.edgesz_edit.Enable = genarronoff{1+handles.filter_bpass_radio.Value};


% --- Executes on button press in do_driftcorr_cc.
function do_driftcorr_cc_Callback(hObject, eventdata, handles)
% hObject    handle to do_driftcorr_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_driftcorr_cc



function edit_driftcorrcc_latsubpix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_driftcorrcc_latsubpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_driftcorrcc_latsubpix as text
%        str2double(get(hObject,'String')) returns contents of edit_driftcorrcc_latsubpix as a double


% --- Executes during object creation, after setting all properties.
function edit_driftcorrcc_latsubpix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_driftcorrcc_latsubpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_driftcorrcc_tempbins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_driftcorrcc_tempbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_driftcorrcc_tempbins as text
%        str2double(get(hObject,'String')) returns contents of edit_driftcorrcc_tempbins as a double


% --- Executes during object creation, after setting all properties.
function edit_driftcorrcc_tempbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_driftcorrcc_tempbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
