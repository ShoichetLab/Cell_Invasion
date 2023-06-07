% This is written by Luke, (TM)
function varargout = App(varargin)
% APP MATLAB code for App.fig
%      APP, by itself, creates a new APP or raises the existing
%      singleton*.
%
%      H = APP returns the handle to a new APP or the handle to
%      the existing singleton*.
%
%      APP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APP.M with the given input arguments.
%
%      APP('Property','Value',   ...) creates a new APP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before App_OpeningFcn gets called.  An
%      unrecognized property name or    invalid value makes property application
%      stop.  All inputs are passed to App_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help App

% Last Modified by GUIDE v2.5 03-Sep-2018 10:07:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @App_OpeningFcn, ...
                   'gui_OutputFcn',  @App_OutputFcn, ...
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


% --- Executes just before App is made visible.
function App_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for App
handles.output = hObject;
handles.FilePath = '';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes App wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = App_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
saveDirectory = get(handles.ImageSaveDirectory, 'String');

exitFunction = VerifyFolder(handles); % verify if all data file folders have been properly filled in
if exitFunction % if VerifyFolder returns true, then exit the 'RunButton_Callback' function
    return
end

oldpointer = get(handles.figure1, 'pointer'); 
set(handles.figure1, 'pointer', 'watch') % change the mouse pointer to a spinning watch symbol when program is running
drawnow;

fileNum = 0;

if get(handles.SurfaceDataRadio, 'Value') == 1.0
    surfaceFilePath = get(handles.SurfaceFilePath, 'String');
    [surfaceArr, errMsg] = ExcelFileProcessing(surfaceFilePath, 'surface file');
    if ~isempty(errMsg) 
        DisplayError(errMsg); % display an error message when surface files aren't arranged properly
        set(handles.figure1, 'pointer', oldpointer);  
        return 
    end
    fileNum = size(surfaceArr, 1);
    set(handles.SurfaceDataToggle, 'Enable', 'on');
    set(handles.CutoffToggle, 'Enable', 'on');
    set(handles.InvasionToggle, 'Enable', 'on');
    set(handles.NonInvadingToggle, 'Enable', 'on');
else % toggle off the following features when surface labeling is not available
    set(handles.SurfaceDataToggle, 'Enable', 'off');
    set(handles.CutoffToggle, 'Enable', 'off');
    set(handles.InvasionToggle, 'Enable', 'off');
    set(handles.NonInvadingToggle, 'Enable', 'off');
end
if get(handles.CellMarkerRadio, 'Value') == 1.0
    cellMarkerFilePath = get(handles.CellMarkerFilePath, 'String');
    [cellMarkerArr, errMsg] = ExcelFileProcessing(cellMarkerFilePath, 'cell marker file');
    if ~isempty(errMsg)
        DisplayError(errMsg);
        set(handles.figure1, 'pointer', oldpointer)
        return
    end
end
if get(handles.InvasionRadio, 'Value') == 1.0
    invasionFilePath = get(handles.InvasionFilePath, 'String');
    [invasionArr, errMsg] = ExcelFileProcessing(invasionFilePath, 'invasion file');
    if ~isempty(errMsg)
        DisplayError(errMsg);
        set(handles.figure1, 'pointer', oldpointer);
        return
    end
    fileNum = size(invasionArr, 1);
end

%% Initialize the arrays
fileName = cell(fileNum + 1, 1);
fileName{1} = 'File Name';
percentInvasion = zeros(fileNum, 1);
avgInvasionDepth = zeros(fileNum, 1);
cellNum = zeros(fileNum, 1);
cutoffThreshold = zeros(fileNum, 1);
invExpPercentage = zeros(fileNum, 1);
nonInvExpPercentage = zeros(fileNum, 1);
Pmax = zeros(fileNum, 1);
Pmin = zeros(fileNum, 1);

rotate3d on
set(handles.NextImageButton, 'Enable', 'on');
set(handles.ManualCutoffRadio, 'Enable', 'off')
set(hObject, 'Enable', 'off');
set(handles.figure1, 'pointer', oldpointer)

%% Only runs when invasion files are not available, and user wants to quantify the surface
if get(handles.InvasionRadio, 'Value') == 0.0
    for i = 1: size(surfaceArr, 1)
        cla('reset');
        surfaceCoord = surfaceArr{i, 1};
        threshold = abs(str2double(get(handles.ThresholdValue, 'String')));
        [f, Pmax_min, Pmin_min] = SurfaceLabeling(surfaceCoord, threshold, handles);
        
        title(surfaceArr{i, 2}, 'Interpreter', 'none');
        if get(handles.SaveRadio, 'Value') == 1.0
            saveFileName = strcat(get(get(gca, 'title'), 'string'), '.fig');
            fullFileName = fullfile(saveDirectory, saveFileName);
            saveas(gcf, fullFileName);
        end
        uiwait;
    end
    return
end

%% When invasion files are available ....
for i = 1: fileNum
    cla('reset');
    fileName{i + 1} = invasionArr{i, 2};
    cellCoord = invasionArr{i, 1};        
    xCellCoord = cellCoord(:, 1);
    yCellCoord = cellCoord(:, 2);
    zCellCoord = cellCoord(:, 3);
    cellNum(i) = size(cellCoord, 1);

    threshold = abs(str2double(get(handles.ThresholdValue, 'String')));

    %% Process the invasion data
    if get(handles.SurfaceDataRadio, 'Value') == 1.0 % if surface labeling file available
        surfaceCoord = surfaceArr{i, 1};
        [f, Pmax_min, Pmin_min] = SurfaceLabeling(surfaceCoord, threshold, handles);
        zSurface = f(xCellCoord, yCellCoord); % calculate the z-position on the gel surface corresponding to the x,y coordinate of the cell
        z_cutoff = zSurface - threshold;
        mask = zCellCoord < z_cutoff;  % all cells that satisfy the mask condition is classified as an invading cell
        invasionDepth = abs(zCellCoord - zSurface); % invasion depth is corrected to a positive value
        invasionDepth = invasionDepth(mask);
        Pmax(i) = Pmax_min;
        Pmin(i) = Pmin_min;
    else % if surface labeling not available, plot cells in 2D            
        offset = max(zCellCoord); 
        zCellCoord = zCellCoord - offset; % zero the top cell
        scatter(xCellCoord, zCellCoord, 12, 'o', 'filled', 'MarkerFaceColor',[0 0 .8]);
        hold on
        plot([min(xCellCoord); max(xCellCoord)], [-threshold; -threshold]);

        if get(handles.ManualCutoffRadio, 'Value') == 1.0 % if user wants to adjust cutoffs for each figure
            uiwait;
            threshold = abs(str2double(get(handles.ThresholdValue, 'String'))); % update threshold value once confirmed
        end

        mask = zCellCoord < -threshold; % reminder: threshold is a positive value, invading when z_position of the cell goes beyond the threshold
        invasionDepth = abs(zCellCoord);
        invasionDepth = invasionDepth(mask); 
    end

    invadingCells = [xCellCoord(mask), yCellCoord(mask), zCellCoord(mask)];
    nonInvadingCells = [xCellCoord(~mask), yCellCoord(~mask), zCellCoord(~mask)];
    avgInvasionDepth(i) = mean(invasionDepth);
    percentInvasion(i) = size(zCellCoord(mask), 1)/ cellNum(i) * 100;
    cutoffThreshold(i) = threshold;

    %% Process the cell marker labeling files
    if get(handles.CellMarkerRadio, 'Value') == 1.0 % if cell marker file available
        markerCoord = cellMarkerArr{i, 1};
        xMarkerCoord = markerCoord(:, 1);
        yMarkerCoord = markerCoord(:, 2);
        zMarkerCoord = markerCoord(:, 3);
        
        if get(handles.SurfaceDataRadio, 'Value') == 1.0 % If surface labeling file is available
            zSurface_exp = f(xMarkerCoord, yMarkerCoord); % calculate the z-position on the gel surface corresponding to the x,y coordinate of the EXPRESSING cell
            zCutoff_exp = zSurface_exp - threshold;
            mask_exp = zMarkerCoord < zCutoff_exp;  % when the cell's coordinate satisfies the 'mask_exp' condition, then it's an invading cell 
            invadingCells_exp = [xMarkerCoord(mask_exp), yMarkerCoord(mask_exp), zMarkerCoord(mask_exp)];
            nonInvadingCells_exp = [xMarkerCoord(~mask_exp), yMarkerCoord(~mask_exp), zMarkerCoord(~mask_exp)];
        else % when surface labeling is not available, plot in 2D
            zMarkerCoord = zMarkerCoord - offset; % shift the marker cells 
            mask_exp = zMarkerCoord < -threshold; % reminder: threshold is a positive value, invading when z_position of the cell is lesser than -threshold
            invadingCells_exp = markerCoord(mask_exp);
            nonInvadingCells_exp = markerCoord(~mask_exp);
        end

        %invadingCells_exp = intersect(markerCoord, invadingCells, 'rows'); % _exp stands for expressing
        %nonInvadingCells_exp= intersect(markerCoord, nonInvadingCells, 'rows');
        invExpPercentage(i) = size(invadingCells_exp, 1) / size(invadingCells, 1) * 100;
        nonInvExpPercentage(i) = size(nonInvadingCells_exp, 1) / size(nonInvadingCells, 1) * 100;
    end

    %% Plotting section
    if get(handles.SurfaceDataRadio, 'Value') == 1.0 % plot invading and non-invading cells in 3D if surface labeling is available
        scatter3(handles.axes, invadingCells(:, 1), invadingCells(:, 2), invadingCells(:, 3)...
            ,12, 'd', 'filled', 'MarkerFaceColor', [0 0 0.8]);
        scatter3(handles.axes, nonInvadingCells(:, 1), nonInvadingCells(:, 2), nonInvadingCells(:, 3)...
            ,12 ,'s', 'filled', 'MarkerFaceColor', [0 0 1]);

        if get(handles.CellMarkerRadio, 'Value') == 1.0 % Highlight the expressing cells in green
            expressingCells = [markerCoord];
            h_expression = scatter3(handles.axes, expressingCells(:, 1), expressingCells(:, 2), expressingCells(:, 3)...
                ,15, 'd', 'filled', 'MarkerFaceColor', [0 0.5 0]);
            uistack(h_expression, 'bottom');
        end 
    else
        if get(handles.CellMarkerRadio, 'Value') == 1.0
            h_expression = scatter(xMarkerCoord, zMarkerCoord, 18, 'o', 'filled', 'MarkerFaceColor',[0.6 0.6 0]);
            uistack(h_expression, 'bottom');                  
        end
    end

    title(fileName{i+1}, 'Interpreter', 'none');
    axis(axis) %fix the axis value to make sure the scale of figure doesn't change when hiding or showing datapoints

    %% Update data on the GUI
    set(handles.CellNumberResult, 'String', num2str(cellNum(i)));
    set(handles.InvasionPercentageResult, 'String', num2str(percentInvasion(i)));
    set(handles.InvasionDepthResult, 'String', num2str(avgInvasionDepth(i)));
    set(handles.CutoffThresholdResult, 'String', num2str(threshold));
    set(handles.nInvExpResult, 'String', num2str(nonInvExpPercentage(i)));
    set(handles.InvExpResult, 'String', num2str(invExpPercentage(i)));
    
    % Save the figure if requested
    if get(handles.SaveRadio, 'Value') == 1.0 
        SaveFigure(handles, saveDirectory);
    end
    uiwait;
end

[resultMatrix] = StoreOutputData(handles, fileName, percentInvasion, avgInvasionDepth, cellNum, cutoffThreshold, Pmin, Pmax, invExpPercentage, nonInvExpPercentage);
handles.myVal = resultMatrix;
guidata(hObject, handles);

set(hObject, 'Enable', 'on');
set(handles.ManualCutoffRadio, 'Enable', 'on')
set(handles.NextImageButton, 'Enable', 'off');
set(handles.SaveButton, 'Enable', 'on');

function NextImageButton_Callback(~, ~, ~)
uiresume;

function [exitFunction] = VerifyFolder(handles)
% Verifies if all data directories on the GUI panel are valid 
exitFunction = false;
if get(handles.SurfaceDataRadio, 'Value') == 1.0
    if exist(get(handles.SurfaceFilePath, 'String'), 'dir') ~= 7
        DisplayError('Invalid Surface File Path')
        exitFunction = true;
    end
end
if get(handles.CellMarkerRadio, 'Value') == 1.0
    if exist(get(handles.CellMarkerFilePath, 'String'), 'dir') ~= 7
        DisplayError('Invalid Cell Marker File Path')
        exitFunction = true;
    end
end
if get(handles.InvasionRadio, 'Value') == 1.0
    if exist(get(handles.InvasionFilePath, 'String'), 'dir') ~= 7
        DisplayError('Invalid Invasion File Path');
        exitFunction = true;
    end
end
if get(handles.SaveRadio, 'Value') == 1.0
    if exist(get(handles.ImageSaveDirectory, 'String'), 'dir') ~= 7
        DisplayError('Invalid Image Path');
        exitFunction = true;
    end
end

function SaveFigure(handles, saveDirectory)
% Save the figure as a standalone '.fig' file, instead of saving the entire GUI interface
% https://www.mathworks.com/matlabcentral/answers/352158-how-to-save-an-axes-within-a-gui-to-a-fig     
fignew = figure('Visible', 'off'); % Invisible figure
newAxes = copyobj(handles.axes, fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
saveFileName = strcat(get(get(gca, 'title'), 'string'), '.fig');
fullFileName = fullfile(saveDirectory, saveFileName);
saveas(fignew, fullFileName);
delete(fignew);

function DisplayError(msg)
errorMsg = errordlg(msg, 'Error');
set(errorMsg, 'WindowStyle', 'modal');
uiwait(errorMsg)

function [finalArr, errMsg] = ExcelFileProcessing(dir_name, fileType)
% Returns an 'errMsg' when the folders are not properly arranged in the
% proper hierarchy 
    batchFile = dir(fullfile(dir_name, '*.xls'));
    if ~isempty(batchFile)
        if size(batchFile, 1) ~= 1 % if the files are not organized correctly (there should only be one batch file in the folder, all individual files should be organized in a subfolder)
           finalArr = {};
           errMsg = ['Files organized incorrectly for ', fileType , ' Note: there should only be one batch file in the folder, and all individual files should be organized in subfolders'];
           return
        else
            batchFilePath = fullfile(batchFile.folder, batchFile.name);
            [batchArr] = BatchFileProcessing(batchFilePath);
        end
    end

    d = dir(dir_name); % https://stackoverflow.com/questions/8748976/list-the-subfolders-in-a-folder-matlab-only-subfolders-not-files/8749121
    isub = [d(:).isdir]; 
    dirList = {d(isub).name}'; % lists out all subdirectories in the folder
    dirList (ismember(dirList,{'.','..'})) = []; % remove '.' and '..'
    if ~isempty(dirList) % process subdirectory files if subdirectory exists
        individualFiles = dir(fullfile(dir_name, dirList{1}, '*.xls'));
        length = size(individualFiles, 1);
        manualArr = cell(length, 2);
        for i = 1: length
           individualFileName = individualFiles(i).name;
           [~, name, ~] = fileparts(individualFileName);

           individualFilePath = fullfile(dir_name, dirList{1}, individualFileName);
           data = xlsread(individualFilePath, 'Position', 'A:C');
           manualArr{i, 1} = data;
           manualArr{i, 2} = name;
        end
    end

    if exist('manualArr', 'var') && exist('batchArr', 'var')
        finalArr = sortrows([batchArr; manualArr], 2); %sort the file based on the file name
        errMsg = '';
    elseif exist('manualArr', 'var') && ~ exist('batchArr', 'var')
        finalArr = sortrows(manualArr, 2);
        errMsg = '';
    elseif ~exist('manualArr', 'var') && exist('batchArr', 'var')
        finalArr = sortrows(batchArr, 2);
        errMsg = '';
    else
        finalArr = {};
        errMsg = ['No valid files in the ', fileType, ' path']; 
    end

function [batchArr] = BatchFileProcessing(filePath) % called by the ExcelFileProcessing function
    [numericData, textData, ~] = xlsread(filePath, 'Position', 'A:M'); % numericData Contains all the numeric data of invasion, ie. x,y,z coordinate of cells

    coordinate = numericData(:, 1:3);
    fileNames = textData(3:end, 13);

    collection = cell(size(coordinate, 1), 4);
    collection(:, 1:3) = num2cell(coordinate);
    collection(:, 4) = fileNames;

    [~,~,X1] = unique(collection(:,4));
    partitionedArr = accumarray(X1, 1:size(collection,1), [], @(row){collection(row,:)});

    batchArr = cell(size(partitionedArr, 1), 2);
    for idx = 1: size(partitionedArr, 1) % reformat the partition array
        batchArr{idx, 1} = cell2mat(partitionedArr{idx}(:, 1:3));
        batchArr{idx, 2} = cell2mat(partitionedArr{idx}(1, 4));
    end 

function [f, Pmax_min, Pmin_min] = SurfaceLabeling(surfaceCoordinates, thresholdCutoff, handles)
x = surfaceCoordinates(:, 1);
y = surfaceCoordinates(:, 2);
z = surfaceCoordinates(:, 3);
    
%     h_originalDataSet = scatter3(handles.axes, x, y, z);
%     hold(handles.axes, 'on')
    
%     exclusion = questdlg('Do you want to exclude certain data?'...
%                 ,'Data Exclusion'...
%                 ,'Yes', 'No', 'No');
%     if strcmp(exclusion, 'Yes')
%         prompt = 'Please enter exclusion criteria, z < ';
%         promptTitle = 'Data Exclusion';
%         definput = {'-800'};
%         exclusionCriteria = inputdlg(prompt, promptTitle, [1, 40], definput);
%         exclusionCriteria =  z < -abs(str2double(exclusionCriteria{1}));
%         
%         f = fit([x, y], z, 'lowess', 'Exclude', exclusionCriteria);
%         f2 = fit([x, y], z - thresholdCutoff, 'lowess', 'Exclude', exclusionCriteria);
%         h_surface = plot(f, [x(~exclusionCriteria),y(~exclusionCriteria)], z(~exclusionCriteria));
%     else
        f = fit([x, y], z, 'poly22');
        f2 = fit([x, y], z - thresholdCutoff, 'poly22');
        
%     end
    
    %delete(h_originalDataSet);
    
h_surface = plot(f, [x,y], z);
hold on
h_surface(1).FaceAlpha = 0.5;
h_surface(1).EdgeColor = 'none';
h_surface(2).MarkerFaceColor = [0.8 0 0];
h_surface(2).MarkerEdgeColor = 'none';
h_surface(2).MarkerSize = 3;

h_cutoffPlane = plot(f2);
h_cutoffPlane.FaceAlpha = 0.5;
h_cutoffPlane.EdgeColor = 'none';
h_cutoffPlane.FaceColor = [0, 0, 0];

[Pmax_min, Pmin_min] = CalculateCurvature(f, x, y, handles);
    
function [Pmax_min, Pmin_min] = CalculateCurvature(f, x, y, handles)
[xmin, ~] = min(x);
[xmax, ~]= max(x);
[ymin, ~] = min(y);
[ymax, ~] = max(y);
xmean = (xmax + xmin) / 2;
ymean = (ymax + ymin) / 2;
xdelta = xmean - xmin;
ydelta = ymean - ymin;

[xmean, ymean, ~] = FindMin(f, xmean, ymean, xdelta, ydelta);

scatter3(handles.axes, xmean, ymean, f(xmean, ymean), 36, 'o', 'filled', 'MarkerFaceColor', [1 0 0]);
x_curvature = linspace(xmean - 20, xmean + 20, 401);
y_curvature = linspace(ymean - 20, ymean + 20, 401);
[X, Y] = meshgrid(x_curvature, y_curvature);
Z = f(X, Y);
[~, ~, Pmax, Pmin] = surfature(X, Y, Z);
Pmax_min = Pmax(201,201); % The value for the principle curvatures at the minimum point of meniscus will be at the center of the Pmax, and Pmin matrix
Pmin_min = Pmin(201,201);

% Utility function used to find the minimum point of the fitted surface,
% called by 'CalculateCurvature()' 
function [x_mean, y_mean, x_delta] = FindMin(func, x_mean, y_mean, x_delta, y_delta) 
xmin = x_mean - x_delta;
xmax = x_mean + x_delta;
ymin = y_mean - y_delta;
ymax = y_mean + y_delta;
x = linspace(xmin, xmax, 700);
y = linspace(ymin, ymax, 700);
[X, Y] = meshgrid(x, y);
Z = func(X, Y);
[row, col] = find(Z == min(Z(:)));
x_mean = X(row(1), col(1));
y_mean = Y(row(1), col(1));
x_delta = (xmax - xmin) / 700;
y_delta = (ymax - ymin) / 700;
while(x_delta > 0.001)
    [x_mean, y_mean, x_delta] = FindMin(func, x_mean, y_mean, x_delta, y_delta);
end

function [resultMatrix] = StoreOutputData(handles, fileName, percentInvasion, avgInvasionDepth, cellNum, cutoffThreshold, Pmin, Pmax, invExpPercentage, nonInvExpPercentage)
percentInvasionArr = [{'Percent Invasion'}; num2cell(percentInvasion)];
avgInvasionDepthArr= [{'Avg Invasion Depth'}; num2cell(avgInvasionDepth)];
cellNumArr = [{'Cell Number'}; num2cell(cellNum)];
cutoffThresholdArr = [{'Cutoff Threshold'}; num2cell(cutoffThreshold)];

resultMatrix = [fileName, percentInvasionArr, avgInvasionDepthArr, cellNumArr, cutoffThresholdArr];
if get(handles.SurfaceDataRadio, 'Value') == 1.0 % If surface labelling was available, also output data about the minimum and maximum curvature
    PminArr = [{'Pmin'}; num2cell(Pmin)];
    PmaxArr = [{'Pmax'}; num2cell(Pmax)];
    resultMatrix = [resultMatrix, PminArr, PmaxArr];
end
if get(handles.CellMarkerRadio, 'Value') == 1.0 % If cell marker labelling was available, output the expressing percentage of given cell markers 
    invExpArr = [{'Expressing Percentage for Invading Cells'}; num2cell(invExpPercentage)];
    nonInvExpArr = [{'Expressing Percentage for Non-Invading Cells'}; num2cell(nonInvExpPercentage)];
    resultMatrix = [resultMatrix, invExpArr, nonInvExpArr];
end

function SaveButton_Callback(~, ~, handles) % Saves the data file when user clicks on the 'Save' button
[file, path]  = uiputfile('Analysis Result.xlsx');
if ~isempty(file) && ~isempty(path)
    writeFileName = fullfile(path, file);   
    resultMatrix = handles.myVal;
    xlswrite(writeFileName, resultMatrix);
end

%% Callback functions for the toggle buttons for the graphics
function SurfaceDataToggle_Callback(hObject, ~, handles)
if get(handles.SurfaceDataRadio, 'Value') == 1
    switch hObject.String
        case 'Hide Surface Data'
            handles.axes.Children(5).Visible = 'off';
            set(hObject, 'String', 'Show Surface Data');
        case 'Show Surface Data'
            handles.axes.Children(5).Visible = 'on';
            set(hObject, 'String', 'Hide Surface Data');
    end
end

function CutoffToggle_Callback(hObject, ~, handles)
if get(handles.SurfaceDataRadio, 'Value') == 1
    switch hObject.String
        case 'Hide Cutoff Plane'
            handles.axes.Children(4).Visible = 'off';
            set(hObject, 'String', 'Show Cutoff Plane');
        case 'Show Cutoff Plane'
            handles.axes.Children(4).Visible = 'on';
            set(hObject, 'String', 'Hide Cutoff Plane');
    end
end

function InvasionToggle_Callback(hObject, ~, handles)
if get(handles.SurfaceDataRadio, 'Value') == 1
    switch hObject.String
        case 'Hide Invading Cells'
            handles.axes.Children(2).Visible = 'off';
            set(hObject, 'String', 'Show Invading Cells');
        case 'Show Invading Cells'
            handles.axes.Children(2).Visible = 'on';
            set(hObject, 'String', 'Hide Invading Cells');
    end
end

function NonInvadingToggle_Callback(hObject, ~, handles)
if get(handles.SurfaceDataRadio, 'Value') == 1
    switch hObject.String
        case 'Hide Non-Invading Cells'
            handles.axes.Children(1).Visible = 'off';
            set(hObject, 'String', 'Show Invading Cells');
        case 'Show Invading Cells'
            handles.axes.Children(1).Visible = 'on';
            set(hObject, 'String', 'Hide Non-Invading Cells');
    end
end

function ExpressionToggle_Callback(hObject, eventdata, handles)
if get(handles.CellMarkerRadio, 'Value') == 1
    switch hObject.String
        case 'Hide Expressing Cells'
            handles.axes.Children(end).Visible = 'off';
            set(hObject, 'String', 'Show Expressing Cells');
        case 'Show Expressing Cells'
            handles.axes.Children(end).Visible = 'on';
            set(hObject, 'String', 'Hide Expressing Cells');
    end
end


%% Callback functions for File Selection
    %% Surface File Buttons
function SurfaceDataRadio_Callback(hObject, ~, handles)
switch get(hObject, 'Value')
    case 0.0
        set(handles.SelectSurfaceFile, 'Enable', 'off');
        set(handles.SurfaceFilePath, 'Enable', 'off');
    case 1.0
        set(handles.SelectSurfaceFile, 'Enable', 'on');
        set(handles.SurfaceFilePath, 'Enable', 'on');
end

function SelectSurfaceFile_Callback(hObject, ~, handles)
guidata(hObject);
surfaceFilePath = uigetdir(handles.FilePath);
if surfaceFilePath ~= 0
    if exist(surfaceFilePath, 'dir') == 7
        set(handles.SurfaceFilePath, 'String', surfaceFilePath);
        handles.FilePath = surfaceFilePath;
        guidata(hObject, handles)
    end
end

    %% Cell Marker File Buttons
function CellMarkerRadio_Callback(hObject, ~, handles)
switch get(hObject, 'Value')
    case 0.0
        set(handles.SelectCellMarkerFile, 'Enable', 'off');
        set(handles.CellMarkerFilePath, 'Enable', 'off');
    case 1.0
        set(handles.SelectCellMarkerFile, 'Enable', 'on');
        set(handles.CellMarkerFilePath, 'Enable', 'on');
end

function SelectCellMarkerFile_Callback(hObject, ~, handles)
guidata(hObject);
cellMarkerFilePath = uigetdir(handles.FilePath);
if cellMarkerFilePath ~= 0 
    if exist(cellMarkerFilePath, 'dir') == 7 
        set(handles.CellMarkerFilePath, 'String', cellMarkerFilePath);
        handles.FilePath = cellMarkerFilePath;
        guidata(hObject, handles)
    end
end

    %% Invasion File Buttons
function InvasionRadio_Callback(hObject, eventdata, handles)
switch get(hObject, 'Value')
    case 0.0
        set(handles.SelectInvasionFile, 'Enable', 'off');
        set(handles.InvasionFilePath, 'Enable', 'off');
    case 1.0
        set(handles.SelectInvasionFile, 'Enable', 'on');
        set(handles.InvasionFilePath, 'Enable', 'on');
end
    
function SelectInvasionFile_Callback(hObject, ~, handles)
guidata(hObject);
invasionFilePath = uigetdir(handles.FilePath);
if invasionFilePath ~= 0
    if exist(invasionFilePath, 'dir') == 7
        set(handles.InvasionFilePath, 'String', invasionFilePath);
        handles.FilePath = invasionFilePath;
        guidata(hObject, handles)
    end
end

    %% Image Saving Buttons
function SaveRadio_Callback(hObject, ~, handles)
switch hObject.Value
    case 1.0
        set(handles.SaveImageButton, 'Enable', 'on');
        set(handles.ImageSaveDirectory, 'Enable', 'on');
    case 0.0
        set(handles.SaveImageButton, 'Enable', 'off');
        set(handles.ImageSaveDirectory, 'Enable', 'off');
end

function SaveImageButton_Callback(hObject, ~, handles)
guidata(hObject);
imageDirectory = uigetdir(handles.FilePath);
if imageDirectory ~= 0 
    if exist(imageDirectory, 'dir') == 7
        set(handles.ImageSaveDirectory, 'String', imageDirectory);
        handles.FilePath = imageDirectory;
        guidata(hObject, handles)
    end
end

function ManualCutoffRadio_Callback(hObject, ~, handles)
switch get(hObject, 'Value') % if manualCutoff radio button is disabled
    case 0.0
        set(handles.ThresholdValue, 'Enable', 'off');
    case 1.0
        set(handles.ThresholdValue, 'Enable', 'on');
end

function ThresholdValue_Callback(hObject, ~, handles)
if get(handles.ManualCutoffRadio, 'Value') == 1.0 && get(handles.SurfaceDataRadio, 'Value') == 0.0
    if ~isempty(handles.axes.Children) 
        if strcmp(handles.axes.Children(1).Type, 'line')
            cutoffLine = handles.axes.Children(1);
            threshold = -abs(str2double(hObject.String));
            cutoffLine.YData = [threshold threshold];
        end
    end
end

%% Matlab Auto-generated functions
function SurfaceFilePath_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CellMarkerFilePath_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InvasionFilePath_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ThresholdValue_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CutoffThresholdResult_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InvasionPercentageResult_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InvasionDepthResult_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CellNumberResult_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ImageSaveDirectory_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InvExpResult_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nInvExpResult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nInvExpResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
