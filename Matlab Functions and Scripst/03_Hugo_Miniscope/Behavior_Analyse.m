
%% Behavior Analyse

% Main call

% Files needed:
% - *.csv files

% - Analyse one or more *.csv files

% - The code relies on the following functions/scripts : -> 

    
% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  02/2024
% Last update: 02/2024

%% Define settings

prompt        = {'Experiment Name';'Experiment (Type (1) to  Open Field / Type (2) to Elevated plus maze):'; ...
    'Plot Figures (Type (1) yes / Type (2) no):';'Frames per second:';'Analysis time (Type (1) to full section / Type time in sec.):';'Arena size -> Width (cm):';'Arena size -> Height (cm):';'Experimental animal (Type (1) to Mouse (C57BL/6) / Type (2) to Rat (Wistar))'};

dlgtitle      = 'Define Settings';
dims         = [1 80];
default_input = {'E03_PT','1','2','30','300','30','30','1'};

input_settings = inputdlg(prompt,dlgtitle,dims,default_input);

%% Organize path to load data and initialize some variables

% Load datafiles (*.csv)
[FileName,PathName] = uigetfile({'*.*'},'MultiSelect', 'on'); % Define file type *.*

% Filename can be organize as a single char or a group char in a cell depending on the number os files selected

% header              -> files information and parameters to analyse
% header.last_frame   -> image of the last frame to take the resolution parameters and plot 

% data                  = raw values / each cell corresponds to data from 1 experiment

header.last_frame = cell(size(FileName));
data = cell(size(FileName));

for ii = 1:length(FileName)
    
    if length(FileName) == 1 % condition for a single file selected 
       header.FilePattern = dir(fullfile(PathName, char(FileName)));
       [~, ~, fExt] = fileparts(FileName);
        
       switch lower(fExt)

       case '.csv' 
       output = readtable([PathName '/' FileName],'Delimiter',',');
       data{1,ii} = table2array(output);
       
       % condition to set analysis time
       if str2num(input_settings{5, 1}) ~= 1
          timeidxstop = str2num(input_settings{4, 1}) * str2num(input_settings{5, 1});
          data{1,ii}(timeidxstop+1:1:end,:) = [];
       end
       
       header.Filename_csv(ii) = FileName{ii};
       fprintf(1, 'Reading %s\n', FileName{ii});

       case '.png' 
       header.last_frame{ii} = imread([PathName '/' FileName]);
       header.Filename_png(ii) = FileName{ii};
       fprintf(1, 'Reading %s\n', FileName{ii});
       
       end

    else      % condition for multiple files selected
        header.FilePattern = dir(fullfile(PathName, char(FileName{ii})));
        [~, ~, fExt] = fileparts(FileName{ii});

        switch lower(fExt)

        case '.csv' 
        output = readtable([PathName '/' FileName{ii}],'Delimiter',',');
        data{1,ii} = table2array(output);
        
         % condition to set analysis time
        if str2num(input_settings{5, 1}) ~= 1
            timeidxstop = str2num(input_settings{4, 1}) * str2num(input_settings{5, 1});
            data{1,ii}(timeidxstop+1:1:end,:) = [];
        end
       
        header.Filename_csv{ii} = FileName{ii};
        fprintf(1, 'Reading %s\n', FileName{ii});

        case '.png' 
        header.last_frame{ii} = rgb2gray(imread([PathName '/' FileName{ii}])); % open image and convert to grey scale
        header.Filename_png{ii} = FileName{ii};
        fprintf(1, 'Reading %s\n', FileName{ii});
        
        end
   end
end 

% Remove empty cells
data = data(~cellfun('isempty',data));
header.last_frame = header.last_frame(~cellfun('isempty',header.last_frame));
header.Filename_csv = header.Filename_csv(~cellfun('isempty',header.Filename_csv));
header.Filename_png = header.Filename_png(~cellfun('isempty',header.Filename_png));

%    val_name = FileName(1:regexp(FileName,'.csv')-1);
%    assignin('base',val_name, data);


%% Set deafult values
set(0,'DefaultFigureWindowStyle','normal')
set(groot,'DefaultFigureColormap',jet)

%%
clear ('fExt','FileName','ii','output','PathName','text_','a', 'aa', 'a_pdf', 'ans', 'ax', 'cb', 'd', 'default_input1', 'dims1', 'dlgtitle1', 'f', 'F', 'f1', 'f10', 'f11', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'factor', 'file', ' file2save', 'h', 'ii', 'mu', 'name', 'newNumberOfCols', 'newNumberOfRows', 'p', 'pd', 'pic', 'plots_rc', 'pos', 'pq', 'prompt1', 'sb', 'scale_size', 'sigma', 'unit', 'v', 'x1', 'x2', 'x_values', 'xq', 'y', 'yq');

%% last update 23/02/2024 - 
% listening: ISIS - So Did We
