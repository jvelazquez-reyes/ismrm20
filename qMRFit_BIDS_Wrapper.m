function qMRFit_BIDS_Wrapper(nii_array, json_array, varargin)

% Supress verbose Octave warnings.
%if moxunit_util_platform_is_octave
%    warning('off','all');
%end

% This env var will be consumed by qMRLab
setenv('ISNEXTFLOW','1');

p = inputParser();

%Input parameters conditions
validNii = @(x) exist(x,'file') && strcmp(x(end-5:end),'nii.gz');
validJsn = @(x) exist(x,'file') && strcmp(x(end-3:end),'json');
validB1factor = @(x) isnumeric(x) && (x > 0 && x <= 1);

%Add OPTIONAL Parameteres
addParameter(p,'mask',[],validNii);
addParameter(p,'b1map',[],validNii);
addParameter(p,'b1factor',[],validB1factor);
addParameter(p,'qmrlab_path',[],@ischar);
addParameter(p,'sid',[],@ischar);
addParameter(p,'containerType',@ischar);
addParameter(p,'containerTag',[],@ischar);
addParameter(p,'description',@ischar);
addParameter(p,'datasetDOI',[],@ischar);
addParameter(p,'datasetURL',[],@ischar);
addParameter(p,'datasetVersion',[],@ischar);

parse(p,varargin{:});

if ~isempty(p.Results.qmrlab_path); qMRdir = p.Results.qmrlab_path; end

try
    disp('=============================');
    qMRLabVer;
catch
    warning('Cant find qMRLab. Adding qMRLab_DIR to the path: ');
    if ~strcmp(qMRdir,'null')
        qmr_init(qMRdir);
    else
        error('Please set qMRLab_DIR parameter in the nextflow.config file.');
    end
    qMRLabVer();
end

protomapper = getMapper(qMR_suffix);

% ================== Instantiate qMRLab object 

eval(['Model=' protomapper.qMRLabModel ';']);

% =================== DATA START 

% ===============================================================
% ROUTE ACTION: MERGE 

if strcmp(protomapper.routeAction,'merge')

% ---------------------------------------------- DATA START  

tmp = load_nii_data(nii_array{1});
sz = size(tmp);
DATA = [];

if ndims(tmp)==2
    
    if strcmp(protomapper.singleton,'1')
        DATA = zeros(sz(1),sz(2),1,length(nii_array));
    else
        DATA = zeros(sz(1),sz(2),length(nii_array));
    end
    
elseif ndims(tmp)==3
    
    DATA = zeros(sz(1),sz(2),sz(3),length(nii_array));
else
    
    error('Data is not a volume or a slice.');
end

for ii=1:length(nii_array)
    
    
    if ndims(tmp)==2   && ~strcmp(protomapper.singleton,'1')
        DATA(:,:,ii) =  double(load_nii_data(nii_array{ii}));
    else
        DATA(:,:,:,ii) =  double(load_nii_data(nii_array{ii}));
        
    end
    
end


data = struct(); 
data.(protomapper.dataFieldName) = DATA;
clear('tmp','DATA'); 
% ---------------------------------------------- DATA END 
if ~isstruct(json_array)
    
    str = cell2struct(json_array,'tmp');
    json_array = [str.tmp];

end


params = setxor('foreach',fieldnames(protomapper.protMap)); 

for ii=1:length(params)
    
    cur_param = cell2mat(params(ii));
    
    if strcmp(protomapper.protMap.(cur_param).order,'col_first')
   
    Model.Prot.(protomapper.protMap.(cur_param).qMRLabProt).Mat(:,ii) = ...
    [json_array(:).(cur_param)].*protomapper.protMap.(cur_param).scale;
    
    else
        
    Model.Prot.(protomapper.protMap.(cur_param).qMRLabProt).Mat(ii,:) = ...
    [json_array(:).(cur_param)].*protomapper.protMap.(cur_param).scale;

    end
end



% ==================================================
elseif strcmp(protomapper.routeAction,'distribute')
    
% ===============================================================
% ROUTE ACTION: DISTRIBUTE 
% ===============================================================

% TODO: For other models like mt_sat, we'll distribute data. 

end

% ===============================================================
% FITTING 
% ===============================================================
             
%disp(['Fitting ' Model.ModelName ' for ' 'sub-' BIDS_subjects{sub_iter} cell2mat(extra) ' ' extras{extra_iter}])
FitResults = FitData(data,Model,0);
             
             
FitResultsSave_nii(FitResults,nii_array{1},pwd);
renameMapsSaveJsons(BIDS,files,BIDS_subjects{sub_iter},curSubDir,protomapper,Model.xnames,cell2mat(extra),extras{extra_iter});
   
             % FIT HERE BUT MANAGE OUTPUT FOLDERS AND OUTPUT NAMES PROPERLY
             % YOU NEED TO CREATE DERIVATIVES FOLDER/SUBJECT/qMRLab 
             %
        
end

function protomapper = getMapper(cur_type)
 % For each type defined in qMRLab, there is a protocol mapper 
 % called <subffix>_BIDSmapper.json under src/common/BIDS_protomaps 
    
        qmrTypeList = dir(fullfile('./BIDS_protomaps', '*.json'));
        tmp = struct2cell(qmrTypeList); % Cell matrix
        qmrTypeList = tmp(1,:); % Cell array, 1st index is "name" field
        % Check if the current type is defined within qMRLab 
        
        if ismember([cur_type '.json'],qmrTypeList)
            
            protomapper = json2struct(['./BIDS_protomaps/' cur_type '.json']);
        else
            protomapper = [];
            disp(['Protocol' cur_type 'not available'])
        end
            

end

function [Model,data] = qMRLabBIDSmapper(datas,metas,protomapper)




end
