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
data = struct(); 

% =================== DATA START 

% ===============================================================
% ROUTE ACTION: MERGE 

if strcmp(protomapper.routeAction,'merge')

% ---- DATA START  

sample = load_nii_data(nii_array{1});
sz = size(sample);

if ndims(sample)==2
    
    if strcmp(protomapper.singleton,'1')
        DATA = zeros(sz(1),sz(2),1,length(nii_array));
    else
        DATA = zeros(sz(1),sz(2),length(nii_array));
    end
    
elseif ndims(sample)==3
    
    DATA = zeros(sz(1),sz(2),sz(3),length(nii_array));
else
    
    error('Data is not a volume or a slice.');
end

for ii=1:length(nii_array)
    
    
    if ndims(sample)==2   && ~strcmp(protomapper.singleton,'1')
        DATA(:,:,ii) =  double(load_nii_data(nii_array{ii}));
    else
        DATA(:,:,:,ii) =  double(load_nii_data(nii_array{ii}));
        
    end
    
end

data.(protomapper.dataFieldName) = DATA;
clear('sample','DATA'); 
% ---------------------------------------------- DATA END 
% if ~isstruct(json_array)
%     
%     str = cell2struct(json_array,'tmp');
%     json_array = [str.tmp];
% 
% end


params = setxor('foreach',fieldnames(protomapper.protMap)); 
qLen = length(nii_array);
for jj=1:qLen
    for ii=1:length(params)
        
        cur_param = cell2mat(params(ii));

        Model.Prot.(protomapper.protMap.(cur_param).qMRLabProt).Mat(jj,ii) = ...
        getfield(json2struct(json_array{jj}),params{ii});
    end
end

%Account for optional inputs and options
if ~isempty(p.Results.mask); data.Mask = double(load_nii_data(p.Results.mask)); end
if ~isempty(p.Results.b1map); data.B1map = double(load_nii_data(p.Results.b1map)); end
if ~isempty(p.Results.b1factor); Model.options.B1correction = p.Results.b1factor; end
if ~isempty(p.Results.sid); SID = p.Results.sid; end

% ==================================================
elseif strcmp(protomapper.routeAction,'distribute')
    
% ===============================================================
% ROUTE ACTION: DISTRIBUTE 
% ===============================================================

input_data = fieldnames(protomapper.dataFieldName);
qLen = length(nii_array);
for ii=1:qLen
    cur_data = cell2mat(input_data{ii});
    data.(cur_data) = nii_array{ii};
end
% ----------- DATA END 
% if ~isstruct(json_array)
%     
%     str = cell2struct(json_array,'tmp');
%     json_array = [str.tmp];
% 
% end

params = setxor('foreach',fieldnames(protomapper.protMap)); 
for jj=1:length(params)
    
    cur_param = cell2mat(params(jj));
        
    for ii=1:length(protomapper.protMap.(cur_param).qMRLabProt)
        
        Model.Prot.(protomapper.protMap.(cur_param).qMRLabProt{ii}).Mat(1,jj) = ...
        getfield(json2struct(json_array{ii}),params{jj});
    end
end

%Account for optional inputs and options
if ~isempty(p.Results.mask); data.Mask = double(load_nii_data(p.Results.mask)); end
if ~isempty(p.Results.b1map); data.B1map = double(load_nii_data(p.Results.b1map)); end
if ~isempty(p.Results.b1factor); Model.options.B1correction = p.Results.b1factor; end
if ~isempty(p.Results.sid); SID = p.Results.sid; end

end

% ===============================================================
% FITTING 
% ===============================================================
             
FitResults = FitData(data,Model,0);

outputs = fieldnames(protomapper.outputMap); 
for ii=1:length(outputs)
    
    cur_output = cell2mat(outputs(ii));
    
    % ==== Weed out spurious values ==== 

    % Zero-out Inf values (caused by masking)
    FitResults.(cur_output)(FitResults.(cur_output)==Inf)=0;
    % Null-out negative values 
    FitResults.(cur_output)(FitResults.(cur_output)<0)=NaN;

end

% ==== Save outputs ==== 
disp('-----------------------------');
disp('Saving fit results...');

FitResultsSave_nii(FitResults,nii_array{1},pwd);

% Save qMRLab object
if ~isempty(SID)
    Model.saveObj([SID protomapper.qMRLabModel '.qmrlab.mat']);
else
    Model.saveObj([protomapper.qMRLabModel '.qmrlab.mat']);
end

% Remove FitResults.mat 
delete('FitResults.mat');

% JSON files for quantitative map(s)
addField = struct();
addField.EstimationReference =  protomapper.estimationPaper;
addField.EstimationAlgorithm =  protomapper.estimationAlgorithm;
addField.BasedOn = [{nii_array},{jsn_array}];

provenance = Model.getProvenance('extra',addField);

for ii=1:length(outputs)
      
    if ~isempty(SID)
        
        % ==== Rename outputs ==== 
        movefile([outputs(ii) '.nii.gz'],[SID '_' outputs(ii) 'map.nii.gz']);
        % ==== Save JSON provenance ==== 
        savejson('',provenance,[pwd filesep SID '_' outputs(ii) 'map.json']);
    else
        movefile([outputs(ii) '.nii.gz'],[outputs(ii) 'map.nii.gz']);
        savejson('',provenance,[pwd filesep outputs(ii) 'map.json']);
    end
   
end

% JSON file for dataset_description
addDescription = struct();
addDescription.Name = 'qMRLab Outputs';
addDescription.BIDSVersion = '1.5.0';
addDescription.DatasetType = 'derivative';
addDescription.GeneratedBy.Name = 'qMRLab';
addDescription.GeneratedBy.Version = qMRLabVer();
addDescription.GeneratedBy.Container.Type = p.Results.containerType;
if ~strcmp(p.Results.containerTag,'null'); addDescription.GeneratedBy.Container.Tag = p.Results.containerTag; end
addDescription.GeneratedBy.Name2 = 'Manual';
addDescription.GeneratedBy.Description = p.Results.description;
if ~isempty(p.Results.datasetDOI); addDescription.SourceDatasets.DOI = p.Results.datasetDOI; end
if ~isempty(p.Results.datasetURL); addDescription.SourceDatasets.URL = p.Results.datasetURL; end
if ~isempty(p.Results.datasetVersion); addDescription.SourceDatasets.Version = p.Results.datasetVersion; end

savejson('',addDescription,[pwd filesep 'dataset_description.json']);

if moxunit_util_platform_is_octave
    warning('on','all');
end            
             

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

function out = json2struct(filename)

tmp = loadjson(filename);

if isstruct(tmp)

    out = tmp;

else

    str = cell2struct(tmp,'tmp');
    out = [str.tmp];

end

end 

function qmr_init(qmrdir)

run([qmrdir filesep 'startup.m']);

end