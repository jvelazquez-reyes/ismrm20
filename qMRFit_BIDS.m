function qMRFit_BIDS(inputDir)

config_file = json2struct([inputDir '/config_file.json']);

try
    BIDS = bids.layout(inputDir);
catch
    error('Cannot load BIDS layout for the provided directory.');
end

% Create derivatives folder 
derivDir = [BIDS.dir filesep 'derivatives' filesep 'qMRLab'];

if ~exist(derivDir, 'dir')
    mkdir(derivDir);
end

BIDS_subjects = bids.query(BIDS,'subjects');

for sub_iter = 1:length(BIDS_subjects)
    
    BIDS_types = bids.query(BIDS,'types','sub',BIDS_subjects(sub_iter));
    BIDS_types = intersect(config_file.Modality{:}, BIDS_types);
    for type_iter = 1:length(BIDS_types)
    
       protomapper = getMapper(BIDS_types{type_iter});
        
       if ~isempty(protomapper)
       % This means that suffix has a correspondance in qMRLab     
        
       % Create subject folder under derivatives 
       curSubDir = [derivDir filesep 'sub-' BIDS_subjects{sub_iter}];
       
       if ~exist(curSubDir, 'dir')
           mkdir(curSubDir);
       end
       
        % If there are additional naming entities, other than those required
        % for a grouping suffix, we need to iterate over them. 
        
        extra = checkExtraEntities(BIDS,protomapper,sub_iter,BIDS_types{type_iter});
       
        if ~isempty(extra) 
        %TODO: Deal with more than one extra entities later. 
        
            extras = unique({BIDS.subjects(sub_iter).anat.(cell2mat(extra))});
        
         for extra_iter = 1:length(extras)   
         
             files = bids.query(BIDS,'data','sub',BIDS_subjects(sub_iter),cell2mat(extra),extras(extra_iter),'type',BIDS_types(type_iter));
             metas = bids.query(BIDS,'metadata','sub',BIDS_subjects(sub_iter),cell2mat(extra),extras(extra_iter),'type',BIDS_types(type_iter));
             
             qMRLabBIDSmapper(files,metas,protomapper, curSubDir);
             
             %disp(['Fitting ' Model.ModelName ' for ' 'sub-' BIDS_subjects{sub_iter} cell2mat(extra) ' ' extras{extra_iter}])

         end
         
        else
            
            files = bids.query(BIDS,'data','sub',BIDS_subjects(sub_iter),'type',BIDS_types(type_iter));
            metas = bids.query(BIDS,'metadata','sub',BIDS_subjects(sub_iter),'type',BIDS_types(type_iter));
            
            qMRLabBIDSmapper(files,metas,protomapper, curSubDir);
            
            %disp(['Fitting ' Model.ModelName ' for ' 'sub-' BIDS_subjects{sub_iter} ]) 
        end
        
             %renameMapsSaveJsons(BIDS,files,BIDS_subjects{sub_iter},curSubDir,protomapper,Model.xnames,cell2mat(extra),extras{extra_iter});
       
       end

    end
        
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
        end
            

end

function extra = checkExtraEntities(BIDS,protomapper,sub_iter,cur_type) 
  
  anat = BIDS.subjects(sub_iter).anat;
  
  for ii =1:length(anat)
      
      if strcmp(anat(ii).type,cur_type)
          
          idx = not(structfun(@isempty, anat(ii)));
          fnames = fieldnames(anat(ii)); 
          extra = setxor(fnames(idx),[protomapper.REQUIREDEntities{:},{'filename'    'ext'    'type'    'sub'}]);
          
      end
      
  end
  
end

function [Model,data] = qMRLabBIDSmapper(datas,metas,protomapper, curSubDir,varargin)

p = inputParser();

%Input parameters conditions
validNii = @(x) exist(x,'file') && strcmp(x(end-5:end),'nii.gz');
validJsn = @(x) exist(x,'file') && strcmp(x(end-3:end),'json');
validB1factor = @(x) isnumeric(x) && (x > 0 && x <= 1);

%Add OPTIONAL Parameteres
addParameter(p,'mask',[],validNii);
addParameter(p,'b1map',[],validNii);
addParameter(p,'b1factor',[],validB1factor);
addParameter(p,'type',[],@ischar);
addParameter(p,'order',[],@isnumeric);
addParameter(p,'dimension',[],@ischar);
addParameter(p,'size',[],@ismatrix);
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


% ================== Instantiate qMRLab object 

eval(['Model=' protomapper.qMRLabModel ';']);

% =================== DATA START 

% ===============================================================
% ROUTE ACTION: MERGE 

if strcmp(protomapper.routeAction,'merge')

% ---- DATA START  

sample = load_nii_data(datas{1});
sz = size(sample);

if ndims(sample)==2
    
    if strcmp(protomapper.singleton,'1')
        DATA = zeros(sz(1),sz(2),1,length(datas));
    else
        DATA = zeros(sz(1),sz(2),length(datas));
    end
    
elseif ndims(sample)==3
    
    DATA = zeros(sz(1),sz(2),sz(3),length(datas));
else
    
    error('Data is not a volume or a slice.');
end

for ii=1:length(datas)
    
    
    if ndims(sample)==2   && ~strcmp(protomapper.singleton,'1')
        DATA(:,:,ii) =  double(load_nii_data(datas{ii}));
    else
        DATA(:,:,:,ii) =  double(load_nii_data(datas{ii}));
        
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

% ==================================================
elseif strcmp(protomapper.routeAction,'distribute')
    
% ===============================================================
% ROUTE ACTION: DISTRIBUTE 
% ===============================================================

input_data = protomapper.dataFieldName;
qLen = length(datas);
for ii=1:qLen
    cur_data = cell2mat(input_data{ii});
    data.(cur_data) = double(load_nii_data(datas{ii}));
end
% ----------- DATA END 

end

%Set Model.Protocol
fields = setxor('foreach',fieldnames(protomapper.protMap));
qLen = length(metas);
for kk=1:length(fields)
cur_field = cell2mat(fields(kk));

if strcmp(protomapper.protMap.(cur_field).fillProtBy, 'files')
    cur_field = cell2mat(fields(kk));
    params = protomapper.protMap.(cur_field).qMRLabProt;
    for jj=1:length(params)
        for ii=1:qLen
            if isfield(metas{ii}, params{jj})
                Model.Prot.(cur_field).Mat(ii,jj) = ...
                    str2double(metas{ii}.(params{jj}));
            else
                disp('No parameter found')
            end
        end
    end
end

if strcmp(protomapper.protMap.(cur_field).fillProtBy, 'parameter')
    cur_field = cell2mat(fields(kk));
    params = protomapper.protMap.(cur_field).qMRLabProt;
    count = 1;
    if kk==1; jj=1; end
    while count < length(params) + 1
        for ii=1:length(params)
            if isfield(metas{jj}, params{ii})
                Model.Prot.(cur_field).Mat(count) = ...
                    str2double(metas{jj}.(params{ii}));
                if ((count ~= length(params) + 1) && (count == length(fieldnames(metas{jj}))))
                    jj = jj + 1;
                end
                count = count + 1;
            else
                disp('No parameter found')
            end 
        end           
    end
end
end

%Account for optional inputs and options
if ~isempty(p.Results.mask); data.Mask = double(load_nii_data(p.Results.mask)); end
if ~isempty(p.Results.b1map); data.B1map = double(load_nii_data(p.Results.b1map)); end
if ~isempty(p.Results.b1factor); Model.options.B1correction = p.Results.b1factor; end
if ~isempty(p.Results.sid); SID = p.Results.sid; end
if ~isempty(p.Results.type); Model.options.Smoothingfilter_Type = p.Results.type; end
if ~isempty(p.Results.order); Model.options.Smoothingfilter_order = p.Results.order; end
if ~isempty(p.Results.dimension); Model.options.Smoothingfilter_dimension = p.Results.dimension; end
if ~isempty(p.Results.size)
    Model.options.Smoothingfilter_sizex = p.Results.size(1);
    Model.options.Smoothingfilter_sizey = p.Results.size(2);
    Model.options.Smoothingfilter_sizez = p.Results.size(3);
end

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
            
        FitResultsSave_nii(FitResults,datas{1},curSubDir);
             
        % Save qMRLab object
        if ~isempty(p.Results.sid)
            Model.saveObj([curSubDir filesep SID '_' protomapper.qMRLabModel '.qmrlab.mat']);
        else
            Model.saveObj([curSubDir filesep protomapper.qMRLabModel '.qmrlab.mat']);
        end
        
        % Remove FitResults.mat
        delete([curSubDir filesep 'FitResults.mat']);
             
        % JSON files for quantitative map(s)
        addField = struct();
        addField.EstimationReference =  protomapper.estimationPaper;
        addField.EstimationAlgorithm =  protomapper.estimationAlgorithm;
        addField.BasedOn = [{datas},{metas{:}}];
             
        provenance = Model.getProvenance('extra',addField);
             
             for ii=1:length(outputs)
                 cur_output = cell2mat(outputs(ii));
                 rename_output = protomapper.outputMap.(cur_output);
                 if ~isempty(p.Results.sid)
                     % ==== Rename outputs ====
                     movefile([curSubDir filesep cur_output '.nii.gz'],[curSubDir filesep SID '_' rename_output '.nii.gz']);
                     % ==== Save JSON provenance ==== 
                     savejson('',provenance,[curSubDir filesep SID '_' rename_output '.json']);
                 else
                     movefile([curSubDir filesep cur_output '.nii.gz'],[curSubDir filesep rename_output '.nii.gz']);
                     savejson('',provenance,[curSubDir filesep rename_output '.json']);
                 end
             end

end

function renameMapsSaveJsons(BIDS,files,curSubName,fileDir,protomapper,xnames,varargin) 

for ii=1:length(xnames)
f = dir(fullfile(fileDir,[xnames{ii} '.nii.gz']));

if ~isempty(f) 
    
    if nargin>6
        newname = ['sub-' curSubName '_' varargin{1} '-' varargin{2} '_' protomapper.outputMap.(xnames{ii})];
    else
        newname = ['sub-' curSubName '_' protomapper.outputMap.(xnames{ii})];

    end

   movefile([fileDir filesep f.name],[fileDir filesep newname '.nii.gz']);
   provenance = getProvenance(BIDS,files,protomapper);
   savejson('',provenance,[fileDir filesep newname '.json']);
    
end
end 

f = dir(fullfile(fileDir,'FitResults.mat'));
if ~isempty(f) 
    newname = newname(1:(max(strfind(newname,'_'))-1));
    movefile([fileDir filesep f.name],[fileDir filesep newname '_FitResults.mat']);
end

end

function FitProvenance = getProvenance(BIDS,files,protomapper)

            FitProvenance = struct();
            FitProvenance.BasedOn = regexprep(files,BIDS.dir(1:max(strfind(BIDS.dir,filesep))),'','ignorecase');
            FitProvenance.EstimationSoftwareName = 'qMRLab';
            FitProvenance.EstimationSoftwareVer  = qMRLabVer;
            FitProvenance.EstimationReference = protomapper.estimationPaper;
            FitProvenance.EstimationAlgorithm = protomapper.estimationAlgorithm;
  
          
            if moxunit_util_platform_is_octave
                
                FitProvenance.EstimationDate = strftime('%Y-%m-%d %H:%M:%S', localtime (time ()));
                [FitProvenance.EstimationSoftwareEnv, FitProvenance.MaxSize, FitProvenance.Endian] = computer;
                FitProvenance.EstimationSoftwareEnvDetails = GetOSDetails();
                FitProvenance.EstimationSoftwareLang = ['Octave ' OCTAVE_VERSION()];
                Fitprovenance.EstimationSoftwareLangDetails = pkg('list');
                
            else 

                FitProvenance.EstimationDate = datetime(now,'ConvertFrom','datenum');
                [FitProvenance.EstimationSoftwareEnv, FitProvenance.MaxSize, FitProvenance.Endian] = computer; 
                FitProvenance.EstimationSoftwareEnvDetails = GetOSDetails();
                FitProvenance.EstimationSoftwareLang = ['Matlab ' version('-release')];
                FitProvenance.EstimationSoftwareLangDetails = ver;

            end
            
end

function details = GetOSDetails
  
    type = computer;

    if moxunit_util_platform_is_octave
        
        if ~isempty(strfind(type,'apple')) % OSX Octave 
            
            [st,out] = unix('cat /etc/os-release');
            
            if ~st
                details = out;
            else
                details = [];
            end
            
        end
        
        if ~isempty(strfind(type,'linux')) % GNU Linux Octave 
            
            [st,out] = unix('system_profiler SPSoftwareDataType');
            
            if ~st
                details = out;
            else
                details = [];
            end
            
        end
        
        if ~isempty(strfind(type,'windows')) % GNU Linux Octave 
            
            [st,out] = system('winver');
            
            if ~st
                details = out;
            else
                details = [];
            end
            
        end
        
    else % MATLAB 
        
        if strncmp(type,'MAC',3)
            
            [st,out] = unix('system_profiler SPSoftwareDataType');
            
            if ~st
                details = rmUserInfoOSX(out);
             
            else
                details = [];
            end
            
        end
        
        if strncmp(type,'GLNX',4)
            
            [st,out] = unix('cat /etc/os-release');
            
            if ~st
                details = out;
            else
                details = [];
            end
            
        end
        
        if ~isempty(strfind(type,'WIN'))
            
            [st,out] = system('winver');
            
            if ~st
                details = out;
            else
                details = [];
            end
            
        end
        
        
    end
end

function out = rmUserInfoOSX(ipt)

    usridx = strfind(ipt,'User Name');
    pcidx = strfind(ipt,'Computer Name');
    nwlines = strfind(ipt,char(10));
   
    if ~isempty(usridx)
    ipt = hideUser(usridx,nwlines,ipt);
    end
    
    if ~isempty(pcidx)
    ipt = hideUser(pcidx,nwlines,ipt);
    end

    out = ipt;

end

function out = hideUser(idx,nwlines,ipt)

    tmp = nwlines - idx;
    tmp = min(tmp(tmp>0));
    interval = idx:idx+min(tmp)-1;

    ipt(interval) = '*';
    ipt(max(interval)+1) = char(10);

    out = ipt;

end