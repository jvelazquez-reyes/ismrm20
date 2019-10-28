function qMRFit_BIDS(inputDir)

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
             [Model,data] = qMRLabBIDSmapper(files,metas,protomapper);
             disp(['Fitting ' Model.ModelName ' for ' 'sub-' BIDS_subjects{sub_iter} cell2mat(extra) ' ' extras{extra_iter}])
             FitResults = FitData(data,Model,0);
             FitResultsSave_nii(FitResults,files{1},curSubDir);
             renameForBIDS(BIDS_subjects{sub_iter},curSubDir,protomapper,Model.xnames,cell2mat(extra),extras{extra_iter});
             
             
             % FIT HERE BUT MANAGE OUTPUT FOLDERS AND OUTPUT NAMES PROPERLY
             % YOU NEED TO CREATE DERIVATIVES FOLDER/SUBJECT/qMRLab 
             %
             
         end
         
        else
            
            % TODO: Put 47-52 into a function and manage here as well. 
        end
       
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
          extra = setxor(fnames(idx),[protomapper.REQUIREDEntities,{'filename'    'ext'    'type'    'sub'}]);
          
      end
      
  end
  
end

function [Model,data] = qMRLabBIDSmapper(datas,metas,protomapper)


% ================== Instantiate qMRLab object 

eval(['Model=' protomapper.qMRLabModel ';']);

% =================== DATA START 

% ===============================================================
% ROUTE ACTION: MERGE 

if strcmp(protomapper.routeAction,'merge')

% ---------------------------------------------- DATA START  

tmp = load_nii_data(datas{1});
sz = size(tmp);
DATA = [];

if ndims(tmp)==2
    
    if strcmp(protomapper.singleton,'1')
        DATA = zeros(sz(1),sz(2),1,length(datas));
    else
        DATA = zeros(sz(1),sz(2),length(datas));
    end
    
elseif ndims(tmp)==3
    
    DATA = zeros(sz(1),sz(2),sz(3),length(datas));
else
    
    error('Data is not a volume or a slice.');
end

for ii=1:length(datas)
    
    
    if ndims(tmp)==2   && ~strcmp(protomapper.singleton,'1')
        DATA(:,:,ii) =  double(load_nii_data(datas{ii}));
    else
        DATA(:,:,:,ii) =  double(load_nii_data(datas{ii}));
        
    end
    
end


data = struct(); 
data.(protomapper.dataFieldName) = DATA;
clear('tmp','DATA'); 
% ---------------------------------------------- DATA END 
if ~isstruct(metas)
    
    str = cell2struct(metas,'tmp');
    metas = [str.tmp];

end


params = setxor('foreach',fieldnames(protomapper.protMap)); 

for ii=1:length(params)
    
    cur_param = cell2mat(params(ii));
    
    if strcmp(protomapper.protMap.(cur_param).order,'col_first')
   
    Model.Prot.(protomapper.protMap.(cur_param).qMRLabProt).Mat(:,ii) = ...
    [metas(:).(cur_param)].*protomapper.protMap.(cur_param).scale;
    
    else
        
    Model.Prot.(protomapper.protMap.(cur_param).qMRLabProt).Mat(ii,:) = ...
    [metas(:).(cur_param)].*protomapper.protMap.(cur_param).scale;

    end
end



% ==================================================
elseif strcmp(protomapper.routeAction,'distribute')
    
% ===============================================================
% ROUTE ACTION: DISTRIBUTE 
% ===============================================================

% TODO: For other models like mt_sat, we'll distribute data. 

end

end

function renameForBIDS(curSubName,fileDir,protomapper,xnames,varargin) 

for ii=1:length(xnames)
f = dir(fullfile(fileDir,[xnames{ii} '.nii.gz']));

if ~isempty(f) 
    
    if nargin>4
        newname = ['sub-' curSubName '_' varargin{1} '-' varargin{2} '_' protomapper.outputMap.(xnames{ii}) '.nii.gz'];
    else
        newname = ['sub-' curSubName '_' protomapper.outputMap.(xnames{ii}) '.nii.gz'];

    end

   movefile([fileDir filesep f.name],[fileDir filesep newname]);
   
end
end 

end