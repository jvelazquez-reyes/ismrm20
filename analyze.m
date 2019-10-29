% Fit data 
qMRFit_BIDS('./data/ds-nist');

rthT1 = load_nii_data('./data/ds-nist/derivatives/qMRLab/sub-01/sub-01_acq-rthawk_T1map.nii.gz');
sieT1 = load_nii_data('./data/ds-nist/derivatives/qMRLab/sub-01/sub-01_acq-siemens_T1map.nii.gz');

% Load labeled masks 
rthMask = load_nii_data('./data/masks/sub-01_acq-rthawk_mask.nii.gz');
sieMask = load_nii_data('./data/masks/sub-01_acq-siemens_mask.nii.gz');

% Generate reference 
refvals = csvread('./System_Phantom/phantom_reference.csv');
refvals = flip(refvals);
refIm = sieMask;
rng('default');

% Confine analysis to 75 voxels per sphere. 
refim = zeros(256,56);
sieTmp = zeros(256,256);
rthTmp = sieTmp;
sieVec = zeros(75,10);
rthVec = sieVec; 
refVec = sieVec; 

for ii=1:10
% Reference image 
idx = find(sieMask(:)==ii); 
idx = idx(1:75);
refIm(idx) = normrnd(refvals(ii,2),refvals(ii,3),[1,75]);
refVec(:,ii) = refIm(idx);
% Siemens T1 (still the same idx) 
sieTmp(idx) = sieT1(idx);
sieVec(:,ii) = sieT1(idx);
% Rth T1 (we need new indexes, masks are different) 
idx = find(rthMask(:)==ii); 
idx = idx(1:75);
rthTmp(idx) = rthT1(idx);
rthVec(:,ii) = rthT1(idx);
end

refIm = double(refIm)./1000;
refVec = double(refVec)./1000;
% Update images to contain speres only.
sieT1 = sieTmp;
rthT1 = rthTmp; 


