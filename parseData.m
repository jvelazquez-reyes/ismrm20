% =====================================================================
% This function is highly customized to parse data collected 
% specifically for this study. This is not a generic BIDS converter. 
%
% Functionality depends on the hardcoded filenames and directory 
% organization.
% 
% Depends on some functions in qMRLab and also loadRthData.m. 
%
% You can download data & repeat this: https://osf.io/4fztd/download/ 
% 
% Author: Agah Karakuzu
% October 2019 
% =====================================================================


function parseData(inputDir)

% inputDir contains two folders: i) rth and ii) sie. 

% For this study, data from RTHawk vfa_t1 sequence, two volumes were exported
% in *.dat format. Images first image: higher FA, second image: lower FA. 

% As the RTHawk sequence was its in early stages of developments, modules 
% to export acquisition parameters were not implemented yet. However, scan
% parameters are known: 

% FlipAngles; 20 and 3 for the first and the second volume, respectively. 
% RepetitionTime: 25ms 
% EchoTime: 4.5ms 
% Resolution: 1X1X5mm 
% Quadratic phase increment: 117 
% Pulse duration: 1.3ms 
% Slab select gradient duration: 3ms 
% FOV: 25.6X25.6 cm in-plane 
% Slab thickness: 25mm 
% Base resolution: 256 X 256 in-plane
% Slab base resolution: 5 
% Spoiler gradient area: 10 cyc/cm 

% I am aware that there are many things hardcoded, but explaning everything
% as much as I can so that the rationale behind operations are clear. 

outputDir = [inputDir filesep 'ds-nist'];

mkdir([outputDir filesep 'sub-01']); %sie
mkdir([outputDir filesep 'sub-01/anat']); % rth

dscr = struct(); 
dscr.Name = 'ISMRM/NIST system phantom scans with Siemens and RTHawk VFA';
dscr.BIDSVersion = '1.2.1';
savejson('',dscr,[outputDir filesep 'dataset_description.json']);

rthData = struct(); 
sieData = struct(); 
 
rthData(1).im = getRthSlice([inputDir '/raw_rth/VFA T12019960020.dat'],3); % FA 3
rthData(1).hdr.FlipAngle = 3; 
rthData(1).hdr.RepetitionTimeExcitation = 25;

rthData(2).im = getRthSlice([inputDir '/raw_rth/VFA T12019960019.dat'],3); % FA 20
rthData(2).hdr.FlipAngle = 20; 
rthData(2).hdr.RepetitionTimeExcitation = 25;

% FA 3 Siemens 
sieData(1).info = dicominfo([inputDir '/raw_sie/T1_FL3D_FLIP3_0004/' ...
'RTHAWKTEST-RETEST.MR.IRM_RECHERCHE_TESTS.0004.0008.2019.10.11.19.32.40.748353.139124355.IMA']);
sieData(1).im = double(dicomread(sieData(1).info));
sieData(1).hdr.FlipAngle = sieData(1).info.FlipAngle; 
% TR in Siemens is given in ms --> s 
sieData(1).hdr.RepetitionTimeExcitation = sieData(1).info.RepetitionTime;  

% FA 20 Siemens 
sieData(2).info = dicominfo([inputDir '/raw_sie/T1_FL3D_FLIP20_0003/' ...
'RTHAWKTEST-RETEST.MR.IRM_RECHERCHE_TESTS.0003.0008.2019.10.11.19.32.40.748353.139123329.IMA']);
sieData(2).im = double(dicomread(sieData(2).info));
sieData(2).hdr.FlipAngle = sieData(2).info.FlipAngle; 
% TR in Siemens is given in ms --> s 
sieData(2).hdr.RepetitionTimeExcitation = sieData(2).info.RepetitionTime; 

% Save siemens data in BIDS-like format 
savejson('',sieData(1).hdr,[outputDir '/sub-01/anat/sub-01_acq-siemens_fa-1_VFA.json']);
nii = make_nii(sieData(1).im, [1,1,5], [0,0,0],64);
save_nii(nii,[outputDir '/sub-01/anat/sub-01_acq-siemens_fa-1_VFA.nii.gz']);

savejson('',sieData(2).hdr,[outputDir '/sub-01/anat/sub-01_acq-siemens_fa-2_VFA.json']);
nii = make_nii(sieData(2).im, [1,1,5], [0,0,0],64);
save_nii(nii,[outputDir '/sub-01/anat/sub-01_acq-siemens_fa-2_VFA.nii.gz']);

% Save rth data in BIDS-like format 

savejson('',rthData(1).hdr,[outputDir '/sub-01/anat/sub-01_acq-rthawk_fa-1_VFA.json']);
nii = make_nii(rthData(1).im, [1,1,5], [0,0,0],64);
save_nii(nii,[outputDir '/sub-01/anat/sub-01_acq-rthawk_fa-1_VFA.nii.gz']);

savejson('',rthData(2).hdr,[outputDir '/sub-01/anat/sub-01_acq-rthawk_fa-2_VFA.json']);
nii = make_nii(rthData(2).im, [1,1,5], [0,0,0],64);
save_nii(nii,[outputDir '/sub-01/anat/sub-01_acq-rthawk_fa-2_VFA.nii.gz']);

disp('==== DONE ===');
disp(['Saved BIDS formatted data under: ' outputDir]); 

end

function slice = getRthSlice(file,sliceNum)

[im,meta] = loadRthData(file); 

% rthImageExport module gives real in the first, imaginary part in the 2nd
im = im(1,:) + 1i*im(2,:);

% Base resolution to reshape data is available in *.dat file. 
vol = reshape(im,[meta.extent0,meta.extent1,meta.extent2]);

% We will rotate these images to bring them to the same orientation with 
% Siemens images. We are only interested in one slice, selecting the 
% 3rd one out of 5 total slices. 
if sliceNum>5 || sliceNum<=0, error('Invalid slice. Available: 1-5'); end
    
slice = vol(:,:,sliceNum);

% Just bringing it to the same orientation with Siemens image. 
slice = imrotate(slice,90); 

end
