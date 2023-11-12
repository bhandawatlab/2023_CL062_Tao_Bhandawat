%% Download Datasets
%% download and save the walking chamber dataset
websave('WalkingChamber.zip','https://www.dropbox.com/scl/fi/jrmppc93i3w6yi8o82bc8/WalkingChamber.zip?rlkey=bke578666jjw4783esqhztsux&dl=1');
movefile WalkingChamber.zip Data.zip
if ~exist('WalkingChamberCodeBase/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','WalkingChamberCodeBase/Data');
disp("WalkingChamber Dataset Downloaded and Extracted");

%% download and save the DMD dataset
websave('DMD.zip','https://www.dropbox.com/scl/fi/7zeqc99z7ucqfs1bmffly/DMD.zip?rlkey=tfkgnafxh5agk133hu3278l40&dl=1');
movefile DMD.zip Data.zip
if ~exist('DMDCodeBase/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','DMDCodeBase/Data')
disp("DMD Dataset Downloaded and Extracted");

%% download and save the connectomics dataset
websave('Connectomics.zip','https://www.dropbox.com/scl/fi/77t2m4smnk0k6cew07fcu/Connectomics.zip?rlkey=zpna11owxi4f6532fgnpio33x&dl=1');
movefile Connectomics.zip Data.zip
if ~exist('Connectomics/RConnectomicsAnalysis/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','Connectomics/RConnectomicsAnalysis/Data');
disp("Connectomics Dataset Downloaded and Extracted");

%% delete zip file
delete Data.zip