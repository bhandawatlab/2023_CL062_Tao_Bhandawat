%% Download Datasets
%% download and save the walking chamber dataset
websave('WalkingChamber.zip','https://www.dropbox.com/scl/fi/6d94o2k25hb2qup10e9ih/WalkingChamber.zip?rlkey=898sc0dxy1xevj5nv5x0inyf1&st=5k81ukyr&dl=1');
movefile WalkingChamber.zip Data.zip
if ~exist('WalkingChamberCodeBase/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','WalkingChamberCodeBase/Data');
disp("WalkingChamber Dataset Downloaded and Extracted");

%% download and save the DMD dataset
websave('DMD.zip','https://www.dropbox.com/scl/fi/6rrish4pt7g91ltfgsx9r/DMD.zip?rlkey=65xvuq917hh3u7hm1cd9tds9b&st=zj6xcm14&dl=1');
movefile DMD.zip Data.zip
if ~exist('DMDCodeBase/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','DMDCodeBase/Data')
disp("DMD Dataset Downloaded and Extracted");

%% download and save the connectomics dataset
websave('Connectomics.zip','https://www.dropbox.com/scl/fi/lue7e5bihmjvo1fe5guv3/Connectomics.zip?rlkey=n8x5t3prownl1l0xs5crq00a7&st=dkkk10je&dl=1');
movefile Connectomics.zip Data.zip
if ~exist('Connectomics/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','Connectomics/Data');
disp("Connectomics Dataset Downloaded and Extracted");

%% delete zip file
delete Data.zip
