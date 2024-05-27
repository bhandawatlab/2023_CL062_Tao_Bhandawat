%% Download Datasets
%% download and save the walking chamber dataset
websave('WalkingChamber.zip','https://www.dropbox.com/scl/fi/2ttw3pfich0b2l16xaeyv/WalkingChamber.zip?rlkey=xt0i2a7d1ej5ojfqdn6e7tbs3&st=ebn5opas&dl=1');
movefile WalkingChamber.zip Data.zip
if ~exist('WalkingChamberCodeBase/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','WalkingChamberCodeBase/Data');
disp("WalkingChamber Dataset Downloaded and Extracted");

%% download and save the DMD dataset
websave('DMD.zip','https://www.dropbox.com/scl/fi/7125p97desuoz2bpsspob/DMD.zip?rlkey=wlrvsiriwqnmcu6ansxus45lw&st=ig83f9p5&dl=1');
movefile DMD.zip Data.zip
if ~exist('DMDCodeBase/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','DMDCodeBase/Data')
disp("DMD Dataset Downloaded and Extracted");

%% download and save the connectomics dataset
websave('Connectomics.zip','https://www.dropbox.com/scl/fi/3xi95nxab0wow59rtm8s2/Connectomics.zip?rlkey=6pf7mocsym6iw0k0qtch67h24&st=ozdelvmv&dl=1');
movefile Connectomics.zip Data.zip
if ~exist('Connectomics/Data', 'dir')
    mkdir(yourFolder)
end
unzip('Data.zip','Connectomics/Data');
disp("Connectomics Dataset Downloaded and Extracted");

%% delete zip file
delete Data.zip