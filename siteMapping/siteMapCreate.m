% siteMapCreate
%
%   siteMapCreate creates site mapping variables necessary for the
%   spike sorting pipeline. This script is mostly just a helpful template
%   that will require extensive editing with data specific to your
%   experiment.

%% Set where to record the data and what name to use

% Set directory to store the site mapping information
siteMapDir = 'C:\Users\Nicholas Timme\Box Sync\NickExperimental\Experiment1\SpikeSorting\Version1\SiteMappings';

% Set the name of the site mapping
siteMapName = 'siteMapAnim15Sept2018';

%% Record a notes variable

% This is a text string that can be used to help identify a site mapping.
% See example.

% Example: siteMapNotes = 'This is the mapping for Nicks animal 15, surgery Sept 2018, 64 channel probe (ASSY-158-P-1, PN 3379), EIB fold on animals left, headstage chip pointed towards animal back.';
siteMapNotes = 'This is the mapping for Nicks animal 15, surgery Sept 2018, 64 channel probe (ASSY-158-P-1, PN 3379), EIB fold on animals left, headstage chip pointed towards animal back.';


%% Record the probe geometry

% This variable describes the physical location in space of the probe
% electrodes. These values should be consistent for all probes of the same
% type, so long as Cambridge doesn't change its probe wiring.

% probeGeo(i,:) = [j,k] means that probe channel i has an x coordinate of j
% and a y coordinate of k (usually in um).

% ASSY-158 (64 channels) (each probe shank sorted separately)
probeGeo = NaN([64,2]);
probeGeo(1,:) = [0,75];
probeGeo(2,:) = [22.5,87.5];
probeGeo(3,:) = [0,175];
probeGeo(4,:) = [0,125];
probeGeo(5,:) = [0,150];
probeGeo(6,:) = [0,0];
probeGeo(7,:) = [0,25];
probeGeo(8,:) = [22.5,12.5];
probeGeo(9,:) = [22.5,37.5];
probeGeo(10,:) = [0,100];
probeGeo(11,:) = [22.5,187.5];
probeGeo(12,:) = [0,50];
probeGeo(13,:) = [22.5,137.5];
probeGeo(14,:) = [22.5,62.5];
probeGeo(15,:) = [22.5,162.5];
probeGeo(16,:) = [22.5,112.5];
probeGeo(17,:) = [22.5,12.5];
probeGeo(18,:) = [0,175];
probeGeo(19,:) = [0,0];
probeGeo(20,:) = [0,125];
probeGeo(21,:) = [0,25];
probeGeo(22,:) = [0,100];
probeGeo(23,:) = [0,75];
probeGeo(24,:) = [0,50];
probeGeo(25,:) = [22.5,187.5];
probeGeo(26,:) = [22.5,162.5];
probeGeo(27,:) = [22.5,87.5];
probeGeo(28,:) = [22.5,137.5];
probeGeo(29,:) = [22.5,112.5];
probeGeo(30,:) = [22.5,62.5];
probeGeo(31,:) = [22.5,37.5];
probeGeo(32,:) = [0,150];
probeGeo(33,:) = [22.5,37.5];
probeGeo(34,:) = [0,150];
probeGeo(35,:) = [22.5,112.5];
probeGeo(36,:) = [22.5,62.5];
probeGeo(37,:) = [22.5,87.5];
probeGeo(38,:) = [22.5,137.5];
probeGeo(39,:) = [22.5,187.5];
probeGeo(40,:) = [22.5,162.5];
probeGeo(41,:) = [0,75];
probeGeo(42,:) = [0,50];
probeGeo(43,:) = [0,25];
probeGeo(44,:) = [0,100];
probeGeo(45,:) = [0,0];
probeGeo(46,:) = [0,125];
probeGeo(47,:) = [22.5,12.5];
probeGeo(48,:) = [0,175];
probeGeo(49,:) = [22.5,162.5];
probeGeo(50,:) = [22.5,112.5];
probeGeo(51,:) = [22.5,137.5];
probeGeo(52,:) = [22.5,62.5];
probeGeo(53,:) = [22.5,187.5];
probeGeo(54,:) = [0,50];
probeGeo(55,:) = [22.5,37.5];
probeGeo(56,:) = [0,100];
probeGeo(57,:) = [0,25];
probeGeo(58,:) = [22.5,12.5];
probeGeo(59,:) = [0,150];
probeGeo(60,:) = [0,0];
probeGeo(61,:) = [0,175];
probeGeo(62,:) = [0,125];
probeGeo(63,:) = [0,75];
probeGeo(64,:) = [22.5,87.5];


%% Record the site mappping

% This variable records the relationship between the open ephys channel 
% numbers and the probe channel numbers. Note, the type of probe, the
% orientation of the EIB in the headcap, and the orientation of the
% headstage when plugged in will affect this. 

% Website for 64 channel Intan Headstages: http://intantech.com/RHD2164_amp_board.html
% Website for 32 channel Intan Headstages: http://intantech.com/RHD2132_RHD2216_amp_accel_board.html
% Note that the website lists channels 0-63 and 0-31. Open ephys uses these
% same channel numbers, except increased by 1 (e.g., channel 0 on the Intan
% diagram goes to channel 1 in open phys).

% openE2probe(i) = j means that open ephys channel i corresponds to probe
% channel j.

% ASSY-158 (64 channels) - fold between EIBs on animal's left, headstage 
% chip toward animal's back
openE2probe = NaN([64,1]);
openE2probe(1) = 16;
openE2probe(2) = 15;
openE2probe(3) = 14;
openE2probe(4) = 13;
openE2probe(5) = 12;
openE2probe(6) = 11;
openE2probe(7) = 10;
openE2probe(8) = 9;
openE2probe(9) = 8;
openE2probe(10) = 7;
openE2probe(11) = 6;
openE2probe(12) = 5;
openE2probe(13) = 4;
openE2probe(14) = 3;
openE2probe(15) = 2;
openE2probe(16) = 1;
openE2probe(17) = 64;
openE2probe(18) = 63;
openE2probe(19) = 62;
openE2probe(20) = 61;
openE2probe(21) = 60;
openE2probe(22) = 59;
openE2probe(23) = 58;
openE2probe(24) = 57;
openE2probe(25) = 56;
openE2probe(26) = 55;
openE2probe(27) = 54;
openE2probe(28) = 53;
openE2probe(29) = 52;
openE2probe(30) = 51;
openE2probe(31) = 50;
openE2probe(32) = 49;
openE2probe(33) = 48;
openE2probe(34) = 47;
openE2probe(35) = 46;
openE2probe(36) = 45;
openE2probe(37) = 44;
openE2probe(38) = 43;
openE2probe(39) = 42;
openE2probe(40) = 41;
openE2probe(41) = 40;
openE2probe(42) = 39;
openE2probe(43) = 38;
openE2probe(44) = 37;
openE2probe(45) = 36;
openE2probe(46) = 35;
openE2probe(47) = 34;
openE2probe(48) = 33;
openE2probe(49) = 32;
openE2probe(50) = 31;
openE2probe(51) = 30;
openE2probe(52) = 29;
openE2probe(53) = 28;
openE2probe(54) = 27;
openE2probe(55) = 26;
openE2probe(56) = 25;
openE2probe(57) = 24;
openE2probe(58) = 23;
openE2probe(59) = 22;
openE2probe(60) = 21;
openE2probe(61) = 20;
openE2probe(62) = 19;
openE2probe(63) = 18;
openE2probe(64) = 17;

%% Record the probe channel shank ID

% This variable identifies to which shank each probe channel belongs. Note
% that the shanks should be numbered 1, 2, 3, and so forth.

% probeChanShankID(i) = j means that probe channel i is on shank j

% ASSY-158 (64 channels) (each probe shank sorted separately)
probeChanShankID = NaN([64,1]);
probeChanShankID(1) = 1;
probeChanShankID(2) = 1;
probeChanShankID(3) = 1;
probeChanShankID(4) = 1;
probeChanShankID(5) = 1;
probeChanShankID(6) = 1;
probeChanShankID(7) = 1;
probeChanShankID(8) = 1;
probeChanShankID(9) = 1;
probeChanShankID(10) = 1;
probeChanShankID(11) = 1;
probeChanShankID(12) = 1;
probeChanShankID(13) = 1;
probeChanShankID(14) = 1;
probeChanShankID(15) = 1;
probeChanShankID(16) = 1;
probeChanShankID(17) = 2;
probeChanShankID(18) = 2;
probeChanShankID(19) = 2;
probeChanShankID(20) = 2;
probeChanShankID(21) = 2;
probeChanShankID(22) = 2;
probeChanShankID(23) = 2;
probeChanShankID(24) = 2;
probeChanShankID(25) = 2;
probeChanShankID(26) = 2;
probeChanShankID(27) = 2;
probeChanShankID(28) = 2;
probeChanShankID(29) = 2;
probeChanShankID(30) = 2;
probeChanShankID(31) = 2;
probeChanShankID(32) = 2;
probeChanShankID(33) = 3;
probeChanShankID(34) = 3;
probeChanShankID(35) = 3;
probeChanShankID(36) = 3;
probeChanShankID(37) = 3;
probeChanShankID(38) = 3;
probeChanShankID(39) = 3;
probeChanShankID(40) = 3;
probeChanShankID(41) = 3;
probeChanShankID(42) = 3;
probeChanShankID(43) = 3;
probeChanShankID(44) = 3;
probeChanShankID(45) = 3;
probeChanShankID(46) = 3;
probeChanShankID(47) = 3;
probeChanShankID(48) = 3;
probeChanShankID(49) = 4;
probeChanShankID(50) = 4;
probeChanShankID(51) = 4;
probeChanShankID(52) = 4;
probeChanShankID(53) = 4;
probeChanShankID(54) = 4;
probeChanShankID(55) = 4;
probeChanShankID(56) = 4;
probeChanShankID(57) = 4;
probeChanShankID(58) = 4;
probeChanShankID(59) = 4;
probeChanShankID(60) = 4;
probeChanShankID(61) = 4;
probeChanShankID(62) = 4;
probeChanShankID(63) = 4;
probeChanShankID(64) = 4;


%% Plot the shank geometry for user verification

nShanks = length(unique(probeChanShankID));
for iShank = 1:nShanks
    figure
    hold on
    chans = find(probeChanShankID == iShank);
    xMin = inf;
    xMax = -inf;
    yMin = inf;
    yMax = -inf;
    for iChan = 1:length(chans)
        text(probeGeo(chans(iChan),1),probeGeo(chans(iChan),2),num2str(chans(iChan)),'HorizontalAlignment','center','VerticalAlignment','middle')
        xMin = min([xMin,probeGeo(chans(iChan),1)]);
        xMax = max([xMax,probeGeo(chans(iChan),1)]);
        yMin = min([yMin,probeGeo(chans(iChan),2)]);
        yMax = max([yMax,probeGeo(chans(iChan),2)]);
    end
    dx = xMax - xMin;
    dy = yMax - yMin;
    xlim([xMin - 0.1*dx,xMax + 0.1*dx])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    title(['Shank ',num2str(iShank)])
end


%% Save the mapping information

save([siteMapDir,filesep,siteMapName],'siteMapNotes','probeGeo','openE2probe','probeChanShankID')
