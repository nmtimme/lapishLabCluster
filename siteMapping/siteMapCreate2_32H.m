% msClusterMapCreateVer1
%
%   msClusterMapCreate creates site mapping variables necessary for the
%   spike sorting pipeline. This script is mostly just a helpful template
%   that will require extensive editing with data specific to your
%   experiment.

%% Set where to record the data and what name to use

% Set directory to store the site mapping information
siteMapDir = 'C:\Users\Nicholas Timme\Box Sync\NickExperimental\Experiment1\SpikeSorting\batch1\SiteMappings';

% Set the name of the site mapping
siteMapName = 'siteMap2_32H';

%% Record a notes variable

% This is a text string that can be used to help identify a site mapping.
% See example.

% Example: siteMapNotes = 'This is the mapping for Nicks animal 15, surgery Sept 2018, 64 channel probe (ASSY-158-P-1, PN 3379), EIB fold on animals left, headstage chip pointed towards animal back.';
siteMapNotes = 'This is a site mapping for an animal with two 32 channel ASSY-158-H4 probes. The probes are implanted with one omnetics print towards animal front and one omnetics print towards animal back. Headstage chip pointed towards animal back. Each probe sorted separately.';



%% Record the probe geometry

% This variable describes the physical location in space of the probe
% electrodes. These values should be consistent for all probes of the same
% type, so long as Cambridge doesn't change its probe wiring.

% probeGeo(i,:) = [j,k] means that probe channel i has an x coordinate of j
% and a y coordinate of k (usually in um).

% 2 ASSY-158-H4 (64 channels) (each probe shank sorted separately, back probe number increased by 32 on the probe sheet)
probeGeo = NaN([64,2]);
probeGeo(1,:) = [0,275];
probeGeo(2,:) = [0,150];
probeGeo(3,:) = [0,375];
probeGeo(4,:) = [0,325];
probeGeo(5,:) = [0,350];
probeGeo(6,:) = [0,0];
probeGeo(7,:) = [0,50];
probeGeo(8,:) = [0,25];
probeGeo(9,:) = [0,75];
probeGeo(10,:) = [0,300];
probeGeo(11,:) = [0,250];
probeGeo(12,:) = [0,100];
probeGeo(13,:) = [0,200];
probeGeo(14,:) = [0,125];
probeGeo(15,:) = [0,225];
probeGeo(16,:) = [0,175];
probeGeo(17,:) = [0,625];
probeGeo(18,:) = [0,575];
probeGeo(19,:) = [0,600];
probeGeo(20,:) = [0,525];
probeGeo(21,:) = [0,650];
probeGeo(22,:) = [0,500];
probeGeo(23,:) = [0,475];
probeGeo(24,:) = [0,700];
probeGeo(25,:) = [0,450];
probeGeo(26,:) = [0,425];
probeGeo(27,:) = [0,750];
probeGeo(28,:) = [0,400];
probeGeo(29,:) = [0,775];
probeGeo(30,:) = [0,725];
probeGeo(31,:) = [0,675];
probeGeo(32,:) = [0,550];
probeGeo(1+32,:) = [0,275];
probeGeo(2+32,:) = [0,150];
probeGeo(3+32,:) = [0,375];
probeGeo(4+32,:) = [0,325];
probeGeo(5+32,:) = [0,350];
probeGeo(6+32,:) = [0,0];
probeGeo(7+32,:) = [0,50];
probeGeo(8+32,:) = [0,25];
probeGeo(9+32,:) = [0,75];
probeGeo(10+32,:) = [0,300];
probeGeo(11+32,:) = [0,250];
probeGeo(12+32,:) = [0,100];
probeGeo(13+32,:) = [0,200];
probeGeo(14+32,:) = [0,125];
probeGeo(15+32,:) = [0,225];
probeGeo(16+32,:) = [0,175];
probeGeo(17+32,:) = [0,625];
probeGeo(18+32,:) = [0,575];
probeGeo(19+32,:) = [0,600];
probeGeo(20+32,:) = [0,525];
probeGeo(21+32,:) = [0,650];
probeGeo(22+32,:) = [0,500];
probeGeo(23+32,:) = [0,475];
probeGeo(24+32,:) = [0,700];
probeGeo(25+32,:) = [0,450];
probeGeo(26+32,:) = [0,425];
probeGeo(27+32,:) = [0,750];
probeGeo(28+32,:) = [0,400];
probeGeo(29+32,:) = [0,775];
probeGeo(30+32,:) = [0,725];
probeGeo(31+32,:) = [0,675];
probeGeo(32+32,:) = [0,550];


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

% ASSY-158-F (64 channels) - fold between EIBs on animal's left, headstage 
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

% ASSY-158-F (64 channels) (each probe shank sorted separately)
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
probeChanShankID(17) = 1;
probeChanShankID(18) = 1;
probeChanShankID(19) = 1;
probeChanShankID(20) = 1;
probeChanShankID(21) = 1;
probeChanShankID(22) = 1;
probeChanShankID(23) = 1;
probeChanShankID(24) = 1;
probeChanShankID(25) = 1;
probeChanShankID(26) = 1;
probeChanShankID(27) = 1;
probeChanShankID(28) = 1;
probeChanShankID(29) = 1;
probeChanShankID(30) = 1;
probeChanShankID(31) = 1;
probeChanShankID(32) = 1;
probeChanShankID(33) = 2;
probeChanShankID(34) = 2;
probeChanShankID(35) = 2;
probeChanShankID(36) = 2;
probeChanShankID(37) = 2;
probeChanShankID(38) = 2;
probeChanShankID(39) = 2;
probeChanShankID(40) = 2;
probeChanShankID(41) = 2;
probeChanShankID(42) = 2;
probeChanShankID(43) = 2;
probeChanShankID(44) = 2;
probeChanShankID(45) = 2;
probeChanShankID(46) = 2;
probeChanShankID(47) = 2;
probeChanShankID(48) = 2;
probeChanShankID(49) = 2;
probeChanShankID(50) = 2;
probeChanShankID(51) = 2;
probeChanShankID(52) = 2;
probeChanShankID(53) = 2;
probeChanShankID(54) = 2;
probeChanShankID(55) = 2;
probeChanShankID(56) = 2;
probeChanShankID(57) = 2;
probeChanShankID(58) = 2;
probeChanShankID(59) = 2;
probeChanShankID(60) = 2;
probeChanShankID(61) = 2;
probeChanShankID(62) = 2;
probeChanShankID(63) = 2;
probeChanShankID(64) = 2;


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
    if dx == 0
        dx = 1;
    end
    dy = yMax - yMin;
    if dy == 0
        dy = 1;
    end
    xlim([xMin - 0.1*dx,xMax + 0.1*dx])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    title(['Shank ',num2str(iShank)])
end


%% Save the mapping information

save([siteMapDir,filesep,siteMapName],'siteMapNotes','probeGeo','openE2probe','probeChanShankID')
