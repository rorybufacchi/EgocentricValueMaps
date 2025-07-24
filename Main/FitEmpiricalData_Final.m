% -----
% USE:
% 
% Run this script to compare various models to empirical data. The models
% are run in 'CreateModelDataToFit_Final.m'
% -----


% First make sure the path is set right though
cFileName       = matlab.desktop.editor.getActiveFilename;
codePath        = cFileName(1:end-29);
addpath(        genpath(codePath));



%% ------------------------------------------------------------------------
%  LOAD Q-values to be fitted to data, for Fig7, and Extended Data Figd 4&5

try
    load('Results\ForFigures\Fig6_7_ModelledData_Precomputed.mat')
catch
    warning('No Precomputed data found. Either Compute results using CreateModelDataToFit_Final.m, or download precompted results from $')
end



%% ------------------------------------------------------------------------
% Create big table d with all empirical data
%
%
% USAGE:
% Each row of d is another set of data
% Each 'd.exp' is fitted with a separate offset, and each 'd.tact' is
% fitted with a separate slope
%
% So, for example, if you want the data in row 1 to be fit with a different 
% offset but the same slope as row 2, you should set:
% d.exp(1) = 1; d.exp(2) = 2;
% d.tact{1} = 'locationA'; d.tact{2} = 'LocationA';
%

% Initialise settings
sBdy.clc.binW    = 5;
sHnd.clc.binW    = 5;
sHndSide.clc.binW= 5;
sHed.clc.binW    = 5;

d = table;

% Put row counter first so it's easier to eyeball things
d.iD(1) = 1;

% Options for how to calculate overall map activity from the set of Q-values
d.psiSettings{1,1} = 'max_then_avg'; % Other options: 'average'; 'average_all'; 'max_above_stationary';
d.psiSettings{1,2} = 'no_norm'; % Other options: %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';
d.psiSettings{1,3} = ''; % This is an additional setting to allow simulating macaque data. In that case it is set to 'avOverRows' (so no need to change here)

% ==============================
% APPROACHING TRUNK 1
cD               = 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
% Set positions of stimuli (rows are euclidian coordinates from bodypart)
d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
                    0,  0,  0,  0,  0,   0 ; ...
                    0,  0,  0,  0,  0,   0  ];
d = AddbinPos(d,sBdy,cD); % approaching body
% Create real values (Taken directly from published figure)
d.realDat{cD} = [0.646 0.755 0.568 0.313 0.028 -0.213] .* 73.39449541;
d.realSTE{cD} = abs([0.474 0.530 0.396 0.038 -0.227 -0.484] .* 73.39449541 - d.realDat{cD})  ;
d.exp(cD)     = 1;
% Define which successor features are part of the egocentric map
% For modelling trunk towards, we take the trunk, head, and hand-by-side
% Average over body-part related q values
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{cD}}}; 

% ==============================
% RECEDING TRUNK 2
cD               = cD + 1;
d(cD,:)          = d(cD-1,:);
d.dir(cD)        = -1;
% Create real values (Taken directly from published figure)
d.realDat{cD}    = [0.183 0.211 0.364 0.479 0.394 0.359] .* 73.39449541;
d.realSTE{cD}    = abs([0.069 -0.039 0.131 0.252 0.207 0.174] .* 73.39449541 - d.realDat{cD})  ;
d.exp(cD)        = 1;
% Define which successor features are part of the egocentric map
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{-1}}};

% ==============================
% APPROACHING HAND ONLY 3
cD               = cD + 1;
d(cD,:)          = d(cD-1,:);
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandSide';
d.tact{cD}       = 'Hand';
d.cmPos{cD}      = [ 5  27, 50, 71, 93 ; ...
                    0,  0,  0,  0,  0, ; ...
                    0,  0,  0,  0,  0, ];
d = AddbinPos(d,sHndSide,cD); % approaching hand on side
% Create real values (Taken directly from published figure)
d.realDat{cD} = [1.044 0.969 0.547 0.285 -0.234] .* 45.62737643;
d.realSTE{cD} = abs([0.831 0.711 0.231 -0.021 -0.516] .* 45.62737643 - d.realDat{cD})  ;
d.exp(cD)     = 2;
% Define which successor features are part of the egocentric map
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% ==============================
% RECEDING HAND ONLY 4
cD            = cD + 1;
d(cD,:)       = d(cD - 1,:);
d.dir(cD)     = -1;
d.realDat{cD} = [0.662 0.624 0.360 0.175 0.208] .* 45.62737643;
d.realSTE{cD} = abs([0.409 0.372 0.078 0.014 0.012] .* 45.62737643 - d.realDat{cD})  ;
d.exp(cD)     = 2;
% Define which successor features are part of the egocentric map
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{-1}}};


% #################################################################
% ==============================
% APPROACHING TRUNK WITH HAND ON TRUNK, TRUNK STIMULATED 5
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
                    0,  0,  0,  0,  0,   0 ; ...
                    0,  0,  0,  0,  0,   0  ];
d = AddbinPos(d,sBdy,cD); % approaching body
% Create real values
d.realDat{cD} = [0.983 0.865 0.751 0.490 0.440 0.046] .* 56.02240896;
d.realSTE{cD} = abs([0.777 0.508 0.404 0.095 0.182 0.334] .* 56.02240896 - d.realDat{cD})  ;
d.exp(cD)     = 3;
% Define which successor features are part of the egocentric map
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, ...
    {'dir' ,{1}}}; 

% ==============================
% APPROACHING HAND WITH HAND ON TRUNK, HAND STIMULATED 6
cD               = cD + 1; handOnTrunkApproachHandStimD = cD;
d(cD,:)          = d(cD - 1,:); % (data is the same as for condition 5 above)
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'Hand';

% ==============================
% RECEDING TRUNK WITH HAND ON TRUNK, TRUNK STIMULATED 7
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.dir(cD)        = -1;
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
% Create real values
d.realDat{cD} = [0.425 0.532 0.579 0.581 0.702 0.362] .* 56.02240896;
d.realSTE{cD} = abs([0.148 0.200 0.219 0.184  0.344 0.012] .* 56.02240896 - d.realDat{cD})  ;
d.exp(cD)     = 3;
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, ...
    {'dir' ,{-1}}}; 

% ==============================
% RECEDING HAND WITH HAND ON TRUNK, HAND STIMULATED 8
cD               = cD + 1;
d(cD,:)          = d(cD-1,:);
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'Hand';

% #################################################################
% ==============================
% APPROACHING BETWEEN HAND AND TRUNK WITH HAND STIMULATED AND HAND AWAY
% FROM CHEST 9
cD              = cD + 1;
d(cD,:)         = d(handOnTrunkApproachHandStimD,:);
% auditory stimulus position offset to the participants' right
d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
                     0,  0,  0,  0,  0,   0 ; ...
                     0,  0,  0,  0,  0,   0  ] + [0 20 0]';
d = AddbinPos(d,sBdy,cD); % approaching body [ but to the side, as defined above]
d.tactloc{cD}    = 'HandSide';
d.tact{cD}       = 'Hand';
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 
% Create real values
d.realDat{cD} = [0.566 0.665 0.046 0.089 -0.303 -0.647] .* 56.11672278;
d.realSTE{cD} = abs([0.415 0.536 -0.143 -0.035 -0.447 -0.839] .* 56.11672278 - d.realDat{cD})  ;
d.exp(cD)     = 4;

% ==============================
% APPROACHING BETWEEN HAND AND TRUNK WITH HAND STIMULATED AND HAND NEAR/ON
% CHEST 10
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'Hand';
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, ...
    {'dir' ,{1}}};
% Create real values
d.realDat{cD} = [0.730 0.694 0.412 -0.006 -0.227 -0.343] .* 56.11672278;
d.realSTE{cD} = abs([0.583 0.538 0.240 -0.182 -0.472 -0.536] .* 56.11672278 - d.realDat{cD})  ;
d.exp(cD)     = 4;

% #################################################################
% ==============================
% APPROACHING FACE WITH FACE STIMULATED 11
cD               = cD + 1;  faceApproachFaceStimD = cD;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'Head';
d.cmPos{cD}      = [ 5, 37, 69, 101, 133, 165, 197 ; ...
                     0,  0,  0,   0,   0,   0,   0 ; ...
                     0,  0,  0,   0,   0,   0,   0];
d = AddbinPos(d,sHed,cD); % approaching head
% Create real values 
d.realDat{cD} = [1.043 0.785 0.566 0.599 0.493 0.543 0.424] .* 63.84676776;
d.realSTE{cD} = abs([0.855 0.629 0.409 0.450 0.331 0.372 0.242] .* 63.84676776 - d.realDat{cD})  ;
d.exp(cD)     = 5;
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% ============================== 12
% APPROACHING BODY WITH FACE STIMULATED
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'Head';
d = AddbinPos(d,sBdy,cD); % approaching body
% Create real values
d.realDat{cD} = [0.38 0.366 0.463 0.314 0.446 0.314 0.394] .* 63.84676776;
d.realSTE{cD} = abs([0.131 0.129 0.252 0.106 0.214 0.084 0.190] .* 63.84676776 - d.realDat{cD}) ;

% ==============================
% APPROACHING BODY WITH BODY STIMULATED 13
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d = AddbinPos(d,sBdy,cD); % approaching body
% Create real values
d.realDat{cD} = [0.959 0.801 0.628 0.433 0.440 0.438 0.358] .* 51.02040816;
d.realSTE{cD} = abs([0.818 0.649 0.467 0.283 0.243 0.259 0.215] .* 51.02040816 - d.realDat{cD})  ;
d.exp(cD)     = 6;

% ==============================
% APPROACHING FACE WITH BODY STIMULATED 14
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d = AddbinPos(d,sHed,cD); % approaching Head
% Create real values
d.realDat{cD} = [0.609 0.539 0.354 0.376 0.189 0.251 0.194] .* 51.02040816;
d.realSTE{cD} = abs([0.351 0.386 0.117 0.120 -0.033 0.006 -0.020] .* 51.02040816 - d.realDat{cD})  ;


% #################################################################
% ==============================
% VISUAL APPROACHING FACE WITH TRUNK STIMULATED 15
cD               = cD + 1;
d(cD,:)          = d(faceApproachFaceStimD ,:); 
% Put position a bit more inbetween head and body for the visual one
d.cmPos{cD}      = [ 5, 37, 69, 101, 133, 165, 197 ; ...
                     0,  0,  0,   0,   0,   0,   0 ; ...
                   -10,-10,-10, -10, -10, -10, -10];
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d = AddbinPos(d,sHed,cD); % approaching Head [ but slightly lower because of cmpos above]
% Create real values
d.realDat{cD} = [1.278 1.000 0.639 0.531 0.380 0.295 -0.073] .* 48.95104895;
d.realSTE{cD} = abs([1.118 0.838 0.363 0.335 0.222 0.085 -0.271] .* 48.95104895 - d.realDat{cD})  ;
d.exp(cD)     = 7;
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% ==============================
% VISUAL APPROACHING FACE WITH FACE STIMULATED 16
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'Head';
% Create real values
d.realDat{cD} = [1.096 0.948 0.092 0.356 0.023 -0.083 -0.114] .* 48.95104895;
d.realSTE{cD} = abs([0.944 0.794 -0.184 0.169 -0.135 -0.300 -0.305] .* 48.95104895 - d.realDat{cD})  ;


% #################################################################
% Create Holmes et al tool data

% ==============================
% NO TOOL for TOOL EXPERIMENT 17
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandCCE';
d.cmPos{cD}      = [5, 45, 70 ; ...
                    0,  0,  0 ; ...
                    0,  0,  0 ];
d = AddbinPos(d,sBdy,cD); % approaching body
% Create real values
d.realDat{cD} = [1.714 0.641 0.342] .* 87.95074758;
% Transform from 95% confidence interval to SE
d.realSTE{cD} = ((abs([2.104 0.972 0.677] .* 87.95074758 - d.realDat{cD})) ./ 47.5) .* 34.1 ;
d.exp(cD)     = 8;
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}}; 


% ==============================
% YES TOOL for TOOL EXPERIMENT 18
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
d.realDat{cD} = [1.587 0.410 0.799] .* 87.95074758;
% Transform from 95% confidence interval to SE
d.realSTE{cD} = ((abs([1.911 0.675 1.071] .* 87.95074758 - d.realDat{cD})) ./ 47.5) .* 34.1 ;
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestPlusTool'}}}; 

% #################################################################
% Create Huijsmans data: spider on a track
% CREATE data table

% ==============================
% BUTTERFLY on track towards hand 19
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandSpider';
d.cmPos{cD}      = [ 5, 20, 35, 50, 65, 80 ; ...
                    10, 10, 10,  0,  0,  0 ; ...
                     0,  0,  0, 0,  0,  0 ] + [0 0 0]';
d = AddbinPos(d,sTrack,cD); % approaching hand on chest
% Create real values
d.realDat{cD} = [1.487 1.830 2.014 1.924 2.754 2.935].* 19.17545542;
d.realSTE{cD} = abs([0.735 1.047 1.256 1.148 1.961 2.172] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestTrack'}} , ...
    {'rew',{1,-1}} }; 

% ==============================
% BUTTERFLY on track away from hand 20
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.cmPos{cD}      = [ 5, 20, 35, 50, 65, 80 ; ...
                    -10, -10, -10,  0,  0,  0 ; ...
                     0,  0,  0, 0,  0,  0 ] + [0 0 0]';
d = AddbinPos(d,sTrack,cD); % approaching hand on chest
% Create real values
d.realDat{cD} = [1.487 1.830 2.014 1.896 2.292 2.427].* 19.17545542;
d.realSTE{cD} = abs([0.735 1.047 1.256 1.092 1.470 1.637] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;


% ==============================
% SPIDER on track towards hand 21
cD               = cD + 1;
d(cD,:)          = d(cD - 2,:);
% % % d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestTrack'}} , ...
% % %     {'rew',{1,-1.25}} }; 
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestTrack'}} , ...
    {'rew',{1,-1}} }; 

% Create real values
d.realDat{cD} = [1.105 1.730 2.129 2.219 3.282 3.164] .* 19.17545542;
d.realSTE{cD} = abs([0.322 0.926 1.391 1.478 2.533 2.27] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;


% ==============================
% SPIDER on track away from hand 22
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.cmPos{cD}      = [ 5, 20, 35, 50, 65, 80 ; ...
                    -10, -10, -10,  0,  0,  0 ; ...
                     0,  0,  0, 0,  0,  0 ] + [0 0 0]';
d = AddbinPos(d,sTrack,cD); % approaching hand on chest
% Create real values
d.realDat{cD} = [1.105 1.730 2.129 1.986 2.608 3.045] .* 19.17545542;
d.realSTE{cD} = abs([0.322 0.926 1.391 1.224 1.797 2.161] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;



% #################################################################
% Create Wamain data: mu power during VR object judgements

% ==============================
% PROTOTYPICAL OBJECT judgements - SCRAMBLED OBJECT judgements 23
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandSide';
d.tact{cD}       = 'HandReachEEG';
d.cmPos{cD}      = [20 87.5 160 ; ...
                     0,  0,   0 ; ...
                     0,  0,   0 ];
d = AddbinPos(d,sHndSide,cD); % approaching hand on chest
% Create real values
tmpD          = [0.926 0.617 0.112; 0.284 0.000 0.201 ; ...
                 0.520 0.481 0.700; 0.666 0.444 0.553];
d.realDat{cD} = (tmpD(1,:) - tmpD(2,:) - ( tmpD(3,:) - tmpD(4,:) )) .* 0.138657793;
tmpErr        = abs([0.494 0.153 -0.260 ; -0.068 -0.363 -0.169; ...
                     0.119 0.094 0.292 ; 0.357 0.109 0.251] - tmpD);
d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.138657793) ./sqrt(17);
d.exp(cD)     = 10;
d.psiSplit{cD} = { {'bodyPart',{'HandBySide'}} , ...
    {'rew',{1}} }; 



% #################################################################
% Create Ronga data: alpha power during tool use

% ==============================
% AFETR COGNITIVE TRAINING: alpha div beta, to baseline because beta is
% unchanged 24
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandToolEEG';
d.cmPos{cD}      = [5, 145 ; ...
                    0,  0 ; ...
                    0,  0 ];
d = AddbinPos(d,sHndToolRake,cD); % approaching hand on chest
% Create real values: divisive normalisation: Account for an OVERALL effect
% of training type on ALL power
tmpD          = [1.748 1.491] ;
d.realDat{cD} = tmpD(1,:) .* 0.343524562  ./ 0.556;
tmpErr        = abs([1.946 1.689] - tmpD)  ;
d.realSTE{cD} = tmpErr .* 0.343524562 ./ 0.556 ;
d.exp(cD)     = 11;
d.psiSplit{cD} = { {'bodyPart',{'HandByChest'}} , ...
    {'rew',{1,-1}} }; 


% ==============================
% AFETR TOOL TRAINING: alpha div beta, to baseline because beta is
% unchanged 25
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values: divisive normalisation: Account for an OVERALL effect
% of training type on ALL power
tmpD          = [2.407 2.403] ;
d.realDat{cD} = tmpD(1,:) .* 0.343524562 ./ 0.661;
tmpErr        = abs([2.933 2.981] - tmpD) ;
d.realSTE{cD} = tmpErr .* 0.343524562 ./ 0.661 ;
d.exp(cD)     = 11;
d.psiSplit{cD} = { {'bodyPart',{'HandByChestPlusTool'}} , ...
    {'rew',{1,-1}} }; 


% #################################################################
% Create Quinlan dPOS data: fMRI moving vs stationary

cD = 25;

% ==============================
% AFETR COGNITIVE TRAINING: alpha div beta, to baseline because beta is
% unchanged 26
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'HeadMRIdPOS';
d.cmPos{cD}      = [15, 38, 84 ; ...
                     0,  0,  0 ; ...
                     0,  0,  0 ];
d = AddbinPos(d,sHed,cD); % approaching hand on chest
% Create real values: divisive normalisation: Account for an OVERALL effect
% of training type on ALL power
tmpD          = [2.237 1.661 1.036] ;
d.realDat{cD} = tmpD(1,:) .* 0.375335121 ;
tmpErr        = abs([3.557 2.675 1.696] - tmpD)  ;
% Transform from 95% confidence interval to SE
d.realSTE{cD} = ((abs([2.104 0.972 0.677] .* 87.95074758 - d.realDat{cD})) ./ 47.5) .* 34.1 ;
d.realSTE{cD} = tmpErr .* 0.375335121 .* (68.2/99.17) ./sqrt(18);
d.exp(cD)     = 12;
d.psiSplit{cD} = { {'bodyPart',{'Head'}} , ...
    {'rew',{1,-1}} }; 


% #################################################################
% Create Holt DIPS PMV data: fMRI approach - withdrawal difference, faces cars spheres
% CREATE data table

% ==============================
% LOOMING - RECEDING FACE, CAR, SPHERE DIPS 27
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'HeadMRIdDIPS';
d.cmPos{cD}      = [20, 150, 150; ...
                     0,   0,   0; ...
                     0,   0,   0];
d = AddbinPos(d,sHed,cD); % approaching hand on chest
% Create real values
tmpD          = [1.866 1.616 1.262; 0.207 1.367 1.315];
d.realDat{cD} = (tmpD(1,:) - tmpD(2,:)) .* 0.049979175;
tmpErr        = abs([2.231 1.911  1.543 ; 0.421  1.456 1.421] - tmpD);
d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.049979175) ;
d.exp(cD)     = 13;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'Head'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 



% ==============================
% LOOMING - RECEDING FACE, CAR, SPHERE PMv 28
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'HeadMRIdDIPS';
d.cmPos{cD}      = [20, 150, 150; ...
                     0,   0,   0; ...
                     0,   0,   0];
d = AddbinPos(d,sHed,cD); % approaching hand on chest
% Create real values
tmpD          = [1.174 0.558 0.661; -0.056 0.548 0.652];
d.realDat{cD} = (tmpD(1,:) - tmpD(2,:)) .* 0.049979175;
tmpErr        = abs([1.546 0.914 1.016; 0.126 0.796 0.935] - tmpD);
d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.049979175) ;
d.exp(cD)     = 13;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'Head'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% #################################################################
% Create Graziano single neuron responses: ARM

% ==============================
% ARM STRAIGHT, ARM RESPONSES 29
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Arm';
d.tact{cD}       = 'Macaque';
d.cmPos{cD}      = [ 20,  20,  20,  20; ...
                    -30,   0,  30,  60; ...
                     10,  10,  10,  10];
d = AddbinPos(d,sArm,cD); % approaching arm
% Create real values
tmpD          = [1.077 1.525 2.593 1.683];
d.realDat{cD} = (tmpD(1,:)) .* 36.954915;
tmpErr        = abs([1.244 1.709 2.652 1.900] - tmpD);
d.realSTE{cD} = tmpErr .* 36.954915 ;
d.exp(cD)     = 14;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'ArmForward'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% ==============================
% ARM LEFT, ARM RESPONSES 30
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
tmpD          = [1.575 2.069 1.682 1.219];
d.realDat{cD} = (tmpD(1,:)) .* 36.954915;
tmpErr        = abs([1.773 2.215 1.825 1.458] - tmpD);
d.realSTE{cD} = tmpErr .* 36.954915 ;
d.exp(cD)     = 14;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'ArmLeft'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% #################################################################
% Create Graziano single neuron responses: HEAD 

% ==============================
% HEAD STRAIGHT, HEAD RESPONSES 31
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'MacaqueHead';
% Stimulus positions
d.cmPos{cD}      = { 5.*[4:13], 5.*[4:15], 5.*[4:16], 5.*[4:15], 5.*[4:13]; ...
                     5.*[4 5 5 6 7 7 8 9 9 10 ], 5.*[ 2 3 3 3 4 4 4 5 5 5 6 6] , 5.*zeros([1 13]) , 5.*-[ 2 3 3 3 4 4 4 5 5 5 6 6 ] , 5.*-[ 4 5 5 6 7 7 8 9 9 10 ]; ... 
                     zeros([1 10]) ,  zeros([1 12]) ,  zeros([1 13]) ,  zeros([1 12]) , zeros([1 10]) };
d = AddbinPos(d,sHed,cD); % approaching arm
% Create real values
tmpD          = [2.187 2.315 4.551 2.742 1.919];
d.realDat{cD} = (tmpD(1,:)) .* 20.55921053;
tmpErr        = abs([2.763 2.635 4.756 2.971 2.184] - tmpD);
d.realSTE{cD} = tmpErr .* 20.55921053 ;
d.exp(cD)     = 15;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'HeadConstr'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 



% ==============================
% HEAD ROTATED, HEAD RESPONSES 32
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
tmpD          = [1.528 1.910 2.816 3.868 2.358];
d.realDat{cD} = (tmpD(1,:)) * 20.55921053;
tmpErr        = abs([1.757 2.146 3.070 4.153 2.580] - tmpD);
d.realSTE{cD} = tmpErr * 20.55921053;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'RotatedHead'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 
% Use the same settings for calculating all psi [just with different
% features]
for iD = 2:size(d,1)
    d.psiSettings(iD,:) = d.psiSettings(1,:);
end

% Change settings for the experiments in which Q-values need to be averaged 
% over rows
d.psiSettings{cD-3,3} = 'avOverRows'; % Macaque single neurons
d.psiSettings{cD-2,3} = 'avOverRows';
d.psiSettings{cD-1,3} = 'avOverRows'; 
d.psiSettings{cD,3}   = 'avOverRows';



% #################################################################
% Create Ferri RT data: Sounds approaching, neutral, negative and positive
% valence


% ==============================
% LOOMING Neutral sound (white noise + Exp 2) 33
cD               = cD + 1;
d(cD,:)          = d(6,:);
d.dir(cD)        = 1;
% % % d.tactloc{cD}    = 'HandSide';
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandSoundValence';
d.cmPos{cD}      = [ 5, 12, 19, 26, 33, 40, 47, 54, 60, 66 ; ...
                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ; ...
                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ];
d = AddbinPos(d,sBdy,cD); % approaching hand with hand on chest
% Create real values
tmpD =  mean([4.220 3.699 2.886 2.886 2.803 2.060 0.734 1.095 1.053 0.393 ; 2.129 2.588 3.150 2.588 2.859 2.657 1.643 1.490 1.456 0.747]);
d.realDat{cD} = tmpD .*  25.54;
% No error bars were reported for this expt, so we copy from other similar
% experiment
similarExp    = 6; % Hand stimulated in the exact position that I assume for this experiment
tmpErr        = nanmean(arrayfun(@(expNum) nanmean(d.realSTE{expNum}), similarExp));
d.realSTE{cD} = tmpErr .* ones(size(d.realDat{cD})) ;
d.exp(cD)     = 16;
% Create psisplit 
d.psiSplit{cD} = { {'bodyPart',{'Trunk','HandByChest'}} , ...
    {'dir',{1}}, {'rew',{1,-1}} }; 
d.psiSettings{cD,4} = {'optimise_reward'}; % Set the reward to be optimised
d.psiSettings{cD,5} = {[2]}; % This is which of the reward magnitudes will be fitted to the data

% ==============================
% LOOMING Negative sound (brown noise + Exp 2) 34
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real value
tmpD = mean( [4.109 3.171 4.076 3.720 2.511 2.310 1.768 1.102 1.261 0.560; 4.004 2.865 3.463 3.178 3.609 2.539 2.504 2.185 1.296 0.574 ]);
d.realDat{cD} = tmpD .*  25.54;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'Trunk','HandByChest'}} , ...
    {'dir',{1}}, {'rew',{1,-1}} }; 
d.psiSettings{cD,4} = {'optimise_reward'}; % Set the reward to be optimised
d.psiSettings{cD,5} = {[2]}; % This is which of the rewards to allow fitting for


% ==============================
% LOOMING Positive sound (Exp 2) 35
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
d.realDat{cD} = [3.296 3.025 2.671 1.997 2.740 2.227 2.428 1.602 0.907 1.011] .*  25.54;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'Trunk','HandByChest'}} , ...
    {'dir',{1}} , {'rew',{1,-1}} }; 
d.psiSettings{cD,4} = {'optimise_reward'};
d.psiSettings{cD,5} = {[2]}; % This is which of the rewards to allow fitting for


% #################################################################
% Create Taffou 2014 RT data: Auditory approaching, negative and positive
% valence. 
% NOTE: Need larger binds because stimuli were at further distance than
% in the other experiments, and moved faster
sBack.clc.binW   = 10;

% ==============================
% APPROACHING BACK, NEUTRAL 36
cD               = 36;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Back';
d.tact{cD}       = 'Back';
tmpDists         = [500:-125:0] .* cos(pi/4) + 20 .* cos(pi/4);
d.cmPos{cD}      = [tmpDists ; ...
                    tmpDists; ...
                    0,  0,  0,  0,  0 ];
d = AddbinPos(d,sBack,cD); % approaching body
% Create real values
d.realDat{cD} = ([2.395 2.121 1.909 1.308 1.082] - 0.346) .* 27.7585;
d.realSTE{cD} = ([2.846 2.534 2.371 1.798 1.511] - 0.346) .* 27.7585 - d.realDat{cD}  ;
% Change it to RT speeding (instead of overall RT)
d.realDat{cD} = 80 - d.realDat{cD};
d.exp(cD)     = 17;
d.psiSplit{cD} = { {'bodyPart',{'Back'}}, ...
    {'dir' ,{1}}, {'rew',{1,-1}}}; 
d.psiSettings{cD,1} = 'max_then_avg';
d.psiSettings{cD,2} = 'no_norm';
d.psiSettings{cD,3} = '';
d.psiSettings{cD,4} = {'optimise_reward'};
d.psiSettings{cD,5} = {[2]}; % This is which of the rewards to allow fitting for

% ==============================
% APPROACHING BACK, THREAT 37
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
d.realDat{cD} = ([2.565 1.676 1.596 1.107 1.036] - 0.346) .* 27.7585;
d.realSTE{cD} = ([3.075 2.162 2.006 1.541 1.463] - 0.346) .* 27.7585 - d.realDat{cD}  ;
% Change it to RT speeding (instead of overall RT)
d.realDat{cD} = 80 - d.realDat{cD};
d.exp(cD)     = 17;
d.psiSplit{cD} = { {'bodyPart',{'Back'}}, ...
    {'dir' ,{1}}, {'rew',{1,-1}}}; 
d.psiSettings{cD,1} = 'max_then_avg'; 
d.psiSettings{cD,2} = 'no_norm';
d.psiSettings{cD,3} = '';
d.psiSettings{cD,4} = {'optimise_reward'};
d.psiSettings{cD,5} = {[2]}; % This is which of the rewards to allow fitting for

% add d-counter to d
for iD = 1:size(d,1)
    d.iD(iD) = iD;
end

% Undo mistake when defining bins:
% Move everything one bin away so that near pos is 5-10 cm instead of 0-5
for iD = 1:size(d,1) - 2 % But don't shift 'the behind the back' data cause it's already correct (so skip last 2 indices)
    if ~iscell(d.cmPos{iD})
        if iD <= 28 || iD >= 33 % don't shift monkey data because it's already correct
            
            d.binPos{iD}(1,:) = d.binPos{iD}(1,:) - 1;
            
            % don't shift if bin position < 9 , because things will go out
            % of bounds, and responses are basically the same that far away
            % from the body part anyway
            d.binPos{iD}(1,d.binPos{iD}(1,:) < 9) = 9;
        end
    end
end



%% ------------------------------------------------------------------------
%  OPTIONAL: normalise empirical data 

% % ==============================
% z score for each default defined experiment
allExps = unique(d.exp);
for iExp = 1:numel(allExps)
    inclD   = find(d.exp == iExp);
    tmpD    = [d.realDat{inclD}];
    tmpMean = mean(tmpD);
    tmpSD = std(tmpD);
    
    % loop through all d that belong to experiment and remove mean
    for cD = inclD'
        d.realDat{cD} = (d.realDat{cD} - tmpMean) ./ tmpSD;
        d.realSTE{cD} = (d.realSTE{cD} ) ./ tmpSD;
    end

    % But also store the means and SDs so it can be converted back for the
    % plot
    % loop through all d that belong to experiment and remove mean
    for cD = inclD'
        d.realDatOrigMean{cD} = tmpMean;
        d.realDatOrigSD{cD}   = tmpSD;
    end
end

% % ==============================
% % Shift to zero mean for each experiment
% allExps = unique(d.exp);
% for iExp = 1:numel(allExps)
%     inclD   = find(d.exp == iExp);
%     tmpD    = [d.realDat{inclD}];
%     tmpMean = mean(tmpD);
%     
%     % loop through all d that belong to experiment and remove mean
%     for cD = inclD'
%         d.realDat{cD} = d.realDat{cD} - tmpMean
%     end
% end


% % ==============================
% % scale between 0 and 1 for each experiment
% allExps = unique(d.exp);
% for iExp = 1:numel(allExps)
%     inclD   = find(d.exp == iExp);
%     tmpD    = [d.realDat{inclD}];
%     tmpMin  = min(tmpD);
%     tmpMax  = max(tmpD);
%     
%     % loop through all d that belong to experiment and remove mean
%     for cD = inclD'
%         d.realDat{cD} = (d.realDat{cD} - tmpMin) ./ (tmpMax - tmpMin);
%         d.realSTE{cD} = (d.realSTE{cD} ) ./ (tmpMax - tmpMin);
%     end
% end


% % ==============================
% % Shift to zero minimum for eachexperiment
% allExps = unique(d.exp);
% for iExp = 1:numel(allExps)
%     inclD   = find(d.exp == iExp);
%     tmpD    = [d.realDat{inclD}];
%     tmpMin  = min(tmpD);
%     
%     % loop through all d that belong to experiment and remove mean
%     for cD = inclD'
%         d.realDat{cD} = d.realDat{cD} - tmpMin;
%     end
% end




%% ------------------------------------------------------------------------
%  Ensure all egocentric maps are set correctly, and create fitting tables
%  for other types of models (e.g. hit probability models, mutlisensory
%  integration, exponential decay etc etc)
% Other options for psi: 

dQ = d;
% ======================= Default Q-values ================================
% (Using rule of Trunk-primacy as described in Paper and Serino et al 2015)
dQ.psiSplit{1}  = { {'bodyPart',{'Trunk','Head','HandByChest'}},            {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{2}  = { {'bodyPart',{'Trunk','Head','HandByChest'}},            {'dir',{-1}} , {'rew',{1,-1}}};
dQ.psiSplit{3}  = { {'bodyPart',{'HandBySide'}},                            {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{4}  = { {'bodyPart',{'HandBySide'}},                            {'dir',{-1}} , {'rew',{1,-1}}};
dQ.psiSplit{5}  = { {'bodyPart',{'Trunk','Head','HandByChest'}},            {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{6}  = { {'bodyPart',{'Trunk','HandByChest'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{7}  = { {'bodyPart',{'Trunk','Head','HandByChest'}},            {'dir',{-1}} , {'rew',{1,-1}}}; 
dQ.psiSplit{8}  = { {'bodyPart',{'Trunk','HandByChest'}},                   {'dir',{-1}} , {'rew',{1,-1}}};
dQ.psiSplit{9}  = { {'bodyPart',{'HandBySide'}},                            {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{10} = { {'bodyPart',{'Trunk','HandByChest'}},                   {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{11} = { {'bodyPart',{'Head'}},                                  {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{12} = { {'bodyPart',{'Head'}},                                  {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{13} = { {'bodyPart',{'Head','Trunk'}},                          {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{14} = { {'bodyPart',{'Head','Trunk'}},                          {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{15} = { {'bodyPart',{'Head','Trunk'}},                          {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{16} = { {'bodyPart',{'Head'}},                                  {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{17} = { {'bodyPart',{'Trunk','HandByChest'}} ,                  {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{18} = { {'bodyPart',{'Trunk','HandByChestPlusTool'}},           {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{19} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,             {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{20} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,             {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{21} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,             {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{22} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,             {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{23} = { {'bodyPart',{'HandBySide','Head','Trunk'}},             {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dQ.psiSplit{24} = { {'bodyPart',{'Trunk','Head','HandByChest'}},            {'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL
dQ.psiSplit{25} = { {'bodyPart',{'Trunk','Head','HandByChestPlusToolRake'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL 
dQ.psiSplit{26} = { {'bodyPart',{'Trunk','Head','HandBySide'}},             {'dir',{1,-1}}  , {'rew',{1,-1}}};  % fMRI moving vs stationary dPOS
dQ.psiSplit{27} = { {'bodyPart',{'Trunk','Head','HandBySide'}},             {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dQ.psiSplit{28} = { {'bodyPart',{'Trunk','Head','HandBySide'}},             {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dQ.psiSplit{29} = { {'bodyPart',{'ArmForward'}},                            {'dir',{ 1}}  ,  {'rew',{1,-1}}}; % Single neuron - arm forward
dQ.psiSplit{30} = { {'bodyPart',{'ArmLeft'}},                               {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dQ.psiSplit{31} = { {'bodyPart',{'HeadConstr'}},                            {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
dQ.psiSplit{32} = { {'bodyPart',{'RotatedHead'}},                           {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side
dQ.psiSplit{33} = { {'bodyPart',{'Trunk','HandByChest'}} ,                  {'dir',{1}}   , {'rew',{1,-1}} }; 
dQ.psiSplit{34} = { {'bodyPart',{'Trunk','HandByChest'}} ,                  {'dir',{1}}   , {'rew',{1,-1}} }; 
dQ.psiSplit{35} = { {'bodyPart',{'Trunk','HandByChest'}} ,                  {'dir',{1}}   , {'rew',{1,-1}} }; 
dQ.psiSplit{36} = { {'bodyPart',{'Back'}},                                  {'dir',{1}}   , {'rew',{1,-1}}};  
dQ.psiSplit{37} = { {'bodyPart',{'Back'}},                                  {'dir',{1}}   , {'rew',{1,-1}}}; 

dQ.psiSettings{19,4} = {'optimise_reward'};
dQ.psiSettings{19,5} = {[2]}; % This is which of the rewards to allow fitting for
dQ.psiSettings{20,4} = {'optimise_reward'};
dQ.psiSettings{20,5} = {[]}; % This is which of the rewards to allow fitting for. Because it doesn't have an index, it should use the same p and the same index as the previous one.
dQ.psiSettings{21,4} = {'optimise_reward'};
dQ.psiSettings{21,5} = {[2]}; % This is which of the rewards to allow fitting for
dQ.psiSettings{22,4} = {'optimise_reward'};
dQ.psiSettings{22,5} = {[]}; % This is which of the rewards to allow fitting for. Because it doesn't have an index, it should use the same p and the same index as the previous one.



dUn = d;
% ======================= Uncertain Q-values ================================
% % Using rule of Trunk-primacy 
dUn.psiSplit{1}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{2}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{-1}} , {'rew',{1,-1}}};
dUn.psiSplit{3}  = { {'bodyPart',{'HandBySideUnc'}},                        {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{4}  = { {'bodyPart',{'HandBySideUnc'}},                        {'dir',{-1}} , {'rew',{1,-1}}};
dUn.psiSplit{5}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{6}  = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}},            {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{7}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{-1}} , {'rew',{1,-1}}}; 
dUn.psiSplit{8}  = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}},            {'dir',{-1}} , {'rew',{1,-1}}};
dUn.psiSplit{9}  = { {'bodyPart',{'HandBySideUnc'}},                        {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{10} = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}},            {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{11} = { {'bodyPart',{'HeadUnc'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{12} = { {'bodyPart',{'HeadUnc'}},                              {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{13} = { {'bodyPart',{'HeadUnc','TrunkUnc'}},                   {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{14} = { {'bodyPart',{'HeadUnc','TrunkUnc'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{15} = { {'bodyPart',{'HeadUnc','TrunkUnc'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{16} = { {'bodyPart',{'HeadUnc'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{17} = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}} ,           {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{18} = { {'bodyPart',{'TrunkUnc','HandByChestPlusToolUnc'}},    {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{19} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,            {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{20} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,            {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{21} = { {'bodyPart',{'TrunkUnc','HandByChestTrack'}} ,         {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{22} = { {'bodyPart',{'TrunkUnc','HandByChestTrack'}} ,         {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{23} = { {'bodyPart',{'HandBySideUnc','HeadUnc','TrunkUnc'}},   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dUn.psiSplit{24} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL
dUn.psiSplit{25} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestPlusToolRakeUnc'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL 
dUn.psiSplit{26} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandBySideUnc'}},   {'dir',{1,-1}}  , {'rew',{1,-1}}};  % fMRI moving vs stationary dPOS
dUn.psiSplit{27} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandBySideUnc'}},   {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dUn.psiSplit{28} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandBySideUnc'}},   {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dUn.psiSplit{29} = { {'bodyPart',{'ArmForwardUnc'}},                        {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dUn.psiSplit{30} = { {'bodyPart',{'ArmLeftUnc'}},                           {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dUn.psiSplit{31} = { {'bodyPart',{'HeadConstrUnc'}},                        {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
dUn.psiSplit{32} = { {'bodyPart',{'RotatedHeadUnc'}},                       {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side
dUn.psiSplit{33} = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}} ,           {'dir',{1}}   , {'rew',{1,-1}} }; 
dUn.psiSplit{34} = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}} ,           {'dir',{1}}   , {'rew',{1,-1}} }; 
dUn.psiSplit{35} = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}} ,           {'dir',{1}}   , {'rew',{1,-1}} }; 
dUn.psiSplit{36} = { {'bodyPart',{'BackUnc'}},                              {'dir',{1}}   , {'rew',{1,-1}}};  
dUn.psiSplit{37} = { {'bodyPart',{'BackUnc'}},                              {'dir',{1}}   , {'rew',{1,-1}}}; 

dUn.psiSettings{19,4} = {'optimise_reward'};
dUn.psiSettings{19,5} = {[2]}; % This is which of the rewards to allow fitting for
dUn.psiSettings{20,4} = {'optimise_reward'};
dUn.psiSettings{20,5} = {[]}; % This is which of the rewards to allow fitting for. Because it doesn't have an index, it should use the same p and the same index as the previous one.
dUn.psiSettings{21,4} = {'optimise_reward'};
dUn.psiSettings{21,5} = {[1 2]}; % This is which of the rewards to allow fitting for
dUn.psiSettings{22,4} = {'optimise_reward'};
dUn.psiSettings{22,5} = {[]}; % This is which of the rewards to allow fitting for. Because it doesn't have an index, it should use the same p and the same index as the previous one.



dSr = d;
% ======================= Uncertain Q-values ================================
% % Using rule of Trunk-primacy 
dSr.psiSplit{1}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{2}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},    {'dir',{-1}} , {'rew',{1,-1}}};
dSr.psiSplit{3}  = { {'bodyPart',{'HandBySideSARSA'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{4}  = { {'bodyPart',{'HandBySideSARSA'}},                              {'dir',{-1}} , {'rew',{1,-1}}};
dSr.psiSplit{5}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{6}  = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}},                {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{7}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},    {'dir',{-1}} , {'rew',{1,-1}}}; 
dSr.psiSplit{8}  = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}},                {'dir',{-1}} , {'rew',{1,-1}}};
dSr.psiSplit{9}  = { {'bodyPart',{'HandBySideSARSA'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{10} = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}},                {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{11} = { {'bodyPart',{'HeadSARSA'}},                                    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{12} = { {'bodyPart',{'HeadSARSA'}},                                    {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{13} = { {'bodyPart',{'HeadSARSA','TrunkSARSA'}},                       {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{14} = { {'bodyPart',{'HeadSARSA','TrunkSARSA'}},                       {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{15} = { {'bodyPart',{'HeadSARSA','TrunkSARSA'}},                       {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{16} = { {'bodyPart',{'HeadSARSA'}},                                    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{17} = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}} ,               {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{18} = { {'bodyPart',{'TrunkSARSA','HandByChestPlusToolSARSA'}},        {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{19} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,                    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{20} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,                    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{21} = { {'bodyPart',{'TrunkSARSA','HandByChestTrack'}} ,               {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{22} = { {'bodyPart',{'TrunkSARSA','HandByChestTrack'}} ,               {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{23} = { {'bodyPart',{'HandBySideSARSA','HeadSARSA','TrunkSARSA'}},     {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dSr.psiSplit{24} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},    {'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL
dSr.psiSplit{25} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestPlusToolRakeSARSA'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL 
dSr.psiSplit{26} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandBySideSARSA'}},     {'dir',{1,-1}}  , {'rew',{1,-1}}};  % fMRI moving vs stationary dPOS
dSr.psiSplit{27} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandBySideSARSA'}},     {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dSr.psiSplit{28} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandBySideSARSA'}},     {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dSr.psiSplit{29} = { {'bodyPart',{'ArmForwardSARSA'}},                              {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dSr.psiSplit{30} = { {'bodyPart',{'ArmLeftSARSA'}},                                 {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dSr.psiSplit{31} = { {'bodyPart',{'HeadConstrSARSA'}},                              {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
dSr.psiSplit{32} = { {'bodyPart',{'RotatedHeadSARSA'}},                             {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side
dSr.psiSplit{33} = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}} ,               {'dir',{1}}   , {'rew',{1,-1}} }; 
dSr.psiSplit{34} = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}} ,               {'dir',{1}}   , {'rew',{1,-1}} }; 
dSr.psiSplit{35} = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}} ,               {'dir',{1}}   , {'rew',{1,-1}} }; 
dSr.psiSplit{36} = { {'bodyPart',{'BackSARSA'}},                                    {'dir',{1}}   , {'rew',{1,-1}}};  
dSr.psiSplit{37} = { {'bodyPart',{'BackSARSA'}},                                    {'dir',{1}}   , {'rew',{1,-1}}}; 

dSr.psiSettings{19,4} = {'optimise_reward'};
dSr.psiSettings{19,5} = {[2]}; % This is which of the rewards to allow fitting for
dSr.psiSettings{20,4} = {'optimise_reward'};
dSr.psiSettings{20,5} = {[]}; % This is which of the rewards to allow fitting for. Because it doesn't have an index, it should use the same p and the same index as the previous one. OK
dSr.psiSettings{21,4} = {'optimise_reward'};
dSr.psiSettings{21,5} = {[1 2]}; % This is which of the rewards to allow fitting for
dSr.psiSettings{22,4} = {'optimise_reward'};
dSr.psiSettings{22,5} = {[]}; % This is which of the rewards to allow fitting for. Because it doesn't have an index, it should use the same p and the same index as the previous one. OK



dHP = d;
% ======================= HITPROB VERSION: ================================
% % Using rule of Trunk-primacy 
dHP.psiSplit{1}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{2}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{-1}} , {'rew',{1}}};
dHP.psiSplit{3}  = { {'bodyPart',{'HandBySideHP'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{4}  = { {'bodyPart',{'HandBySideHP'}},                       {'dir',{-1}} , {'rew',{1}}};
dHP.psiSplit{5}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{6}  = { {'bodyPart',{'TrunkHP','HandByChestHP'}},            {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{7}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{-1}} , {'rew',{1}}}; 
dHP.psiSplit{8}  = { {'bodyPart',{'TrunkHP','HandByChestHP'}},            {'dir',{-1}} , {'rew',{1}}};
dHP.psiSplit{9}  = { {'bodyPart',{'HandBySideHP'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{10} = { {'bodyPart',{'TrunkHP','HandByChestHP'}},            {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{11} = { {'bodyPart',{'HeadHP'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{12} = { {'bodyPart',{'HeadHP'}},                             {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{13} = { {'bodyPart',{'HeadHP','TrunkHP'}},                   {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{14} = { {'bodyPart',{'HeadHP','TrunkHP'}},                   {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{15} = { {'bodyPart',{'HeadHP','TrunkHP'}},                   {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{16} = { {'bodyPart',{'HeadHP'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{17} = { {'bodyPart',{'TrunkHP','HandByChestHP'}} ,           {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{18} = { {'bodyPart',{'TrunkHP','HandByChestPlusToolHP'}},    {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{19} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{20} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{21} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{22} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{23} = { {'bodyPart',{'HandBySideHP','HeadHP','TrunkHP'}} ,   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dHP.psiSplit{24} = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}} ,  {'dir',{1,-1}}  , {'rew',{1}}}; % EEG RONGA TOOL
dHP.psiSplit{25} = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestPlusToolRakeHP'}},{'dir',{1,-1}}  , {'rew',{1}}}; % EEG RONGA TOOL 
dHP.psiSplit{26} = { {'bodyPart',{'TrunkHP','HeadHP','HandBySideHP'}},    {'dir',{1,-1}}  , {'rew',{1}}};  % fMRI moving vs stationary dPOS
dHP.psiSplit{27} = { {'bodyPart',{'TrunkHP','HeadHP','HandBySideHP'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dHP.psiSplit{28} = { {'bodyPart',{'TrunkHP','HeadHP','HandBySideHP'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dHP.psiSplit{29} = { {'bodyPart',{'ArmForwardHP'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dHP.psiSplit{30} = { {'bodyPart',{'ArmLeftHP'}},                          {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dHP.psiSplit{31} = { {'bodyPart',{'HeadConstrHP'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head forward
dHP.psiSplit{32} = { {'bodyPart',{'RotatedHeadHP'}},                      {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head to side
dHP.psiSplit{33} = { {'bodyPart',{'TrunkHP','HandByChestHP'}} ,           {'dir',{1}}   , {'rew',{1}}}; 
dHP.psiSplit{34} = { {'bodyPart',{'TrunkHP','HandByChestHP'}} ,           {'dir',{1}}   , {'rew',{1}}}; 
dHP.psiSplit{35} = { {'bodyPart',{'TrunkHP','HandByChestHP'}} ,           {'dir',{1}}   , {'rew',{1}}}; 
dHP.psiSplit{36} = { {'bodyPart',{'BackHP'}},                             {'dir',{1}}   , {'rew',{1}}};  
dHP.psiSplit{37} = { {'bodyPart',{'BackHP'}},                             {'dir',{1}}   , {'rew',{1}}}; 
% Don't fit reward (because HP can't fit reward)
for iD = 1:size(dHP,1)
    dHP(iD,:).psiSettings{4} = [];
    dHP(iD,:).psiSettings{5} = [];
end

dMI = d;
% ======================= MultiSENSORY INTEGRATION VERSION: ================================
% % Using rule of Trunk-primacy 
dMI.psiSplit{1}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{2}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}} , {'rew',{1}}};
dMI.psiSplit{3}  = { {'bodyPart',{'HandBySideMI'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{4}  = { {'bodyPart',{'HandBySideMI'}},                       {'dir',{1}} , {'rew',{1}}};
dMI.psiSplit{5}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{6}  = { {'bodyPart',{'TrunkMI','HandByChestMI'}},            {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{7}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}} , {'rew',{1}}}; 
dMI.psiSplit{8}  = { {'bodyPart',{'TrunkMI','HandByChestMI'}},            {'dir',{1}} , {'rew',{1}}};
dMI.psiSplit{9}  = { {'bodyPart',{'HandBySideMI'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{10} = { {'bodyPart',{'TrunkMI','HandByChestMI'}},            {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{11} = { {'bodyPart',{'HeadMI'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{12} = { {'bodyPart',{'HeadMI'}},                             {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{13} = { {'bodyPart',{'HeadMI','TrunkMI'}},                   {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{14} = { {'bodyPart',{'HeadMI','TrunkMI'}},                   {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{15} = { {'bodyPart',{'HeadMI','TrunkMI'}},                   {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{16} = { {'bodyPart',{'HeadMI'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{17} = { {'bodyPart',{'HandByChestMI'}} ,                     {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{18} = { {'bodyPart',{'TrunkMI','HandByChestPlusToolMI'}},    {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{19} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{20} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{21} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{22} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{23} = { {'bodyPart',{'HandBySideMI','HeadMI','TrunkMI'}} ,   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dMI.psiSplit{24} = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}} ,  {'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL
dMI.psiSplit{25} = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestPlusToolRakeMI'}},{'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL 
dMI.psiSplit{26} = { {'bodyPart',{'TrunkMI','HeadMI','HandBySideMI'}},    {'dir',{1}}  , {'rew',{1}}};  % fMRI moving vs stationary dPOS
dMI.psiSplit{27} = { {'bodyPart',{'TrunkMI','HeadMI','HandBySideMI'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dMI.psiSplit{28} = { {'bodyPart',{'TrunkMI','HeadMI','HandBySideMI'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dMI.psiSplit{29} = { {'bodyPart',{'ArmForwardMI'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dMI.psiSplit{30} = { {'bodyPart',{'ArmLeftMI'}},                          {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dMI.psiSplit{31} = { {'bodyPart',{'HeadConstrMI'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head forward
dMI.psiSplit{32} = { {'bodyPart',{'RotatedHeadMI'}},                      {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head to side
dMI.psiSplit{33} = { {'bodyPart',{'TrunkMI','HandByChestMI'}} ,           {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{34} = { {'bodyPart',{'TrunkMI','HandByChestMI'}} ,           {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{35} = { {'bodyPart',{'TrunkMI','HandByChestMI'}} ,           {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{36} = { {'bodyPart',{'BackMI'}},                             {'dir',{1}}  , {'rew',{1}}};  
dMI.psiSplit{37} = { {'bodyPart',{'BackMI'}},                             {'dir',{1}}  , {'rew',{1}}}; 
% Don't fit reward (because MI can't fit reward)
for iD = 1:size(dMI,1)
    dMI(iD,:).psiSettings{4} = [];
    dMI(iD,:).psiSettings{5} = [];
end

dDi = d;
% ======================= DISTANCE VERSION: ================================
% % Using rule of Trunk-primacy 
dDi.psiSplit{1}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{2}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}} , {'rew',{1}}};
dDi.psiSplit{3}  = { {'bodyPart',{'HandBySideDist'}},                           {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{4}  = { {'bodyPart',{'HandBySideDist'}},                           {'dir',{1}} , {'rew',{1}}};
dDi.psiSplit{5}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{6}  = { {'bodyPart',{'TrunkDist','HandByChestDist'}},              {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{7}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}} , {'rew',{1}}}; 
dDi.psiSplit{8}  = { {'bodyPart',{'TrunkDist','HandByChestDist'}},              {'dir',{1}} , {'rew',{1}}};
dDi.psiSplit{9}  = { {'bodyPart',{'HandBySideDist'}},                           {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{10} = { {'bodyPart',{'TrunkDist','HandByChestDist'}},              {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{11} = { {'bodyPart',{'HeadDist'}},                                 {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{12} = { {'bodyPart',{'HeadDist'}},                                 {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{13} = { {'bodyPart',{'HeadDist','TrunkDist'}},                     {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{14} = { {'bodyPart',{'HeadDist','TrunkDist'}},                     {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{15} = { {'bodyPart',{'HeadDist','TrunkDist'}},                     {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{16} = { {'bodyPart',{'HeadDist'}},                                 {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{17} = { {'bodyPart',{'HandByChestDist'}} ,                         {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{18} = { {'bodyPart',{'TrunkDist','HandByChestPlusToolDist'}},      {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{19} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{20} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{21} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{22} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{23} = { {'bodyPart',{'HandBySideDist','HeadDist','TrunkDist'}} ,   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dDi.psiSplit{24} = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}} ,  {'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL
dDi.psiSplit{25} = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestPlusToolRakeDist'}},{'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL 
dDi.psiSplit{26} = { {'bodyPart',{'TrunkDist','HeadDist','HandBySideDist'}},    {'dir',{1}}  , {'rew',{1}}};  % fMRI moving vs stationary dPOS
dDi.psiSplit{27} = { {'bodyPart',{'TrunkDist','HeadDist','HandBySideDist'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dDi.psiSplit{28} = { {'bodyPart',{'TrunkDist','HeadDist','HandBySideDist'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dDi.psiSplit{29} = { {'bodyPart',{'ArmForwardDist'}},                           {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dDi.psiSplit{30} = { {'bodyPart',{'ArmLeftDist'}},                              {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dDi.psiSplit{31} = { {'bodyPart',{'HeadConstrDist'}},                           {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head forward
dDi.psiSplit{32} = { {'bodyPart',{'RotatedHeadDist'}},                          {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head to side
dDi.psiSplit{33} = { {'bodyPart',{'TrunkDist','HandByChestDist'}} ,             {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{34} = { {'bodyPart',{'TrunkDist','HandByChestDist'}} ,             {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{35} = { {'bodyPart',{'TrunkDist','HandByChestDist'}} ,             {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{36} = { {'bodyPart',{'BackDist'}},                                 {'dir',{1}}  , {'rew',{1}}};  
dDi.psiSplit{37} = { {'bodyPart',{'BackDist'}},                                 {'dir',{1}}  , {'rew',{1}}}; 
% Don't fit reward (because distance can't fit reward)
for iD = 1:size(dDi,1)
    dDi(iD,:).psiSettings{4} = [];
    dDi(iD,:).psiSettings{5} = [];
end

%% ------------------------------------------------------------------------
% Create 'egocentric maps' and other theoretical constructs
for iD = 1:size(d,1) 
    [dQ.psi{iD},   dQ.posRewPsi{iD},  dQ.negRewPsi{iD},  dQ.psiLines{iD},  dQ.psiLineDescr{iD}]  = MakePsi(allQ, dQ(iD,:));
    [dHP.psi{iD},  dHP.posRewPsi{iD}, dHP.negRewPsi{iD}, dHP.psiLines{iD}, dHP.psiLineDescr{iD}] = MakePsi(allQ, dHP(iD,:));
    [dMI.psi{iD},  dMI.posRewPsi{iD}, dMI.negRewPsi{iD}, dMI.psiLines{iD}, dMI.psiLineDescr{iD}] = MakePsi(allQ, dMI(iD,:));
    [dDi.psi{iD},  dDi.posRewPsi{iD}, dDi.negRewPsi{iD}, dDi.psiLines{iD}, dDi.psiLineDescr{iD}] = MakePsi(allQ, dDi(iD,:));
    [dUn.psi{iD},  dUn.posRewPsi{iD}, dUn.negRewPsi{iD}, dUn.psiLines{iD}, dUn.psiLineDescr{iD}] = MakePsi(allQ, dUn(iD,:));
    [dSr.psi{iD},  dSr.posRewPsi{iD}, dSr.negRewPsi{iD}, dSr.psiLines{iD}, dSr.psiLineDescr{iD}] = MakePsi(allQ, dSr(iD,:));

    % Add Experiment information to psilinedescr
    for iDescr = 1:numel(dQ.psiLineDescr{iD})
        dQ.psiLineDescr{iD}{iDescr}   = [dQ.psiLineDescr{iD}{iDescr}   'Exp == ' num2str(dQ.exp(iD))];
    end
    for iDescr = 1:numel(dHP.psiLineDescr{iD})
        dHP.psiLineDescr{iD}{iDescr} = [dHP.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dHP.exp(iD))];
    end
    for iDescr = 1:numel(dMI.psiLineDescr{iD})
        dMI.psiLineDescr{iD}{iDescr} = [dMI.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dMI.exp(iD))];
    end
    for iDescr = 1:numel(dDi.psiLineDescr{iD})
        dDi.psiLineDescr{iD}{iDescr} = [dDi.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dDi.exp(iD))];
    end
    for iDescr = 1:numel(dUn.psiLineDescr{iD})
        dUn.psiLineDescr{iD}{iDescr} = [dUn.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dUn.exp(iD))];
    end
    for iDescr = 1:numel(dSr.psiLineDescr{iD})
        dSr.psiLineDescr{iD}{iDescr} = [dSr.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dSr.exp(iD))];
    end
end

% Make the distances negative, so that values are maximal near the limbs,
% and positive slopes can fit the data properly
dNegDi = dDi;
dNegDi.psi = cellfun(@(psi) -psi, dDi.psi, 'UniformOutput', false);

%% ========================================================================
% Fit all data: 

% Fit the main model:
[fitRes.Q]      = Fit3Ddat(dQ,'Linear');
% Fit the uncertain model:
[fitRes.QUN]    = Fit3Ddat(dUn,'Linear');
% Fit the SARSA model:
[fitRes.SRS]    = Fit3Ddat(dSr,'Linear');
% Fit the hit-probability model:
[fitRes.HP]     = Fit3Ddat(dHP,'Linear');
% Fit the mutlisensory integration model:
[fitRes.MI]     = Fit3Ddat(dMI,'Linear');
% Fit the exponential model:
[fitRes.EXP]    = Fit3Ddat(dDi,'Exponential');
% Fit the sigmoid model:
[fitRes.SIG]    = Fit3Ddat(dDi,'Sigmoid');
% Fit the linear model:
[fitRes.LIN]    = Fit3Ddat(dNegDi,'Linear');


%% ========================================================================
%  FUNCTIONS

% -------------------------------------------------------------------------
function [psi posRewPsi negRewPsi psiLines lineDescr] = MakePsi(allQ, d)
% Return successor features specific to the dataset

[psiLines lineDescr]   = FindPsiLines(allQ, d.psiSplit{1});

% Add to psi in manner dependent on the settings
psi =[];
for iPsi = 1:numel(psiLines)

    % q values for all actions of a particular [iPsi] bodypart
    if ~iscell(d.cmPos{1})
        % If comparing to some baseline Q-value, don't take the absolute
        if strcmp(d.psiSettings{1,1},'max_above_stationary')
            tmpPsi         = ExtractQ(allQ.qVals(psiLines{iPsi} ),d.binPos{1});
        else
            tmpPsi         = abs(ExtractQ(allQ.qVals(psiLines{iPsi} ),d.binPos{1}));
        end
    end

    switch d.psiSettings{1,3}
        case 'avOverRows'
            cBinPos = d.binPos{1};         

            % If stimulus poisiont is defined as a cell, average over each
            % entry within a give cell
            if iscell(d.cmPos{1})

                % Initialise to NaN so that different lengths can be used
                tmpTmpPsi = nan([max(cellfun(@(x) numel(x), d.binPos{1}),[],'all')  ...
                    size(allQ.qVals{psiLines{iPsi} },1) size(d.cmPos{1},2)]);

                for iTraj = 1:size(d.cmPos{1},2) % loop through trajectories

                    avRows = d.binPos{1}{1,iTraj};

                    for iRow = 1:numel(avRows)
                        cRow = avRows(iRow);

                        cBinPos = [d.binPos{1}{1,iTraj}(iRow) d.binPos{1}{2,iTraj}(iRow) d.binPos{1}{3,iTraj}(iRow)]';

                        % rows within trajectory, then actions, then positions,
                        % i.e. trajectories
                        tmpTmpPsi(iRow,:,iTraj)         = abs(ExtractQ(allQ.qVals(psiLines{iPsi} ),cBinPos));
                    end
                end

            % If no further indication is give, loop over 1m [100cm], stricly along
            % the row dimension
            else
                avRows = min(d.binPos{1}(1,:)) : -1 : (min(d.binPos{1}(1,:)) - 20);
                for iRow = 1:numel(avRows)
                    cRow = avRows(iRow);
                    cBinPos(1,:) = cRow;
                    tmpTmpPsi(iRow,:,:,:,:,:,:)         = abs(ExtractQ(allQ.qVals(psiLines{iPsi} ),cBinPos));
                end
            end

            % Take into account possibility of only 1 possible action
            if size(tmpTmpPsi,2) == 1
                tmpPsi(1,:) = squeeze(nanmean(tmpTmpPsi,1));
            else
                tmpPsi = squeeze(nanmean(tmpTmpPsi,1));
            end
    end
  
    % Placeholder for relative action probability to 'stay still' action
    switch d.psiSettings{1,2}
        case 'rel_still_norm'
            stayPsi = tmpPsi(1,:);
            tmpPsi = tmpPsi - stayPsi;

    end
    switch d.psiSettings{1,1}
        case {'max_above_stationary'}
            tmpPsi = max(tmpPsi,[],1);
        case {'average', 'average_all'}
            tmpPsi = mean(tmpPsi,1);
        case {'max','max_all','max_then_avg'}
            tmpPsi = max(tmpPsi,[],1);
        case 'raw'
    end
    psi = [psi ; tmpPsi];
end

% Make a version that is split by reward sign
rewSign   = allQ(cell2mat(psiLines),:).rew;
posRewPsi = psi(rewSign > 0,:);
negRewPsi = psi(rewSign < 0,:);

switch d.psiSettings{1,1}
    case {'max_above_stationary'}
        posRewPsi = mean(posRewPsi,1);
        negRewPsi = mean(negRewPsi,1);
        psi       = mean( [posRewPsi ; negRewPsi] ,1);
    case {'average_all','max_then_avg'}
        psi       = mean(psi,1);
        posRewPsi = mean(posRewPsi,1);
        negRewPsi = mean(negRewPsi,1);
    case 'max_all'
        psi       = max(psi,[],1);
        posRewPsi = max(posRewPsi,[],1);
        negRewPsi = max(negRewPsi,[],1);
end

% Normalise across successor features --> simulates mutual inhibition, for
% example
switch d.psiSettings{1,2}
    case 'sum_norm'
        psi = psi ./ (sum(psi,1));
    case 'diff_from_mean_norm'
        psi = psi - mean(psi,1);
    case 'diff_from_mean_divmean_norm'
        psi = (psi - mean(psi,1))./ mean(psi,1);
    case 'no_norm'
end
% Adjust for case where everything dividing by zero
psi(isnan(psi)) = 0;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [allLines lineDescr] = FindPsiLines(allQ, psiSplit)
% Finds the lines in allQ of the Psi features that will be included
% Also returns a description of what those lines entail
elements = {};
% Create list of elements: cell array with N vectors to combine
for iSC = 1:numel(psiSplit) % loop through split conditions
    elements = [elements, {psiSplit{iSC}{2}}];
end

combinations        = cell(1, numel(elements)); %set up the varargout result
[combinations{:}]   = ndgrid(elements{:}); % feed each element as separate input
combinations        = cellfun(@(x) x(:), combinations,'uniformoutput',false); %Transform into vectors

allSplits = arrayfun(@(ind) psiSplit{ind}{1},1:numel(psiSplit),'UniformOutput',false);

% Loop through combinations to create logical indexes
allLines = [];
for iSS = 1:numel(combinations{1}) % sub split index

    % Make a Boolean specifying which lines of allQ correspond to this
    % combination of conditions
    inclLines = ones([size(allQ,1) 1]);
    lineDescr{iSS} = '';
    for iSC = 1:numel(allSplits)
        currSpl = allSplits{iSC}; % current split

        % Find the lines in allQ that should be included
        inclLines = inclLines & IsEqual(allQ.(currSpl),{combinations{iSC}{iSS}});
        % Only add to the description if the current split also has
        % alternative options. Otherwise just consider it an effect that
        % Psi naturally takes into account - i.e.e there ren't multiple
        % psis for that effect
        if numel(psiSplit{iSC}{2}) == 1
        else
            lineDescr{iSS} = [lineDescr{iSS} currSpl ' == ' num2str(combinations{iSC}{iSS}) '. '];
        end
    end
    allLines = [allLines; {find( inclLines )}];
end
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [eqBool] = IsEqual(d1,d2)
% Compares any input, whether a string or a number. Takes CELL input. Note
% at least one of the entries must only contain one value

% Convert EVERYTHING to strings inside cells
if ~ iscell(d1)
    d1 = num2cell(d1);
end
if ~ iscell(d2)
    d1 = num2cell(d2);
end
d1 = cellfun(@(subD) num2str(subD), d1, 'UniformOutput', false);
d2 = cellfun(@(subD) num2str(subD), d2, 'UniformOutput', false);

% Make sure they are the same shape
d1 = d1(:);
d2 = d2(:);

% Compare the strings
eqBool = strcmp(d1,d2);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [d] = InsertModelledDat(d,dFitFinal,allDat)
% Takes modelled data and puts it back into the data structure for easier
% plotting
    cP = 0; % current point
    for iD = 1:size(d,1)
        nP = numel(d.realDat{iD}); % number of points
        d.fittedDat{iD} = dFitFinal(cP + [1:nP]);
        cP = cP + nP;
    end
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [psiMat allDat psiLabels] = MakePsiMat(d)
% Makes a successor state matrix that allows faster fitting of the data,
% assuming different multipliers for each tactile location

allTacts = unique(d.tact);
numTacts = numel(allTacts);
numFeats = size(d.psi{1},1);

allDat   = [d.realDat{:}];
numDats  = numel(allDat);

psiMat = nan([numFeats .* numTacts, numDats]);

% Current Row
cRow  = 0;
cRow2 = 0;
doneTacts = [];
for iD = 1:size(d,1)

    cTact = d.tact{iD};

    % Stick the psi in the appropriate tactile zone [the rest is 0]
    tactRow = find(strcmp(allTacts,cTact));
    tactRow = ((tactRow - 1) .* numFeats) + 1;

    % Current data
    cDat    = d.realDat{iD};
    cDLen   = numel(cDat);

    % Re-average or re-max the psi across reward values, in case it's
    % necessary
    rewInd    = find(cellfun(@(splitType) strcmp(splitType{1},'rew') , d.psiSplit{iD}));
    rewMagns  = cell2mat(d.psiSplit{iD}{rewInd}{2});
    nonOneRew = abs(rewMagns) ~= 1;
    if any(nonOneRew) % make new psi
        switch d(iD,:).psiSettings{1}
            case {'average_all','max_then_avg','max_above_stationary'}
                d.psi{iD}  = nanmean(abs(rewMagns(:)) .* [d.posRewPsi{iD}; d.negRewPsi{iD}]);
            case 'max_all'
                d.psi{iD}  = nanmax( abs(rewMagns(:)) .* [d.posRewPsi{iD}; d.negRewPsi{iD}],[],1);
        end
    end
    psiMat(tactRow:tactRow + numFeats - 1, [1:cDLen] + cRow) = d.psi{iD};
    
    % Store the description of the features that the weight refers to
    if ~ismember(cTact, doneTacts)
        for iPsi = 1:numFeats
            doneTacts = [doneTacts {cTact}];
            psiLabels{cRow2 + iPsi} = ['Tact Loc: ' cTact '. ' d.psiLineDescr{iD}{iPsi}];
        end
        cRow2 = cRow2 + iPsi;
    end

    cRow = cRow + cDLen;
end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Calculate fit quality and return best modelled data
function [chiSq, sSqErr, dFit, residuals, d, psiMat] = ErrorFun(psiMat,allDat,d,p,fitType)  % $$$ [sSqErr, dFit, residuals, d, psiMat] = ErrorFun(psiMat,allDat,d,p,fitType)
 
    allExps = unique(d.exp);
    % To use an offset for each experiment
    % Make a list of experiment numbers for each entry in the d table
    expCol  = [];
    iDCol   = [];
    for iD = 1:size(d,1)
        expCol  = [expCol,  d.exp(iD) .* ones(size(d.realDat{iD})) ];
        % Also make a list of ID numbers
        iDCol =   [iDCol,          iD .* ones(size(d.realDat{iD})) ];
    end
    for iExp = 1:numel(allExps)
        expPos = expCol == iExp;
        expOffset(expPos) = p(size(psiMat,1) + iExp);
    end

    switch fitType
        case 'Exponential'
            % Parameter controlling decay rate of the exponential
            psiMat  = exp( (psiMat - p(size(psiMat,1) + iExp + 2) ) .* (p(size(psiMat,1) + iExp + 1))  ); 
            % Also add to the counter to ensure that the retreating stimuli 
            % offset parameter (added below) is in the right place
            iExp    = iExp + 2;
        case 'Sigmoid'
            % Sigmoid fit
            psiMat  = 1 ./ (  1 + exp(-(psiMat -  p(size(psiMat,1) + iExp + 2)     )   ./   (p(size(psiMat,1) + iExp + 1))   )  ) ;
            % Also add to the counter to ensure that the retreating stimuli 
            % offset parameter (added below) is in the right place
            iExp    = iExp + 2;
        case 'Linear'
            % Don't need to add anything for linear cause it's all fit
            % linearly by distance already
    end

    % Add a separate offset for stimuli moving away direction
    if numel(p) > size(psiMat,1) + iExp
        currInd = 0;
        for iD = 1:size(d,1)
            dDat = [d.realDat{iD}];
            if d.dir(iD) == -1
                dirOffset(currInd + [1:numel(dDat)]) = p(size(psiMat,1) + iExp + 1);
            else
                dirOffset(currInd + [1:numel(dDat)]) = 0;
            end
            currInd = numel(dirOffset);
        end
    else
        dirOffset = 0;
    end

    % Fit the reward values, ONLY for those experiments which have the rewardfitting setting turned on
    fitRewards = find(cellfun(@(optSet) strcmp(optSet,'optimise_reward'), d.psiSettings(:,4)));
    cRewNum    = 1;
    if numel(fitRewards) > 0 && strcmp(fitType,'Linear')
        dTmp = d;
        % Update the reward values in dTmp, and calculate new psiMats
        for iFitRewExp = 1:numel(fitRewards)
            cF         = fitRewards(iFitRewExp);

            % Update each of the rewards independently
            % But first find which index refers to rewards (in case I
            % misordered them somewhere)
            rewInd = find(cellfun(@(splitType) strcmp(splitType{1},'rew') , dTmp.psiSplit{cF}));
            % If the to-be-changed-reward index is empty, use the previous
            % intended reward: it means that they share a fitting parameter
            if isempty((dTmp(cF,:).psiSettings{5}{1}))
                tmpI = 1;
                for cInExpRew = (dTmp(lastCF,:).psiSettings{5}{1}) % Loop through the parameters used on the previous condition
                    dTmp(cF,:).psiSplit{1}{rewInd}{2}{cInExpRew} = p(size(psiMat,1) + iExp + cRewNum - numel((dTmp(lastCF,:).psiSettings{5}{1})) + tmpI );
                    tmpI = tmpI + 1;
                end

            % If it's not empty, assing a new parameter
            else
                for cInExpRew = (dTmp(cF,:).psiSettings{5}{1})
                    dTmp(cF,:).psiSplit{1}{rewInd}{2}{cInExpRew} = p(size(psiMat,1) + iExp + cRewNum + 1);
                    cRewNum = cRewNum + 1;
                    lastCF  = cF; % store current CF in case this parameter needs to be re-used
                end
            end
        end

        % Recaltulate the psiMat
        [psiMat, allDat, ~] = MakePsiMat(dTmp);
        d = dTmp;
    end

    % Calculate fit metrics
    dFit = nansum( psiMat .* p(1:size(psiMat,1)) ) + expOffset + dirOffset;
    residuals = (dFit - allDat).^2;
    sSqErr = nansum(residuals);
    dataVar = [d.realSTE{:}].^2;
    chiSq   = nansum(residuals(:)./dataVar(:));

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function d = AddbinPos(d,s,dNums)
for iFitD = dNums

    % If there are multiple entries per 'location' adjust accordingly
    if iscell(d.cmPos{iFitD})

        d.binPos{iFitD} = cell(size(d.cmPos{iFitD}));
        tmpCmPos = d.cmPos{iFitD};

        for iTraj = 1:size(d.cmPos{iFitD},2) % loop through trajectories
            
            % The rows have to be flipped, because small == far
            tmpCmPos{1,iTraj} = -d.cmPos{iFitD}{1,iTraj};
            for iDim = 1:size(d.cmPos{iFitD},1)
                % Convert real positions to voxels
                d.binPos{iFitD}{iDim,iTraj} = round(s.clc.nearPos(iDim) + tmpCmPos{iDim,iTraj} ./ s.clc.binW);
            end
        end
    else

        tmpCmPos = d.cmPos{iFitD};

        % The rows have to be flipped, because small == far
        tmpCmPos(1,:) = -d.cmPos{iFitD}(1,:);
        % Convert real positions to voxels
        d.binPos{iFitD} = round(s.clc.nearPos + tmpCmPos ./ s.clc.binW);
        % But for rows, use FLOOR so that the closest position isn't IN the bodypart
        d.binPos{iFitD}(1,:) = floor(s.clc.nearPos(1,:) + tmpCmPos(1,:) ./ s.clc.binW);
    end
end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Extract appropriate q values using the d.bin positions
function outQ = ExtractQ(allqVals, binPos)
outQ = [];
for iQ = 1:numel(allqVals)

    qVals = allqVals{iQ};
    outQ = [outQ ; cell2mat(arrayfun(@(pos)  qVals(:,binPos(1,pos),binPos(2,pos),binPos(3,pos)),...
    [1:size(binPos,2)], 'UniformOutput', false))];
end
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [ gofScore chiSq pVal1 pVal2 ] = GoFFun( errSq, stdDevSq, k )
%GoFFun Calculates the Goodness of Fit scores

chiSq = nansum(errSq(:)./stdDevSq(:));
gofScore = (chiSq - k) ./sqrt(2.*k);
pVal1=1-chi2cdf(chiSq,k);
pVal2=1-normcdf(gofScore);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [fitRes] = Fit3Ddat(d, fitType)
% Calculates the optimal fit for data andmodel defined in d. 
% Outputs the fitting results


% Initialise parameters
p0 = [ones([numel(unique(d.tact)) 1]); zeros([numel(unique(d.exp)) 1])  ] ; 
% Lower and upper bounds
lb = zeros(size(p0)); lb(1+end-numel(unique(d.exp)):end) = -100; % allow the offsets to be -ve
ub = ones(size(p0)) .* 1000;

% Set optimisation options - keep it simple
A = []; b = []; Aeq = []; beq = []; a = tic;

% Add aditional parameters if necessary for simple curve fitting
switch fitType
    case 'Linear'
    case 'Exponential'
        p0 = [p0;   -.02 ;   -40];
        lb = [lb;    -20 ; -1000];
        ub = [ub;      0 ;  1000];
    case 'Sigmoid'
        p0 = [p0;     -10 ;     -1];
        lb = [lb; -1000 ; -1000];
        ub = [ub;     0 ;  1000];
end

% Extra MULTIPLIER for direction [justified by empirical expectation effec: all
% away-moving conditions seem to have a higher baseline]
p0 = [p0;  -0.2248];
lb = [lb;-1000];
ub = [ub; 1000];

% Add in parameters for each entry into d which requires fitting
% reward values
fitRewards = find(cellfun(@(optSet) strcmp(optSet,'optimise_reward'), d.psiSettings(:,4)));
for iFitRew  = 1:numel(fitRewards)
    cF       = fitRewards(iFitRew);
    rewInd   = find(cellfun(@(splitType) strcmp(splitType{1},'rew') , d.psiSplit{cF}));
    rewSigns = sign(cell2mat(d.psiSplit{cF}{rewInd}{2}));
    % Add a parameter for each of the possible rewards
    for iSepRews = 1 : numel(d.psiSettings{cF,5}{1})

        % Assign possible magnitudes of initial parameters, based on 
        % originally hypothesised reward sign
        if rewSigns(d.psiSettings{cF,5}{1}(iSepRews)) == 1
            p0 = [p0;  1.1]; 
            lb = [lb;  0];
            ub = [ub;  3];
        else
            p0 = [p0;  -1.1]; 
            lb = [lb;  -3];
            ub = [ub;   0];
        end
    end
end

% Put theoretical fitted variables into a matrix for faster fitting. Also
% reformat the data into a matrix
[psiMat allDat psiLabels] = MakePsiMat(d);

% Specify optimised function
FunToOpt = @(p) ErrorFun(psiMat,allDat,d,p,fitType);

% Run optimisation
OPTIONS = optimset('TolCon',1e-10);
[p,chiSq,exitflag,output,lambda,grad,hessian] = fmincon(FunToOpt,p0,A,b,Aeq,beq,lb,ub,[],OPTIONS);
optTime = toc(a)

% Extract optimised data
[chiSqFinal, sSqErrFinal, dFitFinal, residuals, dNew, psiMatNew] = ErrorFun(psiMat,allDat,d,p,fitType);

% Store fitting data in d
d = InsertModelledDat(dNew,dFitFinal,allDat); 

% -------------------
% Make some statistical measure of goodness of fit

% First adjust all the input variance to standard error
dataVar = [d.realSTE{:}].^2;
errSq   = (allDat - dFitFinal).^2;

% Calculate the chi square confidence interval using bootstrapping
chiSqDistr  = arrayfun(@(x) nansum(randsample(errSq(:)./dataVar(:),numel(errSq(:)),1)) ,[1:100000] );
chiSqCI     = prctile(chiSqDistr,[2.5 97.5]);

N   = numel(errSq(:)); 
dof = numel(p);
k   = N - dof;

% Calculate goodness of fit
[ gofScore chiSqCheck pVal1 pVal2 ] = GoFFun( errSq, dataVar, k );

% calculate variance of residuals
resVar  = nansum(errSq(:)./ k );
% calculate log likelihood of data
logLike = -(  N                    .* log(2 .*pi .* resVar) ./ 2 )  - ...
             (1 ./ (2.* resVar))   .* nansum(errSq(:)) ;
% calculate AIC and BIC
AIC     = -2 .* logLike + dof .* 2;
BIC     = -2 .* logLike + dof .* log(N);


f.STATS.f = figure('Position',[20 20 600 300]);

% Plot the fitting results
plot(chi2pdf(0:250 ,k),'LineWidth',2); hold on
yLims = ylim;
plot([chiSqFinal chiSqFinal], yLims,'r','LineWidth',2);
xlim([0 250]);
xlabel('total error (ChiSquare)');
ylabel('probability');
legend('Expected distribution of total error if data is generated by a process like the model', ...
    'Observed total error', 'Location','South');
title(['If P > 0.05, hypothesis that Action Value is a satisfactory explanation of the data ' ...
       'cant be rejected. P = ' num2str(pVal1) ')'])

% Store fitting results
fitRes.p            = p;
fitRes.d            = d;
fitRes.sSqErrFinal  = sSqErrFinal;
fitRes.dFitFinal    = dFitFinal;
fitRes.chiSq        = chiSqFinal;
fitRes.chiSqDistr   = chiSqDistr;
fitRes.chiSqCI      = chiSqCI;
fitRes.gofScore     = gofScore;
fitRes.chiSqCheck   = chiSqCheck;
fitRes.pVal1        = pVal1;
fitRes.pVal2        = pVal2;
fitRes.resVar       = resVar; 
fitRes.logLike      = logLike;
fitRes.AIC          = AIC;
fitRes.BIC          = BIC;

end
% -------------------------------------------------------------------------


