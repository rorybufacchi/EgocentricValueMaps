GW = createGridWorld(9,9)

%% Now, set the initial, terminal and obstacle states.


GW.CurrentState = '[5,9]'; % Put the wasp at position 5,9
% % % GW.TerminalStates = '[5,5]';
% % % GW.ObstacleStates = ["[3,3]";"[3,4]";"[3,5]";"[4,3]"];


%% Next, define the rewards in the reward transition matrix.

nS = numel(GW.States);
nA = numel(GW.Actions);
GW.R = -1*ones(nS,nS,nA);
GW.R(state2idx(GW,"[2,4]"),state2idx(GW,"[4,4]"),:) = 5;
GW.R(:,state2idx(GW,GW.TerminalStates),:) = 10;

% $$$ NOTE: I should make each action that isn't staying still slightly
% sightly punishing