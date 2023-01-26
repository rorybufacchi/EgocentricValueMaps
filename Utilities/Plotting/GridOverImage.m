function [] = GridOverImage(fS,cA)
%GRIDOVERIMAGE puts a grid of lines over an image to make it look like
%separate blocks
%   Detailed explanation goes here

hold on;

xLIMS = xlim(cA); yLIMS = ylim(cA);

if ~isfield(fS,'color')
    fS.color = 'w';
end


xPos=fS.gridXstart:fS.gridXstep:xLIMS(end);
yPos=fS.gridYstart:fS.gridYstep:yLIMS(end);

% for iX=1:length(xPos)
%     lin = vline(xPos(iX)-0.05,'k','LineWidth',1.5);
%     lin = vline(xPos(iX)-0.05,'k','LineWidth',1.5);
% %     uistack(lin, 'top');
% end
% for iY=1:length(yPos)
%     lin = hline(yPos(iY)-0.05,'k','LineWidth',1.5);
%     lin = hline(yPos(iY)-0.05,'k','LineWidth',1.5);
% %     uistack(lin, 'top');
% end

for iX=1:length(xPos)
    lin = vline(xPos(iX),'Color',fS.color,'LineWidth',0.5);
    uistack(lin, 'top');
end
for iY=1:length(yPos)
    lin = hline(yPos(iY),'Color',fS.color,'LineWidth',0.5);
    uistack(lin, 'top');
end

end

