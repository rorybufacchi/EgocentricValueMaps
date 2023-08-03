function [npks pks locs] = RowPeakFind(s,Q)
    Qr = RecenterQNA(s,Q);
    Qline = squeeze(nanmean(Qr(1,:,:,:,:,:),[4]));
    for iLC = 2:size(Qline,2)-1
%         [pks{iLC-1} locs{iLC-1}] = findpeaks(Qline(iLC,:,3));
        [pks{iLC-1} locs{iLC-1}] = findpeaks(nanmean(Qline(iLC,:,:),3));
        npks(iLC-1) = length(locs{iLC-1});
    end
    
end