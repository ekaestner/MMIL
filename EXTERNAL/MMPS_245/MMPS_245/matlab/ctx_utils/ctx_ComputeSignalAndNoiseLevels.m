function [noiselevel,signallevel] = ctx_ComputeSignalAndNoiseLevels(vol)

vol = vol.^2;
thresh = max(vol(:));
for iter = 1:10
  noiselevel = sqrt(mean(vol(find(vol<thresh&vol>0))));
  thresh = 3*noiselevel^2;
%   x = thresh*[0:(10000-1)]/10000; [hc] = hist(vol(find(vol<x(end)&vol>0)),x);
%   figure(1001); bar(x,hc);
%  drawnow; pause(1);
end
signallevel = sqrt(mean(vol(find(vol>25*noiselevel^2))));
% x = 10*signallevel.^2*[0:(10000-1)]/10000; [hc] = hist(vol(find((vol>3*thresh))),x);
% figure(1002); bar(x,hc);
% noiselevel
% signallevel

%keyboard
