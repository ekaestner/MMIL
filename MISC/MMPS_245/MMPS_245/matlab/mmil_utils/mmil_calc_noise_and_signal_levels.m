function [noise,signal] = mmil_calc_noise_and_signal_levels(vol)
%function [noise,signal] = mmil_calc_noise_and_signal_levels(vol)
%
% Purpose: calculate noise and signal levels from image volume
%
% Created:  01/08/15 by Don Hagler
% Last Mod: 01/08/15 by Don Hagler
% 

% based on ctx_ComputeSignalAndNoiseLevels by Anders Dale

noise = []; signal = [];

vol = vol.^2;
thresh = max(vol(:));
for iter = 1:10
  noise = sqrt(mean(vol(find(vol<thresh & vol>0))));
  thresh = 3*noise^2;
end
signal = sqrt(mean(vol(find(vol>25*noise^2))));

