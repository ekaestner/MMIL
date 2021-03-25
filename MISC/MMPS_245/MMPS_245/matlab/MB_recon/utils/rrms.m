function r = rrms(Iref, Irecon)
%
% function r = rrms(Iref, Irecon)
%
% Calculates the relative root mean squared error (RRMS) between two images.
% Formula as in Anja et al. MRM 2008; 59: 382-395.
%
% Inputs
%   Iref   - (Stack of) reference images.
%   Irecon - (Stack of) reconstructed images.
%
% Output
%   r      - The rrms values.
%
% (c) Kangrong Zhu      Stanford University     Nov 2013

r = sqrt(sum(sum(abs(Iref-Irecon).^2,1),2)./(sum(sum(abs(Iref).^2,1),2)+eps));

return