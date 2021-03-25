function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,255) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',255)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,255)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex and colorspec. 
% 

if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

r = hex2dec(hex(:,1:2)); 
g = hex2dec(hex(:,3:4));
b = hex2dec(hex(:,5:6));

if ~exist('range','var')
    range = 1; 
end

if range == 255
    rgb = [r g b]; 
else 
    rgb = [r g b]/255; 
end

end

