function  str = mmil_padstr(o, l, p, d)
%function str = mmil_padstr(o, l, p, [d])
%
% Purpose:
%   pad some object with leading or trailing characters
%
% Input Arguments:
%   o   : initial object to pad; currently supports integers
%   l   : desired length of output string
%   p   : string to pad with
%         {default: 0}
%   d   : direction to pad (left,right, or {left,right})
%         {default: left}
%
% Output Arguments:
%   str : padded string of length l
%
% Examples:
%   padstr(1,5,'x','left')            => xxxx1
%   padstr(1,5,'x','right')          => 1xxxx
%   padstr(1,5,'x',{'left','right'}) => xx1xx
%
% Created:  08/01/07 by Ben Cipollini
% Last Mod: 11/07/11 by Don Hagler
%

if (~exist('o','var')) help(mfilename); end;
if (~exist('l','var')) help(mfilename); end;

% default args dependent on object-type
if (mmil_isint(o))
  o = sprintf('%d',o); % convert to string

  if (~exist('p','var')),   p='0';      end;
  if (~exist('d','var')),   d={'left'}; end;
end;

if (~iscell(d))
  d={d};
end;

ntopad = (l - length(o))/(length(d)*length(p));
str    = '';

if (ismember('left',d))  str    = [str repmat(p, [1 ceil(ntopad)])]; end;
str = [str o];
if (ismember('right',d)) str    = [str repmat(p, [1 floor(ntopad)])]; end;

