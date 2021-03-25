function str = concat_ws(sep, formatstr, varargin)
%
% See MySql's "concat_ws"
%

  str = '';
  
  for i=1:length(varargin)
    
    % Make sure there are no cells in here
    if (iscell(varargin{i}))
      %error('function does not accept cell arguments');
      concat_ws(sep, varargin{i}{:});
    end;
    
    if (i==1)
      str = sprintf(formatstr, varargin{1});
    else
      str = [str sep sprintf(formatstr, varargin{i})];
    end;
  end;