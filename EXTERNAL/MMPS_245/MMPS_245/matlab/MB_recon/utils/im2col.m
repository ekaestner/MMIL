% Copyright (C) 2004 Josep Mones i Teixidor
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%
% -*- texinfo -*-
% @deftypefn {Function File} {@var{B} = } im2col (@var{A}, [@var{m},@var{n}], @var{block_type})
% @deftypefnx {Function File} {@var{B} = } im2col (@var{A}, [@var{m},@var{n}])
% @deftypefnx {Function File} {@var{B} = } im2col (@var{A}, 'indexed', ...)
% Rearranges image blocks into columns.
%
% @code{B=im2col(A, [m, n], blocktype)} rearranges blocks in @var{A}
% into columns in a way that's determined by @var{block_type}, which
% can take the following values:
%
% @table @code
% @item distinct
% Rearranges each distinct @var{m}-by-@var{n} block in image @var{A}
% into a column of @var{B}. Blocks are scanned from left to right and
% the up to bottom in @var{A}, and columns are added to @var{B} from
% left to right. If @var{A}'s size is not multiple @var{m}-by-@var{n}
% it is padded.
% @item sliding
% Rearranges any @var{m}-by-@var{n} sliding block of @var{A} in a
% column of @var{B}, without any padding, so only sliding blocks which
% can be built using a full @var{m}-by-@var{n} neighbourhood are taken.
% In consequence, @var{B} has @var{m}*@var{n} rows and
% (@var{mm}-@var{m}+1)*(@var{nn}-@var{n}+1) columns (where @var{mm} and
% @var{nn} are the size of @var{A}).
%
% This case is thought to be used applying operations on columns of
% @var{B} (for instance using sum(:)), so that result is a
% 1-by-(@var{mm}-@var{m}+1)*(@var{nn}-@var{n}+1) vector, that is what
% the complementary function @code{col2im} expects.
% @end table
%
% @code{B=im2col(A,[m,n])} takes @code{distinct} as a default value for
% @var{block_type}.
%
% @code{B=im2col(A,'indexed',...)} will treat @var{A} as an indexed
% image, so it will pad using 1 if @var{A} is double. All other cases
% (incluing indexed matrices with uint8 and uint16 types and
% non-indexed images) will use 0 as padding value.
%
% Any padding needed in 'distinct' processing will be added at right
% and bottom edges of the image.
%
% @seealso{col2im}
% @end deftypefn
%
% Author:  Josep Mones i Teixidor <jmones@puntbarra.com>

function B = im2col(A, varargin)
  if(nargin<2 || nargin>4)
    usage('B=im2col(B [, ''indexed''], [m,n] [, block_type])');
  end

  % check 'indexed' presence
  indexed=false;
  p=1;
  if(ischar(varargin{1}) && strcmp(varargin{1}, 'indexed'))
    if(nargin<3)
      usage('B=im2col(B [, ''indexed''], [m,n] [, block_type])');
    end
    indexed=true;
    p=p+1;
    if(isa(A,'uint8') || isa(A,'uint16'))
	padval=0;
    else
      padval=1;
    end
  else
    padval=0;
    end

  % check [m,n]
  if(~isvector(varargin{p}))
    error('im2col: expected [m,n] but param is not a vector.');
  end
  if(length(varargin{p})~=2)
    error('im2col: expected [m,n] but param has wrong length.');
  end
  m=varargin{p}(1);
  n=varargin{p}(2);
  p=p+1;

  block_type='sliding';
  if(nargin>p)
    % we have block_type param
    if(~ischar(varargin{p}))
      error('im2col: invalid parameter block_type.');
    end
    block_type=varargin{p};
    p=p+1;
    end

  % if we didn't have 'indexed' but had 4 parameters there's an error
  if(nargin>p)
      usage('B=im2col(B [, ''indexed''], [m,n] [, block_type])');
  end


  % common checks
  if(~ismatrix(A))
    error('im2col: A should be a matrix (or vector).');
  end

  switch(block_type)
    case('distinct')
      % calc needed padding
      sp=mod(-size(A)',[m;n]);

      if(any(sp))
	A=padarray(A,sp,padval,'post');
      end

      % iterate through all blocks
      B=[];
      for i=1:m:size(A,1) % up to bottom
	for j=1:n:size(A,2) % left to right
	  % TODO: check if we can horzcat([],uint8([10;11])) in a
	  % future Octave version > 2.1.58
	  if(isempty(B))
        tmp=A(i:i+m-1,j:j+n-1);
	    B=tmp(:);
      else
        tmp=A(i:i+m-1,j:j+n-1);
	    B=horzcat(B, tmp(:));
	  end
	end
      end

    case('sliding')
      if(indexed)
	disp('WARNING: ''indexed'' has no sense when using sliding.');
      end
      if(m>size(A,1) || n>size(A,2))
	error('im2col: block size can''t be greater than image size in sliding');
      end
      % TODO: check if matlab uses top-down and left-right order
      B=[];
      for j=1:1:size(A,2)-n+1 % left to right
	for i=1:1:size(A,1)-m+1 % up to bottom
	  % TODO: check if we can horzcat([],uint8([10;11])) in a
	  % future Octave version > 2.1.58
	  if(isempty(B))
        tmp=A(i:i+m-1,j:j+n-1);
	    B=tmp(:);
      else
        tmp=A(i:i+m-1,j:j+n-1);
	    B=horzcat(B, tmp(:));
	  end
	end
      end

    otherwise
      error('im2col: invalid block_type.');
  end

  end