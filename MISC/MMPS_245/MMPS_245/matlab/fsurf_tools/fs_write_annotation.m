function fs_write_annotation(fname,vertices,label,ctab)
%function fs_write_annotation(fname,vertices,label,ctab)
%
% Required Input:
%   fname: output file name
%   vertices: vector of 0-based vertex numbers
%   label: vector of annotation (color code numbers)
%     calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0
%   ctab: color table struct
%     ctab.numEntries = number of Entries
%     ctab.orig_tab = name of original ct
%     ctab.struct_names = list of structure names (e.g. central sulcus and so on)
%     ctab.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
%       is b, 4th column is flag, 5th column is resultant integer values
%     calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0
%
% Note: only writes version 2
%
% Created:  10/31/11 by Don Hagler
% Last Mod: 11/29/12 by Don Hagler
%

% based on write_annotation distributed with freesurfer

% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================


fp = fopen(fname, 'w', 'b');

% first write vertices and label

count = fwrite(fp, int32(length(label)), 'int');
if(count~=1)
   error('write_annotation: Writing #vertices/labels not successful');
end

temp = zeros(length(label)*2,1);
temp(1:2:end) = vertices;
temp(2:2:end) = label;
temp = int32(temp);

count = fwrite(fp, int32(temp), 'int');
if(count~=length(temp))
   error('write_annotation: Writing labels/vertices not successful');
end

% write that ctab exists
count = fwrite(fp, int32(1), 'int');
if(count~=1)
   error('write_annotation: unable to write flag that ctab exists');
end

%write version number
count = fwrite(fp, int32(-2), 'int');
if(count~=1)
    error('write_annotation: unable to write version number');
end

%write number of entries
count = fwrite(fp, int32(ctab.numEntries), 'int');
if(count~=1)
    error('write_annotation: unable to write number of entries in ctab');
end

%write original table
orig_tab = [ctab.orig_tab char(0)];
count = fwrite(fp, int32(length(orig_tab)), 'int');
if(count~=1)
    error('write_annotation: unable to write length of ctab source');
end

count = fwrite(fp, orig_tab, 'char');
if(count~=length(orig_tab))
    error('write_annotation: unable to write orig_tab');
end

%write number of entries
count = fwrite(fp, int32(ctab.numEntries), 'int');
if(count~=1)
    error('write_annotation: unable to write number of entries in ctab');
end

%write ct
for i = 1:ctab.numEntries
  count = fwrite(fp, int32(i-1), 'int');
  if(count~=1)
    error('write_annotation: unable to write structure number');
  end

  structure_name = [ctab.struct_names{i} char(0)];
  count = fwrite(fp, int32(length(structure_name)), 'int');
  if(count~=1)
    error('write_annotation: uable to write length of structure name');
  end
  count = fwrite(fp, structure_name, 'char');
  if(count~=length(structure_name))
    error('write_annotation: unable to write structure name');
  end

  count = fwrite(fp, int32(ctab.table(i, 1)), 'int');
  if(count~=1)
    error('write_annotation: unable to write red color'); 
  end

  count = fwrite(fp, int32(ctab.table(i, 2)), 'int');
  if(count~=1)
    error('write_annotation: unable to write blue color'); 
  end

  count = fwrite(fp, int32(ctab.table(i, 3)), 'int');
  if(count~=1)
    error('write_annotation: unable to write green color'); 
  end

  count = fwrite(fp, int32(ctab.table(i, 4)), 'int');
  if(count~=1)
    error('write_annotation: unable to write padded color'); 
  end  
end

fclose(fp);





