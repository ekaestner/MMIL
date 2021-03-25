function [NFV, status] = tessellation_outwards(FV);
%TESSELLATION_OUTWARDS - Ensure that tessellation ordering is "outwards"
% function [NFV, status] = tessellation_outwards(FV);
%
% For a given faces-vertices structure FV, test to see if the majority of
% tesselations are ordered such that the right-hand (counter clockwise) ordering
% yields surface normals outwards from the boundary.
%
% NFV is the new FV structure. If status = 0, then tessellation was already
% ordered outwards, and NFV = FV. If status = 1, then the normals were pointed
% inwards, and the [1 2 3] orderings were reversed to be [1 3 2], reversing the
% normals.
%
% Assumes that the surface represented by FV is a closed manifold, and that the
% majority of triangle normals point in the same direction. This may not be true
% in a highly convoluted surface, but should generally hold for the smoother
% surfaces used in BEMs (inner skull, outer skull, scalp).
%
% This function ASSUMES that all triangles are already consistently ordered in
% the same direction, and does a fast GLOBAL test to see which direction the
% normals are pointed.  If the majority of normals are found to be pointing
% inwards to the center of the vertices, then the ordering of ALL triangles are
% reversed. Individual triangle orderings are not examined. To check the
% ordering of each triangle, use the more extensive tessellation_stats, which
% locally tests every triangle to see if it is consistently ordered relative to
% its neighbors. 
%
% See also TESSELLATION_STATS

%<autobegin> ---------------------- 14-Jun-2004 17:12:40 -----------------------
% --------- Automatically Generated Comments Block Using AUTO_COMMENTS ---------
%
% CATEGORY: Visualization
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\tessellation_stats.m
%
% At Check-in: $Author: Mosher $  $Revision: 4 $  $Date: 6/14/04 3:38p $
%
% This software is part of BrainStorm Toolbox Version 2.0 (Alpha) 14-Jun-2004
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2004 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 14-Jun-2004 17:12:40 -----------------------


% ----------------------------- Script History ---------------------------------
% 27-May-2004 JCM Creation
% ----------------------------- Script History ---------------------------------


[mF,nF] = size(FV.faces);
[mV,nV] = size(FV.vertices);

if nF ~= 3 | nV ~= 3,
   error(sprintf('Fields .vertices and .faces must contain 3 columns'));
end

% map to the output

NFV = FV; 

% Find the center of the tesselations
BoundedCenter = mean(FV.vertices);

% remove the mean from all vertices
FV.vertices = FV.vertices - BoundedCenter(ones(mV,1),:);

% Calculate the face normals of these centered tessellations
[FaceNormal, FaceArea, FaceCenter] =  tessellation_stats(FV,0);

% Each face center and face normal is now in the centered coordinates.
% Easy to see which way these normals mostly face.

% sign of the dot product
Direction = sign(sum(FaceNormal .* FaceCenter)); 

if sum(Direction) < 0,
   % the majority of directions are inwards, reverse the ordering
   status = 1; % reversing the direction
   NFV.faces = NFV.faces(:,[1 3 2]);
else
   status = 0; % no reversal needed
   % the majority are outwards, do nothing
end

