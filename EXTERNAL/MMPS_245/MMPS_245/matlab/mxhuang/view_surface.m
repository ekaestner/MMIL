function [hf,hs,hl] = view_surface(figname,faces,verts,cdata);
%VIEW_SURFACE - Convenient function to consistently plot surfaces
% function [hf,hs,hl] = view_surface(figname,faces,verts,cdata);
% figname is the name of the figure window
% faces is the triangle listing
%       - if faces is a cell array, verts needs to be a cell array of same length
%         alpha transparency is used to visualize the multiple surfaces.
%         Order in array goes from inner to outer surface.
% verts are the corresponding vertices
%       - verts is either an array or a cell array
%       
% cdata is the colordata to use.  If not given, uses random face color
%       - if cdata is a cell array, use on cell per surface for color coding
% hf is the figure handle used
% hs is the handles to the surfaces
% hl is the handles to the lighting

%<autobegin> ---------------------- 26-May-2004 11:34:42 -----------------------
% --------- Automatically Generated Comments Block Using AUTO_COMMENTS ---------
%
% CATEGORY: Visualization
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\windclf.m
%   toolbox\windfind.m
%
% At Check-in: $Author: Mosher $  $Revision: 18 $  $Date: 5/26/04 10:02a $
%
% This software is part of BrainStorm Toolbox Version 2.0 (Alpha) 24-May-2004
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
%<autoend> ------------------------ 26-May-2004 11:34:42 -----------------------

% /---Script Authors-------------------------------------\
% |                                                      |
% |  *** John C. Mosher, Ph.D.                           |
% |  Biophysics Group                                    |
% |  Los Alamos National Laboratory                      |
% |  Los Alamos, New Mexico, USA                         |
% |  mosher@lanl.gov                                     |
% |                                                      |
% |  *** Sylvain Baillet, Ph.D.                          |
% |  Cognitive Neuroscience & Brain Imaging Laboratory   |
% |  CNRS UPR640 - LENA                                  | 
% |  Hopital de la Salpetriere, Paris, France            |
% |  sylvain.baillet@chups.jussieu.fr                    |
% |                                                      |
% \------------------------------------------------------/

% Script History ----------------------------------------------------------------------------------------
%
% SB  10-Mar-2003 : faces  and verts input arguments can now be cell arrays thereby yielding 3D plots 
%                   using alpha transparency for each surface
% SB  03-Jun-2003 : changed axis management and default view point in 3D
% JCM 19-Aug-2003 : updated comments to explain outputs
% SB  21-Oct-2003 : Basic lightning is 'none'
% JCM 11-May-2004 : reset lighting to shiny, CBB would be an external switch
% --------------------------------------------------------------------------------------------------------

if iscell(faces) % Multiple plots requested on same figure window
    if length(faces) ~= length(verts) % sanity check
        errordlg('Faces and Vertices need to be cell arrays of same length', mfilename);
        return
    end
else
    faces = {faces};
    verts = {verts};
    if nargin > 3
        cdata = {cdata};
    end
    
end

if nargin == 4
    if ~iscell(cdata)
        cdata = {cdata};
    end
end
    

for k=1:length(faces) % For each requested surface
    
    if(size(verts{k},2) > 3), % assume transposed
        verts{k} = verts{k}';  % if the assumption is wrong, will crash below anyway
    end
    
    if nargin == 3 
        cdata{k} = repmat(rand(1,3),size(verts{k},1),1);
    end
    
    h = windfind(figname);
    
    figure(h)
    if isempty(h)
        windclf
        hold on
    end
    
    hs(k) = patch('faces',faces{k},'vertices',verts{k},...
        'facevertexcdata',cdata{k},'facecolor','interp','edgecolor','none');
    
    if length(faces)>1
        set(hs(k),'FaceAlpha',k/10);
    end
    
    view(120,30)
    axis equal, axis vis3d
    axis off
    if nargin == 3
        colormap(bone(256))
     end
     
     switch 'shiny'
        case 'dull'
           material dull
           lighting none
        case 'shiny'
           material metal
           lighting phong
     end
     
    if k ==1 | isempty(findobj(h,'type','light')) % avoid accumulating too many light objects in same window
        hl(1) = camlight(-20,30);
        hl(2) = camlight(20,30);
        hl(3) = camlight(-20,-30);
        for i = 1:length(hl),
            set(hl(i),'color',[.8 1 1]/length(hl)/1.2); % mute the intensity of the lights
        end
    end
    
    if(nargout>0),
        hf = h;  % only if the user has output argument
    end
    
end
