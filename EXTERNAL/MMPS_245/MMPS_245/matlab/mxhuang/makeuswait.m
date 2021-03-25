function [] = makeuswait(action)
%MAKEUSWAIT - change the figure pointer according to action
% function [] = makeuswait(action)
% action  = 'start' : pointers turn to 'watch'
% action  = 'stop' : pointers turn to 'arrow'

%<autobegin> ---------------------- 26-May-2004 11:30:51 -----------------------
% --------- Automatically Generated Comments Block Using AUTO_COMMENTS ---------
%
% CATEGORY: GUI and Related
%
% At Check-in: $Author: Mosher $  $Revision: 6 $  $Date: 5/26/04 9:59a $
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
%<autoend> ------------------------ 26-May-2004 11:30:51 -----------------------



% /---Script Authors-------------------------------------\
% |  *** Sylvain Baillet, Ph.D.                          |
% |  Cognitive Neuroscience & Brain Imaging Laboratory   |
% |  CNRS UPR640 - LENA                                  | 
% |  Hopital de la Salpetriere, Paris, France            |
% |  sylvain.baillet@chups.jussieu.fr                    |
% \------------------------------------------------------/
%  
% Script History ---------------------------------------------------------------
% SB  11 Mar-2004 Creation 
% Script History ---------------------------------------------------------------

switch(action)
case 'start'
    set(findobj(0,'type','figure'),'pointer','watch')
case 'stop'
    set(findobj(0,'type','figure'),'pointer','arrow')
otherwise
    return
end
