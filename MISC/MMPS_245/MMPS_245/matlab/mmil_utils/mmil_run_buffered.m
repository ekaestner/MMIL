function mmil_run_buffered( fn, varargin )
%function mmil_groupstats_runanalysis( fn, varargin )
%   Now the catch: MEMORY
%
%   We're going to go into a loop of trying-retrying gathering data and
%   running analyses.  We're going to keep reducing the memory
%   requirements of the analyses until we find an analysis that can be
%   done within the memory requirements of this computer.
%
%   For speed's sake, this slow "discovery" stage of processing can be
%   skipped by entering manually the number of frames to handle at
%   once.
%
%   NOTE:
%       Functions that take advantage of this memory-retry logic MUST
%       return information about which frames they processed.  If an
%       analysis does not return such information, then the memory
%       retry logic cannot work, and so is not run.
%

%   Now, some documentation on memory management:
%
%   frames            - frames to process
%   framebuffersize   - number of frames to process (e.g. max size of
%                       'frames')
%   
%
%

  parms = mmil_args2parms(varargin, [], false);
  
  %   List the groups in the 
  %mmil_logstr('Groups included in surface group average: '); 
  %mmil_logstr('\t(%s)', parms.groups{:});

  % This 'try' block should be enabled for debugging
  %try

    % Execute the code based on the # of output arguments.
    %   0 output args means the function doesn't support
    %   any frame-buffering, so we have to pass all data in.
    %
    %   1 output arg means the function DOES support
    %   frame-bufferings, so we should do that.
    switch (nargout(fn))

      % Execute with no frame-buffering
      case 0
        if (nargin(fn) < 0)
          fn(varargin{:});
        else
          fn(parms);
        end;

      % Execute with frame-buffering.
      case 1

        if (parms.framebuffersize ~= 0)
          framebuffersize    = parms.framebuffersize;
        else
          framebuffersize = length(parms.frames);
        end;

        % Initialize the current frames, making sure
        % that the frames stay within the allowed range.
        frames          = parms.frames(1):min(parms.frames(1)+framebuffersize-1, max(parms.frames));

        %   This is the retry loop, where we retry analyses at
        %   smaller and smaller frame buffer sizes, until we hit 0.
        while (~isempty(frames) || framebuffersize > 0)

          %try
            % determine output log string based on
            % the # of retries
            %
            % First try, no buffering even expected to be needed
            if (parms.framebuffersize==0 && length(frames)==length(parms.frames))
              mmil_logstr(parms, 'Running %s ...', func2str(fn)); 

            % First try, buffering expected
            elseif (parms.framebuffersize==framebuffersize)
              mmil_logstr(parms, 'Running %s, loading %d frame(s) at a time...', func2str(fn), length(frames));

            % Retry after failure
            else
              mmil_logstr(parms, '** re-running %s, loading %d frames at a time **\n', ...
                            func2str(fn), length(frames));
            end;

            % create a new parms structure that
            % also specifies the current frame buffering
            % status
            giparms                     = parms;
            giparms.framebuffersize     = framebuffersize;
            giparms.allframes           = parms.frames;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %   This is the framebuffer loop, where we traverse
            %   the frames to be analyzed until we've walked
            %   through every frame.
            while (~isempty(frames))

              % Process it!!
              giparms.frames      = frames;
              
              % run the function
              if (nargin(fn) < 0)
                framesprocessed   = fn(varargin{:});
              else
                framesprocessed   = fn(giparms);
              end;

              %   Error?
              if (isempty(framesprocessed))
                break;
              end;

              %   Didn't process all the frames; try again!
              if (~isempty(setdiff(frames, framesprocessed)))
                frames  = setdiff(frames, framesprocessed);

              % ???More frames were processed than were told????
              else%if (~isempty(setdiff(framesprocessed,frames)))
                %frames  = setdiff(framesprocessed, frames);

              %   Processed all the frames.  Grab the next batch
              %   next batch.
              %else
                %Determine a string to output
                if (length(frames)==1)
                  mmil_logstr(parms, 'Finished running frame %d', frames(1));
                elseif (isempty(frames))
                  mmil_logstr(parms, 'Finished running all frames.');
                else
                  mmil_logstr(parms, 'Finished running frames %d to %d', frames(1), frames(end));
                end;

                % Get the next set of frames to process, 
                frames = [frames(end)+1:frames(end) + framebuffersize];
                frames = intersect(frames, parms.frames);
              end;

            end; %while processing more frames!

            %   Successful (e.g. no OOM errros), so get outta the retry
            %   memory-retry loop
            break;

          %catch

            %   If we're out of memory, we've got to re-run
            %   the analysis with a lower memory load.
            %
            %   Try re-running this analysis, but by divvying up 
            %   the number of frames in any analysis into smaller 
            %   parcellations
            %if (~isempty(findstr(lower(lasterr), 'out of memory')))
            %  mmil_error(parms, lasterr);
            %end;

            %   Retry with a smaller number of frames per analysis, starting
            %   from scratch.
            framebuffersize     = min(floor(framebuffersize/2), 128);

            %   Resume where we left off
            if (~isempty(frames))
              frames     = [frames(1): min(frames(1)+framebuffersize-1, max(parms.frames))];
            end;
          %end;
        end; % while retry on memory errors

      otherwise
        help('ts_groupstats_runanalysis');
        mmil_error(parms,  'groupaverage analysis functions must have either 0 or 1 outputs.');

    end; % switch on # outputs

%    catch
%      e   = lasterr;
%    end;

  if (isempty(frames) && framebuffersize==0)
    mmil_error(parms, 'Couldn''t run analysis; ran out of memory even on a single frame!');
  end;

  %   Re-throw any error
  if (exist('e', 'var'))
      mmil_error(parms, e);
  end;


