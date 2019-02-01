% DEER resampling into number of points required by the 
% neural network. Syntax:
%
%      deer_trace=deer_resample(deer_trace,out_pts)
% 
% Parameters:
% 
%   deer_trace - input DEER trace, column vector
%
%     out_pts  - number of points in output trace
%
% Output:
% 
%   deer_trace - resampled DEER trace, column vector
% 
% i.kuprov@soton.ac.uk
% s.g.worswick@soton.ac.uk

function deer_trace=deer_resample(deer_trace,out_pts)

% Check consistency
grumble(deer_trace,out_pts);

% Apply a matched Savitzky-Golay filter
frame_length=ceil(numel(deer_trace)/out_pts);
if mod(frame_length,2)==0
    frame_length=frame_length+1;
end
if frame_length>=3
    if exist('sgolayfilt','file') % sgolayfilt is from Signal Processing Toolbox
        deer_trace=sgolayfilt(deer_trace,2,frame_length);
    else
        % TODO: replace with toolbox-free implementation of Savitzky-Golay filter
        deer_trace=sgolayfilt(deer_trace,2,frame_length);
    end
end

% Symmetrise and shift data
deer_trace=[flip(deer_trace); deer_trace];
data_shift=deer_trace(1);
deer_trace=deer_trace-data_shift;

% Run the resampling
if exist('resample','file') % resample is from Signal Processing Toolbox
    deer_trace=resample(deer_trace,(out_pts*2),numel(deer_trace));
else
    deer_trace = interp1(1:numel(deer_trace),deer_trace,1:out_pts*2).';
end

% Desymmetrise data and shift back
deer_trace=deer_trace((end/2+1):end);
deer_trace=deer_trace+data_shift;

end

% Consistency enforcement
function grumble(deer_trace,out_pts)
if (~isnumeric(deer_trace))||(~isreal(deer_trace))||...
   (~iscolumn(deer_trace))
    error('deer_trace must be a real column vector.');
end
if (~isnumeric(out_pts))||(~isreal(out_pts))||...
   (out_pts<1)||(mod(out_pts,1)~=0)
    error('out_pts must be a non-negative real integer.');
end
end

% A notable American commentator, Charles Krauthammer, once 
% explained Rupert Murdoch's success in founding Fox News,
% a cable channel, by pointing out that he had found a niche
% market - half the country.
%
% The Economist

