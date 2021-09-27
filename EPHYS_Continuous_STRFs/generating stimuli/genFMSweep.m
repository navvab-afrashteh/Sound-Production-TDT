function [vSweep varargout] = genFMSweep(nStartFreq, nStopFreq, nLenSec, nFs, varargin)
%genFMSweep - generates a single FM sweep 
%
% Syntax:
%   vSweep = genFMSweep(nStartFreq, nStartFreq, nLenSec, nFs)
%   vSweep = genFMSweep(nStartFreq, nStartFreq, nLenSec, nFs, nRampLen, sRampType)
%   [vSweep vTime] = ...
%
% Input arguments:
%   nStartFreq - The starting frequency in Hz
%   nStopFreq  - The ending frequency in Hz
%   nLenSec    - The sweep length in seconds
%   nFs        - The sampling frequency in Hz
%   nRampLen   - The length of the ramp that is used at the beginning and 
%                ending of the sweep to avoid clicks (default: 0.01s)
%   sRampType  - The ramp type which can be 'lin', 'sin2' or 'non' 
%                (default: 'lin')
%
% Output arguments:
%   vSweep     - A vector with sweep coefficients
%   vTime      - A vector with time points corresponding to the samples in
%                vSweep
%
% See also chirp

% Check number of input arguments
error(nargchk(4,6,nargin))

% Parse input arguments
nRampLen  = 0.01;
sRampType = 'lin';
if nargin > 4
    nRampLen = varargin{1};
end
if nargin > 5
    sRampType = varargin{2};
end

% Create sweep signal
vTime  = 0:1/nFs:nLenSec;
vSweep = chirp(vTime, nStartFreq, nLenSec, nStopFreq)'; % use MATLAB chirp function instead of x = sin(2*pi*f*t) because MATLAB does not compute the frequency correct that way

% Create ramp
nSamples = round(nRampLen / nFs);
switch lower(sRampType)
    case 'lin'
        vRamp = linspace(0,1,nSamples);
        
    case 'sin2'
        vRamp = sin( 0:(pi/2)/nSamples:pi/2 ).^2;
        
    case 'none'
        vRamp = ones(nSamples, 1);
        
    otherwise
        error('Invalid ramp type: %s', sRampType)
end
if isempty(vRamp)
    vRamp = 1;
end
% Apply ramps
vSweep(1:length(vRamp))         = vSweep(1:length(vRamp)) .* vRamp;
vSweep(end-length(vRamp)+1:end) = vSweep(end-length(vRamp)+1:end) .* fliplr(vRamp);

% Assign optional output argument
if nargout
   varargout{1} = vTime; 
end


% ---------------------------------------------------------------------
% Copyright (c) 2010 Arne F. Meyer, Max F.K. Happel, Jan Diepenbrock, 
%                    Frank W. Ohl, Jörn Anemüller
% 
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:
% 
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
% LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% ---------------------------------------------------------------------