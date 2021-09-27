%Script_Generate_FMBanks
%
% This script generates FM-bank stimuli with randomly generated parameters
% as described in:
%
%   [1] Max F. K. Happel, Simon Müller, Jörn Anemüller, and Frank W.
%   Ohl. Predictability of strfs in auditory cortex neurons depends on
%   stimulus class. In INTERSPEECH-2008, page 670, 2008.
% 
%   [2] Arne F. Meyer, Max F. K. Happel, Frank W. Ohl, and Jörn Anemüller.
%   Estimation of spectro-temporal receptive fields based on
%   linear support vector machine classification. BMC Neuroscience,
%   10(Suppl 1):P147, 2009.
%
% See also genFMSweep Script_Generate_FMBanks 

% Author: arne.f.meyer@uni-oldenburg.de
% Date  : 17-Nov-2010 10:17:57

close all
clear
clc


% -----------------------------------------------
% Script settings
% -----------------------------------------------
bSavePlots = true;     % save plots?
bSaveStim  = true;     % save generated stimuli?
sDataDir   = fullfile(fileparts(mfilename('fullpath')), 'Output');
if (bSavePlots || bSaveStim) && ~exist(sDataDir, 'dir')
    mkdir(sDataDir);
end


% -----------------------------------------------
% Generation parameters
% -----------------------------------------------

% Sweep parameters
nMaxFreq  = 16000;
nMinFreq  = 250;
nSweepLen = 0.1;
nRampLen  = 0.005;
sRampType = 'lin';
nGain     = 0.1;

% Stimulus parameters
nSweeps   = 4;
nStimLen  = 10;
nFs       = 44100;


% -----------------------------------------------
% Generate FM-bank stimuli by randomly drawing
% sweep frequencies from an uniform distribution
% -----------------------------------------------
nSweepSamp   = ceil(nSweepLen * nFs) + 1;    % Sweeps start at t=0!
nSweepBlocks = ceil(nStimLen / nSweepLen); 
nSamples     = nSweepSamp * nSweepBlocks;
vStim        = zeros(nSamples,1);
for iBlock=1:nSweepBlocks
    
    % Indices of the current sweep block samples
    vBlockIdx = (iBlock-1)*nSweepSamp+1 : iBlock*nSweepSamp;
    
    % Generate nSweeps sweeps for every sweep block
    for iSweep=1:nSweeps
        
        % Draw frequencies but don't use pure tones!
        nStartFreq = rand() * (nMaxFreq - nMinFreq) - nMinFreq;
        while 1
            nStopFreq  = rand() * (nMaxFreq - nMinFreq) - nMinFreq;
            if nStopFreq ~= nStartFreq
                break;
            end
        end
        
        % Generate sweep
        vStim(vBlockIdx) = vStim(vBlockIdx) + nGain * genFMSweep(nStartFreq, nStopFreq, nSweepLen, nFs, nRampLen, sRampType);
    end
end

if bSaveStim
    wavwrite(vStim, nFs, 16, fullfile(sDataDir, 'Stimulus_FMBanks.wav'));
end


% -----------------------------------------------
% Plotting
% -----------------------------------------------

% Stimulus spectrogram
nWinLen               = 0.01;
nOverlap              = nWinLen / 2;
[mSpec, vFreq, vTime] = spectrogram(vStim, ceil(nWinLen * nFs), ...
    round(nOverlap * nFs), [], nFs);
nFsSpec               = 1 / (vTime(2) - vTime(1));
nPlotLen              = 1;
vPlotIdx              = 1:ceil(nPlotLen * nFsSpec);

figure('Name','Stimulus spectrogram');
imagesc(vTime(vPlotIdx), vFreq, 20*log10(abs(mSpec(:,vPlotIdx))))
axis xy; colorbar;
xlabel('Time / s', 'FontSize', 16), xlim([0 nPlotLen])
ylabel('Frequency / Hz', 'FontSize', 16), ylim([nMinFreq nMaxFreq])
set(gca(), 'FontSize', 14)
if bSavePlots
    saveas(gcf(), fullfile(sDataDir, 'Spectrogram_FMBanks'), 'psc2');
end


% Frequency-wise autocorrelation function
nMaxLagSec = 0.2;
nMaxLag    = ceil(nMaxLagSec * nFsSpec);
mCorr      = zeros(size(mSpec,1), 2*nMaxLag+1);
for iFreq=1:size(mSpec,1)
    mCorr(iFreq,:) = xcorr(abs(mSpec(iFreq,:)), abs(mSpec(iFreq,:)), nMaxLag);
    mCorr(iFreq,:) = mCorr(iFreq,:) / max(mCorr(iFreq,:));
end

figure('Name','Frequency-wise stimulus autocorrelation');
imagesc((-nMaxLag:nMaxLag) / nFsSpec, vFreq, mCorr);
axis xy; colorbar
xlabel('Time Lag / s', 'FontSize', 16), xlim(nMaxLagSec * [-1 1])
ylabel('Frequency / Hz', 'FontSize', 16), ylim([nMinFreq nMaxFreq])
set(gca(), 'FontSize', 14)
if bSavePlots
    saveas(gcf(), fullfile(sDataDir, 'Autocorrelation_FMBanks'), 'psc2');
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
