function Amp_TDT = getTDTAmp(F, FreqOrder, AmpOrder, SpeakerName)
% Get the amplitude corresponding to TDT output.
% F: all different frequencies (in kHz)
% FreqOrder: frequency order for trials (in kHz)
% AmpOrder: amplitude order for trials (in dB)

[calib_V, calib_level, calib_freq, calib_amp] = getCalib(SpeakerName);
calib_freq = calib_freq/1000; % convert to kHz
calib_amp = 10.^(calib_amp/20)*calib_V; % convert to linear and multiplied by calibration voltage level

% calib_freq = [2    2.828 4    5.657 8    11.314 16   22.627 32   45.255 64   90 ]; 
% calib_amp = [0.17   0.21  0.26  0.13   0.15  0.18    0.12   0.05   0.07  0.03    0.053   0.2 ]*calib_V;

Amp = interp1(calib_freq, calib_amp, F);

Amp_TDT = zeros(size(AmpOrder));
for k = 1:length(FreqOrder)
	f = FreqOrder(k);
	[~,idx] = min(abs(F - f));
	a = Amp(idx);
	Amp_TDT(k) = a*10^((AmpOrder(k) - calib_level)/20);
end
n=0;