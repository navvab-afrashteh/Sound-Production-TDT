function Amp_TDT = getTDTAmp(F, FreqOrder, AmpOrder, SpeakerName)
% Get the amplitude corresponding to TDT output.
% F: all different frequencies (in kHz)
% FreqOrder: frequency order for trials (in kHz)
% AmpOrder: amplitude order for trials (in dB)

[calib_V, calib_level, calib_freq, calib_amp] = getCalib(SpeakerName);
calib_freq = calib_freq/1000; % convert to kHz
calib_amp = 10.^(calib_amp/20)*calib_V; % convert to linear and multiplied by calibration voltage level

Amp = interp1(calib_freq, calib_amp, F);

Amp_TDT = zeros(size(AmpOrder));
for r = 1:size(FreqOrder,1)
    for c = 1:size(FreqOrder,2)
        f = FreqOrder(r,c);
        [~,idx] = min(abs(F - f));
        a = Amp(idx);
        Amp_TDT(r,c) = a*10^((AmpOrder(r,c) - calib_level)/20);
    end
end
n=0;