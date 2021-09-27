function Amp_TDT = getTDTAmpClick(AmpClick, SpeakerName)
% Get the amplitude corresponding to TDT output for click stim.
% AmpClick: amplitude order for trials (in dB)

[calib_V, calib_dB] = getCalib_Click(SpeakerName);
dBu = unique(AmpClick);
Vu = interp1(calib_dB, calib_V, dBu);

Amp_TDT = zeros(size(AmpClick));
for k = 1:length(AmpClick)
	amp = AmpClick(k);
	Amp_TDT(k) = Vu(dBu==amp);
end
