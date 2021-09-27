function AmpTDT = genAmpsTDT_PureTones(Freqs,Amp,SpeakerName)

nFreqs = length(Freqs);
AmpTDT = cell(nFreqs,1);
for f = 1:nFreqs
    F = Freqs(f)/1000;
    AmpOrder = Amp{f};
    FreqOrder = F*ones(size(AmpOrder));
    if ~isempty(SpeakerName)
        AmpTDT{f} = getTDTAmp(F, FreqOrder, AmpOrder, SpeakerName);
    else
        AmpTDT{f} = ones(size(AmpOrder));
    end
end
