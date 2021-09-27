function [F, Fspec, FreqOrder, AmpOrder, startDur, Ntrials, NTones] = FI_FATotder(Amp, min_Freq, max_Freq, baseDur, ISI, sweepDur, Mode, Rec_Dur)
% 
% Amp = 80;
% min_Freq = 2;
% max_Freq = 64;
% ISI = 0.250;
% sweepDur = 8;
% Mode = 'backwards Shift';

NTones = floor(sweepDur/ISI);
fmin_oc = log2(min_Freq);
fmax_oc = log2(max_Freq);
freqs_oc = linspace(fmin_oc,fmax_oc,NTones);
FreqOrder = 2.^freqs_oc;
F = FreqOrder;
fpoints = 1/(freqs_oc(2)-freqs_oc(1));
Fspec = [min_Freq, max_Freq,fpoints];

half_idx = round((length(FreqOrder)+1)/2)-1;
switch Mode
    case 'forwards'
        FreqOrder = FreqOrder;
    case 'backwards'
        FreqOrder = fliplr(FreqOrder);
    case 'forwards Shift'
        FreqOrder = [FreqOrder(half_idx:end) FreqOrder(1:half_idx-1)];
    case 'backwards Shift'
        FreqOrder = fliplr(FreqOrder);
        FreqOrder = [FreqOrder(half_idx:end) FreqOrder(1:half_idx-1)];
end

Ntrials = floor((Rec_Dur-baseDur)/(sweepDur+baseDur));
AmpOrder = Amp*ones(1,NTones);

startDur = Rec_Dur - baseDur - Ntrials*(sweepDur+baseDur);
startDur = startDur/2;
% TimeOrder = ISI*ones(1,NTones);
% TimeOrder(1) = TimeOrder(1) + baseDur;
% TimeOrder = repmat(TimeOrder,1,Ntrials);
% TimeOrder(1) = TimeOrder(1) + first_end_dur;
