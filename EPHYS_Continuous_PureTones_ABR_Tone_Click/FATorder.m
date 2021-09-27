function [F, FreqOrder, AmpOrder, PhaseOrder] = FATorder(Fspec, Amp_Array, Ntrials, randomized)

% Fspec : vector of frequencies specifications: [first freq, last freq, frequency per octave],in kHz
% Amp_Array : vector of amplitudes in dB
% ramp : ramp time for the start and the end of sine wave, in milli second
% dur : duration of tones, in milli second
% ISI : time between two consecutive tones, in milli second
% jitter :  jitter for silence period, in mili second
% Ntrials : number of all complex sounds, integer number
% randomized : indicating if the signal be generated in a randomized manner
% or not. 1 or 0
% Edited by Navvab Afrashteh, May 22, 2017

Amp_Array = sort(Amp_Array,'descend');

F = octavefreq(Fspec);
F = [8,16,24,32,40];
nF = length(F);
nA = length(Amp_Array);
FreqOrder = repmat(F(:)',Ntrials*nA,1);
FreqOrder = FreqOrder(:);
AmpOrder = repmat(Amp_Array(:)',Ntrials,1);
AmpOrder = repmat(AmpOrder(:),1,nF);
% AmpOrder = repmat(Amp_Array(:),Ntrials,nF);
PhaseOrder = ones(size(AmpOrder));
PhaseOrder(2:2:end,:) = -1;

AmpOrder = AmpOrder(:);
PhaseOrder = PhaseOrder(:);
if randomized
    p1 = randperm(Ntrials*nF*nA);
    FreqOrder = FreqOrder(p1);
    AmpOrder = AmpOrder(p1);
    PhaseOrder = PhaseOrder(p1);
    p2 = randperm(Ntrials*nF*nA);
    FreqOrder = FreqOrder(p2);
    AmpOrder = AmpOrder(p2);
    PhaseOrder = PhaseOrder(p2);
    p3 = randperm(Ntrials*nF*nA);
    FreqOrder = FreqOrder(p3);
    AmpOrder = AmpOrder(p3);
    PhaseOrder = PhaseOrder(p3);
end
