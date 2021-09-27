function [F,ClickOrder,AmpOrder,PhaseOrder] = ClickAPorder(Amp_Array,Ntrials,randomized)

Amp_Array = sort(Amp_Array,'descend');

F = 1;
nF = length(F);
nA = length(Amp_Array);
ClickOrder = repmat(F(:)',Ntrials*nA,1);
ClickOrder = ClickOrder(:);
AmpOrder = repmat(Amp_Array(:)',Ntrials,1);
AmpOrder = repmat(AmpOrder(:),1,nF);
% AmpOrder = repmat(Amp_Array(:),Ntrials,nF);
PhaseOrder = ones(size(AmpOrder));
PhaseOrder(2:2:end,:) = -1;

AmpOrder = AmpOrder(:);
PhaseOrder = PhaseOrder(:);
if randomized
    p1 = randperm(Ntrials*nF*nA);
    ClickOrder = ClickOrder(p1);
    AmpOrder = AmpOrder(p1);
    PhaseOrder = PhaseOrder(p1);
    p2 = randperm(Ntrials*nF*nA);
    ClickOrder = ClickOrder(p2);
    AmpOrder = AmpOrder(p2);
    PhaseOrder = PhaseOrder(p2);
    p3 = randperm(Ntrials*nF*nA);
    ClickOrder = ClickOrder(p3);
    AmpOrder = AmpOrder(p3);
    PhaseOrder = PhaseOrder(p3);
end