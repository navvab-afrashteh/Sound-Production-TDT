function Stim = applyGate(Stim,Gate)
% apply gate sequence to onset and offset of the stimulus
if size(Stim,1) == 1
    onGate = Gate(:).';
    offGate = fliplr(onGate);
elseif size(Stim,2) == 1
    onGate = Gate(:);
    offGate = flipud(onGate);
end
NGate = length(Gate);
Stim(1:NGate) = Stim(1:NGate) .* onGate;
Stim(end-NGate+1:end) = Stim(end-NGate+1:end) .* offGate;