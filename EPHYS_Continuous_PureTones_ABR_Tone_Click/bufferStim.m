function StimTDT = bufferStim(FreqIndex,Amp_TDT,PhaseOrder,tone,n1,n2,idx1,idx2)

NTones = length(FreqIndex);
f = FreqIndex(n1); a = Amp_TDT(n1); p = PhaseOrder(n1);
firstStim = a*tone(idx1:end,f)*p;

if n1==NTones
    StimTDT = firstStim;
    return;
end
if n2-1>NTones
    k=0;
    for n = n1+1:NTones
        k=k+1;
        f(k) = FreqIndex(n);
        a(k) = Amp_TDT(n);
        p(k) = PhaseOrder(n);
    end
    midStim = a.*tone(:,f).*p;
    StimTDT = [firstStim;midStim(:)];
    return;
end
k=0;
for n = n1+1:n2-1
    k=k+1;
    f(k) = FreqIndex(n);
    a(k) = Amp_TDT(n);
    p(k) = PhaseOrder(n);
end
midStim = a.*tone(:,f).*p;
f = FreqIndex(n2); a = Amp_TDT(n2); p = PhaseOrder(n2);
lastStim = a*tone(1:idx2,f)*p;

StimTDT = [firstStim;midStim(:);lastStim];
