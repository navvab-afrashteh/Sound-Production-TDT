function [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,Nnum,Nden)

idx1 = mod(idx2,Nden)+1;
if idx1 == 1
    n1 = n2+1;
else
    n1 = n2;
end
n2 = n2+ceil(Nnum/Nden);
idx2 = Nnum - (Nden-idx1+1 + (n2-n1-1)*Nden);
if idx2<=0
    n2 = n2-1;
    idx2 = Nnum - (Nden-idx1+1 + (n2-n1-1)*Nden);
end

% p=Nden-idx1+1 + (n2-n1-1)*Nden + idx2;