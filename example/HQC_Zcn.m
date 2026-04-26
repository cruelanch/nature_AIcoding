function [H] = HQC_Zcn(HB,n)
[M,N]=size(HB);
Zc_base=eye(n);
Zc_zero=zeros(n,n);
H=[];
Hr=[];
for i=1:M
    for j=1:N
        sgn=HB(i,j);
        if sgn==-1
            Hr=[Hr Zc_zero];
        else
            Zc_tp=circshift(Zc_base,[0,sgn]);
            Hr=[Hr Zc_tp];
        end
    end
    H=[H
        Hr];
    Hr=[];
end
end

