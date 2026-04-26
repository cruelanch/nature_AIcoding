function [PM] = PMerr(s,M)
m=log2(M);
PM=zeros(1,m/2);
for l=1:m/2
    temp=0;
    for k=0:(1-2^(-l))*sqrt(M)-1
        temp=temp+(-1)^(floor((k*2^(l-1))/sqrt(M)))*(2^(l-1)-floor((k*2^(l-1))/sqrt(M)+0.5))*erfc((2*k+1)*sqrt(3*s/(2*M-2)));
    end
    PM(l)=1/sqrt(M)*temp;
end
end

