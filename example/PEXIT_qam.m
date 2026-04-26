function Iapp = PEXIT_qam(EbN0,rate,M,N,tran_bit,colfind,rowfind,HB,graph,itrnum,qam,Z)
SNR=EbN0+10*log10(rate)+10*log10(log2(qam));
s=10^(SNR/10);
PM=PMerr(s,qam);
Lch=Lch_qam(PM,N,Z,tran_bit);

Ich=zeros(1,N);
for idx=1:N
    Ich(idx)=func_J(2*sqrt((Lch(idx)*2)));
end

Iev=zeros(M,N);Iec=zeros(M,N);Iav=zeros(M,N);Iac=zeros(M,N);
Iapp=zeros(1,N);
for k=1:itrnum
    for i=1:M
        for j=1:N
            if graph(i,j)==1
                temp=0;
                for index=1:length(colfind{j})
                    indexi=colfind{j}(index);
                    if indexi~=i
                        temp=temp+HB(indexi,j)*(arcfunc_J(Iav(indexi,j)))^2;
                    end
                end
                Iev(i,j)=func_J(sqrt(temp+(HB(i,j)-1)*(arcfunc_J(Iav(i,j)))^2+arcfunc_J(Ich(j))^2));
            end
        end
    end
    Iac=Iev;
    for i=1:M
        for j=1:N
            if graph(i,j)==1
                temp=0;
                for index=1:length(rowfind{i})
                    indexj=rowfind{i}(index);
                    if indexj~=j
                        temp=temp+HB(i,indexj)*(arcfunc_J(1-Iac(i,indexj)))^2;
                    end
                end
                Iec(i,j)=1-func_J(sqrt(temp+(HB(i,j)-1)*(arcfunc_J(1-Iac(i,j)))^2));
            end
        end
    end
    Iav=Iec;
    for j=1:N
        temp=0;
        for index=1:length(colfind{j})
            indexi=colfind{j}(index);
            temp=temp+HB(indexi,j)*(arcfunc_J(Iav(indexi,j)))^2;
        end
        Iapp(j)=func_J(sqrt(temp+(arcfunc_J(Ich(j)))^2));
    end
end
end

