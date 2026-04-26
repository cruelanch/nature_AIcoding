function Iapp = layered_PEXIT_qam(EbN0,rate,M,N,tran_bit,colfind,rowfind,HB,graph,itrnum,qam,Z)
SNR=EbN0+10*log10(rate)+10*log10(log2(qam));
s=10^(SNR/10);
PM=PMerr(s,qam);
Lch=Lch_qam(PM,N,Z,tran_bit);

Ich=zeros(1,N);
for idx=1:N
    Ich(idx)=func_J(2*sqrt((Lch(idx)*2)));
end

Iev=zeros(M,N);Iec=zeros(M,N);Iav=zeros(M,N);Iac=zeros(M,N);
Iapp=Ich;
%disp(num2str(Ich))
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
                % num_tmp=-arcfunc_J(Iav(i,j))^2+arcfunc_J(Iapp(j))^2;
                % if num_tmp<0
                %     num_tmp=0;
                % end
                % Iev(i,j)=func_J(sqrt(num_tmp));
            end
        end
        Iac=Iev;
        % disp(num2str(i));
        % disp(num2str(Iac));
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
        Iav=Iec;
        % disp(num2str(i));
        % disp(num2str(Iav));
        % for j=1:N
        %     if graph(i,j)==1
        %         Iapp(j)=func_J(sqrt(arcfunc_J(Iav(i,j))^2+arcfunc_J(Iev(i,j))^2));
        %     end
        % end
        for j=1:N
            temp=0;
            for index=1:length(colfind{j})
                indexi=colfind{j}(index);
                temp=temp+HB(indexi,j)*(arcfunc_J(Iav(indexi,j)))^2;
            end
            Iapp(j)=func_J(sqrt(temp+(arcfunc_J(Ich(j)))^2));
        end
    end
    %disp(num2str(Iapp));
    if min(Iapp(1:N-M))>1-1e-10
        break;
    end
end
end