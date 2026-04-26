function [cj,num_iter]=dec_LNMS(rx,H,max_iteration,Hrow,Hcol,rowfind,f,encdata,alpha)
[rows,cols]=size(H);
Lji=zeros(f,3);
Lji(:,1)=Hrow;
Lji(:,2)=Hcol;
Lij=Lji;
cj=zeros(1,cols);
const=1;
cst=100000;
Lj=rx;
for k=1:max_iteration 
    for i=1:rows
        for m=rowfind(i,1):rowfind(i,2)
            Lji(m,3)=Lj(Lij(m,2))-Lij(m,3);
        end

        % minnum = mink(abs(Lji(rowfind(i,1):rowfind(i,2),3)),2);
        % min1 = minnum(1); min2 = minnum(2);
        % const=1;
        % for m=rowfind(i,1):rowfind(i,2)
        %     const=const*sign(Lji(m,3));
        % end

        for m=rowfind(i,1):rowfind(i,2)
            j=Lij(m,2);
            for n=rowfind(i,1):rowfind(i,2)
                if Lij(n,2)~=j
                    const=const*sign(Lji(n,3));
                    cst=min(abs(Lji(n,3)),cst);
                end
            end
            Lij(m,3)=const*cst*alpha;
            % if abs(Lji(m,3))==min1
            %     Lij(m,3) = min2*const*sign(Lji(m,3))*alpha;
            % else
            %     Lij(m,3) = min1*const*sign(Lji(m,3))*alpha;
            % end
            const=1;
            cst=100000;
            Lj(j)=Lji(m,3)+Lij(m,3);
        end
    end
    cj=Lj<0;
    if isequal(cj(1:cols-rows),encdata(1:cols-rows))
        break;
    end
end
num_iter=k;%µü´ú´ÎÊý
end