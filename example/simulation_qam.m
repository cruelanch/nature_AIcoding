clc;clear;
warning('off', 'all')
qam=4;max_iteration=4;

load("H_SV_out.mat");
Z=64;row=13;
punc_bit_index=[1,2];
[numsum,~]=size(H_SV_out);
EsN0_list=2.8:0.4:4.8;
H_SV_RES=cell(numsum,6);
for i=1:4
    for j=1:numsum
        H_SV_RES{j,i}=H_SV_out{j,i};
    end
end

for num_idx=1:numsum
    rate=22/(22+row-length(punc_bit_index));
    H_p =H_SV_out{num_idx,1};
    H=sparse(logical(HQC_Zcn(H_p,Z)));
    encoder = ldpcEncoderConfig(H);decoder = ldpcDecoderConfig(H,'layered-bp');
    %%
    punc_bit=[];
    for p=1:length(punc_bit_index)
        punc_bit=[punc_bit,(punc_bit_index(p)-1)*Z+1:punc_bit_index(p)*Z];
    end
    short_bit=[];
    [M,N]=size(H);
    %BP_init
    [Hcol,Hrow]=find(H');
    rowweight=sum(H,2);
    rowfind=zeros(M,2);
    rowfind(1,1)=1;
    rowfind(1,2)=rowweight(1);
    for i=2:M
        rowfind(i,1)=rowfind(i-1,2)+1;
        rowfind(i,2)=rowfind(i-1,2)+rowweight(i);
    end
    f=length(Hrow);
    colfind=cell(1,N);
    for j=1:N
        colfind{j}=find(Hcol==j);
    end
    numcode=100000;linenum=1000;
    fer_cc=zeros(1,length(EsN0_list));itr_cc=zeros(1,length(EsN0_list));
    disp('SNR   num     fel     fer     ber     itr');
    SNR_res=[];fer_res=[];
    for EsN0idx=1:length(EsN0_list)
        sigma=(10.^(-EsN0_list(EsN0idx)/10))/2;
        fer1=0;ber1=0;itr1=0;textlength=0;
        for line=1:numcode
            parfor index=1:linenum
                data=randi([0 1],1,22*Z);
                encdata=ldpcEncode(data',encoder)';
                rx=interleave_qam(encdata,punc_bit,short_bit,qam,EsN0_list(EsN0idx),sigma);
                %[dec,num_iter,~]=ldpcDecode(rx',decoder,max_iteration);dec=dec';
                [dec,num_iter]=dec_LNMS(rx,H,max_iteration,Hrow,Hcol,rowfind,f,encdata,0.75);
                ber=mean(xor(dec(1:N-M),data));itr1=itr1+num_iter;ber1=ber1+ber;
                if ber>0
                    fer1=fer1+1;
                end
            end
            fprintf(repmat('\b', 1, textlength));
            dispcontext=[num2str(EsN0_list(EsN0idx)),'     ',num2str(line*linenum),'   ',num2str(fer1),'    ',...
                num2str(fer1/(line*linenum)),'   ',num2str(ber1/(line*linenum)),'    ',num2str(itr1/(line*linenum)),'                          '];
            textlength=length(dispcontext);
            fprintf(dispcontext);
            if fer1>=50
                fprintf('\n');
                fer_cc(EsN0idx)=fer1/(line*linenum);itr_cc(EsN0idx)=itr1/(line*linenum);
                break;
            end
        end
    end
    H_SV_RES{num_idx,5}=fer_cc;
    H_SV_RES{num_idx,6}=itr_cc;
    disp([num2str(num_idx/numsum*100),'%']);
end
save('H_SV_RES');