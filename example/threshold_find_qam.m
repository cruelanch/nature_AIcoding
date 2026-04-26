function [EbN0] = threshold_find_qam(rate,M,N,tran_bit,colfind,rowfind,HB,graph,itrnum,qam,Z,layered_flag)
EbN0_max=20;EbN0_min=0;
while (EbN0_max-EbN0_min)>0.001
    EbN0=(EbN0_max+EbN0_min)/2;
    if layered_flag==1
        Iapp = layered_PEXIT_qam(EbN0,rate,M,N,sort(tran_bit),colfind,rowfind,HB,graph,itrnum,qam,Z);
    else
        Iapp = PEXIT_qam(EbN0,rate,M,N,sort(tran_bit),colfind,rowfind,HB,graph,itrnum,qam,Z);
    end
    if min(Iapp(1:N-M))>1-1e-10
        EbN0_max=EbN0;
    else
        EbN0_min=EbN0;
    end
end
end

