function [Lch] = Lch_qam(PM,N,Z,tran_bit)
Lch=zeros(1,N);

snrM=(erfcinv(2*PM)).^2;
tran_bit_num=length(tran_bit)*Z;
snr_tran_bit=[];
for i=1:length(snrM)
    snr_tran_bit=[snr_tran_bit,snrM(i)*ones(1,ceil(tran_bit_num/length(snrM)))];
end
for i=1:length(tran_bit)
    index_i=tran_bit(i);
    index_min=(i-1)*Z+1;
    index_max=i*Z;
    Lch(index_i)=mean(snr_tran_bit(index_min:index_max));
end
end

