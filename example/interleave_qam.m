function rx_out = interleave_qam(encdata,punc_bit,short_bit,qam,snr,sigma)
tran_bit=encdata;N=length(encdata);
index=1:N;
rx_out=zeros(1,N);
tran_bit([punc_bit,short_bit])=[];
index([punc_bit,short_bit])=[];
m=log2(qam);
if mod(length(tran_bit),m)~=0
    tran_bit_z=[tran_bit,zeros(1,m-mod(length(tran_bit),m))];
else
    tran_bit_z=tran_bit;
end
label=m/2;n=length(tran_bit_z);
tran_M=reshape(tran_bit_z,n/label,label)';
tran_bit_z=reshape(tran_M,1,n);

modulatedSignal = qammod(tran_bit_z', qam,'InputType','bit','UnitAveragePower',1);
noisySignal = awgn(modulatedSignal, snr, 'measured');
demodsig = qamdemod(noisySignal, qam,'OutputType','approxllr','UnitAveragePower',1);
rx=demodsig'/2/sigma;
rx_M=reshape(rx,label,n/label)';
rx=reshape(rx_M,1,n);
if mod(length(tran_bit),m)~=0
    rx_out(index)=rx(1:end-m+mod(length(tran_bit),m));
else
    rx_out(index)=rx;
end
rx_out(short_bit)=10000;
end

