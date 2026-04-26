clear;
load('H_SV_cell.mat');
itrnum=4; qam=4;Z=64; layered_flag=1;row=13;

[num,~]=size(H_SV_cell);
cnt=0;
for i=1:num
    if ~isempty(H_SV_cell{i,1})
        cnt=cnt+1;
    end
end
H_SV_out=cell(cnt+1,4);
cnt=1;
H_SV_out{cnt,1}=genBG1(row,Z);
H_SV_out{cnt,2}=0;
[CBN_r,VBN_r]=size(H_SV_out{cnt,1});
HB = zeros(CBN_r, VBN_r);
HB(H_SV_out{cnt,1}>= 0) = 1;
H_SV_out{cnt,4}=sum(sum(double(HB)));
graph = HB;
rowfind = cell(1, CBN_r);
colfind = cell(1, VBN_r);
for ii = 1:CBN_r
    rowfind{ii} = find(graph(ii, :) == 1);
end
for jj = 1:VBN_r
    colfind{jj} = find(graph(:, jj) == 1);
end
rate = 22/(20 + CBN_r);
tran_bit = 3:(22 + CBN_r);
H_SV_out{cnt,3} = threshold_find_qam(rate, CBN_r, VBN_r, tran_bit, colfind, rowfind, HB, graph, itrnum, qam, Z, layered_flag);
for i=1:num
    if ~isempty(H_SV_cell{i,1})
        cnt=cnt+1;
        H_SV_out{cnt,1}=H_SV_cell{i,1};
        H_SV_out{cnt,2}=H_SV_cell{i,2};
        [CBN_r,VBN_r]=size(H_SV_cell{i,1});
        HB = zeros(CBN_r, VBN_r);
        HB(H_SV_cell{i,1}>= 0) = 1;
        H_SV_out{cnt,4}=sum(sum(double(HB)));
        graph = HB;
        rowfind = cell(1, CBN_r);
        colfind = cell(1, VBN_r);
        for ii = 1:CBN_r
            rowfind{ii} = find(graph(ii, :) == 1);
        end
        for jj = 1:VBN_r
            colfind{jj} = find(graph(:, jj) == 1);
        end
        rate = 22/(20 + CBN_r);
        tran_bit = 3:(22 + CBN_r);
        H_SV_out{cnt,3} = threshold_find_qam(rate, CBN_r, VBN_r, tran_bit, colfind, rowfind, HB, graph, itrnum, qam, Z, layered_flag);
    end
end
seed=7;
lenOfCndd=10;
maxCnsdCyc=32;
parfor i=2:cnt
    H_SVt=H_SV_out{i,1};
    PG=double(H_SVt>=0);
    orgH_DSbs=H_SVt;orgH_DSbs(orgH_DSbs==-1)=-inf;
    orgH_DSbs_tmp=orgH_DSbs(1:4,1:22);
    orgH_DSbs_tmp(orgH_DSbs_tmp>=0)=-inf;
    orgH_DSbs(1:4,1:22)=orgH_DSbs_tmp;
    orgH_DSbs_tmp=orgH_DSbs(5:end,1:26);
    orgH_DSbs_tmp(orgH_DSbs_tmp>=0)=-inf;
    orgH_DSbs(5:end,1:26)=orgH_DSbs_tmp;
    orgH_DSbs=num2cell(orgH_DSbs);
    [bsHCell,~]=consByExCFPEG(PG,[Z],orgH_DSbs,seed,lenOfCndd,maxCnsdCyc);
    H_SVt=cell2mat(bsHCell{1});
    H_SVt(H_SVt==-inf)=-1;
    H_SV_out{i,1}=H_SVt;
end
save('H_SV_out');