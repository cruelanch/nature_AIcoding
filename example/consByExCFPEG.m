function [ bsHCell,tnrGCell ] = consByExCFPEG(PG,typZ,orgH_DSbs,seed,lenOfCndd,maxCnsdCyc)
%consByExCFPEG construction by ExCFPEG algorithm
%   Detailed explanation goes here

% init by computer
H_ExpMat_DSbs=ExCFPEG_getExpMatFrmPG(PG);% get Exponent matrix from a Protograph.
[M_DSbs,N_DSbs]=size(H_ExpMat_DSbs);% number of check nodes and variable nodes in base matrix H
rng(seed);
% get vnDgr_bs,cnDgr_bs
vnDgr_bs=sum(PG);
cnDgr_bs=sum(PG,2);
% get typical multiple
len_typZ=length(typZ);
typMult=typZ(1);
typMult(2:len_typZ)=typZ(2:len_typZ)./typZ(1:len_typZ-1);% typical multiple of typical Z
if(any(fix(typMult)~=typMult))% jude typMult is a integer vector
    error('typMult has non integer elements');
end
% init H_bsCell and tnrGCell
bsHCell=cell(1,len_typZ);% base matrix Cell
tnrGCell=cell(1,len_typZ);
for i=1:len_typZ
    bsHCell{i}=orgH_DSbs;
end

% construct bigger matrix by exponent matrix
matSqn=1;% matrix sequence
while(matSqn<=len_typZ)
    curZ=typZ(matSqn);
    M=M_DSbs*curZ;
    H_DSbs=bsHCell{matSqn};
    H_bs=getMQCFromDSMQC(H_DSbs,curZ);
    tnrG=getTGfromMQCM(H_bs,curZ);
    % get dpthOfTE, maxTau and maxD
    dpthOfTE=2*curZ*M_DSbs;% depth of tree expending when count minimum e-cycle set. its max value equal to 2*Z_table*M_bs
    % init remain edge number matirx
    rmnEdgeNumInPG=PG;% remain edge number in PG
    for i=1:M_DSbs
        for j=1:N_DSbs
            len_Hbs=length(H_DSbs{i,j});
            if(len_Hbs==1)
                if(H_DSbs{i,j}~=-inf)
                    rmnEdgeNumInPG(i,j)=rmnEdgeNumInPG(i,j)-1;
                    if(rmnEdgeNumInPG(i,j)<0)
                        error('rmnEdgeNumInPG(i,j)<0');
                    end
                end
            else
                rmnEdgeNumInPG(i,j)=rmnEdgeNumInPG(i,j)-len_Hbs;
                if(rmnEdgeNumInPG(i,j)<0)
                    error('rmnEdgeNumInPG(i,j)<0');
                end
            end
        end
    end
    % other init work
    ccnDgr_bs=tnrG.rowDist(1:curZ:end);% column weight distribution
    cvnDgr_bs=tnrG.colDist(1:curZ:end);% current VN degree base matrix
    % begin construction
    i=1;
    while(i<=N_DSbs)
        j=cvnDgr_bs(i)+1;
        while(j<=vnDgr_bs(i))
            %clc;
            % fprintf('current Z = %d\n\n',curZ);
            % fprintf('current variable nodes i = %d\n\n',i);
            % fprintf('current edge j = %d\n\n',j);
            if(j==1)
                % construct first edge
                [rs_e1,shift_e1]=ExCFPEG_firEdge(rmnEdgeNumInPG,i,H_ExpMat_DSbs,curZ,typZ);% row sequence of edge 1, shift parameter of edge 1
                H_DSbs_tmp=zeros(M_DSbs,N_DSbs)-1;
                H_DSbs_tmp(rs_e1,i)=shift_e1;
                tnrG=addCPMtoTnrG_v2(H_DSbs_tmp,tnrG);% add circular permutation matrix to bitwise matrix
                if(H_DSbs{rs_e1,i}(end)==-inf)
                    H_DSbs{rs_e1,i}=shift_e1;
                else
                    H_DSbs{rs_e1,i}=union(H_DSbs{rs_e1,i},shift_e1);
                end
                rmnEdgeNumInPG(rs_e1,i)=rmnEdgeNumInPG(rs_e1,i)-1;
                ccnDgr_bs(rs_e1)=ccnDgr_bs(rs_e1)+1;
            else
                % construct other edge
                CN_W=cntMCEVforQCEMPEGinT_v3(tnrG,(i-1)*curZ+1,vnDgr_bs(i),dpthOfTE);
                [isRchOfEC_bs,CN_W_tmp]=PGEMPEG_isRchOfEC(CN_W, tnrG, rmnEdgeNumInPG,i);% is reached of each check node in base matrix
                if(all(isRchOfEC_bs))
                    % estimation
                    [slctCyc,crspEMD,crspCA,crspRN,crspSP]=...
                        ExCFPEG_estOfEC(CN_W,tnrG,H_DSbs,i,lenOfCndd,rmnEdgeNumInPG,typZ,H_ExpMat_DSbs);% estimation of each check node
                    if(all(slctCyc(1:lenOfCndd)>maxCnsdCyc))
                        slctCyc2=slctCyc;
                        crspCA2=crspCA;% correspond cycle amount
                        crspEMD2=crspEMD;
                        crspRN2=crspRN;
                        crspSP2=crspSP;
                    else
                        % additional enumeration
                        [slctCyc2,crspEMD2,crspCA2,crspRN2,crspSP2]=...
                            PGCFPEG_addEnu(crspRN,crspSP,lenOfCndd,dpthOfTE,i,tnrG,vnDgr_bs);
                    end
                    % generate edge
                    k=1;
                    H_DSbs_tmp=zeros(M_DSbs,N_DSbs)-1;
                    H_DSbs_tmp(crspRN2(k),i)=crspSP2(k);
                    tnrG=addCPMtoTnrG_v2(H_DSbs_tmp,tnrG);
                    if(H_DSbs{crspRN2(k),i}(end)==-inf)
                        H_DSbs{crspRN2(k),i}=crspSP2(k);
                    else
                        H_DSbs{crspRN2(k),i}=union(H_DSbs{crspRN2(k),i},crspSP2(k));
                    end
                    rmnEdgeNumInPG(crspRN2(k),i)=rmnEdgeNumInPG(crspRN2(k),i)-1;
                    ccnDgr_bs(crspRN2(k))=ccnDgr_bs(crspRN2(k))+1;
                else
                    % generate unreached edge
                    [rs_unRchE,shift_unRchE]=ExCFPEG_unRchE(isRchOfEC_bs,H_ExpMat_DSbs,i,curZ,typZ);
                    H_DSbs_tmp=zeros(M_DSbs,N_DSbs)-1;
                    H_DSbs_tmp(rs_unRchE,i)=shift_unRchE;
                    tnrG=addCPMtoTnrG_v2(H_DSbs_tmp,tnrG);% add circular permutation matrix to bitwise matrix
                    if(H_DSbs{rs_unRchE,i}(end)==-inf)
                        H_DSbs{rs_unRchE,i}=shift_unRchE;
                    else
                        H_DSbs{rs_unRchE,i}=union(H_DSbs{rs_unRchE,i},shift_unRchE);
                    end
                    rmnEdgeNumInPG(rs_unRchE,i)=rmnEdgeNumInPG(rs_unRchE,i)-1;
                    ccnDgr_bs(rs_unRchE)=ccnDgr_bs(rs_unRchE)+1;
                end
            end
            % control while(j<=vn_dgr_bsH(i))
            j=j+1;
            cvnDgr_bs(i)=cvnDgr_bs(i)+1;
        end
        % control while(i<=N_bsH)
        i=i+1;
    end
    % control the while loop
    H_ExpMat_DSbs=H_DSbs;
    bsHCell{matSqn}=H_DSbs;
    tnrGCell{matSqn}=tnrG;
    matSqn=matSqn+1;
end


end

function [ H_DSbs ] = ExCFPEG_getExpMatFrmPG( PG )
%EXEMPEG3_GETEXPMATFRMPG get exponent matrix from protograph
%   PG: protograph with single or multiple edge
% 
%   H_bs: base matrix with QC structure


[M_DSbs,N_DSbs]=size(PG);
% H_bs=num2cell(PG-1);
H_DSbs=num2cell(-inf*ones(M_DSbs,N_DSbs));
for i=1:M_DSbs
    for j=1:N_DSbs
        if(PG(i,j)==1)
            H_DSbs{i,j}=0;
        elseif(PG(i,j)>=2)
            H_DSbs{i,j}=zeros(1,PG(i,j));
        end
    end
end

end

function [ H_bs ] = getMQCFromDSMQC( H_DSbs, Z )
%GETMQCFROMDSMQCMAT get MQC from double shift MQC matrix
%   H_DSbs: double shift base matrix
% 
%   H_bs: base matrix
% 
%   Z: lifting factor

[M_bs,N_bs]=size(H_DSbs);
H_bs=cell(M_bs,N_bs);
for i=1:M_bs
    for j=1:N_bs
        if(H_DSbs{i,j}(1)==-inf)
            H_bs{i,j}=-1;
        else
            H_bs{i,j}=mod(H_DSbs{i,j},Z);
            H_bs(i,j)=shortOfDSMQCMat(H_bs(i,j));
        end
        
    end
end
end

function [ H_out_bs ] = shortOfDSMQCMat( H_bs )
%SHORTOFDSMQCMAT Summary of this function goes here
%   Detailed explanation goes here

[M_bs,N_bs]=size(H_bs);
H_out_bs=cell(M_bs,N_bs);
for i=1:M_bs
    for j=1:N_bs
        uni_H_bs=unique(H_bs{i,j});
        histc_H_bs=histc(H_bs{i,j},uni_H_bs);
        histc_H_bs_mod2=mod(histc_H_bs,2);
        H_out_bs{i,j}=uni_H_bs(logical(histc_H_bs_mod2));
        len_H_out_bs=length(H_out_bs{i,j});
        if(len_H_out_bs~=1)
            H_out_bs{i,j}(H_out_bs{i,j}==-inf)=[];
        end
    end
end
end

function [ tnrG ] = getTGfromMQCM( H_MQC,Z )
%GETTGFROMQCM get Tanner Graph from Multiple-QC-matrix
%   H_QC is multiple QC matrix, including shift parameter, and is right
%   shift, and is a cell.
% 
%   Z is expanded multiple
% 
%   tnrG is a tanner Graph, containing M(scalar), N(scalar), Z(scalar),
%   VN(struction), CN(struction), colDist(vector), rowDist(vector)

tnrG.Z=Z;
[M_QC,N_QC]=size(H_MQC);
M=M_QC*Z;
N=N_QC*Z;
tnrG.M=M;
tnrG.N=N;
VN(N).nbr=[];
CN(M).nbr=[];
colDist_QC=zeros(1,N_QC);
rowDist_QC=zeros(1,M_QC);
for i=1:N_QC
    num_QCnbr_VN=0;
    QCnbr_VN=zeros(1,M_QC);
    sp_CN=zeros(1,M_QC);
    for j=1:M_QC
        for k=1:length(H_MQC{j,i})
            if(H_MQC{j,i}(k)~=-1)
                num_QCnbr_VN=num_QCnbr_VN+1;
                QCnbr_VN(num_QCnbr_VN)=j;
                sp_CN(num_QCnbr_VN)=H_MQC{j,i}(k);
                colDist_QC(i)=colDist_QC(i)+1;
            end
        end
    end
%     QCnbr_VN=find(H_MQC(:,i)~=-1);% neighbor of variable nodes
%     [~,~,sp_CN]=find(H_MQC(:,i)+1);% shift parameter of each check node
%     sp_CN=sp_CN-1;
    for j=1:num_QCnbr_VN
        position_CN_tmp=mod(-sp_CN(j)+Z,Z);
        position_CN_tmp2=(QCnbr_VN(j)-1)*Z+1;
        position_CN=position_CN_tmp2+position_CN_tmp;
        VN((i-1)*Z+1).nbr=[VN((i-1)*Z+1).nbr,position_CN];
        for k=2:Z
            position_CN_tmp=mod(position_CN_tmp+1,Z);
            position_CN=position_CN_tmp2+position_CN_tmp;
            VN((i-1)*Z+k).nbr=[VN((i-1)*Z+k).nbr,position_CN];
        end
    end
end
for i=1:M_QC
    num_QCnbr_CN=0;
    QCnbr_CN=zeros(1,N_QC);
    sp_VN=zeros(1,N_QC);
    for j=1:N_QC
        for k=1:length(H_MQC{i,j})
            if(H_MQC{i,j}(k)~=-1)
                num_QCnbr_CN=num_QCnbr_CN+1;
                QCnbr_CN(num_QCnbr_CN)=j;
                sp_VN(num_QCnbr_CN)=H_MQC{i,j}(k);
                rowDist_QC(i)=rowDist_QC(i)+1;
            end
        end
    end
%     QCnbr_CN=find(H_MQC(i,:)~=-1);% neighbor of variable nodes
%     [~,~,sp_VN]=find(H_MQC(i,:)+1);% shift parameter of each check node
%     sp_VN=sp_VN-1;
    for j=1:num_QCnbr_CN
        position_VN_tmp=mod(sp_VN(j)+Z,Z);
        position_VN_tmp2=(QCnbr_CN(j)-1)*Z+1;
        position_CN=position_VN_tmp2+position_VN_tmp;
        CN((i-1)*Z+1).nbr=[CN((i-1)*Z+1).nbr,position_CN];
        for k=2:Z
            position_VN_tmp=mod(position_VN_tmp+1,Z);
            position_CN=position_VN_tmp2+position_VN_tmp;
            CN((i-1)*Z+k).nbr=[CN((i-1)*Z+k).nbr,position_CN];
        end
    end
end
tnrG.VN=VN;
tnrG.CN=CN;
% colDist_QC=sum(H_MQC~=-1);% column weight distribution
for i=1:N_QC
    colDist((i-1)*Z+1:i*Z)=colDist_QC(i)*ones(1,Z);
end
% rowDist_QC=sum(H_MQC~=-1,2);% column weight distribution
for i=1:M_QC
    rowDist((i-1)*Z+1:i*Z)=rowDist_QC(i)*ones(1,Z);
end
tnrG.colDist=colDist;
tnrG.rowDist=rowDist;
end

function [ rs_e1,shift_e1 ] = ExCFPEG_firEdge( rmnEdgeNumInPG,i,H_ExpMat_bs,curZ,typZ)
%ExCFPEG_firEdge exponent matrix expanding by CFPEG. This is the first edge step
% modified form ExEMPEG3_firEdge to add that random choose shift_e1 
%   rmnEdgeNumInPG(matrix):remain edge number in Protograph
% 
%   rs_e1: row sequence of edge 1
% 
%   shift_e1: shift parameter of edge 1
% 
%   i: current variable node sequence
% 
%   H_ExpMat_bs: the base matrix of exponent matrix
% 
%   curZ: current Z
% 
%   typZ: typical Z

% get rs_e1
availPstn=find(rmnEdgeNumInPG(:,i)~=0);% available position
ind_rand_e1=round(rand*(length(availPstn)-1))+1;
rs_e1=availPstn(ind_rand_e1);% row sequence number of edge 1
% get shift_e1
ind_typZ=find(typZ==curZ);
if(ind_typZ==1)
    multOftypZ=typZ(ind_typZ);
    prvOftypZ=1;
else
    multOftypZ=typZ(ind_typZ)/typZ(ind_typZ-1);% multiple
    prvOftypZ=typZ(ind_typZ-1);% previous of typical Z
end
avlValVec=[];
for i1=0:multOftypZ-1
    avlValVec=union(avlValVec,H_ExpMat_bs{rs_e1,i}+i1*prvOftypZ);
end
ind_rand_avlValVec=round(rand*(length(avlValVec)-1))+1;
shift_e1=avlValVec(ind_rand_avlValVec);

end

function [ tnrG ] = addCPMtoTnrG_v2( H_QC_add, tnrG )
%ADDCPMTOTNRG add circular permutation matrix to tanner graph
% new feature in version 2
% 1.tnrG.VN((c_Hqc(i)-1)*Z+1).nbr=union(tnrG.VN((c_Hqc(i)-1)*Z+1).nbr,p_CN);
% 2.tnrG.VN((c_Hqc(i)-1)*Z+k).nbr=union(tnrG.VN((c_Hqc(i)-1)*Z+k).nbr,p_CN);
% 3.tnrG.CN((r_Hqc(i)-1)*Z+1).nbr=union(tnrG.CN((r_Hqc(i)-1)*Z+1).nbr,p_VN);
% 4.tnrG.CN((r_Hqc(i)-1)*Z+k).nbr=union(tnrG.CN((r_Hqc(i)-1)*Z+k).nbr,p_VN);
% 
% 
% 
%   H_QC_add is a QC matrix, including shift parameter. it has M_bs row and
%   N_bs column.
% 
%   tnrG is a tanner graph, including Z,M,N,VN,CN,colDist,rowDist

[r_Hqc,c_Hqc,v_Hqc]=find(H_QC_add+1);% row, column, value
v_Hqc=v_Hqc-1;
Z=tnrG.Z;
for i=1:length(r_Hqc)
    % modify the neighbor of variable node in Tanner Graph
    p_CN_tmp=(r_Hqc(i)-1)*Z+1;% position of check node
    sp_CN=mod(-v_Hqc(i)+Z,Z);% shift paramter of position_CN_tmp
    p_CN=p_CN_tmp+sp_CN;
    tf_p=ismember(p_CN,tnrG.VN((c_Hqc(i)-1)*Z+1).nbr);% tf means trueFalse, loc means location
    if(tf_p)
        continue
    else
        tnrG.VN((c_Hqc(i)-1)*Z+1).nbr=union(tnrG.VN((c_Hqc(i)-1)*Z+1).nbr,p_CN);
        tnrG.colDist((c_Hqc(i)-1)*Z+1)=tnrG.colDist((c_Hqc(i)-1)*Z+1)+1;
    end
    for k=2:Z
        sp_CN=mod(sp_CN+1,Z);
        p_CN=p_CN_tmp+sp_CN;
        tnrG.VN((c_Hqc(i)-1)*Z+k).nbr=union(tnrG.VN((c_Hqc(i)-1)*Z+k).nbr,p_CN);
        tnrG.colDist((c_Hqc(i)-1)*Z+k)=tnrG.colDist((c_Hqc(i)-1)*Z+k)+1;
    end
    % modify the neighbor of check node in Tanner Graph
    p_VN_tmp=(c_Hqc(i)-1)*Z+1;% position of variable node
    sp_VN=mod(v_Hqc(i)+Z,Z);% shift paramter of position_VN_tmp
    p_VN=p_VN_tmp+sp_VN;
    tnrG.CN((r_Hqc(i)-1)*Z+1).nbr=union(tnrG.CN((r_Hqc(i)-1)*Z+1).nbr,p_VN);
    tnrG.rowDist((r_Hqc(i)-1)*Z+1)=tnrG.rowDist((r_Hqc(i)-1)*Z+1)+1;
    for k=2:Z
        sp_VN=mod(sp_VN+1,Z);
        p_VN=p_VN_tmp+sp_VN;
        tnrG.CN((r_Hqc(i)-1)*Z+k).nbr=union(tnrG.CN((r_Hqc(i)-1)*Z+k).nbr,p_VN);
        tnrG.rowDist((r_Hqc(i)-1)*Z+k)=tnrG.rowDist((r_Hqc(i)-1)*Z+k)+1;
    end
end
end

function CN_W = cntMCEVforQCEMPEGinT_v3( tnrG,i_tG,dv_i, maxD )
% count e-minimum cycle of each variable at the the position i_QC of H_QC for
% QCEMPEG construction in Tanner Graph. writen by Wu XiaoNing 2015/04/19.
% 
% Reference:
% refer to the appendix A of the paper 'Constructing LDPC Codes by Error 
% Minimization Progressive Edge Growth', but not the same as it.there are
% many difference from it. Also refer to 'Selective avoidance of cycles in
% irregular LDPC code construction' for tree spanning.
% 
% difference form version 2
% 1.colDstr(i_tG)=dv_i;
% 2.if(Emd_new<VN(T_v_aprx(i)).Emd&&~(Emd_new==1&&length(S)<colDstr(T_v_aprx(i))))
% 
% tnrG is a Tanner Graph. It contain M,N,Z,VN,CN
% 
% i_tG is the i-th variable nodes in Tanner Graph
% 
% dv_i is the i-th variable node degree
% 
% CN_W is a structure including W with the same mean of the algorithm of
% appendix A. cn.W=[d,/tau,m;d1,/tau1,m1;...]. W is a matrix.
% 
% maxD is max tree expanding depth 

% maxD=10;
% init the variale
M=tnrG.M;
N=tnrG.N;
colDstr=tnrG.colDist;% column weight distribution
colDstr(i_tG)=dv_i;
D=0;% the depth of the tree spanned from variable_N, start from 0
T_v=i_tG;% tiers of active variable at depth D,i_ptl = i_partial
VN=tnrG.VN;
CN=tnrG.CN;
VN(N).Emd=inf;
CN(M).Emd=inf;
CN_W(M).W=zeros(1,3);
CN_W(M).W_ind=1;% record the position of row which should be writen in the next step in the W. every row of W means a set of triplets
for i=1:M
    CN(i).Emd=inf;
    CN_W(i).W=zeros(1,3);
    CN_W(i).W_ind=1;
end
for i=1:N
    VN(i).Emd=inf;
end
VN(i_tG).Mult=1;% multiplicity of variable node N
VN(i_tG).prnt=[];% parent of variable node N
VN(i_tG).Emd=dv_i;% the reason why "+1" is that the "+1" means the added one edge

% count the cycle formed by i_ptl variable in the H_ptl
while(D<maxD&&~isempty(T_v))
    D=D+1;
    % expand tier D of check nodes:
    T_c=[];% tiers of active check at depth D
    for i=1:length(T_v)% how many variable nodes in T_v
        T_c_tmp=setdiff(VN(T_v(i)).nbr,VN(T_v(i)).prnt);
        for j=1:length(T_c_tmp)% how many check nodes in T_c_tmp
            if(VN(T_v(i)).Emd<CN(T_c_tmp(j)).Emd)
                CN(T_c_tmp(j)).Emd=VN(T_v(i)).Emd;
                CN(T_c_tmp(j)).prnt=T_v(i);
                if(CN_W(T_c_tmp(j)).W_ind==1||D~=CN_W(T_c_tmp(j)).W(CN_W(T_c_tmp(j)).W_ind-1,1))
                    %add new (d,/tau,m) to W
                    CN_W(T_c_tmp(j)).W(CN_W(T_c_tmp(j)).W_ind,:)=[D,CN(T_c_tmp(j)).Emd,VN(T_v(i)).Mult];
                    CN_W(T_c_tmp(j)).W_ind=CN_W(T_c_tmp(j)).W_ind+1;
                else
                    %cover current value
                    CN_W(T_c_tmp(j)).W(CN_W(T_c_tmp(j)).W_ind-1,:)=[D,CN(T_c_tmp(j)).Emd,VN(T_v(i)).Mult];
                end
                T_c=union(T_c,T_c_tmp(j));%
            elseif((VN(T_v(i)).Emd==CN(T_c_tmp(j)).Emd) && any(T_c==T_c_tmp(j)))
                CN(T_c_tmp(j)).prnt=union(CN(T_c_tmp(j)).prnt,T_v(i));
                %no need to judge, just add 'm' to origin 
                CN_W(T_c_tmp(j)).W(CN_W(T_c_tmp(j)).W_ind-1,:)=CN_W(T_c_tmp(j)).W(CN_W(T_c_tmp(j)).W_ind-1,:)+[0,0,VN(T_v(i)).Mult];
            end
        end
    end
    D=D+1;
    T_v=[];% 
    T_v_aprx=[];
    for i=1:N
        VN(i).preEmd=VN(i).Emd;% previous Emd means Emd^{D-2}
        VN(i).preMult=VN(i).Mult;
        VN(i).prnt=[];
        VN(i).Mult=0;
    end
    for j=1:length(T_c)
        T_v_tmp=setdiff(CN(T_c(j)).nbr,CN(T_c(j)).prnt);
        for i=1:length(T_v_tmp)
            VN(T_v_tmp(i)).prnt=union(VN(T_v_tmp(i)).prnt,T_c(j));
            T_v_aprx=union(T_v_aprx,T_v_tmp(i));
        end
    end
    for i=1:length(T_v_aprx)
        P=VN(T_v_aprx(i)).prnt;
        prnt_P=[];
        for j=1:length(P);
            prnt_P=union(CN(P(j)).prnt,prnt_P);
        end
        for j=1:length(prnt_P)
            S=intersect(P,VN(prnt_P(j)).nbr);
            Emd_new=VN(prnt_P(j)).preEmd+colDstr(T_v_aprx(i))-2*length(S);
            if(Emd_new<VN(T_v_aprx(i)).Emd&&~(Emd_new==1&&length(S)<colDstr(T_v_aprx(i))))% avoid tau==1
                VN(T_v_aprx(i)).Emd=Emd_new;
                VN(T_v_aprx(i)).Mult=VN(prnt_P(j)).preMult;
                VN(T_v_aprx(i)).prnt=S;
                T_v=union(T_v,T_v_aprx(i));
            elseif(Emd_new==VN(T_v_aprx(i)).Emd && any(T_v==T_v_aprx(i)))
                VN(T_v_aprx(i)).Mult=VN(T_v_aprx(i)).Mult+VN(prnt_P(j)).preMult;
                VN(T_v_aprx(i)).prnt=union(VN(T_v_aprx(i)).prnt,S);
            end
        end
    end
end
end

function [ isRchOfEC_bs,CN_W_tmp ] = PGEMPEG_isRchOfEC( CN_W, tnrG, rmnEdgeNumInPG,i )
%QCEMPEG_ISRCHOFEC Whether each Check Node is reached by tree expanding
%   CN_W(structure):including W with the same mean of the algorithm of
%   appendix A. cn.W=[d,/tau,m;d1,/tau1,m1;...]. W is a matrix.
% 
%   isRchOfEC(0/1): Whether each Check Node is reached
% 
%   rmnEdgeNumInPG: remain edge number in Protograph
% 
%   i: current variable node sequence.

M=tnrG.M;
Z=tnrG.Z;
M_bs=M/Z;
for k=1:M
    CN_W_tmp(1,(k-1)*3+1:(k-1)*3+3)=CN_W(k).W(1,1:3);
end
isRchOfEC_bs=zeros(1,M_bs);% is zeros of each check node
% isZofEC_bsH=zeros(1,M_bs);
for k=1:M_bs
    k1=(k-1)*Z+1:(k-1)*Z+1+Z-1;
    if(any(CN_W_tmp(3*(k1-1)+1)))
        %exist non-zeros
        isRchOfEC_bs(k)=1;
    end
end
isRchOfEC_bs(rmnEdgeNumInPG(:,i)==0)=1;
end

function [ rs_unRchE,shift_unRchE ] = ExCFPEG_unRchE( isRchOfEC_bs,H_ExpMat_bs,i,curZ,typZ )
%ExEMPEG3_unRchE Exponent empeg algorithm 3
% modified from ExEMPEG_unRchE to get random shift value
% 1.
% 
%   isRchOfEC_bs: is reached of each check node
% 
%   rs_unRchE: row sequence of unReached edge
% 
%   shift_unRchE: shift parameter of unreached edge
% 
%   H_ExpMat_bs: the base matrix of exponent matrix
% 
%   i: i-th column

% get rs_unRchE
fndRch=find(isRchOfEC_bs==0);
ind_rand=round(rand*(length(fndRch)-1))+1;
rs_unRchE=fndRch(ind_rand);
% get shift_unRchE
ind_typZ=find(typZ==curZ);
if(ind_typZ==1)
    multOftypZ=typZ(ind_typZ);
    prvOftypZ=1;
else
    multOftypZ=typZ(ind_typZ)/typZ(ind_typZ-1);% multiple
    prvOftypZ=typZ(ind_typZ-1);% previous of typical Z
end
avlValVec=[];
for i1=0:multOftypZ-1
    avlValVec=union(avlValVec,H_ExpMat_bs{rs_unRchE,i}+i1*prvOftypZ);
end
ind_rand_avlValVec=round(rand*(length(avlValVec)-1))+1;
shift_unRchE=avlValVec(ind_rand_avlValVec);


end

function [ slctCyc,crspEMD,crspCA,crspRN,crspSP ] = ...
    ExCFPEG_estOfEC( CN_W,tnrG,H_DSbs,i,lenOfCndd,rmnEdgeNumInPG,typZ,H_ExpMat_DSbs )
%ExCFPEG_estOfEC exponent matrix cfpeg. estimation of each check node
% modified by PGCFPEG_estOfEC_V2 to support ExCFPEG
% 
%   CN_W:CN_W is a structure including W with the same mean of the algorithm of
%        appendix A. cn.W=[d,/tau,m;d1,/tau1,m1;...]. W is a matrix.
% 
%   lenOfLBP: length of lower bound of P_B^LB(j,p+1)
% 
%   uB:upper bound
% 
%   LBofP_B:Lower bound of P_B
% 
%   minLBofP_j:row sequence of minimum Lower bound of P
% 
%   minLBofP_p:shift parameter of minimum Lower bound of P
% 
%   minLBofP:minimum Lower bound of P
% 
%   tnrG:tnrG is a Tanner Graph. It contain M,N,Z,VN,CN,colDist,rowDist
% 
%   H_bs:base matrix
% 
%   i: current variable sequence.
% 
%   ccnDgr_bs: current check node degree in base matrix
% 
%   maxCND: maximum check node degree
% 
%   Z_table:Z of the table
% 
%   rmnEdgeNumInPG: remain edge number in protograph
% 
%   maxNumOfsd: max number of s-d. start from 0. prefer to 4.
M=tnrG.M;
Z=tnrG.Z;
M_bs=M/Z;
minCycOfEachSP=zeros(M_bs,Z);% minimum cycle of each shift parameter
minCycAmntOfESP=zeros(M_bs,Z);% minimum cycle amount of each shift parameter
EMDCrspMCOfESP=zeros(M_bs,Z);% EMD correspondding to minimum cycle of each shift parameter
slctCyc=zeros(1,lenOfCndd);% the lenOfLBP minimum value of Lower bound of P_B
crspEMD=zeros(1,lenOfCndd);
crspCA=inf*ones(1,lenOfCndd);
crspRN=zeros(1,lenOfCndd)-1;% according to equation(14)
crspSP=zeros(1,lenOfCndd)-1;% according to equation(14)
ind_typZ=find(typZ==Z);
if(ind_typZ==1)
    multOftypZ=typZ(ind_typZ);
    prvOftypZ=1;
else
    multOftypZ=typZ(ind_typZ)/typZ(ind_typZ-1);% multiple
    prvOftypZ=typZ(ind_typZ-1);% previous of typical Z
end
for k=1:M_bs
    if(any(H_DSbs{k,i}~=-inf))
        if(rmnEdgeNumInPG(k,i)==0)
            continue;
        end
        % find available Exponent shift parameter(avlExpSP)
        len_H_DSbs_ki=length(H_DSbs{k,i});
        avlExpSP=H_ExpMat_DSbs{k,i};
        for j=1:len_H_DSbs_ki
            for j1=0:multOftypZ-1
                ExpSP_tmp=H_DSbs{k,i}(j)-j1*prvOftypZ;
                if(ismember(ExpSP_tmp,H_ExpMat_DSbs{k,i}))
                    ind_avlExpSP=find(avlExpSP==ExpSP_tmp,1,'first');
                    avlExpSP(ind_avlExpSP)=[];
                    break;
                end
            end
        end
        % find avlValVec
        avlValVec=[];
        for i1=0:multOftypZ-1
            avlValVec=union(avlValVec,avlExpSP+i1*prvOftypZ);
        end
        for p=0:Z-1 %each shift parameter
            if(any(H_DSbs{k,i}==p))
                continue;
            elseif(~ismember(p,avlValVec))
                continue;
            elseif(any(H_DSbs{k,i}-p==Z/2)||any(H_DSbs{k,i}-p==-Z/2))% avoid cycle 4 in one position with multiple shift parameter
                minCycOfEachSP(k,p+1)=2;
                minCycAmntOfESP(k,p+1)=Z/2;% maximum
                EMDCrspMCOfESP(k,p+1)=0;
            else
                minCycOfEachSP(k,p+1)=inf;
                for q=0:Z-1
                    k1=(k-1)*Z+q+1;
                    i1=1;
                    if(CN_W(k1).W(i1,1)==0)
                        continue;
                    elseif(CN_W(k1).W(i1,1)==1)
                       continue;
                    else
                       d=(CN_W(k1).W(i1,1)+1)/2;
                       tau=CN_W(k1).W(i1,2)-2;
                       d1=d*Z/gcd(q+p,Z);%correct
                       tau1=tau*Z/gcd(q+p,Z);%correct
                       if(d1<minCycOfEachSP(k,p+1))
                           minCycOfEachSP(k,p+1)=d1;
                           minCycAmntOfESP(k,p+1)=CN_W(k1).W(i1,3);
                           EMDCrspMCOfESP(k,p+1)=tau1;
                       end
                    end
                end
            end
            
            % sort
            for i2=1:lenOfCndd
                if(minCycOfEachSP(k,p+1)<slctCyc(i2))
                    continue;
                elseif(minCycOfEachSP(k,p+1)>slctCyc(i2))
                    slctCyc(i2+1:lenOfCndd)=slctCyc(i2:lenOfCndd-1);
                    slctCyc(i2)=minCycOfEachSP(k,p+1);
                    crspEMD(i2+1:lenOfCndd)=crspEMD(i2:lenOfCndd-1);
                    crspEMD(i2)=EMDCrspMCOfESP(k,p+1);
                    crspCA(i2+1:lenOfCndd)=crspCA(i2:lenOfCndd-1);
                    crspCA(i2)=minCycAmntOfESP(k,p+1);
                    crspRN(i2+1:lenOfCndd)=crspRN(i2:lenOfCndd-1);
                    crspRN(i2)=k;
                    crspSP(i2+1:lenOfCndd)=crspSP(i2:lenOfCndd-1);
                    crspSP(i2)=p;
                    break;
                else
                    if(EMDCrspMCOfESP(k,p+1)<crspEMD(i2))
                        continue;
                    elseif(EMDCrspMCOfESP(k,p+1)>crspEMD(i2))
                        slctCyc(i2+1:lenOfCndd)=slctCyc(i2:lenOfCndd-1);
                        slctCyc(i2)=minCycOfEachSP(k,p+1);
                        crspEMD(i2+1:lenOfCndd)=crspEMD(i2:lenOfCndd-1);
                        crspEMD(i2)=EMDCrspMCOfESP(k,p+1);
                        crspCA(i2+1:lenOfCndd)=crspCA(i2:lenOfCndd-1);
                        crspCA(i2)=minCycAmntOfESP(k,p+1);
                        crspRN(i2+1:lenOfCndd)=crspRN(i2:lenOfCndd-1);
                        crspRN(i2)=k;
                        crspSP(i2+1:lenOfCndd)=crspSP(i2:lenOfCndd-1);
                        crspSP(i2)=p;
                        break;
                    else
                        if(minCycAmntOfESP(k,p+1)>=crspCA(i2))
                            continue;
                        elseif(minCycAmntOfESP(k,p+1)<crspCA(i2))
                            slctCyc(i2+1:lenOfCndd)=slctCyc(i2:lenOfCndd-1);
                            slctCyc(i2)=minCycOfEachSP(k,p+1);
                            crspEMD(i2+1:lenOfCndd)=crspEMD(i2:lenOfCndd-1);
                            crspEMD(i2)=EMDCrspMCOfESP(k,p+1);
                            crspCA(i2+1:lenOfCndd)=crspCA(i2:lenOfCndd-1);
                            crspCA(i2)=minCycAmntOfESP(k,p+1);
                            crspRN(i2+1:lenOfCndd)=crspRN(i2:lenOfCndd-1);
                            crspRN(i2)=k;
                            crspSP(i2+1:lenOfCndd)=crspSP(i2:lenOfCndd-1);
                            crspSP(i2)=p;
                            break;
                        end
                    end
                end
            end
        end
    else
        if(rmnEdgeNumInPG(k,i)==0)
            continue;
        end
        % find avlValVec
        avlValVec=[];
        for i1=0:multOftypZ-1
            avlValVec=union(avlValVec,H_ExpMat_DSbs{k,i}+i1*prvOftypZ);
        end
        for p=0:Z-1
            if(~ismember(p,avlValVec))
                continue;
            end
            minCycOfEachSP(k,p+1)=inf;
            for q=0:Z-1
                k1=(k-1)*Z+q+1;
                i1=1;
                if(CN_W(k1).W(i1,1)==0)
                   continue;
                else
                   d=(CN_W(k1).W(i1,1)+1)/2;
                   tau=CN_W(k1).W(i1,2)-2;
                   d1=d*Z/gcd(q+p,Z);%correct
                   tau1=tau*Z/gcd(q+p,Z);%correct
                   if(d1<minCycOfEachSP(k,p+1))
                       minCycOfEachSP(k,p+1)=d1;
                       minCycAmntOfESP(k,p+1)=CN_W(k1).W(i1,3);
                       EMDCrspMCOfESP(k,p+1)=tau1;
                   end
                end
            end
            % sort LBofP_B
            for i2=1:lenOfCndd
                if(minCycOfEachSP(k,p+1)<slctCyc(i2))
                    continue;
                elseif(minCycOfEachSP(k,p+1)>slctCyc(i2))
                    slctCyc(i2+1:lenOfCndd)=slctCyc(i2:lenOfCndd-1);
                    slctCyc(i2)=minCycOfEachSP(k,p+1);
                    crspEMD(i2+1:lenOfCndd)=crspEMD(i2:lenOfCndd-1);
                    crspEMD(i2)=EMDCrspMCOfESP(k,p+1);
                    crspCA(i2+1:lenOfCndd)=crspCA(i2:lenOfCndd-1);
                    crspCA(i2)=minCycAmntOfESP(k,p+1);
                    crspRN(i2+1:lenOfCndd)=crspRN(i2:lenOfCndd-1);
                    crspRN(i2)=k;
                    crspSP(i2+1:lenOfCndd)=crspSP(i2:lenOfCndd-1);
                    crspSP(i2)=p;
                    break;
                else
                    if(EMDCrspMCOfESP(k,p+1)<crspEMD(i2))
                        continue;
                    elseif(EMDCrspMCOfESP(k,p+1)>crspEMD(i2))
                        slctCyc(i2+1:lenOfCndd)=slctCyc(i2:lenOfCndd-1);
                        slctCyc(i2)=minCycOfEachSP(k,p+1);
                        crspEMD(i2+1:lenOfCndd)=crspEMD(i2:lenOfCndd-1);
                        crspEMD(i2)=EMDCrspMCOfESP(k,p+1);
                        crspCA(i2+1:lenOfCndd)=crspCA(i2:lenOfCndd-1);
                        crspCA(i2)=minCycAmntOfESP(k,p+1);
                        crspRN(i2+1:lenOfCndd)=crspRN(i2:lenOfCndd-1);
                        crspRN(i2)=k;
                        crspSP(i2+1:lenOfCndd)=crspSP(i2:lenOfCndd-1);
                        crspSP(i2)=p;
                        break;
                    else
                        if(minCycAmntOfESP(k,p+1)>=crspCA(i2))
                            continue;
                        elseif(minCycAmntOfESP(k,p+1)<crspCA(i2))
                            slctCyc(i2+1:lenOfCndd)=slctCyc(i2:lenOfCndd-1);
                            slctCyc(i2)=minCycOfEachSP(k,p+1);
                            crspEMD(i2+1:lenOfCndd)=crspEMD(i2:lenOfCndd-1);
                            crspEMD(i2)=EMDCrspMCOfESP(k,p+1);
                            crspCA(i2+1:lenOfCndd)=crspCA(i2:lenOfCndd-1);
                            crspCA(i2)=minCycAmntOfESP(k,p+1);
                            crspRN(i2+1:lenOfCndd)=crspRN(i2:lenOfCndd-1);
                            crspRN(i2)=k;
                            crspSP(i2+1:lenOfCndd)=crspSP(i2:lenOfCndd-1);
                            crspSP(i2)=p;
                            break;
                        end
                    end
                end
            end
        end
    end
end

end

function [slctCyc2,crspEMD2,crspCA2,crspRN2,crspSP2] = ...
    PGCFPEG_addEnu( crspRN,crspSP,lenOfCndd,dpthOfTE,i,tnrG,vnDgr_bs )
%PGEMPEGCF_addEnu addtional enumeration in PGEMPEG cycle first algorithm
% 
%   minLBofP_j;row sequence of minimum Lower bound of P. 
% 
%   minLBofP_p:shift parameter of minimum Lower bound of P
% 
%   minLBofP:minimum Lower bound of P
% 
%   uB:upper bound
% 
%   lenOfLBP:length of lower bound of P_B^LB(j,p+1)

M_bs=tnrG.M/tnrG.Z;
N_bs=tnrG.N/tnrG.Z;
Z=tnrG.Z;
slctCyc2=zeros(1,lenOfCndd);% the lenOfLBP minimum value of Lower bound of P_B
crspEMD2=zeros(1,lenOfCndd);
crspCA2=inf*ones(1,lenOfCndd);
crspRN2=zeros(1,lenOfCndd)-1;% according to equation(14)
crspSP2=zeros(1,lenOfCndd)-1;% according to equation(14)
for k=1:lenOfCndd% the last one isn't the smallest lenOfLBP, so just to lenOfLBP-1
    if(crspRN(k)==-1)
        continue;
    end
    H_base_tmp=zeros(M_bs,N_bs)-1;
    H_base_tmp(crspRN(k),i)=crspSP(k);
    tnrG_tmp=addCPMtoTnrG_v2(H_base_tmp,tnrG);% versiton 2 new
    i1=1;
    i2=(i-1)*Z+i1;%the position of variable nodes in H
    j2=(crspRN(k)-1)*Z+1+mod(Z-crspSP(k)+i1-1,Z);%the position of check nodes in H
    % version 2 new
    e_dlt.pOfCN=j2;%edge_delete
    e_dlt.pOfVN=i2;%edge_delete
    tnrG_tmp2=delOneEfromTG(tnrG_tmp,e_dlt);
    CN_W=cntMCEVforQCEMPEGinT_v3(tnrG_tmp2,i2,vnDgr_bs(i),dpthOfTE);
    i3=1;
    if(CN_W(j2).W(i3,1)==0)
        minCyc=inf;
        minCycAmnt=1;
        EMDCrspMC=inf;
    else
        d=(CN_W(j2).W(i3,1)+1)/2;
        tau=CN_W(j2).W(i3,2)-2;
        minCyc=d;
        minCycAmnt=CN_W(j2).W(i3,3);
        EMDCrspMC=tau;
    end
    
    % fprintf('Cycle %d, EMD %d, CycleAmount %d',minCyc,EMDCrspMC,minCycAmnt);
    % fprintf('\n')
    
    for i2=1:lenOfCndd
        if(minCyc<slctCyc2(i2))
            continue;
        elseif(minCyc>slctCyc2(i2))
            slctCyc2(i2+1:lenOfCndd)=slctCyc2(i2:lenOfCndd-1);
            slctCyc2(i2)=minCyc;
            crspEMD2(i2+1:lenOfCndd)=crspEMD2(i2:lenOfCndd-1);
            crspEMD2(i2)=EMDCrspMC;
            crspCA2(i2+1:lenOfCndd)=crspCA2(i2:lenOfCndd-1);
            crspCA2(i2)=minCycAmnt;
            crspRN2(i2+1:lenOfCndd)=crspRN2(i2:lenOfCndd-1);
            crspRN2(i2)=crspRN(k);
            crspSP2(i2+1:lenOfCndd)=crspSP2(i2:lenOfCndd-1);
            crspSP2(i2)=crspSP(k);
            break;
        else
            if(EMDCrspMC<crspEMD2(i2))
                continue;
            elseif(EMDCrspMC>crspEMD2(i2))
                slctCyc2(i2+1:lenOfCndd)=slctCyc2(i2:lenOfCndd-1);
                slctCyc2(i2)=minCyc;
                crspEMD2(i2+1:lenOfCndd)=crspEMD2(i2:lenOfCndd-1);
                crspEMD2(i2)=EMDCrspMC;
                crspCA2(i2+1:lenOfCndd)=crspCA2(i2:lenOfCndd-1);
                crspCA2(i2)=minCycAmnt;
                crspRN2(i2+1:lenOfCndd)=crspRN2(i2:lenOfCndd-1);
                crspRN2(i2)=crspRN(k);
                crspSP2(i2+1:lenOfCndd)=crspSP2(i2:lenOfCndd-1);
                crspSP2(i2)=crspSP(k);
                break;
            else
                if(minCycAmnt>=crspCA2(i2))
                    continue;
                elseif(minCycAmnt<crspCA2(i2))
                    slctCyc2(i2+1:lenOfCndd)=slctCyc2(i2:lenOfCndd-1);
                    slctCyc2(i2)=minCyc;
                    crspEMD2(i2+1:lenOfCndd)=crspEMD2(i2:lenOfCndd-1);
                    crspEMD2(i2)=EMDCrspMC;
                    crspCA2(i2+1:lenOfCndd)=crspCA2(i2:lenOfCndd-1);
                    crspCA2(i2)=minCycAmnt;
                    crspRN2(i2+1:lenOfCndd)=crspRN2(i2:lenOfCndd-1);
                    crspRN2(i2)=crspRN(k);
                    crspSP2(i2+1:lenOfCndd)=crspSP2(i2:lenOfCndd-1);
                    crspSP2(i2)=crspSP(k);
                    break;
                end
            end
        end
    end
end


end

function [ tnrG ] = delOneEfromTG( tnrG,edge )
%DELONEEFROMTG delete one edge from tanner graph
%   tnrG is a tanner graph, including Z,M,N,VN,CN,colDist,rowDist
% 
%   edge is an edge in Tanner graph. it is a structure, including
%   pOfVN,pOfCN,

if(ismember(edge.pOfCN,tnrG.VN(edge.pOfVN).nbr)&&ismember(edge.pOfVN,tnrG.CN(edge.pOfCN).nbr))
    tnrG.VN(edge.pOfVN).nbr=setdiff(tnrG.VN(edge.pOfVN).nbr,edge.pOfCN);
    tnrG.CN(edge.pOfCN).nbr=setdiff(tnrG.CN(edge.pOfCN).nbr,edge.pOfVN);
    tnrG.colDist(edge.pOfVN)=tnrG.colDist(edge.pOfVN)-1;
    tnrG.rowDist(edge.pOfCN)=tnrG.rowDist(edge.pOfCN)-1;
else
    warning('matlab:delOneEfromTG','no edge is deleted');
end



end


































