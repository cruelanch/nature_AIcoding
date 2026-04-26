clear; clc; close all;
%diary('GAoutput_multi_region2.txt');  % 日志文件（记录所有区域优化过程）

%% ========================== 1. 基础参数配置（用户可按需修改）==========================
itrnum = 4;          % 迭代相关基础参数（原逻辑保留）
qam = 4;             % QAM调制阶数
layered_flag = 1;    % 分层标志
row = 13;            % 基础行数
Z = 384;             % 基础列数
rowlist = [4, 13];
global_assess=1;

% 正交约束参数
orthrow = 0;         % 正交行块大小
orth = 0;           % 正交约束权重
relufactor = 1;      % 性能惩罚系数
alpha = 0.01;

% 遗传算法核心参数（每个区域共用此参数）
Population = 100;    % 种群数量
G = 100;             % 单个区域的迭代次数（每个区域都迭代G次）

%% ========================== 2. 自定义n个独立优化区域（核心修改：用户在此定义区域）==========================
regions = [
    struct('startRow',1, 'endRow',4, 'startCol',1, 'endCol',22)
    struct('startRow',5, 'endRow',13, 'startCol',1, 'endCol',26)
];
n_regions = length(regions);  % 区域总数

H_SV_cell=cell(G*n_regions,2);

%% ========================== 3. 初始化基础矩阵 & 区域合法性检查==========================
% 生成初始块矩阵（原逻辑保留）
BL_tem = genBG1(row, Z);

[CBN, VBN] = size(BL_tem);  % 获取矩阵维度（行：CBN，列：VBN）
fprintf('初始矩阵维度：%d行 × %d列\n', CBN, VBN);

% 全局最优矩阵初始化（第一个区域以原始矩阵为起点，后续区域以上一个区域的最优结果为起点）
global_BL_opt = BL_tem;  
global_min_gamma = inf;  % 全局最小gamma值（所有区域优化完成后的最终性能）

%% ========================== 4. 串行优化每个区域（核心逻辑：逐个区域迭代G次）==========================
for region_idx = 1:n_regions
    fprintf('==================================== 开始优化第%d个区域 ====================================\n', region_idx);
    current_region = regions(region_idx);  % 当前待优化区域
    fprintf('当前区域范围：行[%d-%d]，列[%d-%d]\n', ...
        current_region.startRow, current_region.endRow, current_region.startCol, current_region.endCol);
    
    % -------------------------- 4.1 初始化当前区域的种群（以上一个区域的最优矩阵为起点）--------------------------
    t = cell(3*Population, 1);  % 种群存储（3倍种群：原始+交叉+变异）
    fits = zeros(3*Population, 1);  % 适应度值
    fal = zeros(3*Population, 1);   % 正交约束惩罚值
    nump = zeros(3*Population, 1);
    fal_val = 100;                  % 初始惩罚值
    nump_val = 100;

    rowlist_tmp = rowlist; 
    if global_assess==0 || current_region.startRow == 1
        rowlist_tmp(rowlist_tmp>current_region.endRow)=[];
    end
    rowlist_tmp(rowlist_tmp<current_region.endRow)=[];

    EbN0_re = zeros(1, length(rowlist_tmp));
    for ridx = 1:length(rowlist_tmp)
        % 计算当前区域优化前的基准性能（EbN0_re）
        row = rowlist_tmp(ridx);
        VBN_r = 22+row; CBN_r = row;
        HB_base = zeros(CBN_r, VBN_r); 
        HB_base(global_BL_opt(1:CBN_r,1:VBN_r) >= 0) = 1; 
        graph_base = HB_base;
        % 预计算行/列非零索引（加速后续性能评估）
        rowfind_base = cell(1, CBN_r); 
        colfind_base = cell(1, VBN_r);
        for ii = 1:CBN_r
            rowfind_base{ii} = find(graph_base(ii, :) == 1);
        end
        for jj = 1:VBN_r
            colfind_base{jj} = find(graph_base(:, jj) == 1);
        end
        rate = 22/(20 + row);
        tran_bit = 3:(22 + row);
        EbN0_re(ridx) = threshold_find_qam(rate, CBN_r, VBN_r, tran_bit, colfind_base, rowfind_base, HB_base, graph_base, itrnum, qam, Z, layered_flag);
    end
    gamma_re = sum(EbN0_re);
    fprintf('第%d个区域优化前的基准gamma值：%.4f\n', region_idx, gamma_re);
    
    % 初始化种群：第1个个体为上一区域的最优矩阵，其余个体仅在当前区域内随机调整
    for i = 1:Population
        if i == 1
            t{i} = global_BL_opt;  % 继承上一区域的最优结果
        else
            temp_BL = global_BL_opt;  % 以全局最优矩阵为基础
            
            % 仅在当前区域内随机交换（-1和非-1）
            L = randi([2, 10]);  % 交换次数（原逻辑保留）
            for l = 1:L
                % 循环直到找到符合条件的两个点（不同位置、类型相反）
                % while 1
                %     m1 = randi([current_region.startRow, current_region.endRow]);
                %     n1 = randi([current_region.startCol, current_region.endCol]);
                %     m2 = randi([current_region.startRow, current_region.endRow]);
                %     n2 = randi([current_region.startCol, current_region.endCol]);
                %     % 条件：位置不同 + 一个是-1、一个是非-1
                %     if (m1 ~= m2 || n1 ~= n2) && ...
                %        ((temp_BL(m1,n1) == -1 && temp_BL(m2,n2) >= 0) || ...
                %         (temp_BL(m1,n1) >= 0 && temp_BL(m2,n2) == -1))
                %         break;
                %     end
                % end
                % % 交换两点值
                % temp_Z = temp_BL(m1, n1);
                % temp_BL(m1, n1) = temp_BL(m2, n2);
                % temp_BL(m2, n2) = temp_Z;
                m1 = randi([current_region.startRow, current_region.endRow]);
                n1 = randi([current_region.startCol, current_region.endCol]);
                if temp_BL(m1, n1)>=0
                    temp_BL(m1, n1)=-1;
                else
                    temp_BL(m1, n1)=randi([0,Z-1],1,1);
                end
            end
            t{i} = temp_BL;
        end
    end
    
    % -------------------------- 4.2 当前区域的遗传算法迭代（G次）--------------------------
    current_min_gamma = inf;  % 当前区域的最小gamma值
    current_min_fits = inf;   % 当前区域的最小适应度值
    current_BL_opt = global_BL_opt;  % 当前区域的最优矩阵（初始为上一区域结果）
    
    for g = 1:G
        % -------------------------- 4.2.1 并行计算种群适应度（原逻辑保留，仅评估当前矩阵）--------------------------
        parfor i = 1:length(t)
            if isempty(t{i})
                fits(i) = 1000;  % 空矩阵赋大适应度（惩罚）
            else
                EbN0 = zeros(1, length(rowlist_tmp));
                for ridx = 1:length(rowlist_tmp)
                    row = rowlist_tmp(ridx);
                    % 1. 构建HB矩阵（-1→0，非-1→1）
                    CBN_r=row;VBN_r=22+row;
                    HB = zeros(CBN_r, VBN_r);
                    HB(t{i}(1:CBN_r,1:VBN_r)>= 0) = 1;
                    graph = HB;
                    
                    % 2. 计算行/列非零索引
                    rowfind = cell(1, CBN_r);
                    colfind = cell(1, VBN_r);
                    for ii = 1:CBN_r
                        rowfind{ii} = find(graph(ii, :) == 1);
                    end
                    for jj = 1:VBN_r
                        colfind{jj} = find(graph(:, jj) == 1);
                    end
                    
                    % 3. 计算EbN0（性能指标）
                    rate = 22/(20 + row);
                    tran_bit = 3:(22 + row);
                    EbN0(ridx) = threshold_find_qam(rate, CBN_r, VBN_r, tran_bit, colfind, rowfind, HB, graph, itrnum, qam, Z, layered_flag);
                end
                nump(i) = sum(sum(double(HB)));

                % 4. 计算正交约束惩罚值（fal）
                fal(i) = 0;
                if current_region.startRow > 4 && orthrow>1
                    % 分块计算正交约束（原逻辑保留）
                    for ridx_block = 1:floor((current_region.endRow-current_region.startRow+1)/orthrow)
                        row_start = current_region.startRow+orthrow * ridx_block - orthrow;
                        row_end = current_region.startRow+orthrow * ridx_block-1;
                        matrix_tmp = HB(row_start:row_end, :);
                        matrix_sum = sum(matrix_tmp);
                        matrix_sum(matrix_sum == 0) = 1;  % 避免除以0
                        fal(i) = fal(i) + sum(matrix_sum - 1);
                    end
                    % 处理剩余行（若有）
                    if current_region.startRow + orthrow * ridx_block < current_region.endRow
                        row_start = current_region.startRow+orthrow * ridx_block;
                        row_end = current_region.endRow;
                        matrix_tmp = HB(row_start:row_end, :);
                        matrix_sum = sum(matrix_tmp);
                        matrix_sum(matrix_sum == 0) = 1;
                        fal(i) = fal(i) + sum(matrix_sum - 1);
                    end
                end
                
                % 5. 计算适应度（性能惩罚 + 正交约束惩罚）
                EbN0_delta = EbN0 - EbN0_re;
                EbN0_delta(EbN0_delta > 0) = relufactor * EbN0_delta(EbN0_delta > 0);  % 性能下降额外惩罚
                fits(i) = sum(EbN0_delta) + orth * fal(i) + alpha * nump(i);
            end
        end
        
        % -------------------------- 4.2.2 更新当前区域的最优解--------------------------
        for i = 1:length(t)
            if fits(i) < current_min_fits
                % 更新最小适应度
                current_min_fits = fits(i);
                % 重新计算当前个体的gamma值（排除惩罚项，仅性能指标）
                EbN0 = zeros(1, length(rowlist_tmp));
                for ridx = 1:length(rowlist_tmp)
                    row = rowlist_tmp(ridx);
                    % 1. 构建HB矩阵（-1→0，非-1→1）
                    CBN_r=row;VBN_r=22+row;
                    HB = zeros(CBN_r, VBN_r);
                    HB(t{i}(1:CBN_r,1:VBN_r)>= 0) = 1;
                    graph = HB;
                    % 2. 计算行/列非零索引
                    rowfind = cell(1, CBN_r);
                    colfind = cell(1, VBN_r);
                    for ii = 1:CBN_r
                        rowfind{ii} = find(graph(ii, :) == 1);
                    end
                    for jj = 1:VBN_r
                        colfind{jj} = find(graph(:, jj) == 1);
                    end
                    % 3. 计算EbN0（性能指标）
                    rate = 22/(20 + row);
                    tran_bit = 3:(22 + row);
                    EbN0(ridx) = threshold_find_qam(rate, CBN_r, VBN_r, tran_bit, colfind, rowfind, HB, graph, itrnum, qam, Z, layered_flag);
                end
                current_min_gamma = sum(EbN0);
                fal_val = fal(i);
                nump_val = nump(i);
                current_BL_opt = t{i};  % 更新当前区域的最优矩阵
                H_SV_cell{g+(region_idx-1)*G,1}=current_BL_opt;
                H_SV_cell{g+(region_idx-1)*G,2}=g+(region_idx-1)*G;
            end
        end
        
        % -------------------------- 4.2.3 输出当前迭代结果--------------------------
        fprintf('第%d个区域：迭代G=%d，当前最小gamma=%.4f，正交惩罚fal=%.2f， 元素数量nump=%d\n', ...
            region_idx, g, current_min_gamma, fal_val, nump_val);
        if mod(g, 50) == 0  % 每50次迭代输出一次当前最优矩阵（避免日志冗余）
            fprintf('第%d个区域迭代G=%d的最优矩阵：\n', region_idx, g);
            disp('H_SV = [ '); disp(num2str(current_BL_opt)); disp(' ];');
        end
        
        % -------------------------- 4.2.4 选择操作（原逻辑保留）--------------------------
        if g ~= 1  % 第一次迭代不选择（种群未更新）
            t1 = t;
            % 选择适应度最小的前(Population/2 + 1)个个体（精英保留）
            [~, popidx] = mink(fits, Population/2 + 1);
            % 复制精英个体到新种群
            for i = 1:Population/2
                t{2*i - 1} = t1{popidx(1)};    % 最优个体复制两次
                t{2*i} = t1{popidx(i + 1)};    % 其余精英个体复制一次
            end
            t(Population + 1:end) = [];  % 截断种群到原大小
        end
        
        % -------------------------- 4.2.5 交叉操作（仅在当前区域内进行）--------------------------
        t2 = t(1:Population);  % 父代种群
        for i = 1:floor(Population/2)
            temp1 = t{2*i};      % 父代1
            temp2 = t{2*i - 1};  % 父代2
            
            % 1. 提取当前区域的HB矩阵（仅关注0/1分布差异）
            HB1 = temp1(current_region.startRow:current_region.endRow, current_region.startCol:current_region.endCol) >= 0;
            HB2 = temp2(current_region.startRow:current_region.endRow, current_region.startCol:current_region.endCol) >= 0;
            crossHB = xor(HB1, HB2);  % 差异位置（仅这些位置需要交叉）
            [rows_cross, cols_cross] = find(crossHB);
            
            % 2. 映射差异位置到全局矩阵索引
            rows_cross = rows_cross + current_region.startRow - 1;
            cols_cross = cols_cross + current_region.startCol - 1;
            N_cross = length(rows_cross);  % 差异位置数量
            
            % 3. 交叉逻辑（原逻辑保留，仅操作当前区域差异位置）
            if N_cross > 0
                L_cross = randi(N_cross);  % 随机选择交叉位置数量
                selectedIdx = randperm(N_cross, L_cross);  % 选中的交叉位置
                remainingIdx = setdiff(1:N_cross, selectedIdx);  % 未选中的交叉位置
                
                % 交换父代1和父代2的选中位置
                selectedCoords = [rows_cross(selectedIdx), cols_cross(selectedIdx)];
                remainingCoords = [rows_cross(remainingIdx), cols_cross(remainingIdx)];
                
                t2{2*i - 1}(selectedCoords) = temp2(selectedCoords);
                t2{2*i - 1}(remainingCoords) = temp1(remainingCoords);
                t2{2*i}(selectedCoords) = temp1(selectedCoords);
                t2{2*i}(remainingCoords) = temp2(remainingCoords);
            else
                % 无差异时直接复制父代
                t2{2*i - 1} = temp2;
                t2{2*i} = temp1;
            end
        end
        t(Population + 1:2*Population) = t2;  % 交叉后代加入种群
        
        % -------------------------- 4.2.6 变异操作（仅在当前区域内进行）--------------------------
        t3 = t(1:Population);  % 待变异种群
        for i = 1:Population
            temp_BL = t{i};
            
            % 仅在当前区域内随机交换（-1和非-1，同初始化逻辑）
            L_mut = randi([2, 10]);  % 变异交换次数
            for l = 1:L_mut
                % while 1
                %     m1 = randi([current_region.startRow, current_region.endRow]);
                %     n1 = randi([current_region.startCol, current_region.endCol]);
                %     m2 = randi([current_region.startRow, current_region.endRow]);
                %     n2 = randi([current_region.startCol, current_region.endCol]);
                %     % 条件：位置不同 + 类型相反
                %     if (m1 ~= m2 || n1 ~= n2) && ...
                %        ((temp_BL(m1,n1) == -1 && temp_BL(m2,n2) >= 0) || ...
                %         (temp_BL(m1,n1) >= 0 && temp_BL(m2,n2) == -1))
                %         break;
                %     end
                % end
                % % 交换两点值（变异核心）
                % temp_Z = temp_BL(m1, n1);
                % temp_BL(m1, n1) = temp_BL(m2, n2);
                % temp_BL(m2, n2) = temp_Z;
                m1 = randi([current_region.startRow, current_region.endRow]);
                n1 = randi([current_region.startCol, current_region.endCol]);
                if temp_BL(m1, n1)>=0
                    temp_BL(m1, n1)=-1;
                else
                    temp_BL(m1, n1)=randi([0,Z-1],1,1);
                end
            end
            t3{i} = temp_BL;
        end
        t(2*Population + 1:3*Population) = t3;  % 变异后代加入种群
    end
    
    % -------------------------- 4.3 当前区域优化完成：更新全局最优--------------------------
    global_BL_opt = current_BL_opt;  % 当前区域的最优矩阵作为下一个区域的初始值
    global_min_gamma = current_min_gamma;  % 更新全局最小gamma
    fprintf('\n==================================== 第%d个区域优化完成 ====================================\n', region_idx);
    fprintf('第%d个区域优化后的最小gamma值：%.4f\n', region_idx, current_min_gamma);
    fprintf('第%d个区域优化后的最优矩阵：\n', region_idx);
    disp('H_SV = [ '); disp(num2str(global_BL_opt)); disp(' ];');
    fprintf('===========================================================================================\n\n');
end

%% ========================== 5. 所有区域优化完成：输出最终结果==========================
fprintf('所有%d个区域优化完成！最终全局最小gamma值：%.4f\n', n_regions, global_min_gamma);
fprintf('最终优化矩阵（H_SV）：\n');
disp('H_SV = [ '); disp(num2str(global_BL_opt)); disp(' ];');
save('H_SV_cell.mat','H_SV_cell');