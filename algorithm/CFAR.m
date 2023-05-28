function [res, det] = CFAR(data, disable_range)
    if nargin < 2
        disable_range = 64;
    end
    data(1:disable_range, :) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CA-CFAR %%%%%%%%%%%%%%%%%%%%%%%%%
%     data = abs(data);
    det = [];
    det_ = []; % 
    if length(size(data)) >= 3
        fprintf("2D data is required!");
        return;
    end
    M=size(data,1);
    R=16;%参考单元数
    P=5; % 保护单元数
    L_slipper=R;%滑窗长度
    L_move=1;%滑窗间隔
    ratio = 2.0;
    remove_ratio = 0.5;
    L_num=floor((M-L_slipper)/L_move);%滑窗次数
    res = zeros(M - L_slipper, size(data, 2));
    res_ = zeros(M - L_slipper, size(data, 2)); % 预先分配数据大小
    for jj = 1:size(data,2)
        for ii = 1:L_num
            start_idx = (ii-1)*L_move+1;
            middle_idx_1 = start_idx+floor(R/2)-P;
            middle_idx_2 = start_idx + floor(R/2) + P;
            end_idx = start_idx + R;
            index = [start_idx:middle_idx_1, middle_idx_2:end_idx];
            Z = sum(abs(data(index, jj))) / (R - P * 2);
            val = abs(data(ii+floor(R/2), jj));
            if Z * ratio < val
                res_(ii,jj) = data(ii+floor(R/2), jj);
%                 det_ = [det_; ii + R, jj, data(ii+floor(R/2)), abs(data(ii+floor(R/2)))];
%                 if localmaxValDetect(data(ii+floor(R/2)), data(start_idx:start_idx+R))
%                 end
            end
        end
    end
    
    %% 去除反射强度过低的目标
    max_det = max(max(abs(res_)));
    for ii = 1:size(res_, 1)
        for jj = 1:size(res_, 2)
            if abs(res_(ii, jj)) >= max_det * remove_ratio
                det = [det; ii, jj, res_(ii, jj)];
                res(ii, jj) = res_(ii, jj);
            end
        end
    end
    
end

function res = localmaxValDetect(val, win)
    res = 0;
    if val == max(win) res = 1; end
end