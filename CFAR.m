function [res, det] = CFAR(data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CA-CFAR %%%%%%%%%%%%%%%%%%%%%%%%%
%     data = abs(data);
    res = [];
    det = []; % 
    if length(size(data)) >= 3
        fprintf("2D data is required!");
        return;
    end
    M=size(data,1);
    R=16;%�ο���Ԫ��
    P=5; % ������Ԫ��
    L_slipper=R;%��������
    L_move=1;%�������
    L_num=floor((M-L_slipper)/L_move);%��������
    res = zeros(M - L_slipper, size(data, 2)); % Ԥ�ȷ������ݴ�С
    for jj = 1:size(data,2)
        for ii = 1:L_num
            start_idx = (ii-1)*L_move+1;
            Z = abs(sum(data(start_idx:start_idx+R/2-P, jj))) / (R - 2 * P) + ...
                abs(sum(data(start_idx + R/2 + P:start_idx+L_slipper, jj))) / (R - 2 * P);
            if Z < abs(data(ii+R/2, jj))
                res(ii,jj) = data(ii+R/2, jj);
                if localmaxValDetect(data(ii+R/2), data(start_idx:start_idx+R))
                    det = [det; ii + R, jj, data(ii+R/2)];
                end
            end
        end
    end
end

function res = localmaxValDetect(val, win)
    res = 0;
    if val == max(win) res = 1; end
end