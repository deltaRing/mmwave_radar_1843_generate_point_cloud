function [res, detect_res] = AngleCFAR(angle_range_map)
    % ������
    armap = angle_range_map;
    detect_res = [];
    res = [];
    if length(size(armap)) ~= 2
        fprintf("ARMap ��Ҫ��ά�ĽǶ�_������ͼ");
        return;
    end
    
    M=size(armap,2);
    R=32;%�ο���Ԫ��
    P=12; % ������Ԫ��
    L_slipper=R;%��������
    L_move=1;%�������
    L_num=floor((M-L_slipper)/L_move);%��������
    res = zeros(size(armap, 1), M - L_slipper); % Ԥ�ȷ������ݴ�С
    Z_map = [];
    d_map = [];
    
    for jj = 1:size(armap,1)
        for ii = 1:L_num
            start_idx = (ii-1)*L_move+1;
            Z = abs(sum(armap(jj, start_idx:start_idx+R/2-P))) / (R - 2 * P) + ...
                abs(sum(armap(jj, start_idx + R/2 + P:start_idx+L_slipper))) / (R - 2 * P);
            Z_map(jj, ii) = Z;
            d_map(jj, ii) = abs(armap(jj, ii+R/2));
            if Z < abs(armap(jj, ii+R/2))
                res(ii,jj) = armap(jj, ii+R/2);
                if localmaxValDetect(armap(jj, ii+R/2), armap(jj, start_idx:start_idx+R))
                    detect_res = [detect_res; jj, ii + R, armap(jj, ii+R/2)];
                end
            end
        end
    end
end


function res = localmaxValDetect(val, win)
    res = 0;
    if val == max(win) res = 1; end
end