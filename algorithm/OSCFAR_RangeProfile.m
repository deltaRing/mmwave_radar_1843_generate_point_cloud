% ��Ծ������OS-CFAR�㷨
% ���룺������ RangeProfile
% ���룺cell_num ����Ԫֵ
% ���룺gain ����ֵ
% ���룺L����Ϳɽ��ܵ�Ԫ
% ���룺H����߿ɽ��ܵ�Ԫ
% �����result�������
function result = OSCFAR_RangeProfile(rangeprofile, cell_num, gain, L, H)
    if nargin < 2
        cell_num = 32;
        gain = 0.075; % ������õ�
        L = 8;
        H = 24;
    end
    
    result = []; % �����ʾ
    for tt = 1:size(rangeprofile, 2)
        for rr = 1:size(rangeprofile, 1) 
            if rr + cell_num > size(rangeprofile, 1)
                break; % ��������Ƿ񳬱꣬����һ�������ϵ�OS-CFAR
            end
            test_index       = rr + cell_num / 2;                 % ������� 
            test_units       = rangeprofile(test_index, tt);      % ���Ե�Ԫ
            back_ground_unit = [rangeprofile(rr:test_index - 1, tt) ...
                rangeprofile(test_index+1:rr + cell_num, tt)];    % ������Ԫ
            sorted_units     = sort(abs(back_ground_unit));       % �����Ԫ
            selected_units   = sorted_units(L:H);                 % ѡ���Ԫ
            sum_units        = sum(selected_units);               % ���
            gain_units       = gain * sum_units;                  % ��������
            if abs(test_units) >= gain_units
                result(rr, tt) = test_units;
            else
                result(rr, tt) = 0;
            end
        end
    end
    
end