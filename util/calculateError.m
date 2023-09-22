% ����·�����
% ����1��·����ʼλ��
% ����2��·������λ��
% ����3����⵽Ŀ��λ��
% ���1��ÿ��Ŀ��λ�õ�·����ƫ��
function derivation = calculateError(startLoc, EndLoc, targets)
PathLine = EndLoc - startLoc;        % ����б��
K = atan2(PathLine(2), PathLine(1)); % ����б��
B = startLoc(2) - K * startLoc(1);   % ����ƫ����
derivation = [];                     % ��¼ƫ��ֵ

for tt = 1:size(targets, 1)
    ttt = targets(tt, :);
    if abs(K - pi / 2) < 1e-5 % ��ֱ��
        startX = startLoc(1);
        endY1  = startLoc(2);
        endY2  = EndLoc(2);
        if ttt(2) > endY1 && ttt(2) < endY2
            derivation = [derivation; abs(startX - ttt(1))];
        else
            if ttt(2) < endY1
                derivation = [derivation; norm(ttt - startLoc)];
            elseif ttt(2) > endY2
                derivation = [derivation; norm(ttt - EndLoc)];
            end
        end
    elseif mod(K, pi) < 1e-5 % ˮƽ��
        startY = startLoc(2);
        endX1  = startLoc(1);
        endX2  = EndLoc(1);
        if ttt(1) < endX1 && ttt(1) > endX2
            derivation = [derivation; abs(startY - ttt(2))];
        else
            if ttt(1) > endX1
                derivation = [derivation; norm(ttt - startLoc)];
            elseif ttt(1) < endX2
                derivation = [derivation; norm(ttt - EndLoc)];
            end
        end
    else % ���������
        PathLine1 = startLoc - ttt;
        PathLine2 = EndLoc - ttt;
        K1 = atan2(PathLine1(2), PathLine1(1)); % б��1
        K2 = atan2(PathLine2(2), PathLine2(1)); % б��2

        if K1 <= K && K1 >= -1 / K && K2 <= K && K2 >= -1 / K
            dis = abs(K * ttt(1) + ttt(2) + B) / sqrt(K^2 + 1);
            derivation = [derivation; dis];
        else
            disFromStartLoc = norm(ttt - startLoc);
            disFromEndLoc   = norm(ttt - EndLoc);
            % �鿴�ĸ���������
            if disFromStartLoc > disFromEndLoc
                derivation = [derivation; disFromEndLoc];
            else
                derivation = [derivation; disFromStartLoc];
            end
        end
    end
end

end

