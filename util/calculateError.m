% 计算路径误差
% 输入1：路径起始位置
% 输入2：路径结束位置
% 输入3：检测到目标位置
% 输出1：每个目标位置到路径的偏离
function derivation = calculateError(startLoc, EndLoc, targets)
PathLine = EndLoc - startLoc;        % 计算斜率
K = atan2(PathLine(2), PathLine(1)); % 计算斜率
B = startLoc(2) - K * startLoc(1);   % 计算偏置项
derivation = [];                     % 记录偏离值

for tt = 1:size(targets, 1)
    ttt = targets(tt, :);
    if abs(K - pi / 2) < 1e-5 % 垂直线
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
    elseif mod(K, pi) < 1e-5 % 水平线
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
    else % 其他情况下
        PathLine1 = startLoc - ttt;
        PathLine2 = EndLoc - ttt;
        K1 = atan2(PathLine1(2), PathLine1(1)); % 斜率1
        K2 = atan2(PathLine2(2), PathLine2(1)); % 斜率2

        if K1 <= K && K1 >= -1 / K && K2 <= K && K2 >= -1 / K
            dis = abs(K * ttt(1) + ttt(2) + B) / sqrt(K^2 + 1);
            derivation = [derivation; dis];
        else
            disFromStartLoc = norm(ttt - startLoc);
            disFromEndLoc   = norm(ttt - EndLoc);
            % 查看哪个距离点最近
            if disFromStartLoc > disFromEndLoc
                derivation = [derivation; disFromEndLoc];
            else
                derivation = [derivation; disFromStartLoc];
            end
        end
    end
end

end

