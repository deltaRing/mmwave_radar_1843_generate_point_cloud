% 初始化EKF状态
% 输入1：X X坐标
% 输入2：Y Y坐标
% 输入3：Z Z坐标
% 输入4：Xn 原来的状态 
% 输入5：Pn 原来的协方差
% 可选输入1：EKF2 使用2维 还是 3维的EKF
% 输出1：Xn 状态
% 输出2：Pn 协方差
function [Xn, Pn] = init_EKF(X, Y, Z, Xn, Pn, EKF2)
    if nargin == 5
        EKF2 = 0;
    end
    initNum = length(Pn);
    for ii = 1:length(X)
        if ~EKF2
            Xn = [Xn; X(ii) 0.0 Y(ii) 0.0 Z(ii) 0.0];
            Pn{initNum + ii} = eye(6);
        else
            Xn = [Xn; X(ii) 0.0 Y(ii) 0.0];
            Pn{initNum + ii} = eye(4);
        end
    end
end