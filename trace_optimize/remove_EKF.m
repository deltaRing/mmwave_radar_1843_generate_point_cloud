% 删除EKF状态
% 输入1：Xn 所有的状态
% 输入2：Pn 所有的协方差
% 输入3：removeIndex 删除的索引
function [Xn, Pn] = remove_EKF(Xn, Pn, removeIndex)
    Xn(removeIndex, :) = [];
    Pn(removeIndex)    = [];
end