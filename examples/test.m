clc
clear
close all

% 设置随机数种子以确保结果的可复现性
rng(123);

% 生成随机的残差协方差矩阵 R
R = randn(5,5);
R = R'*R;  % 使其为对称正定矩阵

% 生成随机的雅可比矩阵 J
J = randn(5,3);

% 计算左式
M1 = (J' / R * J) \ (J' / R);

% 计算右式
M2 = pinv(J);

% 显示两个表达式的结果
disp('Result of the left-hand side:');
disp(M1);
disp('Result of the right-hand side:');
disp(M2);
