%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com Date:
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-07-01 18:35:58
%FilePath: \exercise\exercise9_simple.m
%Description: exercise: 只取一个lambda计算
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
% tic
currentDir = pwd;
parentDir = fileparts(currentDir);
folderToAdd = fullfile(parentDir, '函数');
addpath(folderToAdd)
subStreamNumberDefault = 2132;
run('InitScript.m');
addpath(genpath('D:\Users\xushe\Documents\GitHub\ssm_tools\'))
addpath(genpath('F:\git\ssm_tools\'))
addpath(genpath('C:\Users\shengyix\Documents\GitHub\ssm_tools'))
addpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\2课程学习\卡尔曼滤波\从加速度积分至位移\Ma")
%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 6; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;




%% plot

if fig_bool == ON
    figureIdx = figureIdx + 1;
    row_idx = floor((figureIdx - 1) / num_figs_in_row);
    col_idx = mod((figureIdx - 1), num_figs_in_row);
    figPos_temp = figPos;
    figPos_temp(1:2) = figPos(1:2) + [col_idx * (figPos(3) + gap_between_images(1)) row_idx * (figPos(4) + gap_between_images(2))];
    hFigure = figure('Position', figPos_temp);
    hAxes = axes(hFigure);
    hLineObj = plot(t, P);
    legend("measure", "filtered")
    set(hLineObj, 'LineWidth', lineWidthThin);
    title(['\lambda = ', num2str(lambda)])

end


