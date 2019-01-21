tic;
close all; clear; clc;

%% Constants

GROUP_NAME = 'stakkastakka';
IMAGE_NAME = 'mulino';
IMAGE_FORMAT = 'bmp';
w_FILE_NAME = 'stakkastakka.mat';
ORIGINAL_IM_NAME = [IMAGE_NAME '.' IMAGE_FORMAT];
w_IM_NAME = [IMAGE_NAME '_' GROUP_NAME '.' IMAGE_FORMAT];
N_OF_w=1000;
%Parameter
step=8;
l=step-1;
rng(65);

%% Read Images
I = imread(ORIGINAL_IM_NAME,IMAGE_FORMAT);
I_m = imread(w_IM_NAME,IMAGE_FORMAT);

%% DWT
[LL_o,~,~,~] = dwt2(I,'haar');
[LL_m,~,~,~] = dwt2(I_m,'haar');

%% DCT block extraction from I
D_o=zeros(size(LL_o)/step);
i_b=1;
j_b=1;

for i=1:step:(length(LL_o))-l
    for j=1:step:(length(LL_o))-l
        tmp=dct2(LL_o(i:i+step-1,j:j+step-1));
        D_o(i_b,j_b)=tmp(step,step);
        j_b=j_b+1;
    end
    j_b=1;
    i_b=i_b+1;
end

 %%  DCT block extraction from I_m
D_m=zeros(size(LL_m)/step);
i_b=1;
j_b=1;

for i=1:step:(length(LL_m))-l
    for j=1:step:(length(LL_m))-l
        tmp=dct2(LL_m(i:i+step-1,j:j+step-1));
        D_m(i_b,j_b)=tmp(step,step);
        j_b=j_b+1;
    end
    j_b=1;
    i_b=i_b+1;
end


%% Block subtraction
D = D_m - D_o;

%% Watemark retrieving
%extraceted watermark
w_retr = D > 0;
w_retr = double(w_retr);

% Arnold inverse transform
w_retr = iarnold(w_retr,1000);

%reshape watermark extracted
w_retr = reshape(w_retr,1,length(w_retr)*length(w_retr));

%% Threshold computation
% load original watermark
load(w_FILE_NAME);

%reshape
w = reshape(w,1,length(w_retr));

% Creation of random watermark
w_rand = round(rand(N_OF_w,length(w_retr)));

% threshold
sim = zeros(1,N_OF_w);

for i=1:N_OF_w-1
    sim(i) = abs((w_rand(i,:)*w_retr') / sqrt((w_rand(i,:) * w_rand(i,:)')));
end

sim(N_OF_w) = abs(w_retr*w' / sqrt(w*w'));



sim_ordered = sort(sim,'descend');
t = sim_ordered(2);
T = t+0.1*t;

fprintf('Threshold value = %f\n',T);

toc;