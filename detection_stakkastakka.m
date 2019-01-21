function [output_1,output_2] = detection_stakkastakka(input_1,input_2,input_3)


%% Constants
ORIGINAL_IM_NAME = input_1;
w_IM_NAME = input_2;
ATTACKED_IM_NAME = input_3;
IMAGE_FORMAT = 'bmp';
T=13.585565;
%Parameter
step=8;
l=step-1;
ARNOLD_ITER=3000;


%% Read Images
I = imread(ORIGINAL_IM_NAME,IMAGE_FORMAT);
I_m = imread(w_IM_NAME,IMAGE_FORMAT);
I_a =imread(ATTACKED_IM_NAME,IMAGE_FORMAT);

%% DWT
[LL_o,LH_o,HL_o,HH_o] = dwt2(I,'haar');
[LL_m,LH_m,HL_m,HH_m] = dwt2(I_m,'haar');
[LL_a,LH_a,HL_a,HH_a] = dwt2(I_a,'haar');

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

%% DCT block extraction from I_a
D_a=zeros(size(LL_a)/step);
i_b=1;
j_b=1;

for i=1:step:(length(LL_a))-l
    for j=1:step:(length(LL_a))-l
        tmp=dct2(LL_a(i:i+step-1,j:j+step-1));
        D_a(i_b,j_b)=tmp(step,step);
        j_b=j_b+1;
    end
    j_b=1;
    i_b=i_b+1;
end

%% Block subtraction
D_m = D_m - D_o;
D_a = D_a - D_o;

%% Watemark retrieving
%original watermark
w = D_m > 0;
w=double(w);

%attacked watermark
w_retr = D_a > 0;
w_retr=double(w_retr);

%% Arnold inverse transform
try
w= iarnold(w,ARNOLD_ITER);
w_retr= iarnold(w_retr,ARNOLD_ITER);
catch   
end
%% Decision
w=reshape(w,1,1024);
w_retr=reshape(w_retr,1,1024);

sim = abs(w_retr*w' / sqrt(w_retr*w_retr'));

presence = sim > T;
output_1 = presence;

wpsnr = WPSNR(uint8(I_m),uint8(I_a));
output_2 = wpsnr;


%% FUNCTION
% Arnold inverse tranform
function [ out ] = iarnold( in, iter )
    if (~ismatrix(in))
        error('Oly two dimensions allowed');
    end
    [m, n] = size(in);
    if (m ~= n)
        error(['Arnold Transform is defined only for squares. ' ...
        'Please complete empty rows or columns to make the square.']);
    end
    out = zeros(m);
    n = n - 1;
    for j=1:iter
        for y=0:n
            for x=0:n
                p = [ 2 -1 ; -1 1 ] * [ x ; y ];
                out(mod(p(2), m)+1, mod(p(1), m)+1) = in(y+1, x+1);
            end
        end
        in = out;
    end
end

end

