tic;
close all; clear; clc;
warning('off')

%% Constants
GROUP_NAME = 'stakkastakka';
IMAGE_NAME = 'mulino';
IMAGE_FORMAT = 'bmp';
w_FILE_NAME = 'stakkastakka.mat';
SHOW_IMG=1;
%% parameter ...
STEP=8;
l=STEP-1;
k1 = 32;
k2 = 1.6;
ARNOLD_ITER=3000;

%% Read Image
%load image double datatype
I = imread([IMAGE_NAME '.' IMAGE_FORMAT],IMAGE_FORMAT);
[dimx,dimy] = size(I);
%Load watermark
load(w_FILE_NAME);

%% DWT
[LL_c,LH_c,HL_c,HH_c] = dwt2(I,'haar');

%% Arnold tranform

w=arnold(w,ARNOLD_ITER);

%% DCT block
D=zeros(size(LL_c)/STEP);
i_b=1;
j_b=1;

for i=1:STEP:((length(LL_c))-(STEP-1))
    for j=1:STEP:((length(LL_c))-(STEP-1))
        tmp=dct2(LL_c(i:i+STEP-1,j:j+STEP-1));
        D(i_b,j_b)=tmp(STEP,STEP);
        j_b=j_b+1;
    end
    j_b=1;
    i_b=i_b+1;
end
%% Watermark embedding
for i=1:length(w)
    for j=1:length(w)
        if w(i,j)==1
            D(i,j) = D(i,j) + k1 * k2;
        else
            D(i,j) = D(i,j) - k1;
        end
    end
end
%% ...
i_b=1;
j_b=1;

for i=1:STEP:((length(LL_c))-(STEP-1))
    for j=1:STEP:((length(LL_c))-(STEP-1))
        tmp=dct2(LL_c(i:i+STEP-1,j:j+STEP-1));
        tmp(STEP,STEP)=D(i_b,j_b);
        LL_star(i:i+STEP-1,j:j+STEP-1)=idct2(tmp);
        j_b=j_b+1;
    end
    j_b=1;
    i_b=i_b+1;
end
%% Create Watermarked image
I_m = uint8(idwt2(LL_star, LH_c, HL_c, HH_c, 'haar'));

%% Save watermarked image 
cd 'C:\Users\Alessandro\Documents\ALE\Scuola\UNI\Anno V\MDS\Battle';
imwrite(uint8(I_m),[IMAGE_NAME '_' GROUP_NAME '.bmp'],'bmp'); 

%% computation WPSNR
wpsnr_result = WPSNR(I_m,I);
fprintf('WPSNR = + %f dB\n',wpsnr_result)
toc;

%% Visualize image
if SHOW_IMG
    figure();
    imshow(I,[])
    figure();
    imshow(I_m,[])
    figure();
    imshow((I_m-I),[])
end

%% Embedding quality feedback
% print feedback criteria to evaluate embedding quality
PSNR(uint8(I_m),uint8(I));

fprintf('SIM = %f\n',abs(reshape(double(I_m),1,dimx*dimy)*reshape(double(I),1,dimx*dimy)' / sqrt(reshape(double(I_m),1,dimx*dimy)*reshape(double(I_m),1,dimx*dimy)')))

if ( wpsnr_result < 50 )
    fprintf('1 PUNTO - %d dB\n',wpsnr_result)
elseif ( 50 <= wpsnr_result ) && ( wpsnr_result < 54 )
    fprintf('2 PUNTI - %d dB\n',wpsnr_result)
elseif ( 54 <= wpsnr_result ) && ( wpsnr_result < 58 )
    fprintf('3 PUNTI - %d dB\n',wpsnr_result)
elseif ( 58 <= wpsnr_result ) && ( wpsnr_result < 62 )
    fprintf('4 PUNTI - %d dB\n',wpsnr_result)
elseif ( 62 <= wpsnr_result ) && ( wpsnr_result < 65)
    fprintf('5 PUNTI - %d dB\n',wpsnr_result)
elseif ( wpsnr_result > 65)
    fprintf('6 PUNTI - %d dB\n',wpsnr_result)
    if ( wpsnr_result > 68)
        fprintf('VALUE TOO HIGH - %d dB\n',wpsnr_result)
    end
end

%% FUNCRION
% Arnold Function
function [ out ] = arnold( in, iter )
    if (ndims(in) ~= 2)
        error('Oly two dimensions allowed');
    end
    [m n] = size(in);
    if (m ~= n)
        error(['Arnold Transform is defined only for squares. ' ...
        'Please complete empty rows or columns to make the square.']);
    end
    out = zeros(m);
    n = n - 1;
    for j=1:iter
        for y=0:n
            for x=0:n
                p = [ 1 1 ; 1 2 ] * [ x ; y ];
                out(mod(p(2), m)+1, mod(p(1), m)+1) = in(y+1, x+1);
            end
        end
        in = out;
    end
end

