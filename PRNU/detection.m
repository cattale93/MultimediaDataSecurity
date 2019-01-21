function [map] = detection(Img,finger_print)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    Noisex = NoiseExtractFromImage(Img,2);
    Noisex = WienerInDFT(Noisex,std2(Noisex));
    Ix = double(rgb2gray(Img));
    C = crosscorr(Noisex,Ix .* finger_print)
    map = PCE(C);
end

