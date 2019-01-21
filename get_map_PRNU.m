function [map] = get_map_PRNU(image_name)
    
    %% Inizialisation
    %Camera fingerprints loading
    load('fingerprints.mat');
    block_size = 25;
    step_size = 5;
    x = 1;
    y = 1;
    imx = imread(image_name);
    [dimx,dimy,~] = size(imx);
    
    %Manage RGBA Images
    [~,~,colour] = size(imx);
    if colour == 4
        imx = imx(:,:,1:3);
    end
    
    %% Camera identification
    Noisex = NoiseExtractFromImage(imx,2);
    Noisex = WienerInDFT(Noisex,std2(Noisex));
    Ix = double(rgb2gray(imx));

    C(:,:,1) = crosscorr(Noisex,Ix .* fingerprints.c1);
    C(:,:,2) = crosscorr(Noisex,Ix .* fingerprints.c2);
    C(:,:,3) = crosscorr(Noisex,Ix .* fingerprints.c3);
    C(:,:,4) = crosscorr(Noisex,Ix .* fingerprints.c4);

    camera(1) = PCE(C(:,:,1));
    camera(2) = PCE(C(:,:,2));
    camera(3) = PCE(C(:,:,3));
    camera(4) = PCE(C(:,:,4));
    
    %Decide which camera took the photo
    camera = find(camera==max(camera));
    
    
    %Recompute specific camera parameters
    Noisex = NoiseExtractFromImage(imx,2);
    Noisex = WienerInDFT(Noisex,std2(Noisex));
    
    %% Tampered area localization
    for i=1:step_size:dimx-block_size
    
        for j=1:step_size:dimy-block_size
            C = crosscorr(Noisex(i:block_size+i-1,j:block_size+j-1,:),Ix(i:block_size+i-1,j:block_size+j-1,:) .* fingerprints.(['c' num2str(camera)])(i:block_size+i-1,j:block_size+j-1));
            tmp = PCE(C);
            map(x,y) = tmp;
            y=y+1;
        end
        
        y=1;
        x=x+1;
    end
    
    %% Output map filtering phase
    
    map=im2bw(map,graythresh(map));
    map=not(map);
    map=bwareaopen(map,500);
    map = imgaussfilt(double(map),3);
    map=im2bw(map,graythresh(map));
    map = imfill(map,'holes');
end