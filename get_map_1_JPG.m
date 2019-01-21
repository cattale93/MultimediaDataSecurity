function [Output_map]=get_map_1_JPG(image_path)

    %% Inizialization
    block_size = 100;
    step_size = 10;
    x = 1;
    y = 1;
    
    %% Image pre-processing
    im_c =imread(image_path); 
    im_c = RGB2YCbCr(im_c);
    dct_c = bdct(im_c(:,:,1)-128);
   
    [dimx,dimy] = size(dct_c);

    %% Map extraction
    for i=1:step_size:dimx-block_size+step_size

        for j=1:step_size:dimy-block_size+1
            
            if (1496<block_size+i-1)
                i=i-4;
            end
            logL0_c(x,y) = L0_test(dct_c(i:block_size+i-1,j:block_size+j-1,:),2:10,8*pi);
            y=y+1;
        end
        y=1;
        x=x+1;
    end
    
    %% Post processing 
    temp = im2bw(logL0_c,graythresh(logL0_c));
    temp = imgaussfilt(double(temp),2);
    temp = im2bw(temp,graythresh(temp));
    Output_map = bwareaopen(temp,1000);

end


