
close all; clear; clc;

%% INITIALIZATION
gt=imread('C:\Users\Alessandro\Documents\MATLAB\2Competition\dev-dataset\dev-dataset-maps\dev_0030.bmp');
% Question window settings
opts.Interpreter = 'tex';
opts.Default = 'No';
quest = ['Before execute code be sure to be in the projet folder ',...
         '\bf"stakkastakka"\rm. ',...
         'If you are ready you will be asked to provide the dataset \bfimage\rm path.',...
         '                                                       Are you ready?'];
   
answer = questdlg(quest,'Info',...
                  'Yes','No',opts);

if ~(strcmp('Yes',answer))
    fprintf('User not ready\n')
    return;
end

%add folder project to the path
addpath(genpath(pwd));
PROJECT_FOLDER = pwd;

RESULT_FOLDER = [PROJECT_FOLDER '\DEMO_RESULTS'];

%% DATASET PATH

% asks for dataset path 
DATASET_PATH = uigetdir();

MAP_GT_PATH = uigetdir();

% add dataset path
addpath(genpath(DATASET_PATH));

addpath(genpath(MAP_GT_PATH));

cd (DATASET_PATH)
image_name = dir;
image_name = {image_name.name};

cd (PROJECT_FOLDER)

cd (MAP_GT_PATH)
map_name = dir;
map_name = {map_name.name};

cd (PROJECT_FOLDER)

%% IMAGE ANALISYS

N_of_detect_strategies = 5;
low_thresh = 1;
high_thresh = 90;
map_DIMx = 300;
map_DIMy = 400;
mean_reg_low = 15;
mean_reg_high = mean_reg_low * 6;
active_output=1;
map = zeros(300,400,5);
alg_name = {'map_1_JPG','map_2_JPG','map_BRC','map_CFA','map_PRNU'};

first_image=1;



for image_index = 803: 803
%% Map creation
    try
        if ~(contains(image_name{image_index},'jpg'))
            if active_output
                fprintf('Running single JPG detection\n');
            end
            tmp = imresize(get_map_1_JPG([DATASET_PATH '\' image_name{image_index}]),[map_DIMx,map_DIMy]);
            map(:,:,1) = tmp;
        else
            map(:,:,1) = NaN; 
        end
    catch
        map(:,:,1) = NaN;   
    end
    
    try
        if ~(contains(image_name{image_index},'jpg'))
             map(:,:,2) = NaN; 
        else
            if active_output
                fprintf('Running double JPG detection\n')
            end
            tmp = imresize(get_map_2_JPG([DATASET_PATH '\' image_name{image_index}]),[map_DIMx,map_DIMy]);
            map(:,:,2) = tmp;
        end
    catch
        map(:,:,2) = NaN;   
    end
    
    try
        if ~(contains(image_name{image_index},'jpg'))
            if active_output
                fprintf('Running BRC detection\n')
            end
            tmp = imresize(get_map_BRC([DATASET_PATH '\' image_name{image_index}]),[map_DIMx,map_DIMy]);
            map(:,:,3) = tmp;
        else
            map(:,:,3) = NaN; 
        end
    catch
        map(:,:,3) = NaN;   
    end
    
    try
        if active_output
           fprintf('Running CFA detection\n')
        end
        tmp = imresize(get_map_CFA([DATASET_PATH '\' image_name{image_index}]),[map_DIMx,map_DIMy]);
        map(:,:,4) = tmp;
    catch
        map(:,:,4) = NaN;   
    end
    
    try
        if active_output
            fprintf('Running PRNU detection\n')
        end
        tmp = imresize(get_map_PRNU([DATASET_PATH '\' image_name{image_index}]),[map_DIMx,map_DIMy]);
        map(:,:,5) = tmp;
    catch
        map(:,:,5) = NaN;   
    end

%% Map Analysis
    im_tf_v = zeros(map_DIMx,map_DIMy);
    im_tf_o = zeros(map_DIMx,map_DIMy);
    failing = 0;
    
    border_length = ones(5,802) * 10000000;
    mean_m = ones(5,802) * 10000000;
    mean_v = ones(5,802) * 10000000;
    mean_o = ones(5,802) * 10000000;
    F = ones(5,802) * 10000000;
    F_C = ones(1,802) * 10000000;
    map_gt = imread([MAP_GT_PATH '\' map_name{image_index}]);
    for map_index = 1 : N_of_detect_strategies
        %% MAP SELECTION - AREA DIMENSION BASED
               
        output_map = map(:,:,map_index);
        tmp = imresize(output_map,[1500,2000]);
        tmp = im2bw(tmp,graythresh(tmp));
        % computation of F_measure, if not needed comment it
        
        F(map_index,image_index) = f_measure(map_gt,tmp); 
        
        if active_output
            fprintf('\n----------------');
            fprintf('Analysing %s \n',alg_name{map_index});
        end
        
        N_tampered_pixels = length(find(output_map==1));
        N_pixels_total = map_DIMx * map_DIMy;
        pixel_ratio=(N_tampered_pixels/N_pixels_total)*100;
        
        if active_output
           fprintf('Pixel ratio %d \n',pixel_ratio);  
        end
        
        if (isnan(output_map))
            failing = failing + 1;
            if active_output
                fprintf('NaN \n');
            end
        elseif ((pixel_ratio < low_thresh) || (pixel_ratio > high_thresh))
            failing = failing + 1;
        else
            
        %% MAP SELECTION - HIGH FREQUENCY BASED   
        [Gmag,~] = imgradient(output_map);
        border_length(map_index,image_index) = length(find(Gmag>0.01));
        
        %fft computation and image procesing
        for j=1:map_DIMy
            im_tf_v(:,j) = fft(output_map(:,j));
        end
        
        im_tras = im_tf_v';
        tmp = mean(im_tras);
        mean_v(map_index,image_index) = abs(mean(tmp(mean_reg_low:mean_reg_high)));
        
        
        for j=1:map_DIMx
            im_tf_o(j,:) = fft(output_map(j,:));
        end
        
        im_tras = im_tf_o';
        tmp = mean(im_tras);
         mean_o(map_index,image_index) = abs(mean(tmp(mean_reg_low:mean_reg_high)));
        
        im_tf_m = fft(output_map);
        tmp = im_tf_m((mean_reg_low:mean_reg_high),(mean_reg_low:mean_reg_high));
        mean_m(map_index,image_index) = abs(mean(tmp(:)));
        
        imwrite(output_map,[RESULT_FOLDER,'\trials\', num2str(image_index-2),alg_name{map_index}, '.bmp'],'bmp');
        
        end
    end
    
%     fprintf('Border_length:\n');
%     fprintf('%.3f | %.3f | %.3f | %.3f | %.3f \n',border_length(1),border_length(2),border_length(3),border_length(4),border_length(5));
%     fprintf('Mean_m:\n');
%     fprintf('%.3f | %.3f | %.3f | %.3f | %.3f \n',mean_m(1),mean_m(2),mean_m(3),mean_m(4),mean_m(5));
%     fprintf('Mean_o:\n');
%     fprintf('%.3f | %.3f | %.3f | %.3f | %.3f \n',mean_o(1),mean_o(2),mean_o(3),mean_o(4),mean_o(5));
%     fprintf('Mean_v:\n');
%     fprintf('%.3f | %.3f | %.3f | %.3f | %.3f \n',mean_v(1),mean_v(2),mean_v(3),mean_v(4),mean_v(5));
    
%     %% Save best Output Map
%     fprintf('\n-------------------\nFailed alg = %d\n',failing);
%     output_map = ones(1500,2000);
%     if (failing==5)
%         output_map = ones(1500,2000);
%         output_map(500:999,500:1299) = zeros(500,800);
%     else
%         
%         best_m = find(mean_m == min(mean_m));
%         if (length(best_m) > 1) 
%             if ~isempty(find(best_m == 2))
%                 best_m = 2;
%             elseif ~isempty(find(best_m == 3))
%                 best_m = 3;
%             elseif ~isempty(find(best_m == 4))
%                 best_m = 4;
%             elseif ~isempty(find(best_m == 5))
%                 best_m = 5;
%             else
%                 best_m = 1; 
%             end
%         end
%         fprintf('Best m %d\n',best_m);
%         
%         best_o = find(mean_o == min(mean_o));
%         if (length(best_o) > 1) 
%              if ~isempty(find(best_o == 2))
%                 best_o = 2;
%             elseif ~isempty(find(best_o == 3))
%                 best_o = 3;
%             elseif ~isempty(find(best_o == 4))
%                 best_o = 4;
%             elseif ~isempty(find(best_o == 5))
%                 best_o = 5;
%             else
%                 best_o = 1; 
%             end
%         end
%         fprintf('Best o %d\n',best_o);
%    
%         best_v = find(mean_v == min(mean_v));
%         if (length(best_v) > 1) 
%              if ~isempty(find(best_v == 2))
%                 best_v = 2;
%             elseif ~isempty(find(best_v == 3))
%                 best_v = 3;
%             elseif ~isempty(find(best_v == 4))
%                 best_v = 4;
%             elseif ~isempty(find(best_v == 5))
%                 best_v = 5;
%             else
%                 best_v = 1; 
%             end
%         end
%         fprintf('Best v %d\n',best_v);
%         
%         if (best_m == best_o)
%             best = best_m;
%         elseif (best_m == best_v) 
%             best = best_v;     
%         elseif (best_o == best_v)
%             best = best_o;  
%         else
%             if ((best_m >= 2) && (best_m >= 4)) || ((best_o >= 2) && (best_o >= 4)) || ((best_v >= 2) && (best_v >=4))
%                 best = find(border_length(2:4) == min(border_length(2:4)),1)+1;
%             else
%                 best = find([border_length(1),border_length(5)] == min([border_length(1),border_length(5)]),1);
%                 if (best == 2)
%                     best = 5;
%                 end
%             end
%         end
%         
%        if ((best == 5) && ((contains(image_name{image_index},'jpg'))) && ((best_o==2||(best_m==2)||(best_v==2))))
%             best = 5;
%        elseif (best == 5)
%            diff(:,1) = mean_m(1:4) - mean_m(5);
%            diff(:,2) = mean_o(1:4) - mean_o(5);
%            diff(:,3) = mean_v(1:4) - mean_v(5);
% 
%            a=min(diff(1));
%            b=min(diff(2));
%            c=min(diff(3));
%        
%            if ((a+b<3) || (a+c<3) || (b+c<3 ))
%                if (best_v~=5)
%                    best=best_v;
%                elseif (best_o~=5)
%                    best=best_o;
%                elseif (best_m~=5)
%                    best=best_m;
%                else
%                    best=5;
%                end
%            else
%                best = 5;
%            end
%         end
%         
%         output_map = imresize(map(:,:,best),[1500,2000]);
%         output_map = im2bw(output_map,graythresh(output_map));
%         %imwrite(output_map,[RESULT_FOLDER '\map_',num2str(image_index-2), '.bmp'],'bmp');
%     end
%     imwrite(output_map,[RESULT_FOLDER '\map_',num2str(image_index-2), '.bmp'],'bmp');
%     F_C(1,image_index) = F(best,image_index);
% end


fprintf('\n-------------------\nFailed alg = %d\n',failing);
    output_map = ones(1500,2000);
    if (failing==5)
        output_map = ones(1500,2000);
        output_map(500:999,500:1299) = zeros(500,800);
    else  
        
        best_m = find(mean_m(:,image_index) == min(mean_m(:,image_index)));
        if (length(best_m) > 1) 
            if ~isempty(find(best_m == 2))
                best_m = 2;
            elseif ~isempty(find(best_m == 3))
                best_m = 3;
            elseif ~isempty(find(best_m == 4))
                best_m = 4;
            elseif ~isempty(find(best_m == 5))
                best_m = 5;
            else
                best_m = 1; 
            end
        end
        %fprintf('Best m %d\n',best_m);
        
        best_o = find(mean_o(:,image_index) == min(mean_o(:,image_index)));
        if (length(best_o) > 1) 
             if ~isempty(find(best_o == 2))
                best_o = 2;
            elseif ~isempty(find(best_o == 3))
                best_o = 3;
            elseif ~isempty(find(best_o == 4))
                best_o = 4;
            elseif ~isempty(find(best_o == 5))
                best_o = 5;
            else
                best_o = 1; 
            end
        end
        %fprintf('Best o %d\n',best_o);
   
        best_v = find(mean_v(:,image_index) == min(mean_v(:,image_index)));
        if (length(best_v) > 1) 
             if ~isempty(find(best_v == 2))
                best_v = 2;
            elseif ~isempty(find(best_v == 3))
                best_v = 3;
            elseif ~isempty(find(best_v == 4))
                best_v = 4;
            elseif ~isempty(find(best_v == 5))
                best_v = 5;
            else
                best_v = 1; 
            end
        end
        %fprintf('Best v %d\n',best_v);
        
        if (best_m == best_o)
            best = best_m;
        elseif (best_m == best_v) 
            best = best_v; 
        elseif (best_o == best_v)
            best = best_o; 
        else
            if ((best_m >= 2) && (best_m >= 4)) || ((best_o >= 2) && (best_o >= 4)) || ((best_v >= 2) && (best_v >=4))
                best = find(border_length(2:4,image_index) == min(border_length(2:4,image_index)),1)+1;
                
            else
                best = find([border_length(1,image_index),border_length(5,image_index)] == min([border_length(1,image_index),border_length(5,image_index)]),1);
               
                if (best == 2)
                    k=k+1
                    best = 5;
                end
            end
        end
   if ((best == 5) && ((contains(image_name{image_index},'jpg'))) && ((best_o==2||(best_m==2)||(best_v==2))))
       
      best=2;
   
   elseif (best == 5)
       diff(1:4,1) = mean_m(1:4,image_index)-mean_m(5,image_index);
       diff(1:4,2) = mean_o(1:4,image_index)-mean_o(5,image_index);
       diff(1:4,3) = mean_v(1:4,image_index)-mean_v(5,image_index);
       
       a=min(diff(1));
       b=min(diff(2));
       c=min(diff(3));
       
       if ((a+b<3) || (a+c<3) || (b+c<3 ))
           if (best_v~=5)
               best=best_v;
           elseif (best_o~=5)
               best=best_o;
           elseif (best_m~=5)
               best=best_m;
           else
               best=5;
           end
       else
           best = 5;
       end
   end
    
       output_map = imresize(map(:,:,best),[1500,2000]);
       output_map = im2bw(output_map,graythresh(output_map));
       %imwrite(output_map,[RESULT_FOLDER '\map_',num2str(image_index-2), '.bmp'],'bmp');
    end

    imwrite(output_map,[RESULT_FOLDER '\map_',num2str(image_index-2), '.bmp'],'bmp');
    F_C(1,image_index) = F(best,image_index);
 end
    