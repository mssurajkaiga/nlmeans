%DIRECTORY and FILE PARAMETERS
extension = '.png';
scene_name = 'street_scene';
denoised_online = strcat(scene_name,'\denoised_online', extension);
denoised = strcat(scene_name,'\denoised', extension);
noisy = strcat(scene_name, '\noisy', extension);
original = strcat(scene_name, '\original', extension);
%ANALYSIS PARAMETERS
scale = 10;

input1 = imread(denoised_online);
input2 = imread(denoised);
input3 = imread(noisy);
input4 = imread(original);
diff1 = input2 - input1;
diff2 = input3 - input1;
diff3 = input3 - input2;
diff4 = input4 - input2;
diff5 = input4 - input1;

% PLOT STUFF
ROWS = 3;
COLUMNS = 3;

figure(1);
subplot(ROWS,COLUMNS,1);
image(input4);
title('Reference image');

subplot(ROWS,COLUMNS,2);
image(input3);
title('Noisy image input');

subplot(ROWS,COLUMNS,3);
image(input2);
title('CPU Denoised output');

subplot(ROWS,COLUMNS,4);
image(input1);
title('Online Denoised output');

subplot(ROWS,COLUMNS,5);
image(diff1*scale);
title('Difference between cpu denoiser and online denoiser results');
colorbar;

subplot(ROWS,COLUMNS,6);
image(diff2*scale);
title('Difference between noisy image and online denoised image');
colorbar;

subplot(ROWS,COLUMNS,7);
image(diff3*scale);
title('Difference between noisy image and cpu denoised image');
colorbar;

subplot(ROWS,COLUMNS,8);
image(diff4*scale);
title('Difference between original image and cpu denoised image');
colorbar;

subplot(ROWS,COLUMNS,9);
image(diff5*scale);
title('Difference between original image and online denoised image');
colorbar;