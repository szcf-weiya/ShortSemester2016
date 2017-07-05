lena512 = imread('lena512.bmp');
lena_all_noise = imread('lena_all_noise.bmp');
lena_guassian_0_02 = imread('lena_gaussian_0_02.bmp');
lena_saltpepper_05 = imread('lena_saltpepper_05.bmp');

lena512 = im2double(lena512);
lena_all_noise = im2double(lena_all_noise);
lena_guassian_0_02 = im2double(lena_guassian_0_02);
lena_saltpepper_05 = im2double(lena_saltpepper_05);

subplot(2,2,1);
imshow(lena512);
title('Ô­Í¼');
subplot(2,2,2);
imshow(lena_saltpepper_05);
title('½·ÑÎÔëÉù');
subplot(2,2,3);
img = myMidFilter(lena_saltpepper_05,lena512);
imshow(img);
title('ÖĞÖµÂË²¨');
subplot(2,2,4);
imshow