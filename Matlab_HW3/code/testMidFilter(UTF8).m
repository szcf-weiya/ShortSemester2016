subplot(2,2,1);
imshow(lena512);
title('`原始图像`');
subplot(2,2,2);
imshow(lena_saltpepper_05);
title('`椒盐噪声`');
% K1 = filter2(fspecial('average',3), lena_saltpepper_05)/255;
myimg = myMidFilter(lena_saltpepper_05, lena512);
subplot(2,2,3);
imshow(myimg);
title('`中值滤波`(myself)')
matimg = medfilt2(lena_saltpepper_05);
subplot(2,2,4);
imshow(matimg);
title('`中值滤波`(Matlab)');
disp('medfilt2 error');
norm(matimg - lena512)