subplot(2,3,1);
imshow(lena512);
title('原始图像');
subplot(2,3,2);
imshow(lena_guassian_0_02);
title('高斯噪声');
% K1 = filter2(fspecial('average',3), lena_saltpepper_05)/255;
myimg = myMidFilter(lena_guassian_0_02, lena512);
subplot(2,3,3);
imshow(myimg);
title('中值滤波')
%matimg = medfilt2(lena_guassian_0_02);

subplot(2,3,4);
avg = myAverageFilter(lena_guassian_0_02,lena512);
imshow(avg);
norm(avg-lena512)
title('均值滤波');

%imshow(matimg);
%title('中值滤波(Matlab)');
%disp('medfilt2 error');
%norm(matimg - lena512)

%img = diffusionAnisotropic(lena_guassian_0_02, 'function', 'complex', 'sigmaGauss', 0.01, 'time', 400 , 'maxIter', 1000);
%imshow(img);
%title('0.01,371');
%norm(img-lena512)

%img = diffusionAnisotropic(lena_guassian_0_02, 'function', 'complex', 'sigmaGauss', 0.01, 'time', 750 , 'maxIter', 1000);
%imshow(img);
%norm(img-lena512)
%title('0.01,271');
%imgAnisotropic = diffusionAnisotropic(lena_guassian_0_02, 'function', 'tuckey', 'sigma', 0.004, 'time', 2 , 'maxIter', 300);
%imshow(imgAnisotropic);
%norm(imgAnisotropic-lena512)
%subplot(2,2,3);
%0.002
%imgAnisotropic = diffusionAnisotropic(lena_guassian_0_02, 'function', 'tuckey', 'sigma', 0.00010656, 'time', 2 , 'maxIter', 300);
%imshow(imgAnisotropic);
%norm(imgAnisotropic-lena512)
subplot(2,3,5);
imgAnisotropic = diffusionAnisotropic(lena_guassian_0_02, 'function', 'tuckey', 'sigma', 0.00010655, 'time', 2 , 'maxIter', 300);
imshow(imgAnisotropic);
title('AF滤波(sigma=0.0001)');
norm(imgAnisotropic-lena512)
subplot(2,3,6);
imgAnisotropic = diffusionAnisotropic(lena_guassian_0_02, 'function', 'perona', 'sigma', 5.1, 'time', 1, 'maxIter', 300);
norm(imgAnisotropic-lena512)
imshow(imgAnisotropic);
title('AF滤波(sigma=5.1)');