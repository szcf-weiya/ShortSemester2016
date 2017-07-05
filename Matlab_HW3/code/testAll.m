subplot(3,3,1);
imshow(lena512);
title('ԭʼͼ��');
subplot(3,3,2);
imshow(lena_all_noise);
title('ȫ����');
% K1 = filter2(fspecial('average',3), lena_saltpepper_05)/255;
myimg = myMidFilter(lena_all_noise, lena512);
disp('------------------Med----------------------------');
subplot(3,3,3);
imshow(myimg);
['��ֵ�˲�',num2str(norm(myimg-lena512))]
title('��ֵ�˲�');
%matimg = medfilt2(lena_guassian_0_02);
disp('--------------------Average-----------------------------');
subplot(3,3,4);
avg = myAverageFilter(lena_all_noise,lena512);
imshow(avg);
['��ֵ�˲�',num2str(norm(avg-lena512))]
title('��ֵ�˲�');
disp('---------------------wiener------------------------------');
subplot(3,3,5);
[err, wiener] = mywiener(lena_all_noise,lena512,5);
message = [num2str(err),'wiener'];
disp(message);
imshow(wiener);
title('ά���˲�');
disp('--------------------ord---------------------------------');
subplot(3,3,6);
[err, ord] = myordfilt2(lena_all_noise,lena512, 5,5,5);
err
imshow(ord);
title('ͳ�ƶ�ά�˲�');
disp('------------------gaussian--------------------------');
subplot(3,3,7);
h = fspecial('gaussian',3,0.02);
gaussianimg = imfilter(lena_all_noise,h);
imshow(gaussianimg);
norm(gaussianimg-lena512)
title('��˹�˲�(Matlab)');
disp('-----------------Anisotropic1-------------------');
subplot(3,3,8);
imgAnisotropic = diffusionAnisotropic(lena_all_noise, 'function', 'tuckey', 'sigma', 0.00010655, 'time', 4 , 'maxIter', 300);
imshow(imgAnisotropic);
title('AF�˲�(sigma=0.0001)');
norm(imgAnisotropic-lena512)
disp('----------------Anisotropicl--------------------');
subplot(3,3,9);
%imgAnisotropic = diffusionAnisotropic(lena_all_noise, 'function', 'perona', 'sigma', 5.1, 'time', 1, 'maxIter', 300);
imgAnisotropic = diffusionAnisotropic(lena_all_noise, 'function', 'tuckey', 'sigma', 5.1, 'time', 1, 'maxIter', 300);
norm(imgAnisotropic-lena512)
imshow(imgAnisotropic);
title('AF�˲�(sigma=5.1)');

subplot(2,2,1);
imshow(lena512);
title('ԭʼͼ��');
subplot(2,2,2);
imshow(lena_all_noise);
title('ȫ����');
subplot(2,2,3);
imshow(img3);
title('��ֵ�˲�+��ֵ�˲�');
subplot(2,2,4)
imshow(imgAnisotropic1);
title('��ֵ�˲�+AF�˲�');

