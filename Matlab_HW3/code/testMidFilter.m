subplot(2,2,1);
imshow(lena512);
title('ԭʼͼ��');
subplot(2,2,2);
imshow(lena_saltpepper_05);
title('��������');
% K1 = filter2(fspecial('average',3), lena_saltpepper_05)/255;
myimg = myMidFilter(lena_saltpepper_05, lena512);
subplot(2,2,3);
imshow(myimg);
title('��ֵ�˲�(myself)')
matimg = medfilt2(lena_saltpepper_05);
subplot(2,2,4);
imshow(matimg);
title('��ֵ�˲�(Matlab)');
disp('medfilt2 error');
norm(matimg - lena512)