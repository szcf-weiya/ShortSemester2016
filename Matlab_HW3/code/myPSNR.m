function PSNR = myPSNR(origin, img)
[h,w] = size(origin);
MES = sum(sum((origin - img).^2))/(h*w);
PSNR = 20 * log10(255/sqrt(MES));
end