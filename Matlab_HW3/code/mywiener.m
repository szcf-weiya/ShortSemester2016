function [err, img1] = mywiener(img0, origin, k)
img1 = wiener2(img0, [k k]);
err = norm(origin - img1);
end