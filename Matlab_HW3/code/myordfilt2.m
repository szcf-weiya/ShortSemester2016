function [err, img] = myordfilt2(img0,origin, order,m,n)
img = ordfilt2(img0,order,ones(m,n));
err = norm(img - origin,2);
end