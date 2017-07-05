function paraTime = myAnisotropic(origin,img0)
    N = 300;
    err = zeros(1,N);
    for i = 30:N
        img = diffusionAnisotropic(img0, 'function', 'complex', 'sigmaGauss', 0.2, 'time', i , 'maxIter', 1000);
        err(i) = norm(img - origin);
        disp('...');
    end
    [~,index] = min(err(1,30:N));
    
    paraTime = index;
    
end