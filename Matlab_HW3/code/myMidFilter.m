function img = myMidFilter(img0,origin)
N = size(img0, 1);
err = 100000; % a large number
threshold = 1e-3;
img = zeros(N, N);
maxIter = 10000;
errstop = 0.001;
message = ['iteration ', 'norm(img-img0) ', 'norm(img-origin)'];
disp(message);
for iter = 1:maxIter
    for i = 1:N
        for j = 1:N
            if i==1 && j==1
                tmp = sort([img0(i,j), img0(i,j+1), img0(i+1, j)]);                
                img(i, j) = tmp(2);
            elseif i==1 && j==N
                tmp = sort([img0(i,j), img0(i,j-1), img0(i+1, j)]);
                img(i, j) = tmp(2);                
            elseif i==N && j==1
                tmp = sort([img0(i, j), img0(i,j+1), img(i-1, j)]);
                img(i, j) = tmp(2);                
            elseif i==N && j==N
                tmp = sort([img0(i,j), img0(i,j-1), img0(i-1,j)]);
                img(i, j) = tmp(2);                
            elseif i==1 && j~=1 && j~=N
                tmp = sort([img0(i,j), img0(i,j+1), img0(i+1,j), img0(i,j-1)]);
                img(i,j) = 0.5 * (tmp(2) + tmp(3));
            elseif i==N && j~=1 && j~=N
                tmp = sort([img0(i,j), img0(i,j+1), img0(i-1,j), img0(i, j-1)]);
                img(i,j) = 0.5 * (tmp(2) + tmp(3));
            elseif j==1 && i~=1 && i~=N
                tmp = sort([img0(i,j), img0(i,j+1), img0(i-1,j), img0(i, j+1)]);
                img(i,j) = 0.5 * (tmp(2) + tmp(3));
            elseif j==N && i~=1 && i~=N
                tmp = sort([img0(i,j), img0(i,j-1), img0(i-1,j), img0(i+1, j)]);
                img(i,j) = 0.5 * (tmp(2) + tmp(3));
            else
                tmp = sort([img0(i,j), img0(i,j-1), img0(i-1,j), img0(i+1, j), img0(i, j+1)]);
                img(i,j) = tmp(3);
            end
            if abs(img(i,j)-img0(i,j)) > threshold
                img0(i,j) = img(i,j);
            end
        end
    end
    if norm(img0 - img) < errstop || norm(img0 - img) > err
        img = img0;
        message = [iter, norm(img-img0),norm(img - origin)];
        disp(message);
        message = ['After ',num2str(iter),' iterations reach err ',num2str(norm(origin - img0))];
        disp(message);
        break;
    end
    message = [iter, norm(img-img0),norm(img - origin)];
    disp(message);
    err = norm(img0 - img);
end
if iter == maxIter
    message = ['After ',num2str(iter),' iterations reach err',num2str((norm(origin - img))),'but it can be continue when you change the max iterations'];
    disp(message);
end
end 
