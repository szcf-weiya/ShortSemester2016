function img = myAverageFilter(img0, origin)
N = size(img0, 1);
err = 100000; % a large number
threshold = 1e-3;
img = zeros(N, N);
maxIter = 10000;
errstop = 0.001;
for iter = 1:maxIter
    for i = 1:N
        for j = 1:N
            if i==1 && j==1
                img(i,j) = (img0(i,j)+img0(i,j+1)+img0(i+1,j))/3;
            elseif i==1 && j==N
                img(i,j) = (img0(i,j)+img0(i,j-1)+img0(i+1,j))/3;
            elseif i==N && j==1
                img(i,j) = (img0(i,j)+img0(i,j+1)+img0(i-1,j))/3;
            elseif i==N && j==N
                img(i,j) = (img0(i,j)+img0(i,j-1)+img0(i-1,j))/3;
            elseif i==1 && j~=1 && j~=N
                img(i,j) = (img0(i,j)+img0(i,j+1)+img0(i+1,j)+img0(i,j-1))/4;
            elseif i==N && j~=1 && j~=N
                img(i,j) = (img0(i,j)+img0(i,j+1)+img0(i-1,j)+img0(i,j-1))/4;
            elseif j==1 && i~=1 && i~=N
                img(i,j) = (img0(i,j)+img0(i,j+1)+img0(i-1,j)+img0(i+1,j))/4;
            elseif j==N && i~=1 && i~=N
                img(i,j) = (img0(i,j)+img0(i,j-1)+img0(i-1,j)+img0(i+1,j))/4;
            else
                img(i,j) = (img0(i,j)+img0(i,j-1)+img0(i-1,j)+img0(i+1,j)+img0(i,j+1))/5;
            end
            if abs(img(i,j)-img0(i,j)) > threshold
                img(i,j) = img0(i,j);
            end
        end
    end
    message = [iter, norm(img-img0)];
    disp(message);
    if norm(img0 - img) < errstop || norm(img0 - img) >= err
        img = img0;
        message = ['After ',num2str(iter),' iterations reach err ',num2str(norm(origin - img0))];
        disp(message);
        break;
    end
    err = norm(img0 - img);
end
if iter == maxIter
    message = ['After ',num2str(iter),' iterations reach err',num2str((norm(origin - img))),'but it can be continue when you change the max iterations'];
    disp(message);
end
end