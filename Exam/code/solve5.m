function d = solve5(n)
c = rand(n);
d = zeros(n);
 for i = 1:n
    for j = 1:n
        if i==1 && j==1
            d(i,j) = (c(i,j)+c(i,j+1)+c(i+1,j))/5;
        elseif i==1 && j==n
            d(i,j) = (c(i,j)+c(i,j-1)+c(i+1,j))/5;
        elseif i==n && j==1
            d(i,j) = (c(i,j)+c(i,j+1)+c(i-1,j))/5;
        elseif i==n && j==n
            d(i,j) = (c(i,j)+c(i,j-1)+c(i-1,j))/5;
        elseif i==1 && j~=1 && j~=n
            d(i,j) = (c(i,j)+c(i,j+1)+c(i+1,j)+c(i,j-1))/5;
        elseif i==n && j~=1 && j~=n
            d(i,j) = (c(i,j)+c(i,j+1)+c(i-1,j)+c(i,j-1))/5;
        elseif j==1 && i~=1 && i~=n
            d(i,j) = (c(i,j)+c(i,j+1)+c(i-1,j)+c(i+1,j))/5;
        elseif j==n && i~=1 && i~=n
            d(i,j) = (c(i,j)+c(i,j-1)+c(i-1,j)+c(i+1,j))/5;
        else
            d(i,j) = (c(i,j)+c(i,j-1)+c(i-1,j)+c(i+1,j)+c(i,j+1))/5;
        end
    end
 end
end