x = 5:0.001:7;
x1 = 5:0.001:6;
x2 = 6:0.001:7;
y1 = cos(6*abs(x1));
x2 = x2(2:end);
y2 = sin(6.^x2);
plot(x1,y1,x2,y2);
title('f(x)����ͼ��');
xlabel('x');
ylabel('y');