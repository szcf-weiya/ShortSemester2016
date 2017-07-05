subplot(1,2,1);
x = 1:0.001:4;
y1 = sqrt(x.^2 + 6) + x;
plot(x, y1,'r:');
xlabel('x');
ylabel('y');
title('y = $$\sqrt{(x^2+6)}$$+x','interpreter','latex');
axis equal
subplot(1,2,2);
y2 = sin(x.^2);
plot(x, y2,'g-.');
title('y=$$sin(x^2)$$','interpreter','latex');
xlabel('x');
ylabel('y');
axis square