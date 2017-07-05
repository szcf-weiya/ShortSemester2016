syms x
eqn = x^3+7*x+1 == 0;
solx = solve(eqn,x);
x = (-10:0.001:10);
y = x.^3+7*x+1;
plot(x,y)
S_vpa = vpa(solx);
plot(x,y);
hold on;
line([-10,10],[0,0]);
xlabel('x');
ylabel('y');
