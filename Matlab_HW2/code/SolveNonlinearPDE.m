function [err, step] = SolveNonlinearPDE
global N;
N = 10;
global h;
global x;
global y;
h = 1/(N+1);
maxLoop = 100;
errStop = 1.0e-6;

u = ones(N*N, 1);
exu = zeros(N*N, 1);

x = h * (1:N);
y = h * (1:N);

rowIdx = zeros(1,N*N+1);
F = zeros(N*N, 1);
for lp = 1:maxLoop
%% Assign value 
%% Matrix A
%% the first row
i = 1;
j = 1;
rowIdx(1) = 1;
colIdx(1) = 1;
entries(1) = p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) +p4func(u,i,j);
F(1,1) = h * RightF(x(1),y(1));
for i = 2:N-1
	% (i,j) (i-1,j) (i+1,j) (i,j+1)
	rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
	colIdx = [colIdx, i+(j-1)*N, i-1+(j-1)*N, i+1+(j-1)*N, i + j*N];
	entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) +p4func(u,i,j), -1*p1func(u, i, j) ,-1*p2func(u,i,j),-1*p4func(u, i,j) ];
	F(i+(j-1)*N) = h*RightF(x(i),y(j));
end
% (i,j) (i-1,j) (i,j+1)
rowIdx(N) = size(colIdx, 2) + 1;
colIdx = [colIdx, i+(j-1)*N, i-1+(j-1)*N, i+j*N];
entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j), -1*p1func(u, i, j), -1*p4func(u, i,j)]; 
F(i+(j-1)*N) = h*RightF(x(i), y(j));
for j = 2:N-1
	for i = 1:N
		if i == 1				
			% (i,j) (i,j-1) (i+1,j) (i,j+1)			
			rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
			colIdx = [colIdx, i+(j-1)*N, i+(j-2)*N, i+1+(j-1)*N, i+j*N ];
			entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j),-1*p3func(u,i,j), -1*p2func(u,i,j),-1*p4func(u,i,j)];
            F(i+(j-1)*N) = h*RightF(x(i),y(j));
		elseif i == N				
			% (i,j) (i,j-1), (i-1, j) (i, j+1)
			rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
			colIdx = [colIdx, i+(j-1)*N, i+(j-2)*N, i-1+(j-1)*N,i+j*N ];
			entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j),-1*p3func(u,i,j), -1*p1func(u,i,j), -1*p4func(u,i,j)];
            F(i+(j-1)*N) = h*RightF(x(i),y(j)) + p2func(u,i,j)*(u(N+(j-1)*N)+1/h*y(j));
		else	
			% (i,j) (i,j-1) (i-1,j) (i+1,j) (i,j+1)
			rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
			colIdx = [colIdx, i+(j-1)*N, i+(j-2)*N, i-1+(j-1)*N, i+1+(j-1)*N ,i+j*N];
			entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j),-1*p3func(u,i,j), -1*p1func(u,i,j), -1*p2func(u,i,j),-1*p4func(u,i,j)];
            F(i+(j-1)*N) = h*RightF(x(i),y(j));
		end
	end
end
% (i,j) (i,j-1) (i+1,j) 
j = N;
i = 1;
rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
colIdx = [colIdx, i+(j-1)*N, i+(j-2)*N, i+1+(j-1)*N];
entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j), -1*p2func(u,i,j), -1*p4func(u,i,j)];
F(i+(j-1)*N) = h*RightF(x(i),y(j)) + p4func(u,i,j)*(u(i+(N-1)*j)+1/h*x(i));
for i = 2:N-1
	rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
	colIdx = [colIdx, i+(j-1)*N,i+(j-2)*N,i-1+(j-1)*N,i+1+(j-1)*N];
	entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j),-1*p3func(u,i,j), -1*p1func(u,i,j), -1*p2func(u,i,j)];
	F(i+(j-1)*N) = h*RightF(x(i),y(j)) + p4func(u,i,j)*(u(i+(N-1)*j)+1/h*x(i));
end
i = N;
rowIdx(i+(j-1)*N) = size(colIdx, 2) + 1;
colIdx = [colIdx, i + (j-1)*N];
entries = [entries, p1func(u,i,j) + p2func(u,i,j) + p3func(u,i,j) + p4func(u,i,j)];
F(i+(j-1)*N) = h*RightF(x(i),y(j)) + p4func(u,i,j)*(u(i+(N-1)*j)+1/h*x(i)) + p2func(u,i,j)*(u(i+(j-1)*N)+1/h*y(j));
rowIdx(N*N+1) = N*N;
%%
u0 = ones(N*N,1);

	%%matspA = mysp2matsp(rowIdx, colIdx, entries);
    myfull = mysparse2full(rowIdx,colIdx,entries);
	u = myfull \ F;
	err = norm(abs(u - exu));
	if err < errStop
		break;
	end
	u0 = u;
end

%p1 = p1func(u, i, j);%p_{i-1/2,j}
%p2 = p2func(u, i, j);%p_{i+1/2,j}
%p3 = p3func(u, i, j);%p_{i,j-1/2}
%p4 = p4func(u, i, j);%p_{i,j+1/2}


end
function vp1 = p1func(u, m, n)	
    global N;
    global h;
    global x;
    global y;
	if m == 1||n == N
        vpx = (u(m+(n-1)*N))^2;
		vpy1 = u(1+(n-1)*N) - 0.25*u(2+(n-1)*N);
		%vpx = (u(m+(n-1)*N))^2;
		%vpy1 = u(1+n*N) - 0.25*u(2+n*N);
	else
		vpx = (u(m+(n-1)*N) - u(m-1 + (n-1)*N))^2;
		vpy1 = u(m-1+n*N) + 0.25*(u(m+n*N)-u(m-2+n*N));
	end
	if n == 1 || m==1
		vpy2 = 0;
    elseif m==1 && n==2
		vpy2 = u(m+(n-1)*N) + 0.25*(u(m+(n-1)*N)-u(m+(n-1)*N));
    elseif m == 2 && n ==2
        vpy2 = u(m-1+(n-2)*N) + 0.25*(u(m+(n-2)*N)-u(m-1+(n-2)*N));
    else
        vpy2 = u(m-1+(n-2)*N) + 0.25*(u(m+(n-2)*N)-u(m-2+(n-2)*N));
	end
	vpy = 0.25*(vpy1 - vpy2)^2;
	vp1 = 1.0/(vpx + vpy);
end
function vp2 = p2func(u, m, n)
	global N;
	global h;
    global x;
    global y;
    if m < N
        vpx = (u(m+1+(n-1)*N) - u(m + (n-1)*N))^2;
    else
        vpx = (u(m+(n-1)*N)+1/h*y(n) - u(m + (n-1)*N))^2;
    end
	if m < N && n < N
		vpy1 = u(m+n*N) + 0.25*(u(m+1+n*N)-u(m-1+n*N));
    elseif m < N && n == N 
		vpy1 = u(m+(n-1)*N)+1/h*x(m) + 0.25*(u(m+1+(n-1)*N)+1/h*x(m) - u(m-1+(n-1)*N)-1/h*x(m));
    elseif m == N && n < N
        vpy1 = u(m+n*N) + 0.25*(u(m+n*N)+1/h*y(n) - u(m-1+n*N));
    else
        vpy1 = u(m+(n-1)*N)+1/h*x(m) + 0.25*(u(m+(n-1)*N)+1/h*y(n) + u(m-1+(n-1)*N));
    end
	if n == 1
		vpy2 = 0;
    elseif m == N
        vpy2 = u(m+(n-2)*N) + 0.25*(u(m+(n-2)*N)+1/h*y(n)-u(m-1+(n-2)*N));
    else
        if m ==1
            vpy2 = 0;
        else
            vpy2 = u(m+(n-2)*N) + 0.25*(u(m+1+(n-2)*N)-u(m-1+(n-2)*N));
        end
	end
	vpy = 0.25*(vpy1 - vpy2)^2;
	vp2 = 1.0/(vpx + vpy);
end
function vp3 = p3func(u, m, n)
	global N;
    global x;
    global y;
    global h;
	if n == 1
		vpx1 = u(m+1+(n-1)*N) - 0.25*u(m+1+(2-1)*N);
    elseif n == N
        vpx1 = u(m+(n-1)*N) - 0.25*(u(m + (n-1)*N)+1/h*x(m) +1/h*y(n)- u(m+1 + (n-2)*N));
    elseif m == N
        vpx1 = u(m+(n-1)*N)+1/h*y(n) - 0.25*(u(m + (n-1)*N)+1/h*x(m)+1/h*y(n) - u(m + (n-2)*N)+1/h*y(n));
    else
		vpx1 = u(m+1+(n-1)*N) - 0.25*(u(m+1 + n*N) - u(m+1 + (n-2)*N));
	end
	if m == 1||n == 1
		vpx2 = 0;
    elseif n == N
        vpx2 = u(m-1+(n-1)*N) - 0.25*(u(m-1 + (n-1)*N) + 1/h*x(m) - u(m-1 + (n-2)*N));
    else
		vpx2 = u(m-1+(n-1)*N) - 0.25*(u(m-1 + n*N) - u(m-1 + (n-2)*N));
	end
	
	vpx = 0.25*(vpx1 - vpx2)^2;
	if n == 1
		vpy = 0;
	else
		vpy = (u(m+(n-1)*N) - u(m+(n-2)*N))^2;
	end
	vp3 = 1.0/(vpx + vpy);
end
function vp4 = p4func(u, m, n)
    global h;
	global N;
    global x;
    global y;
	if n == 1
		vpx1 = u(m+1+(n-1)*N) + 0.25*(u(m+1+n*N));
		vpx2 = 0;
    elseif n == N || m == N
        vpx1 = u(m+(n-1)*N)+1/h*y(n) + 0.25*(u(m + (n-1)*N)+1/h*x(m) - u(m + (n-2)*N));
        vpx2 = u(m-1+(n-1)*N) + 0.25*(u(m-1 + (n-1)*N) - u(m-1 + (n-2)*N));
    else
		vpx1 = u(m+1+(n-1)*N) + 0.25*(u(m + n*N)+1/h*y(n) - u(m+1 + (n-2)*N));
        if m ==1 && n==1
            vpx2 = 0;
        else
            if m == 1 && n==2
                vpx2 = 0;
            else                
                vpx2 = u(m-1+(n-1)*N) + 0.25*(u(m-1 + n*N) - u(m-1 + (n-2)*N));
            end
        end
	end
	vpx = 0.25*(vpx1 - vpx2)^2;
    if n == N
        vpy = (u(m+(n-1)*N)+1/h*x(m) - u(m+(n-1)*N))^2;
    else
        vpy = (u(m+n*N) - u(m+(n-1)*N))^2;
    end
	vp4 = 1.0/(vpx + vpy);
end
function vboundary = boundaryfunc(m, n)
	if m == 1 || n == 1
		vboundary = 0;
	elseif m == N
		vboundary = 1;
	elseif n == N
		vboundary = 1;
	else 
		vboundary = 0;
	end
end
function v = RightF(x, y)
v = 2*x*y*(x^2+y^2)^(-3/2);
end

