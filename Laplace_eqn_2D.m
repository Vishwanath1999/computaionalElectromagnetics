

%%
%Specifying parameters
lx      = 10;        %Length in x-direction
ly      = 10;        %Length in y-direction
nx      = 600;       %Number of steps in space(x)
ny      = 600;       %Number of steps in space(y)       
niter   = 1000;     %Number of iterations 
dx      = lx/(nx-1); %Width of space step(x)
dy      = ly/(ny-1); %Width of space step(y)
x       = 0:dx:lx;   %Range of x(0,10) specifying the grid pts
y       = 0:dy:ly;   %Range of y(0,10) specifying the grid pts

%%
%Initial Conditions
p       = zeros(ny,nx);      %Pre-allocating p
pn      = zeros(ny,nx);      %Pre-allocating pn

%Boundary conditions (Dirichlet conditions)     
p(:,1)   =10;
p(:,nx)  =0;
p(1,:)   =0;                   
p(ny,:)  =0; 



%Explicit iterative CD scheme in space (5-point difference)
i = 2:nx-1;
j = 2:ny-1;
for it = 1:niter
    pn = p;
    p(i,j) = ((dy^2*(pn(i+1,j)+pn(i-1,j)))+...
              (dx^2*(pn(i,j+1)+pn(i,j-1))))...
              /(2*(dx^2+dy^2));
end

%Plotting the solution
surf(x,y,p,'EdgeColor','none');       
shading interp
%shading faceted
title({'2-D Laplace''s equation';...
      ['{\itNumber of iterations} = ',num2str(it)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (P) \rightarrow')
view(0,90);

