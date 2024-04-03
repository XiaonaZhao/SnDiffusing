% Suraj Shankar (2024). Diffusion in 1D and 2D
% (https://www.mathworks.com/matlabcentral/fileexchange/38088-diffusion-in-1d-and-2d),
% MATLAB Central File Exchange. Retrieved April 3, 2024. Simulating the 2-D

% Simulating the 2-D Diffusion equation by the Finite Difference
...Method 
% Numerical scheme used is a first order upwind in time and a second
...order central difference in space (Implicit and Explicit)

Mask = imread('H:\MATLAB\Mask3.tif');

%%
%Specifying parameters
nx=size(Mask,1);                           %Number of steps in space(x)
ny=size(Mask,2);                           %Number of steps in space(y)       
nt=500;                           %Number of time steps, frames 
dt=1;                         %Width of each time step
wid = 200; % s3, 200; s1, 32;
dx=wid/(nx-1);                     %Width of space step(x)
dy=wid/(ny-1);                     %Width of space step(y)
x=0:dx:wid;                        %Range of x(0,2) and specifying the grid points
y=0:dy:wid;                        %Range of y(0,2) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
vis=1.24;                         %Diffusion coefficient/viscocity
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)
%%
%Initial Conditions
% for i=1:nx
%     for j=1:ny
%         if (((0.4*wid)<=y(j))&&(y(j)<=(0.6*wid))&&((0.4*wid)<=x(i))&&(x(i)<=(0.6*wid)))
%             u(i,j)=2;
%         else
%             u(i,j)=0;
%         end
%     end
% end



u = double(Mask)'/255;

% for i=1:nx
%     for j=1:ny
%         if (((0.42*wid)<=y(j))&&(y(j)<=(0.58*wid))&&((0.42*wid)<=x(i))&&(x(i)<=(0.58*wid)))
%             u(i,j)=0;
%         end
%     end
% end
u0 = u;

%%
%B.C vector
bc=zeros(nx-2,ny-2);
bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs
%bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.Cs
%bc(:,1)=-UnS/dy; bc(:,nx-2)=UnN/dy;  %Neumann B.Cs
%B.Cs at the corners:
bc(1,1)=UW/dx^2+US/dy^2; bc(nx-2,1)=UE/dx^2+US/dy^2;
bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
bc=vis*dt*bc;
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
%Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;
%%
%Calculating the field variable for each time step
saveRoute = 'H:\MATLAB\Result_SnSe2\simulation3';
i=2:nx-1;
j=2:ny-1;
ucell = cell(nt, 1);
speed = 1/200;
for it=0:nt
    % it = 0;
    un=u;
    z = (u-u0)';
    z(z<0) = 0;
    ucell{it+1, 1} = z;
    % h=surf(x,y,z,'EdgeColor','none');       %plotting the field variable
    % shading interp
    % axis ([0 wid 0 wid 0 1])
    % title({['\rm2-D Diffusion with \itD\rm = ',num2str(vis/1e11),' m^2s^-^1'];...
    %     ['\rmtime (\itt\rm) = ',num2str(it*0.02),' s']})
    % xlabel('Spatial co-ordinate (x) pixel \rightarrow')
    % ylabel('{\leftarrow} Spatial co-ordinate (y) pixel')
    % zlabel('Optical intensity of mass transport (z) \rightarrow')
    % drawnow; 
    % refreshdata(h);
    h = imshow(z, 'DisplayRange',[], 'InitialMagnification', 'fit');
    % set(gca, 'CLim', [0 0.4]);
    colormap fire % parula
    hb = colorbar;
    set(get(hb,'title'),'string','Component');
    scalebar;
    axis on
    xlabel('Spatial co-ordinate (x) pixel \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y) pixel')
    title({['\rm2-D Diffusion with \itD\rm = ',num2str(vis/1e9),' m^2s^-^1'];...
        ['\rmtime (\itt\rm) = ',num2str(it*dt*speed, '%.3f'),' s']})
    % pause(0.1)
    drawnow; 
    refreshdata(h);
    
    figPath = [saveRoute '\S3_' num2str(it, '%03d')];
    saveas(h, figPath, 'jpg')

    pause(0.1)
    %Uncomment as necessary
    %Implicit method:
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    U=D\U;
    U=reshape(U,nx-2,ny-2);
    u(2:nx-1,2:ny-1)=U;
    %Boundary conditions
    %Dirichlet:
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    %Neumann:
    %u(1,:)=u(2,:)-UnW*dx;
    %u(nx,:)=u(nx-1,:)+UnE*dx;
    %u(:,1)=u(:,2)-UnS*dy;
    %u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
    %Explicit method:
    %{
    u(i,j)=un(i,j)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(dy*dy));
    %Boundary conditions
    %Dirichlet:
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    %Neumann:
    %u(1,:)=u(2,:)-UnW*dx;
    %u(nx,:)=u(nx-1,:)+UnE*dx;
    %u(:,1)=u(:,2)-UnS*dy;
    %u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
end