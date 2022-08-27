% Lid Driven Cavity by Stream Function Vorticity Method
clear all; 
clc;
Nx=129;
Ny=129;
Re=100;
visc=0.01;
dx=1.0/(Nx-1);
dy=1.0/(Ny-1);
x=zeros(Nx,Ny); y=zeros(Nx,Ny);
for i=1:Nx
for j=1:Ny
x(i,j)=dx*(i-1);
y(i,j)=dy*(j-1);
end
end
L=x(Nx,1);
H=y(1,Ny);
end
U=(Re*visc)/L;
dt=min(((visc*dx)/(L*(U^2))),((visc*dy)/(H*(U^2))));
%Max Iterations
time=0.0;
MaxStep=5000;
MaxIt=500;
%SOR
beta=1.5;
MaxErr=0.001;
sf=zeros(Nx,Ny); vt=zeros(Nx,Ny); vto=zeros(Nx,Ny);
for istep=1:MaxStep
for iter=1:MaxIt
vto=sf;
for i=2:Nx-1
for j=2:Ny-1
sf(i,j)=0.25*beta*(sf(i+1,j)+sf(i-1,j)+sf(i,j+1)+sf(i,j-1)+h*h*vt(i,j))+(1.0-beta)*sf(i,j);
end
end
Err=0.0;
for i=1:Nx
for j=1:Ny 
% check error
Err=Err+abs(vto(i,j)-sf(i,j)) ;
end
end
if Err <= MaxErr
break,

end
end
vt(2:Nx-1,1)=-2.0*sf(2:Nx-1,2)/(h*h); % vorticity on bottom wall
vt(2:Nx-1,Ny)=-2.0*sf(2:Nx-1,Ny-1)/(h*h)-2.0*U/h; % top wall
vt(1,2:Ny-1)=-2.0*sf(2,2:Ny-1)/(h*h); % right wall
vt(Nx,2:Ny-1)=-2.0*sf(Nx-1,2:Ny-1)/(h*h); % left wall
vto=vt;
for i=2:Nx-1
for j=2:Ny-1
vt(i,j)=vt(i,j)+dt*(-0.25*((sf(i,j+1)-sf(i,j-1))*...
(vto(i+1,j)-vto(i-1,j))-(sf(i+1,j)-sf(i-1,j))*...
(vto(i,j+1)-vto(i,j-1)))/(h*h)...
+(1/Re)*(vto(i+1,j)+vto(i-1,j)+vto(i,j+1)+...
vto(i,j-1)-4.0*vto(i,j))/(h^2) );
end
end
time=time+dt;
pause(0.01)
end
%Calculating the velocities
u=zeros(Nx,Ny);
v=zeros(Nx,Ny); 
for i=2:Nx-1
for j=2:Ny-1
u(i,j)=(sf(i,j+1)-sf(i,j-1))/(2*h);
v(i,j)=(sf(i+1,j)-sf(i-1,j))/(2*h);
u(:,Ny)=U;
end
end
u;
v;
%Pressure
p=zeros(Nx,Ny);
rho=1;
for it=1:200
    % Pressure Boundary Conditions
        p(1,:)=p(2,:);
        p(Nx,:)=p(Nx-1,:);
        p(:,1)=p(:,2);
        p(:,Ny)=p(:,Ny-1);
 
    for j=2:Ny-1
        for i=2:Nx-1

           p(i,j)=0.25*((p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))+((v(i+1,j)-v(i,j))^2 ...
                     +(u(i,j+1)-u(i,j))^2+(2*(u(i+1,j)-u(i,j))*(v(i,j+1)-v(i,j)))));
        end
    end
end


%Post-Processing
figure()
contourf(x,y,p);
axis('square');
title('Pressure');
figure()
contourf(x,y,sf);
axis('square');
title('Stream Function');
figure()
contourf(x,y,vt);
axis('square');
title('Vorticity');
figure()
contourf(x,y,u),axis('square');title('x-component of velocity');
figure()
contourf(x,y,v),axis('square');title('y-component of velocity');
