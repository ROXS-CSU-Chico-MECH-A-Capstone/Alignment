clear variables
clc; clf; close all;


%% BUILD SPIRAL
spacing=10;
Rs=20;
revs=Rs/spacing
k=spacing/(2*pi);                       %scale
theta=linspace(0,revs*2*pi,300); %how many points
r=k*theta;                   %spiral
x=r.*cos(theta);             %convert to x cartesian
y=r.*sin(theta);             %convert to y cartesian

%% 
xo=5.6;
yo=-7.9;
Intensity=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)


        R = 10*sqrt((x(i)+xo).^2 + (y(j)+yo).^2) + eps;
        Intensity(i,j) = sin(R)./R;
        
  
    end
end
%%
[X,Y] = meshgrid(-20:.1:20);
R = 10*sqrt((X+xo).^2 + (Y+yo).^2) + eps;
Z = sin(R)./R;
%%
figure(1)
subplot(2,1,1)
plot(x,y,'r.') 
axis equal
grid on


subplot(2,1,2)
surf(X,Y,Z)
shading interp

figure(2)
surf(x,y,Intensity)
xlabel('X Translation')
ylabel('Y Translation')
zlabel('Intensity')
axis equal
