clear variables
clc; clf; close all;


%% BUILD SPIRAL
spacing=.5;
Rs=5;
revs=Rs/spacing;
k=spacing/(2*pi);                       %scale
theta=linspace(0,revs*2*pi,1777); %how many points
r=k*theta;                   %spiral
x=r.*cos(theta);             %convert to x cartesian
y=r.*sin(theta);             %convert to y cartesian

%% 
xo=.9;
yo=-.7;
Intensity=zeros(length(x),1);
for i=1:length(x)   


        R = sqrt((x(i)-xo).^2 + (y(i)-yo).^2) + eps;
        Intensity(i) = sin(R)./R;
        
  
end
%%
[X,Y] = meshgrid(-20:.1:20);
R = sqrt((X+xo).^2 + (Y+yo).^2) + eps;
Z = sin(R)./R;
%%
figure(1)
subplot(2,1,1)
plot(x,y,'r.') 
axis equal
grid on
axis off

subplot(2,1,2)
hold on
surf(X,Y,Z)
axis off
plot(x,y,'r.')

shading interp
%%
figure(2)
hold on
plot3(x,y,Intensity,'r.')
plot3(x,y,Intensity)
axis off
colormap('jet(200)')
xlabel('X Translation')
ylabel('Y Translation')
zlabel('Intensity')

%axis equal
shading interp

%figure(3)
%surf(x,y,Intensity)
%%
theta=linspace(0,revs*2*pi,3000); %how many points
r=k*theta;                   %spiral
XX=r.*cos(theta);             %convert to x cartesian
YY=r.*sin(theta);             %convert to y cartesian

figure(4)
hold on
plot(x,y,'r.') 
plot(XX,YY) 
axis equal
grid on
axis off
%%
figure(5)
hold on
surf(X,Y,Z)
axis off
plot(x,y,'r.')


shading interp

%%
figure(2)
hold on
plot3(x,y,Intensity,'r.')
plot3(x,y,Intensity)
%axis off
colormap('jet(200)')
xlabel('X Translation')
ylabel('Y Translation')
zlabel('Intensity')

%axis equal
shading interp

[pks, locs]=findpeaks(Intensity)
for a=1:length(pks)
    figure(2)
    plot3(x(locs(a)),y(locs(a)),pks(a),'b.');
    plot3(x(locs(a)),y(locs(a)),zeros(length(locs),1),'b+');
    
end
    
[M,I] = max(Intensity);
a=x(I);
b=y(I);


rv=sqrt(x.^2+y.^2)./(k*2*pi);
revsV=floor(rv);
revsV=revsV+1;
a=0;
I_Revs=nan(length(x),revs+1);
for a=1:length(revsV)   
    I_Revs(a,revsV(a))=Intensity(a);
    X_Revs(a,revsV(a))=x(a);
    Y_Revs(a,revsV(a))=y(a);
end
[M,I] = max(I_Revs);

I=I';
a=0;
for a=1:length(I)
    figure(2)
    plot3(x(I(a)),y(I(a)),M(a),'k+');
    
end
figure(10)
a=0;
for a=1:width(X_Revs)
    hold on
plot(X_Revs(:,a),Y_Revs(:,a));
end