clear variables
clc; clf; close all;

%% BUILD SPIRAL
spacing=.2;
Rs=5;                               %max radius
revs=Rs/spacing;                    %calc number of revolutions
k=spacing/(2*pi);                       
theta=linspace(0,revs*2*pi,1777); %how many points and angle coord
r=k*theta;                       %radius coord
x=r.*cos(theta);                 %convert to x cartesian
y=r.*sin(theta);                 %convert to y cartesian

%% build sample function
xo=.9;
yo=-.7;
Intensity=zeros(length(x),1);
for i=1:length(x)   
        R = sqrt((x(i)-xo).^2 + (y(i)-yo).^2) + eps;
        Intensity(i) = sin(R)./R;
end
%% plot surface of function
[X,Y] = meshgrid(-20:.1:20);
R = sqrt((X+xo).^2 + (Y+yo).^2) + eps;
Z = sin(R)./R;
%% Sampling Spiral and Spiral on Function
figure(1)
subplot(2,1,1)
plot(x,y,'r.') 
axis equal
grid on
axis off
title('Sampling Spiral')

subplot(2,1,2)
hold on
surf(X,Y,Z)
axis off
plot(x,y,'r.')
title('Function and Sampling Spiral')
shading interp
%% Nice Spiral Mapping for Presentaions
figure(2)
hold on
plot3(x,y,Intensity,'r.')
plot3(x,y,Intensity)
axis off
colormap('jet(200)')
title('Spiral Mapping')
xlabel('X Translation')
ylabel('Y Translation')
zlabel('Intensity')

%axis equal
shading interp

% 
% %% Nice Spiral Graphing for Presentaions
% theta=linspace(0,revs*2*pi,3000); %how many points
% r=k*theta;                   %spiral
% XX=r.*cos(theta);             %convert to x cartesian
% YY=r.*sin(theta);             %convert to y cartesian
% 
% figure(4)
% hold on
% plot(x,y,'r.') 
% plot(XX,YY) 
% axis equal
% grid on
% axis off
% %% Function with Spiral Overlay for Presentaions
% figure(5)
% hold on
% surf(X,Y,Z)
% axis off
% plot(x,y,'r.')
% 
% 
% shading interp

%% Map of Function from Spiral
figure(2) 
hold on
plot3(x,y,Intensity,'r.')   %dots
plot3(x,y,Intensity)        %lines

colormap('jet(200)')
shading interp

xlabel('X Translation')
ylabel('Y Translation')
zlabel('Intensity')
%axis off
%axis equal

%%
[pks, locs]=findpeaks(Intensity); %find peaks (values where dz/dxdy=0) and locations 
avg=mean(pks); 
for a=1:length(locs)
    if pks(a) >= avg            %check if peak value is over average to filter out false peaks
        Peaks(a,1)=x(locs(a)); %If true pass values into matrix Peaks
        Peaks(a,2)=y(locs(a));
        Peaks(a,3)=Intensity(locs(a));
    else
        Peaks(a,1)=nan;         %If false pass nan placeholders
        Peaks(a,2)=nan;
        Peaks(a,3)=nan;

    end
end
Peaks=rmmissing(Peaks); %remove nans

%% Peaks in different rev numbers

for a=1:length(pks)%plot points of global peaks on the 
    figure(2)
    plot3(x(locs(a)),y(locs(a)),pks(a),'b.');
    plot3(x(locs(a)),y(locs(a)),zeros(length(locs),1),'b+');
end
    
rv=sqrt(x.^2+y.^2)./(k*2*pi); %Calc distance to center for each point/pitch for rev number
revsV=floor(rv);    %round rev number down           
revsV=revsV+1; %rev1 is now 1 not zero

I_Revs=nan(length(x),revs+1);%preallocate for speed
for a=1:length(revsV)  %build matrixes where column number corresponds to rev number
    I_Revs(a,revsV(a))=Intensity(a);
    X_Revs(a,revsV(a))=x(a);
    Y_Revs(a,revsV(a))=y(a);
end


[M,I] = max(I_Revs);

hold on
for a=1:length(I)
    figure(2)
    plot3(x(I(a)),y(I(a)),M(a),'k+'); 
end

fit=[Peaks(:,1),Peaks(:,2)];
figure(10)

for a=1:width(X_Revs)
hold on
plot(X_Revs(:,a),Y_Revs(:,a));
end

%% linear fit
[Int,Pos]=max(Intensity);
a1=x(Pos);
b1=y(Pos);


M_Guess=-1;
Xpos=Peaks(:,1);
FUN=@(S,Xpos)S(1).*(Xpos+a1)-b1;

Slope=lsqcurvefit(FUN,M_Guess,Peaks(:,1),Peaks(:,2));
figure(2)
hold on
plot(Peaks(:,1),FUN(Slope,Peaks(:,1)))
plot(Peaks(:,1),Peaks(:,2),'b.')
plot(a1,b1,'r+')
%% find direction
Step1=.001;

x1n=(a1:-Step1:a1-(Step1*4));
x1p=(a1:Step1:a1+(Step1*4));

y1n=FUN(Slope,x1n);
y1p=FUN(Slope,x1p);

Int1n=GetIntensity(x1n,y1n);
Int1p=GetIntensity(x1p,y1p);

if mean(Int1p) > mean(Int1n) %if dI/dx is positive step=+
    Step1=Step1;
    Int1=(Int1p);
    X1=x1p(end)+Step1;
    Linear1M(:,1)=x1p;
    Linear1M(:,2)=y1p;
elseif mean(Int1p) < mean(Int1n)
    Step1=-Step1;
    Int1=(Int1n);
    X1=x1n(end)+Step1;
    XX1=[x1n X1]
else
    print('Error1')
end

%%
a=4;

while Int1(a)>=.99999*mean([Int1(a-1) Int1(a-2) Int1(a-3)])
a=a+1;
Y1=FUN(Slope,X1); %To do:make vectors of displacements maybe : maybe growing vector
Int1(a)=GetIntensity(X1,Y1);
X1=X1+Step1;


end
x1=Peaks(:,1);
y1=Slope*(x1+a1)-b1;
plot(x1,y1)
linear1=[x1,y1];