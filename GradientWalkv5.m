clear variables
clc; clf; close all;


cycles=1; %cyles of Spiral searching
%spiral offsets
x0=0; 
y0=0;

c=jet(5); %color code for later
%%
for q=1:cycles
points=2000*q;
spacing=.2^q;
Radius=5/q^q;                               %max radius
[S,revs,k]=Build_Spiral(spacing,Radius,points,x0,y0);
x=S(:,1);
y=S(:,2);

%% build sample function
xo=3.24; %x offset for intensity function 
yo=-3.5; %y offset for intensity function 
Intensity=GetIntensity(x,y,xo,yo);

%% plot surface of function 
[X,Y] = meshgrid(-20:.1:20);
Z = GetIntensity(X,Y,xo,yo);
%% Sampling Spiral and Spiral on Function
%GUI1(x,y,X,Y,Z)
%% Nice Spiral Mapping for Presentaions
%GUI2(x,y,Intensity)
%GUI3(revs,k,x,y,X,Y,Z) %nice figs for presentations


%% Finding which peaks to run linear fit on
[pks, locs]=findpeaks(Intensity); %find peaks (values where dz/dxdy=0) and locations 
avg=mean(pks);  %average of magnitude of all detected peaks
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


%% linear fit
[Int,Pos]=max(Intensity); %find maximun detected point of intensity
a1=x(Pos);      %x coord of max
b1=y(Pos);      %y coord of max


M_Guess=-1; %initial guess of slope
Xpos=Peaks(:,1); %x positions of peak values
FUN=@(S,Xpos)S(1).*(Xpos+a1)-b1; %defining slope function of peak values
Slope=lsqcurvefit(FUN,M_Guess,Peaks(:,1),Peaks(:,2)); %find slope thru R^2 fit

%figure(2)
hold on
%plot(Peaks(:,1),FUN(Slope,Peaks(:,1)))
%plot(Peaks(:,1),Peaks(:,2),'b.')
%plot(a1,b1,'r+')
%% find direction
Step1=.001; %linear step size

x1n=(a1:-Step1:a1-(Step1*4)); %look in negative x direction
x1p=(a1:Step1:a1+(Step1*4));  %look in positive x direction

y1n=FUN(Slope,x1n); %find -x y values
y1p=FUN(Slope,x1p); %find +x y values

Int1n=GetIntensity(x1n,y1n,xo,yo); %find intensities on linear movements
Int1p=GetIntensity(x1p,y1p,xo,yo);
Linear1M=[]; %clear previous run

if mean(Int1p) > mean(Int1n) %if dI/dx is positive step=+
    Step1=Step1; %define positive step size
    Int1=(Int1p); %use pos x Int data
    X1=x1p(end)+Step1;
    Linear1M(:,1)=x1p;
    Linear1M(:,2)=y1p;
    Linear1M(:,3)=(Int1p);
elseif mean(Int1p) < mean(Int1n)
    Step1=-Step1;
    Int1=(Int1n);
    X1=x1n(end)+Step1;
    Linear1M(:,1)=x1n;
    Linear1M(:,2)=y1n;
    Linear1M(:,3)=(Int1n);
else
    print('Error1')
end

%%
a=length(x1n);

%while Linear1M(a,3)>=1*mean([Linear1M(a-1,3) Linear1M(a-2,3) Linear1M(a-3,3)]) %continue if current value is above threshold
while Linear1M(a,3)>=Linear1M(a-1,3) %continue if current value is above threshold
a=a+1;
Linear1M(a,1)=Linear1M(a-1,1)+Step1; %create new x value
Linear1M(a,2)=FUN(Slope,Linear1M(a,1)); %create new y value
Linear1M(a,3)=GetIntensity(Linear1M(a,1),Linear1M(a,2),xo,yo);%create new I value
end

x1=Peaks(:,1); %call peak x coords
y1=FUN(Slope,x1); %solve for y coords from function
%plot(x1,y1,'g-') %plot linear function
%%


[pksL1, locsL1]=findpeaks(Linear1M(:,3)); %Find peak value in intensity on linear movement
Spiral2=[Linear1M(locsL1,1),Linear1M(locsL1,2)]; %coord of peak on linear move
x0=Linear1M(locsL1,1);
y0=Linear1M(locsL1,2);
I0=Linear1M(locsL1,3);

figure(11) %Plot searching algorythym
hold on
axis equal
title('Searching Algorythm')
%GUI2(x,y,Intensity)

LColor=c(q, :); %get linear color 

%plot3(x,y,Intensity,'r-') %Plot Spiral
plot3(x,y,Intensity,'-','color',LColor,'DisplayName','Spiral') %Plot Spiral
plot3(a1,b1,Int,'*','color',LColor,'LineWidth',1,'DisplayName','Spiral Peak') %plot spiral peak

plot3(Linear1M(1,1),Linear1M(1,2),Linear1M(1,3),'o','color',LColor,'LineWidth',2,'DisplayName','Linear Start') %Linear Start
plot3(Linear1M(:,1),Linear1M(:,2),Linear1M(:,3),'color',LColor,'LineWidth',1,'DisplayName','Linear Move') %Plot linear move
plot3(x0,y0,I0,'+','color',LColor,'LineWidth',2,'DisplayName','Linear Peak') %plot linear peak/Next Spiral Start

%plot3([x1n x1p],[y1n y1p],[Int1n Int1p],'k--','Linewidth',2,'DisplayName','Linear Fit')
plot3(Peaks(:,1),Peaks(:,2),Peaks(:,3),'x','color',LColor,'Linewidth',1,'DisplayName','Linear Fit')
%L=GUITestLegend(cycles)
%legend(L)
%pause(1)

%pause(1)
end
%L=GUITestLegend(cycles)
legend()
disp('Success!')

