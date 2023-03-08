function z=GetIntensity(x,y,x0,y0) %get intensity on coords x and y with resective offsets xo and yo
R = sqrt((x-x0).^2 + (y-y0).^2) + eps;
z = sin(R)./R;
end