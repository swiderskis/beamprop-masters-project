%{
Program that calculates the required inputs in order to produce a scanning
beam at the output facet of a photonic lantern.
%}

clear;
clc;

%---
%Set optical fibre & Gaussian beam properties.

load tm.mat

r = 1*10^(-6); %core radius
sigma = r/sqrt(2); %Gaussian beam standard deviation
core_number = sqrt(size(tm(1,:),2));
r_mmf = 3*10^(-6) * (core_number - 1); %MMF width
r_mmf = r_mmf * sqrt(2); %MMF diagonal distance
wavelength = 500*10^(-9); %light wavelength

%---
%Set size and resolution of target frame.

N = sqrt(size(tm(:,1),1)) - 1; %number of samples taken in each axis -1
L = 10*r*(N/64); %x & y side length of sample space in metres
ds = L/N; %sample interval in metres
x = -L/2:ds:L/2; %array of x coordinates
y = x; %array of y coordinates
[X,Y] = meshgrid(x,y);

a = 0;
input = zeros(core_number^2,1225);
input(1,:) = 1;
step_size = 10;
scan = zeros(385);

for jx = 3:step_size:N+1-2
    for jy = 3:step_size:N+1-2
        if abs(X(jx,jy)) < r_mmf/sqrt(2) && abs(Y(jx,jy)) < r_mmf/sqrt(2)
            a = a + 1;
            disp(['Frame ',num2str(a)]);
            
            for j = 2:core_number^2
                for m = 1:4
                    input(j,a) = exp(1i * (m-1) * pi/2);
                    output = reshape(tm*input(:,a),N+1,N+1);
                    intensity(m) = abs(output(jx,jy))^2;
                    intensity_ft = fft(intensity);
                    angles = angle(intensity_ft);
                end

                input(j,a) = exp(1i * angles(4));
            end
        end
    end
end

save('input.mat','input');