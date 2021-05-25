%{
Program that calculates the TM of a square array of single mode fibres that
propagate light into a square multimode fibre.
Requires gif.m file to generate gifs, along with beamprop.m and heatmap.m
to calculate and plot wavefronts.
This program only supports an even number of cores along each side of the
fibre array.
%}

clear;
clc;

%---
%Set optical fibre properties.

r = 1*10^(-6); %core radius
n_core = 1.4; %core refractive index
n_cladding = 1.39; %cladding refractive index
n_diff = n_core - n_cladding; %difference between refractive indices
n_delta = n_diff/n_core; %relative refractive number difference
core_number = 10; %number of cores along each side of the array
core_distance = 8*10^(-6); %distance between COMs of each core horizontally 
                         %and vertically
core_value = core_number/2 - 0.5; %a value calculated from the core number,
                                  %used in various parts of the program 
r_mmf = 3*10^(-6) * (core_number - 1); %MMF width

%---
%Set Gaussian beam properties.

wavelength = 500*10^(-9); %light wavelength
sigma = r/sqrt(2); %Gaussian beam standard deviation
wavelength_eff = wavelength/n_cladding; %effective wavelength for angular
                                        %spectrum calculations in cladding
kmag = 2*pi/wavelength; %beam wave vector magnitude
kmag_eff = 2*pi/wavelength_eff; %beam wave vector effective magnitude

%---
%Set initial wavefront array, coordinates of points in both real and
%Fourier space, and create arrays required for calculation of phase mask.

L = 40*r; %x & y side length of sample space in metres
N = 256; %number of samples taken in each axis minus 1

while (core_number - 1) * core_distance > L - core_distance
    L = L + 40*r;
    N = N + 256;
end
%extends the frame to ensure it can accomodate all cores with sufficient
%zero padding

ds = L/N; %sample interval in metres
x = -L/2:ds:L/2; %array of x coordinates
y = x; %array of y coordinates
[X,Y] = meshgrid(x,y);

kx = -pi/ds:2*pi/(N*ds):pi/ds;
ky = kx;
[KX,KY] = meshgrid(kx,ky);

kz = sqrt(kmag_eff^2 - KX.^2 - KY.^2); %array of z components of k for
                                       %plane waves when FT is taken
                                   
dz = 2.5*wavelength; %propagation distance for each interval

index = 0;
index_array = zeros(1,core_number^2);

for jx = -core_value:-0.5
    for jy = -core_value:jx
        
        index = index + 1;
        
        r = 1*10^(-6); %core radius
        core_distance = 8*10^(-6); %distance of outer cores COM from origin
        
        wavefront = exp(-(((X + jx*core_distance).^2 + ...
            (Y + jy*core_distance).^2) / (2*sigma^2)));
        
        n = zeros(N+1);
        
        for core_x = -core_value:core_value
            for core_y = -core_value:core_value
                n_temp = n_core*sqrt(1-2*(n_delta*(sqrt((X + core_x*...
                    core_distance).^2 + (Y + core_y*core_distance).^2)/...
                    r).^8)) - n_cladding;
                n_temp((X + core_x*core_distance).^2 + ...
                    (Y + core_y*core_distance).^2 >= r^2) = 0;

                n = n + n_temp;
            end
        end
        
        %---
        %Calculate propagation across a number of wavelengths using angular 
        %spectrum method, and applies a further phase mask based on if
        %light propagated through core or cladding. This is repeated j 
        %times.

        frame = 0;
        
        for m = 1:50
            disp(['Frame ',num2str(frame+1),' (',num2str(index),')'])
            frame = frame + 1;
            
            wavefront = beamprop(wavefront,kz,dz,n,kmag); %performs BPM on 
                                                          %wf
        end
        %loop that allows the light to settle into the fibre
        
        bend = pi/180; %bend of fibre in radians
        dx = bend*dz; %shift of fibre COM in each frame due to bend
        steps = 185; %no of steps the cores are brought together
        dr = (1/2)*(r/steps); %reduction in core radius each step
        
        for m = 1:steps
            core_distance = core_distance - dx;
            r = r - dr;
            
            disp(['Frame ',num2str(frame+1),' (',num2str(index),')'])
            frame = frame + 1;
            
            n = zeros(N+1);
            
            for core_x = -core_value:core_value
                for core_y = -core_value:core_value
                    n_temp = n_core*sqrt(1-2*(n_delta*(sqrt((X + core_x*...
                        core_distance).^2 + (Y + core_y*...
                        core_distance).^2)/r).^8)) - n_cladding;
                    n_temp((X + core_x*core_distance).^2 + ...
                        (Y + core_y*core_distance).^2 >= r^2) = 0;

                    n = n + n_temp;
                end
            end
            
            wavefront = beamprop(wavefront,kz,dz,n,kmag); %performs BPM on 
                                                          %wf   
        end
        %loop that brings the cores together
        
        %STEP INDEX
        n = zeros(N+1) + n_core - n_cladding;
        n(abs(X) > r_mmf) = 0;
        n(abs(Y) > r_mmf) = 0;
        
        %GRADE INDEX
        %{
        r_mmf = r_mmf * sqrt(2); %MMF diagonal distance
        
        n = zeros(N+1);
        n = n_core*sqrt(1-2*(n_delta*((abs(X)+abs(Y))/r_mmf).^8)) - ...
            n_cladding;
        n(abs(X) + abs(Y) > r_mmf) = 0;
        n = imrotate(n,45,'bicubic','crop');
        %}
        %plots the phase mask of the square MMF
        
        for m = 1:500
            disp(['Frame ',num2str(frame+1),' (',num2str(index),')'])
            frame = frame + 1;
            
            wavefront = beamprop(wavefront,kz,dz,n,kmag); %performs BPM on 
                                                          %wf
        end
        %loop that propagates the light through the MMF
        
        tm(:,index) = reshape(wavefront,[],1);
        
        if jx ~= jy
            index_array(index) = 1; %indexes if a TM column corresponds to 
                                    %a fibre that does not sit on a
                                    %diagonal of the core array
        end
    end
end

temp_index = 0;

%---
%Go through every measured fibre output not on a diagonal, and transpose
%its field to get the TM of one quarter of the fibre array.

for m = 1:index
    input = zeros(index,1);
    
    if index_array(m) == 1
        input(m) = 1;
        temp_index = temp_index + 1;
        output = reshape(tm*input,N+1,N+1);
        output = output';
        tm_temp(:,temp_index) = reshape(output,[],1);
    end
end

tm = horzcat(tm,tm_temp);
clear tm_temp;

%---
%Obtain the rest of the TM by rotating the known output fields of one
%quarter of the fibre array.
core_quarter = (core_number/2)^2;

input = zeros(core_quarter,1);

for m = 1:core_quarter
    input = zeros(core_quarter,1);
    input(m) = 1;
    output = reshape(tm*input,N+1,N+1);
    tm_temp_1(:,m) = reshape(flipud(output),[],1);
    tm_temp_2(:,m) = reshape(rot90(output,2),[],1);
    tm_temp_3(:,m) = reshape(fliplr(output),[],1);
end

tm = horzcat(tm,tm_temp_1);
clear tm_temp_1;
tm = horzcat(tm,tm_temp_2);
clear tm_temp_2;
tm = horzcat(tm,tm_temp_3);
clear tm_temp_3;

%---
%Trim the TM to reduce the size of the saved file.

for j = 1:core_number^2
    input = zeros(core_number^2,1);
    input(j) = 1;
    
    output = tm*input;
    output = reshape(output,N+1,N+1);
    output_new = output((N/4)+1:(3*N/4)+1,(N/4)+1:(3*N/4)+1);
    
    tm_new(:,j) = reshape(output_new,[],1);
end

clear tm;
tm = tm_new;
clear tm_new;

save('tm10step2.mat','tm');