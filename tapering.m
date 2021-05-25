%{
Program which calculates the phase masks for the tapering region of a
photonic lantern.
%}

%---
%Set up photonic lantern properties.

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

dz = 2.5*wavelength; %propagation distance for each interval
taper = pi/(180*core_value); %bend of fibre in radians
dx = taper*dz; %shift of fibre COM in each frame due to bend
steps = round(185*core_value); %no of steps the cores are brought 
                               %together
dr = (9/10)*(r/steps); %reduction in core radius each step
bend = 2*10^(-2);

%---
%Calculate phase masks.

frame = 0;

for m = 1:steps
    core_distance = core_distance - dx;
    r = r - dr;

    disp(['Frame ',num2str(frame+1)])
    frame = frame + 1;

    if m == 1
        n_taper = zeros(N+1,N+1,steps);
    end

    for core_x = -core_value:core_value
        for core_y = -core_value:core_value
            n_temp = n_core*sqrt(1-2*(n_delta*(sqrt((X + ...
                core_x*core_distance).^2 + (Y + core_y*...
                core_distance).^2)/r).^8)) - n_cladding;
            n_temp((X + core_x*core_distance).^2 + ...
                (Y + core_y*core_distance).^2 >= r^2) = 0;

            n_taper(:,:,m) = n_taper(:,:,m) + n_temp;
        end
    end
    
    n_taper(:,:,m) = n_taper(:,:,m) + n_cladding;
    %n_taper(:,:,m) = n_taper(:,:,m).*(1+Y/bend) - n_cladding;
end

save('n_taper.mat','n_taper');