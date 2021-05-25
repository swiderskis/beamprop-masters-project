%{
Program that scans an image using a 10x10 photonic lantern.
Requires heatmap.m.
%}

clear;
clc;

%---
%Set optical fibre & Gaussian beam properties.

load tm10long.mat
load inputlong3.mat

image = imread('penguin2cropped.png');
image = imresize(image,[385 385]);
image = rgb2gray(image);
image = flipud(image);
image = cast(image,'double');
image = image/256;
image(1:24,:) = 0;
image(:,1:24) = 0;
image(362:385,:) = 0;
image(:,362:385) = 0;

r = 1*10^(-6); %core radius
sigma = r/sqrt(2); %Gaussian beam standard deviation
core_number = sqrt(size(tm(1,:),2));
core_value = core_number/2 - 0.5;
r_mmf = 3*10^(-6) * (core_number - 1); %MMF width
wavelength = 500*10^(-9); %light wavelength

%---
%Set size and resolution of target frame.

N = sqrt(size(tm(:,1),1)) - 1; %number of samples taken in each axis -1
L = 10*r*(N/64); %x & y side length of sample space in metres
ds = L/N; %sample interval in metres
x = -L/2:ds:L/2; %array of x coordinates
y = x; %array of y coordinates
[X,Y] = meshgrid(x,y);

%{
%---
%Apply a random phase offset to each core, ranging from -offset to offset.
max_offset = 8*pi/8;
for m = 1:100
    rand_offset = -max_offset + 2*max_offset*rand;
    input(m,:) = input(m,:) * exp(1i*rand_offset);
end
%}
%input = guidestar(tm,input,N);

a = 0;
step_size = 10;
scan = zeros(N+1);

%---
%Apply the scanning beam input (calculated from scanning_beam.m) to the
%loaded TM input to give a scanning beam output.
for jx = 3:step_size:N+1-2
    for jy = 3:step_size:N+1-2
        if abs(X(jx,jy)) < r_mmf && abs(Y(jx,jy)) < r_mmf
            a = a + 1;
            %disp(['Frame ',num2str(a)]);

            output = reshape(tm*input(:,a),N+1,N+1);

            scan(jx-(step_size/2):jx+(step_size/2), ...
                jy-(step_size/2):jy+(step_size/2)) = ...
                sum(image.*abs(output).^2,'all');

            %heatmap(x,y,abs(output).^2);
            %colormap('hot');
            %frame = getframe(gcf);
            %writeVideo(v,frame);

        end
    end
end

scan = scan/max(scan,[],'all');

%---
%Plot the scanned image and find the 2D correlation coefficient.
heatmap(x,y,scan);
xlim([-2.7,2.7]*10^(-5));
ylim([-2.7,2.7]*10^(-5));
set(gca, 'Visible', 'off')

hold on
plot([-2.4; -1.4]*10^(-5), [-2; -2]*10^(-5), 'r-', 'LineWidth', 5)
hold off
text(-1.9*10^(-5), -1.8*10^(-5), '$$ 10\mu m $$', 'HorizontalAlignment',...
    'center','Interpreter','latex','color','red','FontSize',24)

image_cropped = image(25:361,25:361);
scan_cropped = scan(25:361,25:361);
R = corr2(image_cropped,scan_cropped)
saveas(gcf,'penguin2scan.png');