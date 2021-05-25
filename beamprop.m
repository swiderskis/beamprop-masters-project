%Performs the beam propagation method onto a given wavefront.
function wavefront = beamprop(wavefront,kz,dz,n,kmag)
    %Obtains FT of wavefront
    wavefront = ifftshift(wavefront);
    wavefrontft = fft2(wavefront);
    wavefrontft = fftshift(wavefrontft);

    %Applies phase shift to wavefront FT, calculates inverse FT
    wavefrontft = wavefrontft .* exp(1i * kz * dz);
    wavefrontft = ifftshift(wavefrontft);
    wavefront = ifft2(wavefrontft);
    wavefront = fftshift(wavefront);
    
    wavefront = wavefront .* exp(1i .* n .* kmag .* dz); %applies phase
                                                         %mask
                                                         
end