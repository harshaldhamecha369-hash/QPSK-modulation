N = 1000;
fs = 1;
time = linspace(0,1,N);
x = sin(2*pi*fs*time);


bits_list = [4,6,8,10,12,14];
quantisation_error_rms = zeros(length(bits_list));
for k = 1:length(bits_list)
    bits = bits_list(k);
    levels = 2^bits;
    Vpi = 1;
    resolution = 2*Vpi/(levels - 1);
    quantisation = round(x/resolution)*resolution;
    quantisation_error = x - quantisation;
    quantisation_error_rms(k) = rms(quantisation_error);

    figure('Position', [100,100,600,600])
    plot(time,quantisation)
    hold on 
    plot(time,x)
    title(sprintf('DAC Output (%d bits)',bits))
    ylim([-1.6 1.6])
    yticks(-1.6:0.2:1.6) 
    grid("on")

    figure('Position', [100,100,600,600])
    plot(time,quantisation_error)
    title(sprintf('Quantization Error (%d bits)',bits))
    ylim([-1.6 1.6])
    yticks(-1.6:0.2:1.6)
    grid("on")

end

figure('position', [100 100 600 600])
plot(bits_list,quantisation_error_rms,'LineWidth',2)
grid on
xlabel('DAC Resolution (bits)')
ylabel('RMS Quantization Error')
title('Quantization Error vs DAC Resolution')