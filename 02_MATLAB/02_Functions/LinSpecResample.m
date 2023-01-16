function xn_down = LinSpecResample(xn_0, numptsnew)
N = length(xn_0);

switch mod(N,2)
    case 0
        XN  = fft(xn_0(:));
    if numptsnew < N
        XN  = XN([1:N-numptsnew/2 N/2+1+numptsnew/2:end]);
    else 
        XN  = [XN(1:N/2-1); zeros(N-numptsnew,1); XN(N/2+1:end)];
    end
    case 1
end

xn_down = ifft(XN);

end

