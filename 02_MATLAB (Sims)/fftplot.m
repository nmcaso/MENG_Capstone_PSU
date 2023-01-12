function ff = fftplot(linear_spec_1d, varargin)

    plotoptions = varargin;
    linear_spec_1d = linear_spec_1d(:);

    if size(linear_spec_1d,2) ~= 1
        error("Must be a 1-Dimensional Array")
    end

    siglength = length(linear_spec_1d);
    switch mod(siglength,2)
        case 0 
            belowfs = linear_spec_1d(2:siglength/2);
            abovefs = linear_spec_1d(siglength/2+1:end);
        case 1
            belowfs = linear_spec_1d(2:(siglength+1)/2);
            abovefs = linear_spec_1d((siglength+1)/2+1:end);
    end

    if belowfs ~= conj(abovefs)
        linear_spec_1d = fftshift(linear_spec_1d);
    end

    u = zeros(size(linear_spec_1d));
    v = real(linear_spec_1d);
    w = imag(linear_spec_1d);

    x = (0:siglength-1)';
    y = zeros(size(u));
    z = y;    
   
    plotobj = plot3(x,v,w, plotoptions{:}); hold on;
    plot2 = quiver3(x,y,z,u,v,w,rms(abs(linear_spec_1d)), '-c','LineWidth',1.3);
    grid on;
    xlabel("\omega");
    ylabel("Real Axis");
    zlabel("Imaginary Axis")
    daspect([1 1 1])

    ff = gcf;

end