function forwardMatrix = TriangleInterp(obj,forwardMatrix,SensIndex)
    %TriangleInterp will take the nearest 3 points on the grid to
    %   interpolate the weight coefficients for the corresponding
    %   points and sensor times.

    switch max(obj.nnDistance)
        case obj.nnDistance(1) %bottom right triangle
            x(1) = obj.xNeighbors(1);
            x(2) = obj.xNeighbors(2);
            x(3) = obj.xNeighbors(2);
            y(1) = obj.yNeighbors(2);
            y(2) = obj.yNeighbors(1);
            y(3) = obj.yNeighbors(2);
            z(1) = obj.nnIndex(2);
            z(2) = obj.nnIndex(3);
            z(3) = obj.nnIndex(4);
            
            triCase = 1;
            coeffs = CoeffCalc(x,y,obj.x0,obj.y0,R,triCase);
        case obj.nnDistance(2) %top right triangle
            x(1) = obj.xNeighbors(1);
            x(2) = obj.xNeighbors(2);
            x(3) = obj.xNeighbors(2);
            y(1) = obj.yNeighbors(1);
            y(2) = obj.yNeighbors(1);
            y(3) = obj.yNeighbors(2);
            z(1) = obj.nnIndex(1);
            z(2) = obj.nnIndex(3);
            z(3) = obj.nnIndex(4);
            
            triCase = 2;
            coeffs = CoeffCalc(x,y,obj.x0,obj.y0,R, triCase);
        case obj.nnDistance(3) %bottom left triangle
            x(1) = obj.xNeighbors(1);
            x(2) = obj.xNeighbors(1);
            x(3) = obj.xNeighbors(2);
            y(1) = obj.yNeighbors(1);
            y(2) = obj.yNeighbors(2);
            y(3) = obj.yNeighbors(2);
            z(1) = obj.nnIndex(1);
            z(2) = obj.nnIndex(2);
            z(3) = obj.nnIndex(4);
            
            triCase = 3;
            coeffs = CoeffCalc(x,y,obj.x0,obj.y0,R,triCase);
        case obj.nnDistance(4) %top left triangle
            x(1) = obj.xNeighbors(1);
            x(2) = obj.xNeighbors(1);
            x(3) = obj.xNeighbors(2);
            y(1) = obj.yNeighbors(1);
            y(2) = obj.yNeighbors(2);
            y(3) = obj.yNeighbors(1);
            z(1) = obj.nnIndex(1);
            z(2) = obj.nnIndex(2);
            z(3) = obj.nnIndex(3);
            
            triCase = 4;
            coeffs = CoeffCalc(x,y,obj.x0,obj.y0,R,triCase);
    end

    forwardMatrix(SensIndex,z(1)) = coeffs(1);
    forwardMatrix(sensIndex,z(2)) = coeffs(2);
    forwardMatrix(sensIndex,z(3)) = coeffs(3);
end