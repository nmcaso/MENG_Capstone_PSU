classdef Neighbors
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0;
        y0;
        xNeighbors;
        yNeighbors;
        nnDistance;
        nnIndex;
        farthest;
    end
    
    methods
        function obj = Neighbors(R,theta,image,sensor)
            %Neighbors Construct an instance of this class
            %   Find the nearest neighbors to the interpolated point along
            %   the arc at radius R and angle theta with respect to the
            %   sensor location.
            obj.x0 = sensor(1);
            obj.y0 = sensor(2);
            xT = R*cos(theta)+obj.x0;
            yT = R*sin(theta)+obj.y0;
            obj.xNeighbors = image.xArr(abs(image.xArr-xT)<image.res);
            obj.yNeighbors = image.yArr(abs(image.yArr-yT)<image.res);
            if length(obj.xNeighbors)==length(obj.yNeighbors) && length(obj.xNeighbors)==2
                obj.nnDistance(1) = hypot((xNeighbors(1)-xT),(yNeighbors(1)-yT)); %top-left
                obj.nnDistance(2) = hypot((xNeighbors(1)-xT),(yNeighbors(2)-yT)); %bottom-left
                obj.nnDistance(3) = hypot((xNeighbors(2)-xT),(yNeighbors(1)-yT)); %top-right
                obj.nnDistance(4) = hypot((xNeighbors(2)-xT),(yNeighbors(2)-yT)); %bottom-right
                
                %pixel index for each neighbor
                nnXindex = find(abs(xT-setup.x_pixMat)<=res);
                nnYindex = find(abs(yT-setup.y_pixMat)<=res);
                obj.nnIndex = intersect(nnXindex,nnYindex);
            else
                %outside imaging region
            end
            obj.farthest = find(obj.nnDistance==max(obj.nnDistance));
        end
    end
end

