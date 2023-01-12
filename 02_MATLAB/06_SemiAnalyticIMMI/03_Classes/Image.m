classdef Image
    %IMAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xArr
        yArr
        xMat
        yMat
        res
        reconstruction
        totalPix
    end
    
    methods
        function obj = Image(xMin,xMax,yMin,yMax,resolution)
            %IMAGE Construct an instance of this class
            %   Detailed explanation goes here
            obj.xArr = xMin:resolution:xMax;
            obj.yArr = yMin:resolution:yMax;
            obj.res = resolution;
            obj.xMat = repmat(obj.xArr,length(obj.yArr),1);
            obj.yMat = repmat(obj.yArr,length(obj.xArr),1)';
            obj.reconstruction = zeros(size(obj.yMat));
            obj.totalPix = size(obj.yMat,1)*size(obj.yMat,2);
        end
        
        function img = FinalShape(obj,recon)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.reconstruction = reshape(recon,size(obj.reconstruction,1),size(obj.reconstruction,2));
            img = obj.reconstruction;
        end
    end
end

