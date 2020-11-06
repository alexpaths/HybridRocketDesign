% Alexander Athougies
% Senior Project
% 
% Table Lookup Function (Linear)
% 
% Inputs
% 1) X = Vector of X's
% 2) Y = Vector of Y's
% 3) pt = data point ( X(1) < pt < X(end) )
%
%    X, Y must have same dimensions
% 
% Outputs
% 1) data = value of curve at point
% 

function data = table(X,Y,pt)
    
    index = find(X < pt);
    index2 = find(X > pt);
    
    if isempty(index)
        data = Y(end);
    elseif isempty(index2)
        data = Y(1);
    else
        x2 = X(index2(1));
        x1 = X(index(end));
        
        y2 = Y(index2(1));
        y1 = Y(index(end));
        
        m = (y2 - y1) / (x2 - x1);
        b = y2 - m * x2;
        
        data = m*pt + b;
    end
    
end
