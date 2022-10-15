function [result] = traprule(x,y)
%% Trapezoidal Rule Calculator (traprule): Numeric Integration with Trapezoidal Rule
% Takes inputs for x (independent variable) and y (dependent variable)
% Outputs one value for numeric integration
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 08/31/21

%% Define Variable to Store Summation

sum = 0;

%% Run Numeric Integration

for i = 1:length(x)-1
    
    %Define delta variable to reduce clutter
    
    delta = (x(i+1) - x(i)) / 2;
    
    %Individual (non-composite) value for Trapezoidal Rule
    trapInd = delta * (y(i) + y(i+1));
    
    %Sum individual components to make Composite Trapezoidal Rule
    sum = sum + trapInd;
end


%% Return Integration Result

result = sum;

end
