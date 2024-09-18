% Author: Reagan Gakstatter / reg0052@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw01

classdef hw01
    methods (Static)

        function p1()
            % This function only contains comments. Fill the following table. Do not write any code here.
            % :return: no returns
        
            % Write your result and explanation for each command here.
            % 
            % commands         |  results          | explanations
            % -----------------|-------------------|-----------------------------------
            % eps              | 2.2204e-16         | Epsilon (`eps`) is the smallest difference between 1 and the next largest floating-point number.
            % realmax          | 1.7977e+308        | The largest representable floating-point number.
            % realmin          | 2.2251e-308        | The smallest positive normalized floating-point number.
            % 1 + eps - 1      | 2.2204e-16         | Due to floating-point precision, this results in `eps`, as `1 + eps` is slightly greater than 1.
            % 1 + eps/2 - 1    | 0                  | `eps/2` is too small to affect the value of `1`, so the result is 0.
            % realmin/1e10     | 2.2251e-318        | Dividing the smallest normalized number by `1e10` still results in a representable value.
            % realmin/1e16     | 0                  | Dividing `realmin` by `1e16` results in underflow, producing zero.
            % realmax*10       | Inf                | Multiplying `realmax` by 10 exceeds the largest representable number, resulting in `Inf` (infinity).
        end
        
function s_n = p2(n, choice)
           % This function computes the Archimedes' method for pi.
           % :param n: the number of sides of the polygon
           % :param choice: 1 or 2, the formula to use
           % :return: s_n, the approximation of pi using Archimedes' method.
           % Tabulate the error of |s_n - pi| for n = 0, 1, 2, ..., 15 and choice = 1 and 2.
           % for both choices of formulas.
           % n     | choice 1 | choice 2
           % ------|----------|----------
           % 0     | 0.322509 | 0.322509
           % 1     | 0.073798 | 0.073798
           % 2     | 0.018067 | 0.018067
           % 3     | 0.004494 | 0.004494
           % 4     | 0.001122 | 0.001122
           % 5     | 0.000280 | 0.000280
           % 6     | 0.000070 | 0.000070
           % 7     | 0.000018 | 0.000018
           % 8     | 0.000004 | 0.000004
           % 9     | 0.000001 | 0.000001
           % 10    | 0.000000 | 0.000000
           % 11    | 0.000000 | 0.000000
           % 12    | 0.000000 | 0.000000
           % 13    | 0.000000 | 0.000000
           % 14    | 0.000000 | 0.000000
           % 15    | 0.000001 | 0.000000

           % Explanation of the results (why there is a difference between the two choices):
           %
           % Formula 1 tends to converge slower relative to Formula 2 because as the number of sides
           % increases, the approximation improves. While Formula 2,
           % converges quickly to pi as n increases. The errors for both methods begin to decrease as the number of sides increase.
           % Formula 2 does reach a closer approximation to pi more quickly.
           %
      
       
           % Write your code here
           if choice == 1
               % Use the 1st formula
               p0 = (1/sqrt(3));
               P = zeros(1,n+1);
               P(1) = p0;
               for i = 1:n
                   P(i+1) = (sqrt(1+P(i)^2)-1)/P(i);
               end
               s_n = 2^n * 6 * P(n+1); % Write your code here
               error = abs(s_n-pi);
           else
               % Use the 2nd formula
               p0 = 1/sqrt(3);
               P = zeros(1,n+1);
               P(1) = p0;
               for i = 1:n
                   P(i+1) = P(i)/(1+sqrt(1+P(i)^2));
               end
               s_n = 2^n * 6 * P(n+1); % Write your code here
               error = abs(s_n-pi);
           end
       end



        function s = p3(a)
            % This function computes the Kahan summation algorithm.
            % :param a: a vector of numbers
            % :return: summation of the vector a using Kahan summation algorithm
            
            s = a(1);  % Initialize sum
            e = 0;     % Initialize the error term

            for j = 2:length(a)
                y = a(j) - e;         % Compensate for error
                temp = s + y;         % Perform summation
                e = (temp - s) - y;   % Update error term
                s = temp;             % Update sum
            end
        end

      function p4(a)
            % This function tests the performance of Kahan summation algorithm against native sum.
            % :param a: a vector of numbers in double precision.
            % :return: no returns
        
            % Test this function with a = rand(n, 1) with various size n multiple times. 
            % Summarize your findings below.
            %
            % Findings:
            % Kahan summation in single precision produces more accurate results than the naive
            % summation, specifically for large arrays. This difference becomes greater as the size
            % of the array increases. While the built-in sum function is faster, it accumulates
            % more floating-point errors relative to Kahan summation.
            % Kahan summation reduces the rounding error, especially when the array consists
            % of numbers with large variations in magnitude.
            %
        
            % Convert a to single precision
            single_a = single(a); 
            
            % Kahan sum of 'a' under double precision (regarded as truth)
            s = hw01.p3(a); 
            
            % Kahan sum of 'single_a' under single precision
            single_Kahan_s = hw01.p3(single_a); 
            
            % Naive sum of 'single_a' under single precision
            single_naive_s = sum(single_a); 
            
            % Display the error for both summation methods in single precision
            disp(['Error of naive sum under single precision: ', num2str(single_naive_s - s)]);
            disp(['Error of Kahan sum under single precision: ', num2str(single_Kahan_s - s)]);
        end

    end
end
