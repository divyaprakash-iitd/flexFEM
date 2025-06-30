function [C] = cofactor(A)
% cofactor: Calculates the co-factors of each element of a matrix
% [C] = cofactor(A):
%   Calculates the co-factors (signed minors) of each element of a matrix
%   https://en.wikipedia.org/wiki/Minor_(linear_algebra)
% Example:
%   A = magic(4);
%   C = cofactor(A);
% input: 
%   A       = A square matrix
% output:
%   C       = Matrix containing co-factors of each element of A
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 June 2022

    C = zeros(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            AA = A;
            AA(i,:) = [];
            AA(:,j) = [];
        
            C(i,j) = (-1)^(i+j)*det(AA);
        end
    end
end