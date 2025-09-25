function sliced = sliceDim(array, dim, idx)
%SLICEDIM Extract a slice along a specific dimension and remove that dimension
%   SLICED = SLICEDIM(ARRAY, DIM, IDX) extracts a slice from ARRAY at
%   index IDX along dimension DIM, and returns an array with one fewer
%   dimension than the original array.
%   If the result has only one dimension, a dimension of size 1 is
%   added at the end to ensure the result is at least 2D.
%
% Inputs:
%   ARRAY - Input multidimensional array
%   DIM   - Dimension along which to slice (1, 2, 3, ...)
%   IDX   - Index value along that dimension
%
% Output:
%   SLICED - Extracted slice with one fewer dimension, guaranteed to be at least 2D

% Get the original array size
    array_size = size(array);
    
    % Create indexing structure
    idx_cell = repmat({':'}, 1, ndims(array));
    idx_cell{dim} = idx;
    
    % Extract the slice
    sliced = array(idx_cell{:});
    
    % Create new target dimensions
    new_size = array_size;
    new_size(dim) = [];
    
    % Ensure the result is at least 2D
    if length(new_size) == 1
        new_size(2) = 1;  % Add a dimension of size 1 at the end to ensure 2D array
    end
    
    % Reshape array to remove the specified dimension, ensuring at least 2D
    sliced = reshape(sliced, new_size);
end