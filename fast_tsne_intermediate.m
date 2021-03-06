function mappedX = fast_tsne(X, no_dims, initial_dims, perplexity, theta, random_seed, init_Y, save_interval)
%FAST_TSNE Runs the C++ implementation of Barnes-Hut t-SNE
%
%   mappedX = fast_tsne(X, no_dims, initial_dims, perplexity, theta, random_seed, init_Y, save_interval)
%
% Runs the C++ implementation of Barnes-Hut-SNE. The high-dimensional 
% datapoints are specified in the N x D matrix X. The dimensionality of the 
% datapoints is reduced to initial_dims dimensions using PCA (default = 50)
% before t-SNE is performed. Next, t-SNE reduces the points to no_dims
% dimensions. The perplexity of the input similarities may be specified
% through the perplexity variable (default = 30). The variable theta sets
% the trade-off parameter between speed and accuracy: theta = 0 corresponds
% to standard, slow t-SNE, while theta = 1 makes very crude approximations.
% Appropriate values for theta are between 0.1 and 0.7 (default = 0.5).
% The function returns the two-dimensional data points in mappedX. The
% parameter random_seed can be used to set the random seed for this run
% (default = 0). The variable init_Y can be used to specify an initial map;
% the input matrix should have size N x no_dims. The variable save_interval
% can used to set the number of iterations between two intermediate saves
% of the map during learning (default = [], which means no intermediate
% saving is performed).
%
%
% NOTE: The function is designed to run on large (N > 5000) data sets. It
% may give poor performance on very small data sets (it is better to use a
% standard t-SNE implementation on such data).


% Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the Delft University of Technology.
% 4. Neither the name of the Delft University of Technology nor the names of 
%    its contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
% EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
% OF SUCH DAMAGE.


    if nargin < 2
        error('Function should have at least two inputs.');
    end
    if ~exist('initial_dims', 'var') || isempty(initial_dims)
        initial_dims = 50;
    end
    if ~exist('perplexity', 'var') || isempty(perplexity)
        perplexity = 30;
    end
    if ~exist('theta', 'var') || isempty(theta)
        theta = 0.5;
    end
    if ~exist('random_seed', 'var') || isempty(random_seed)
        random_seed = 0;
    end
    if ~exist('init_Y', 'var') || isempty(init_Y)
        init_Y = [];
    else
        if (size(init_Y, 1) ~= size(X, 1)) || (size(init_Y, 2) ~= no_dims)
            error('Size of initial map does not match other inputs.');
        end
    end
    if ~exist('save_interval', 'var') || isempty(save_interval)
        save_interval = -1;
    end
    
    % Perform the initial dimensionality reduction using PCA
    X = double(X);
    X = bsxfun(@minus, X, mean(X, 1));
%     covX = X' * X;
%     [M, lambda] = eig(covX);
%     [~, ind] = sort(diag(lambda), 'descend');
%     if initial_dims > size(M, 2)
%         initial_dims = size(M, 2);
%     end
% 	M = M(:,ind(1:initial_dims));
%     X = X * M;
%     clear covX M lambda

if ~isempty(initial_dims)
    [~,~,M] = svds(X,initial_dims);
    %         covX = X' * X;
    %         [M, lambda] = eig(covX);
    %         [~, ind] = sort(diag(lambda), 'descend');
    %         if initial_dims > size(M, 2)
    %             initial_dims = size(M, 2);
    %         end
    %         M = M(:,ind(1:initial_dims));
    X = X * M;
    clear covX M lambda
end
    
    % Run the fast diffusion SNE implementation
    write_data(X, no_dims, theta, perplexity, random_seed, init_Y, save_interval);
    tic, system('./bh_tsne_intermediate'); toc
    [mappedX, landmarks, costs] = read_data('result.dat');   
    landmarks = landmarks + 1;                 % correct for Matlab indexing
    delete('data.dat');
end


% Writes the datafile for the fast t-SNE implementation
function write_data(X, no_dims, theta, perplexity, random_seed, init_Y, save_interval)
    [n, d] = size(X);
    h = fopen('data.dat', 'wb');
	fwrite(h, n, 'integer*4');
	fwrite(h, d, 'integer*4');
    fwrite(h, theta, 'double');
    fwrite(h, perplexity, 'double');
	fwrite(h, no_dims, 'integer*4');
    fwrite(h, X', 'double');
    fwrite(h, int32(round(random_seed)), 'integer*4');
    fwrite(h, save_interval, 'integer*4');
    if ~isempty(init_Y)
        fwrite(h, init_Y', 'double');
    end
	fclose(h);
end

