function [ ltrBs, times, mycputimes, Fnorms ] = plotsDlingam
% Generates scatterplots to visualize the estimation results
% ltrBs: Numbers of non-zero elements in Bp when it is ordered based on the estimated ordering
% times: actual elapsed times
% mycputimes: CPU times
% Fnorms: Frobenius norms
%
% Before using this function, download KernelICA codes from
% "http://www.di.ens.fr/~fbach/", extract the files,
% and add a path to "kernel-ica1_2" containing the extracted files.
%
% Version: 0.9
% Shohei Shimizu (14 Dec 2010)
% based on the LiNGAM code package and Yasuhiro Sogawa and Takanori Inazumi 's codes




% Set the randseed to the same each time
randseed = 0;
fprintf('Using randseed: %d\n',randseed);
rand('seed',randseed);
randn('seed',randseed);

% Set the fonts properly
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');

% Clear figure
figure(1);
clf; drawnow;

% Choose an estimation method:
method = 'DirectLiNGAM'
% method = 'ICA-LiNGAM'

% Choose sparse or full networks:
data = 'sparse'
% data = 'full'

% How many tests to run for each dimensionality and sample size?
ntests = 5

% amount of prior knowledge (range: 0.0 - 1.0) only for DirectLiNGAM
known_rate = 0.0


% Iterate through different data dimensionalities
itotal = 0;
idims = 0;

% testdims = [10 20 50 100] % number of variables
testdims = [10 ] % number of variables

for dims = testdims,
    
    idims = idims+1;
    
    % Iterate through different dataset sizes
    isamples = 0;
    %     testsamples = [500 1000 2000] % sample size
    testsamples = [500 1000] % sample size
    
    for samples = testsamples,
        
        isamples = isamples+1;
        itotal = itotal+1;
        
        % --- Initialize variables to store data to be plotted ---
        
        scatterp = [];
        ltrB = [];
        time = [];
        mycputime = [];
        Fnorm = [];
        
        % Do a number of tests with these parameters
        for itests = 1:ntests,
            itests
            
            % --- Network and data generation ---
            
            switch data
                case 'sparse'
                    
                    disp('sparse networks');
                    
                    % First create the network following Kalisch and Bulman
                    % (JMLR2007)
                    B = zeros( dims );
                    
                    % number of neighbors
                    if rand(1) > 1/2
                        Nneig = 2;
                    else
                        Nneig = 5;
                    end
                    
                    if Nneig/(dims-1) >= 1 || Nneig/(dims-1) <= 0
                        disp('Change the number of neighers');
                    end
                    
                    for irowB = 1: dims
                        for jcolB = 1: irowB
                            B( irowB, jcolB ) = binornd( 1, Nneig/(dims-1), 1, 1);
                        end
                    end
                    
                    % Set the diagonal elements to be zeros
                    for irowB = 1: dims
                        B( irowB, irowB ) = 0;
                    end
                    
                case 'full'
                    
                    disp('fully-connected networks');
                    
                    % Generate full networks
                    % [min max] standard deviation owing to parents
                    parminmax = [0.5 1.5]; % Dummy
                    
                    % [min max] standard deviation owing to disturbance
                    errminmax = [0.5 1.5]; % Dummy
                    
                    % N of parents
                    indegree = Inf;
                    B = randnetbalanced( dims, indegree, parminmax, errminmax );
                    
            end
            
            % Second, determine the values of the path coefficients
            % following Silva et al. (JMLR2006), i.e., following U(0.5,1.5) or
            % U(-1.5,-0.5) and the variances of external influences
            % following U(1,3)
            
            for irowB = 1: dims
                for jcolB = 1: irowB
                    if B( irowB, jcolB ) == 1
                        B( irowB, jcolB ) = sign(rand( 1 )-1/2) * ( rand( 1 ) + 0.5 );
                    end
                end
            end
            
            disturbancevar = ( 2 * rand( dims, 1 ) + 1 );
            disturbancestd = disturbancevar.^(1/2);
            
            % constants, giving non-zero means
            c = 2*randn(dims,1);
            
            % This generates the disturbance variables, which are mutually
            % independent, and non-gaussian
            S = randn(dims,samples);
            
            dist_name = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r']; % Bach and Jordan (JMLR, 2002)
            dist_num = randi(size(dist_name,2),dims,1); % Randomly select distributions
            
            for i=1:dims
                S(i,:) = usr_distrib(dist_name(dist_num(i,1)),'rnd',samples);
                fprintf('[%c]',dist_name(dist_num(i,1)));
            end
            fprintf('\n');
            disp('Data generated.');
            
            % This normalizes the disturbance variables to have the
            % appropriate scales
            S = S./((sqrt(mean((S').^2)')./disturbancestd)*ones(1,samples));
            
            % Now we generate the data one component at a time
            Xorig = zeros(dims,samples);
            for i=1:dims,
                Xorig(i,:) = B(i,:)*Xorig + S(i,:) + c(i);
            end
            
            % Select a random permutation because we do not assume that
            % we know the correct ordering of the variables
            p = randperm(dims);
            
            % Permute the rows of the data matrix, to give us the
            % observed data
            X = Xorig(p,:);
            
            % Permute the rows and columns of the original generating
            % matrix B so that they correspond to the actual data
            Bp = B(p,p);
            
            % Permute the generating disturbance stds so that they
            % correspond to the actual data
            %             disturbancestdp = disturbancestd(p);
            
            % Permute the generating constants so that they correspond to
            % the actual data
            %             cp = c(p);
            
            
            % --- Call LiNGAM to do the actual estimation ---
            
            tstart = cputime; % start: to measure CPT times
            tic; % start: to measure actual elapsed times
            switch method
                
                case 'DirectLiNGAM'
                    
                    % Generate prior knowledge matrix
                    A = inv(eye(dims,dims)-B);
                    M = double(tril(A ~= 0, -1));
                    mask = (rand(dims,dims) > known_rate) & ~eye(dims,dims);
                    M(mask) = -1;
                    % Permute M
                    Mp = M(p,p);
                    
                    [Best, stde, ci, k] = Dlingam(X, 'pk', Mp);
                    
                case 'ICA-LiNGAM'
                    
                    [Best , stde, ci, k] = estimate(X);
                    
            end
            tstopcputime = cputime-tstart;  % end: to measure CPT times
            tstop = toc; % end: to measure actual elapsed times
            
            % --- Collect the data for plotting ---
            
            % Gather everything for one scatter plot
            ltrB = [ ltrB; sum( sum( triu( Bp( k, k ), 1 ) ~= 0 ) )]; % Number of non-zero elements in Bp when it is ordered based on the estimated ordering. See Shimizu et al. (UAI09)
            time = [ time; tstop]; % actual elapsed times
            mycputime = [ mycputime; tstopcputime]; % CPU times
            Fnorm = [ Fnorm; norm( Best - Bp,'fro' ) ]; % Frobenius norms between true B (Bp) and estimated B (Best)
            
            if dims < 20,
                scatterp = [scatterp; [Bp(:) Best(:)]];
            else
                % Take a random selection of the points so that
                % they sum up to 1000.
                nselection = floor(1000/ntests);
                indexperm = randperm(dims * dims);
                selection = indexperm(1:nselection);
                
                orig = Bp(:);
                est = Best(:);
                scatterp = [scatterp; [orig(selection) est(selection)]];
            end
        end
        
        % Create scatterplot for this combination of dims/samples
        figure(1);
        subplot(length(testdims),length(testsamples),itotal);
        mv = 3;
        figure(1); plot(scatterp(:,1),scatterp(:,2),'ko', ...
            [-mv mv],[-mv mv],'r-');
        axis([-mv mv -mv mv]);
        drawnow;
        
        % Median number of errors
        ltrBs( idims, isamples ) = median( ltrB );
        times( idims, isamples ) = median( time );
        mycputimes( idims, isamples ) = median( mycputime );
        Fnorms( idims, isamples ) = median( Fnorm );
        
    end
    
end

% Show some parameters to remember under what setting the simulation was
% conducted
testdims
testsamples

switch data
    case 'sparse'
        disp('sparse networks');
    case 'full'
        disp('fully-connected networks');
end

switch method
    case 'DirectLiNGAM'
        disp('DirectLiNGAM');
    case 'ICA-LiNGAM'
        disp('ICA-LiNGAM');
end

known_rate % Amount of prior knowledge

% Show some statistics to see how well the estimation worked:
ltrBs % Numbers of non-zero elements in Bp when it is ordered based on the estimated ordering
times % actual elapsed times
mycputimes % CPU times
Fnorms % Frobenius norms



