function pivots = locateStarCntrs (SmplStrtgy, numStars, numDim, sliceSize, StarFile, funPath, seedNumber)

% =============================================================== %
% Description:
% This function generates "Star Centers" for star-based sampling strategy
% implemented in VARS TOOL framework using various sampling techniques including:
% - Crude Random Sampling (RND)
% - Latin Hypercube Sampling (LHS) (default)
% - Symmetric Latin Hypercube Sampling (SymLHS)
% - Progressive Latin Hypercube Sampling (PLHS)
% - Sobol Sequences (Sobolseq)
% - Halton Sequences (Halton)
%
% =============================================================== %
% Inputs:
%          SmplStrtgy ==========>>> sampling strategy i.e. 'RND', 'LHS', 'PLHS', 'Sobolseq', or 'Halton
%             numStars ==========>>> number of sample points
%              numDim  ==========>>> dimension i.e. number of input parameters/factors
%           sliceSize  ==========>>> number of sample points in each slice for PLHS/reporting frequency
%             Starfile ==========>>> star centers file name (optional)
%             funPath  ==========>>> folder address of star centers if applicable (optional)
%           seedNumber ==========>>> seed number for randomization sampling strategies (optional)
%
%  Outputs:
%           pivots  [numStars x numDim)]  ==========>>>  star centers


% =============================================================== %
% References:
% (1) R. Sheikholeslami & S. Razavi (2017):
% Progressive Latin Hypercube Sampling: An Efficient Approach for Robust
% Sampling-based Analysis of Environmental Models,Env. Model. Software.  93, 109-126
% (2) S. Razavi, & H.V. Gupta (2016):
% A new framework for comprehensive, robust, and efficient global sensitivity analysis:
% 1. Theory. Water Resources Research, 52,423-439.
% (3) S. Razavi & H.V. Gupta (2016):
% A new framework for comprehensive, robust, and efficient global sensitivity analysis:
% 2. Application, Water Resources Research, 52, 440–455.
%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami & Saman Razavi (2017)
% Global Institute for Water Security, School of Environment
% and Sustainability, University of Saskatchewan
% E-mail: razi.sheikholeslami@usask.ca; saman.razavi@usask.ca
% ----------------------------------------------------------------

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isscalar(numStars); error('''numStars'' must be a scalar'); end
if numStars<=1; error('''numStars'' must be larger than 1' ); end
if abs(numStars-round(numStars)); error('''n'' must be integer'); end

if ~isscalar(numDim); error('''numDim'' must be a scalar'); end
if numDim<=1; error('''numDim'' must be larger than 1' ); end
if abs(numDim-round(numDim)); error('''numDim'' must be integer'); end

if isempty ( seedNumber ) == false
    rand('state', seedNumber); %automatic randomization
end

%<><><><><><><><><><><><><><><><><> MAIN CODE <><><><><><><><><><><><><><>

if isempty ( StarFile ) == false
    currFolder = pwd;
    cd (funPath);
    pivots = dlmread(StarFile);
    %	    TEMP = importdata(StarFile); % Razi added
    %   pivots = TEMP;					% Razi added
    cd(currFolder);
else
    switch  SmplStrtgy
        case 'RND'
            pivots = rand(numStars, numDim); %uniformly distributed random numbers
        case 'LHS'
            pivots = lhsdesign(numStars, numDim,'criterion','maximin', 'iterations',100 );
            % generates optimal LHS based on 'maxmin' criterion.
            % by default the number of iteration is 100.
        case 'SymLHS'
            pivots = SymLHS( numStars , numDim, 'maxmin', 10 );
        case 'PLHS'
            numSlices = ceil(numStars/sliceSize);
            if numSlices > 1
                smpSize = numSlices * sliceSize;
                p = PLHSdesign(smpSize, numDim, numSlices, 10, 'maxmin');
                pivots = p(1:numStars,:);
            else
                pivots = lhsdesign(numStars, numDim,'criterion','maximin', 'iterations',100 );
            end
            % generates optimal PLHS based on 'maxmin' criterion.
            % by default the number of iteration is 10.
        case 'Sobolseq'
            p = sobolset(numDim,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
            p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
            pivots = p(1:numStars,:);
        case 'Halton'
            p = haltonset(numDim,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
            p = scramble(p,'RR2'); %PR2 scramling strategy. For no scrambling disable this line.
            pivots = p(1:numStars,:);
        otherwise
            pivots = lhsdesign(numStars, numDim,'criterion','maximin', 'iterations',50 );
    end
end
end