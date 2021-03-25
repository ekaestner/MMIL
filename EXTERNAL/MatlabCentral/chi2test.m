function [chi2, critical] = chi2test (data, n, alpha, dist, a1, a2, a3);
%CHI2TEST   Chi-square Test for Continuous Distributions.
%   [A, B] = CHI2TEST(DATA, N, ALPHA, DIST, X, Y, Z) returns 
%   the chi-square statistic for the samples contained in the 
%   row vector DATA. 
%
%   N specifies the number of equal-probability class intervals
%   for the test. ALPHA is the confidence level parameter used 
%   to find the critical chi-square value. 
%
%   DIST is a string containing the probability distribution 
%   that we are testing against.  See the staitsctics toolbox 
%   for supported distributions - 'exp', 'gam', 'unif' are 
%   some of them.
%
%   X, Y, and Z specify the estimated parameters for the 
%   selected DIST.  Some distributions require only one of 
%   these parameters, and the order that these parameters are 
%   provided follows the values given to the cummulative 
%   distribution functions UNIFCDF, GAMCDF, EXPCDF, and others.
%   
%   A is the computed chi-square statistic, and B is the 
%   critical tabulated value at the degrees of freedom.  The 
%   degree of freedom is the number of intervals minus the 
%   number of estimated parameters. 
%   
%   In general, if A is less than B, the H0 hypothesis 
%   that DATA follows the DIST distribution is accepted. 
%
%   An attempt to fit some data with the uniform distribution
%   on the interval from 1.5 to 2.9. The test fails, since A > B: 
%
%   [a, b] = chi2test (data, 10, 0.05, 'unif', 1.5, 2.9)
%   a =
%      38.7500
%   b =
%      14.0671
%   
%   See also MLE, CHI2INV, CHI2STAT, HIST, CDF, ICDF, PDF

%   Copyright 2004 Leonardo Salomone, Carleton University, Ottawa, Canada

% check input 
if nargin < 4, error('Not enough input arguments'); end

% check if number of bins complies to suggested interval range
nsamples = length(data); 

if nsamples < 20,
    error('Sample data too small, chi-square test not recommended');    
elseif nsamples < 50,
    if n < 5,               error('Number of intervals too small'); end
    if n > 10,              error('Number of intervals too large'); end  
elseif nsamples < 100,
    if n < 10,              error('Number of intervals too small'); end
    if n > 20,              error('Number of intervals too large'); end
else    
    if n < sqrt(nsamples),  error('Number of intervals too small'); end    
    if n > nsamples/5,      error('Number of intervals too large'); end    
end;     

% create functions for the bin probabilities and for the inverse cdf for the parameters given
switch (7 - nargin)
    case 2
        prob = inline(sprintf('cdf(''%s'', b, %10.10g) - cdf(''%s'', a, %10.10g)', dist, a1, dist, a1), 'a', 'b'); 
        invcdf = inline(sprintf('icdf(''%s'', x, %10.10g)', dist, a1), 'x');            
    case 1
        prob = inline(sprintf('cdf(''%s'', b, %10.10g, %10.10g) - cdf(''%s'', a, %10.10g, %10.10g)', dist, a1, a2, dist, a1, a2), 'a', 'b'); 
        invcdf = inline(sprintf('icdf(''%s'', x, %10.10g, %10.10g)', dist, a1, a2), 'x');            
    case 0 
        prob = inline(sprintf('cdf(''%s'', b, %10.10g, %10.10g, %10.10g) - cdf(''%s'', a, %10.10g, %10.10g, %10.10g)', dist, a1, a2, a3, dist, a1, a2, a3), 'a', 'b'); 
        invcdf = inline(sprintf('icdf(''%s'', x, %10.10g, %10.10g, %10.10g)', dist, a1, a2, a3), 'x');            
    otherwise
        return;
end;

% find out the bin edges, for equal probabilities of continous distributions, using the inverse CDF
pi = (1/n) .* [0:n]; 
intvls = invcdf(pi);

% ensure that last item is infinity for exp distribution
switch (dist)
case 'exp'
    intvls(end) = inf;
end;

% find bin counts for intervals 
o_freq = histc(data, intvls);

% remove the last bin, it only reflects exact counts at the last item (see histc)
o_freq = o_freq(1:end-1);

% find expected bin probabilities, they are all the same
e_freq = prob(intvls(1),intvls(2)) .* ones(1,n);

% multiply by number of samples for expected frequency
e_freq = length(data) .* e_freq;

% find the chi2 statistic for each interval 
chi2bins = ((o_freq - e_freq).^2)./e_freq;

% sum up the statistics
chi2 = sum(chi2bins);

% find degrees of freedom
df = n - (nargin - 3); 

% find critical value
critical = chi2inv(1-alpha, df);