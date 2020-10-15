function pt = tstat3( v, tp, stat )
% TSTAT3 computes one of three t_statistics: one-sided t-probability, or
%        two-sided t-probability, or the inverse t-statistic in any single 
%        call.  It does not take vectors as arguments.  
% 
%   INPUT ARGUMENTS: 
%       v       Degrees-of freedom (integer)
% 
%       tp      T-statistic: to produce t-probability
%                   OR
%               Probability (alpha): to produce t-statistic
% 
%       stat    Desired test — 
%                   'one'   One-tailed t-probability
%                   'two'   Two-tailed t-probability
%                   'inv'   Inverse t-test
% 
%   OUTPUT ARGUMENT:
%       pt      T-probability OR T-statistic, as requested in ‘stat’
% 
%   USE:
%       FIND ONE-TAILED PROBABILITY GIVEN t-STATISTIC & DEGREES-OF-FREEDOM
%           p = tstat3(v, t, 'one')
% 
%       FIND TWO-TAILED PROBABILITY GIVEN t-STATISTIC & DEGREES-OF-FREEDOM
%           p = tstat3(v, t, 'two')
% 
%       FIND ONE-TAILED t-STATISTIC GIVEN PROBABILITY & DEGREES-OF-FREEDOM
%           t = tstat3(v, p, 'inv')
%       
%  
%  
% Star Strider — 2016 01 24 — 

% Copyright (c) 2016, Star Strider
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% % % % T-DISTRIBUTIONS — 
% Variables: 
% t: t-statistic
% v: degrees of freedom

tdist2T = @(t, v) 1 - betainc(v ./ (v + t .* t), v * 0.5, 0.5);                              % 2-tailed t-distribution
tdist1T = @(t, v) 0.5 * (1 + tdist2T(t,v));                                          % 1-tailed t-distribution

% This calculates the inverse t-distribution (parameters given the
%   probability ‘alpha’ and degrees of freedom ‘v’: 
t_inv = @(alpha,v) fzero(@(tval) (max(alpha, (1 - alpha)) - tdist1T(tval, v)), 5);  % T-Statistic Given Probability ‘alpha’ & Degrees-Of-Freedom ‘v’

statcell = {'one' 'two' 'inv'};                                                 % Available Options
nc = cellfun(@(x)~isempty(x), regexp(statcell, stat));                          % Logical Match Array
n = find(nc);                                                                   % Convert ‘nc’ To Integer

if isempty(n)                                                                   % Error Check ‘if’ Block
    error('                    —> The third argument must be either ''one'', ''two'', or ''inv''.')
elseif (n == 3) && ((tp < 0) || (tp > 1))
    error('                    —> The probability for ''inv'' must be between 0 and 1.')
elseif (isempty(v) || (v <= 0))
    error('                    —> The degrees-of-freedom (''v'') must be > 0.')
end

switch n                                                                        % Calculate Requested Statistics
    case 1
        pt = tdist1T(tp, v);
    case 2
        pt = tdist2T(tp, v);
    case 3
        pt = t_inv(tp, v);
    otherwise
        pt = NaN;
end

end
% ———————————————————————————  END: tstat3.m  ————————————————————————————

