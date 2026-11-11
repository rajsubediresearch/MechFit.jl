function [objfunction, scoringCurve] = quantdiffObjectiveValue(fitcurve, observations, method1, alpha, d)
%QUANTDIFFOBJECTIVEVALUE Score a supplied curve without another ODE solve.
% All six supported formulas are retained from parameterSearchODE, including:
%   * 0.001 replacement for exact zeros (also for LS and SAD);
%   * the existing gamma-term treatment of nonpositive observations;
%   * the existing parameter-independent likelihood constants being omitted.
% fitcurve itself is not modified. Invalid accepted scores raise an error;
% parameterSearchODE maps such errors to Inf for scalar optimizer calls only
% (fit_model then bounds it above every valid objective value).

if ~isvector(fitcurve) || ~isvector(observations) || ...
        isempty(fitcurve) || numel(fitcurve) ~= numel(observations) || ...
        ~isreal(fitcurve) || ~isreal(observations) || ...
        any(~isfinite(fitcurve(:))) || any(~isfinite(observations(:)))
    error('QuantDiffForecast:InvalidScoreInput', ...
        'A score requires matching finite real observation and curve vectors.');
end
yfit = fitcurve(:);
ydata = observations(:);
if sum(yfit) == 0
    error('QuantDiffForecast:ZeroCurve', ...
        'The original objective assigns a penalty to a zero-sum curve.');
end
yfit(yfit == 0) = 0.001;
yg = max(ydata,0);

switch method1

        case 0  %Least squares

            objfunction=sum((ydata-yfit).^2);


        case 1 % MLE for Poisson distribution (negative log-likelihood)

            objfunction=-sum(ydata.*log(yfit)-yfit);


        case 3  % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean;

            % Vectorized log-likelihood using the identity
            %   sum_{j=0}^{y-1} log(j+m) = gammaln(y+m) - gammaln(m),
            % with m = yfit/alpha.
            m = yfit./alpha;
            sum1 = sum( gammaln(yg+m) - gammaln(m) ...
                        + ydata.*log(alpha) ...
                        - (ydata + m).*log(1+alpha) );

            objfunction=-sum1;

        case 4
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^2;

            % Vectorized via the same gammaln identity, with constant m = 1/alpha.
            m = 1/alpha;
            sum1 = sum( gammaln(yg+m) - gammaln(m) ...
                        + ydata.*log(alpha.*yfit) ...
                        - (ydata + m).*log(1+alpha.*yfit) );

            objfunction=-sum1;

        case 5
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^d;

            % Vectorized via the same gammaln identity, with
            % m = (1/alpha)*yfit.^(2-d).
            m = (1./alpha).*yfit.^(2-d);
            sum1 = sum( gammaln(yg+m) - gammaln(m) ...
                        + ydata.*log(alpha.*yfit.^(d-1)) ...
                        - (ydata + m).*log(1+alpha.*yfit.^(d-1)) );

            objfunction=-sum1;

        case 6

           objfunction=sum(abs(ydata-yfit));


        otherwise
            error('QuantDiffForecast:UnsupportedMethod', ...
                'Supported estimation methods are 0, 1, 3, 4, 5 and 6.');

    end

if ~isscalar(objfunction) || ~isreal(objfunction) || ~isfinite(objfunction)
    error('QuantDiffForecast:InvalidObjective', ...
        'The returned parameters do not give a finite real objective.');
end
if nargout > 1
    scoringCurve = yfit;
end
end
