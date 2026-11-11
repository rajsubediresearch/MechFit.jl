function curves=AddErrorStructure(yi,M,dist1,factor1,d)

%yi is cumulative curve

yi = yi(:);

curves = zeros(numel(yi),M);

increments = diff(yi);

for realization = 1:M

    yirData = zeros(numel(yi),1);
    yirData(1) = yi(1);

    % Existing distribution-specific calculations
    % Use increments(t-1) instead of repeatedly calculating:
    % yi(t)-yi(t-1)

    switch dist1

        case 0

            %Normal distribution
            for t=2:length(yi)
                lambda=normrnd(yi(t)-yi(t-1),factor1);
                yirData(t,1)=lambda;
            end

        case 1

            %Poisson dist
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                yirData(t,1)=poissrnd(lambda,1,1);
            end

        case 2

            % Negative binomial dist with VAR=factor*mean1;

            eps=0.001;
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                if mean1==0
                    mean1=eps;
                end

                var1=mean1*factor1;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);

            end

            % Quasi Poisson distribution

            %             eps=0.001;
            %             for t=2:length(yi)
            %                 lambda=abs(yi(t)-yi(t-1)); % Mean of the Poisson distribution
            %                 phi = factor1;    % Dispersion parameter (> 1 for overdispersion)
            %
            %                 if lambda==0
            %                     lambda=eps;
            %                 end
            %
            %                 % Step 1: Generate Gamma distributed random effects
            %                 % The shape parameter of the Gamma distribution is 1/phi
            %                 % The scale parameter is lambda*phi
            %                 gamma_shape = 1 / phi;
            %                 gamma_scale = lambda * phi;
            %
            %                 gamma_samples = gamrnd(gamma_shape, gamma_scale, 1, 1);
            %
            %                 % Step 2: Generate Poisson samples with the adjusted mean
            %                 yirData(t,1) = poissrnd(gamma_samples);
            %
            %             end


        case 3
            % Negative binomial dist with parameter VAR= MEAN + alpha*MEAN
            eps1=0.001;
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                if mean1==0     % guard against zero mean: mean1=0 gives var1=0,
                    mean1=eps1; % then p1=0/0=NaN and nbinrnd(NaN,NaN)=NaN
                end             % would enter the replicate data.
                var1=mean1+mean1*factor1;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end

        case 4
            % Negative binomial dist with parameter VAR= MEAN +
            % alpha*MEAN^2

            eps1=0.001;
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                if mean1==0     % guard against zero mean: mean1=0 gives var1=0,
                    mean1=eps1; % then p1=0/0=NaN and nbinrnd(NaN,NaN)=NaN
                end             % would enter the replicate data.
                var1=mean1+factor1*mean1^2;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end

        case 5
            % Negative binomial dist with parameter VAR= MEAN +
            % alpha*MEAN^d

            eps1=0.001;
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                if mean1==0     % guard against zero mean: mean1=0 gives var1=0,
                    mean1=eps1; % then p1=0/0=NaN and nbinrnd(NaN,NaN)=NaN
                end             % would enter the replicate data.
                var1=mean1+factor1*mean1^d;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end

        case 6 % Laplace distribution

            for t=2:length(yi)

                % Step 1: Generate uniform random numbers between -0.5 and 0.5
                U = rand(1, 1) - 0.5;

                % Step 2: Apply the inverse CDF of the Laplace distribution
                lambda = abs(yi(t) - yi(t-1));  % Location parameter (e.g., difference between observations)
                b = factor1;             % Scale parameter

                % Calculate the Laplace-distributed random variable using the inverse CDF

                yirData(t,1)=lambda - b * sign(U) .* log(1 - 2 * abs(U));
            end

    end


    curves(:,realization) = yirData;
end
