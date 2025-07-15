% nakaRush.m
% written by Brandon Nanfito
% function for fitting the naka rushton function to empirical data

% INPUTS
% 1. r = vector of responses for each corresponding contrast value in c
% 2. c = vector of contrast values for each corresponding response in r
%           value should be 0 >= c <= 100
% 3. p = plot the results? 1:yes, 0:no

% OUTPUTS
% 1. nk = structure with fields:
%           1. fit = function that takes in vector of contrast values and
%           outputs the predicted response from the model fit.
%           2. rMax = model parameter for the asymptotic value the response
%           saturates to at high contrast
%           3. c50 = model parameter for contrast value at which response
%           reaches half of its maximum
%           4. n = model parameter controling the steepness of the
%           inflection point in the curve
%           5. resnorm = the sum of squared residuals (measure of error)
%           6. residuals = the difference (error) between the model fit and
%           empirical data points
%           7. aic = akaike information criterion
%           8. bic = bayesian information criterion

function [nk] = nakaRush(r,c,p)

    %initial guesses at model parameters
    rMax = max(r); 
    c50 = 50;
    n = 1;
    x0 = [n c50 rMax]; %initial guesses for model parameters
    lb = [0 0 0]; %lower bound on model parameters
    ub = [inf 100 inf]; %upper bound on model parameters

    %define error function that computes the residuals between observed
    %response (r) and the model predictions given a set of parameters
    %(pars) for the empirically evaluated contrast values (c)
    err = @(pars) (pars(1)*( (c.^pars(3))./((c.^pars(3))+(pars(2)^pars(3))) )) - r;
    [params,resnorm,residuals] = lsqnonlin(err,x0,lb,ub);
    rMax = params(1);
    c50 = params(2);
    n = params(3);

    fit = @(cont) rMax*( (cont.^n)./((cont.^n)+(c50^n)) );
        
    %calculate goodness of fit metrics
    nC = length(c);
    nP = length(params);
    rms = sqrt(resnorm/nC);
    aic = (nC*log(rms))+(2*nP);
    bic = (nC*log(rms))+(log(nC)*nP);

    %create output structure
    nk.fit = fit;
    nk.rMax = rMax;
    nk.c50 = c50;
    nk.n = n;
    nk.resnorm = resnorm;
    nk.residuals = residuals;
    nk.aic = aic;
    nk.bic = bic;

    if p == 1
        figure; hold on
        pl(1) = plot(c,r,'bo','LineWidth',2); %plot empirical data points
        x = 1:100;
        pl(2) = plot(x,nk.fit(x),'k','LineWidth',2); %plot model fit
        pl(3) = yline(rMax,'r--'); %plot asymptotic saturation point of response at high contrasts
        pl(4) = plot([0 c50],[rMax/2 rMax/2],'r'); %plot the c50 value and its corresponding response value
        plot([c50 c50],[0 rMax/2],'r')
        xlim([0 100])
        legend(pl,{'empirical data','model fit','response saturation','c50/half max response'},'Location','southeast')
    end

end