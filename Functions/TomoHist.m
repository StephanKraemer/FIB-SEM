function [ival,nval,varargout] = TomoHist(IMS, varargin)
%  calculate histogram of image, image series or volume
% -----------------------------------------------------------------------
%
%  SYNTAX   TomoHist(IMS);
%           [ival, nval] = TomoHist(IMS)
%           [ival, nval] = TomoHist(IMS,hpar)
%           [ival,nval,fit] = TomoHist(IMS,hpar,fit)
%
%  INPUT    IMS     data
%           hpar    (o) MATLAB structure with parameters for non-default
%                   settings
%                   .nbin    alternative step size for histogram
%                   .mima    values of lowest and highest bin
%                   .YNplot  'yes' plot dat
%                   .YNomit0 'yes' omit 0-values
%           fitinfo  MATLAB structure with fitting data
%                   .ng      number of Gaussians to fit  
%
%  OUTPUT   ival    intensity values at center of bins
%           nval    number of pixels per intensitya
% -----------------------------------------------------------------------


%% Check input and set variables

% default values
nbin     = 256;
mima     = MiMa(IMS);
YNplot   = 'no';
YNomit0  = 'no';

if nargin>1
    hpar = varargin{1};
    
    if isfield(hpar,'nbin');    nbin    = hpar.nbin;    end
    if isfield(hpar,'mima');    mima    = hpar.mima;    end
    if isfield(hpar,'YNplot');  YNplot  = hpar.YNplot;  end
    if isfield(hpar,'YNomit0'); YNomit0 = hpar.YNomit0; end
end


binWidth  = ( mima(2)-mima(1) )/ nbin;



%% Calcualate histogram

% NOTE  Matlab recommended to replace 'histc' with 'histcounts'

% PROC  - calcualte edge locations based on desired bin width and min/max
%         range
%       - ival represents the center of the bin



if strcmp(YNomit0,'yes')
    M = IMS~=0;
    [nval,edges] = histcounts( IMS(M), 'BinWidth',binWidth, 'BinLimits',mima );
else
    [nval,edges] =  histcounts( IMS, 'BinWidth',binWidth, 'BinLimits',mima );
end

% center of bins
ival      = edges(1:end-1) + binWidth/2;



%% Fit Guassian

if nargin>2

    fitinfo = varargin{2};
    
    ftstring = sprintf('gauss%u',fitinfo.ng);
    ft = fittype(ftstring);
    
    % take out all 0-valued bins
    nnull = find(nval~=0);
    
    % prepare variables for fit
    [ivnn,nvnn] = prepareCurveData( ival(nnull), nval(nnull) );
    
    % add start values if given
    if isfield(fitinfo,'StartPoint')
        cf = fit(ivnn,nvnn,ft,'StartPoint',fitinfo.StartPoint);
    else
        cf = fit(ivnn,nvnn,ft);
    end
    
    yf = feval(cf,ivnn);

    fitinfo.result = cf;
    fitinfo.ival = ivnn;
    fitinfo.nval = nvnn;
    fitinfo.yf = yf;
    fitinfo.coeff = coeffvalues(cf);
    
    varargout{1} = fitinfo;
end



%% Plot result

if nargout==0 || strcmp(YNplot,'yes')
    
    figure('Name','Histogram')
    stairs(edges(1:end-1),nval,'LineWidth',1.5)
    
    if exist('fitinfo','var')
        hold on
        % envelope
        plot(ivnn,yf,'r-','LineWidth',1.5)
        
        %individual Gaussians
        coeff = fitinfo.coeff;
        whos coeff

        for i=1:fitinfo.ng
            i0 = 3*(i-1);
            a = squeeze(coeff(i0+1));
            b = squeeze(coeff(i0+2));
            c = squeeze(coeff(i0+3));
            ygauss = a*exp(-((ivnn-b)/c).^2);
            plot(ivnn,ygauss,'r-');
        end
        
        hold off
        legend('histogram','fit')
    end
    xlabel('Intensity','FontSize',14)
    ylabel('Number of pixels','FontSize',14)
    
end


