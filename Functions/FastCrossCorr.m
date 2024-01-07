function ccOut = FastCrossCorr(I,IR,ccIn,varargin)
%  Calculate cross correlation of two images using fft
% -----------------------------------------------------------------------
%
%  Determine displacement of image relative to reference image from
%  location of maximum in cross correlation. Apply polynomial surface
%  fit to achieve potential sub-pixel resolution. 
%
%  SYNTAX   ccOut = FastCrossCorr(I,IR,ccIn)
%           ccOut = FastCrossCorr(I,IR,ccIn,plotfit)
% 
%  INPUT    I       current image
%           IR      reference image: same size, could be another
%                   image of a tilt or focal series or template for
%                   template recognition
%
%           ccIn    MATLAB structure with relevant parameters
%           .mode   use of the CrossCorrelation
%                   'template'  Template matching. 
%                               Taper only image with width w of
%                               ~>=2*diameter of template.
%                   'shift'     Shift between two images.
%                               Taper both images.
%           .w      taper witdth (taper if w>0)
%           .YNfit  'yes': apply sub-pixel surface fit
%           .fit    structure with fit parameters
%            .deg   degree of polynomial (2 or 4]
%            .sx2, .sy2     half width of area cropped around
%                           pixel-resolved cc center
%            .dpf   sub-pixel resolution of first fit [0.1,0.01]
%            .plotfitnum        1 or 2
%
%           plotfit (o) 'yes' plot fit, residuals and contour plot
%           
%
%  OUTPUT   ccOut   updated MATLAB structure
%           .d      shift between images d=(dx,dy)
%           .df1    sub-pixel resolution shift from fit 1
%           .df2    sub-pixel resolution shift from fit 2 
%           .df2conf   confidence interval of position from fit2
%           .CC     cross correlation
%           .fit1,2     results of fit1 and fit2
%;;
%  NOTES    Fit 1 uses a general polynomial fit with 6 coefficients for
%           degree 2 and 15 coefficients for degree 4. The position of the
%           maximum is calculated by interpolating on a finer grid.

%           Fit 2 uses some of the resulting coefficients as starting
%           values for a fit of a paraboloid centered around the max. This
%           fit depends sensitively on the starting value, hence fit 1
%           which behaves more benignly. Fit 2 contains the position of the
%           maximum and a corresponding estimate of its error.
%
% -----------------------------------------------------------------------

% %% Dummy input from workspace
% 
% % images
% I  = double(IMS(134+1:134+310,24+1:24+680,49));
% 
% % IR  = double(IMS(134+1:134+310,24+1:24+680,49));
% IR = double(IMS(134+1:134+310,24+1:24+680,48));
% 
% 
% % Gaussian filter
% H  = fspecial('gaussian',10,3);
% I  = imfilter(I,H,'replicate');
% IR = imfilter(IR,H,'replicate');
% 
% % gradient
% [DX,DY] = gradient(I);
% I = sqrt(DX.^2 + DY.^2);
% [DX,DY] = gradient(IR);
% IR = sqrt(DX.^2 + DY.^2);
% 
% 
% % parameters
% ccIn.mode = 'shift';
% ccIn.YNfit = 'yes';
% 
% ccIn.w = 0;
% 
% ccIn.fit.deg = 4;
% ccIn.fit.sx2 = 2;
% ccIn.fit.sy2 = 3;
% ccIn.fit.dpf = 0.01;
% 
% plotfit = 'yes';
% ccIn.plotfitnum = 2;



%% FastCrossCorr: prepare data

[sx,sy]   = size(I);
[sx2,sy2] = size(IR);

if (sx2 ~= sx) || (sy2 ~= sy)
    error('[Error] FastCrossCorr: Images do not have the same size');
end

if ccIn.w > 0
    I = Taper(double(I),ccIn.w);
    
    if strcmpi(ccIn.mode,'shift')
        IR = Taper(double(IR),ccIn.w);
    else
        IR = double(IR);
    end
else
    I = double(I);
    IR = double(IR);
end

% copy structure
ccOut = ccIn;

% plot
if nargin>3
    plotfit = varargin{1};
else
    plotfit = 'no';
end


%%  FastCrossCorr: FFT based correlation

%  2nd image needs to be rotated by 90deg
CC = fftshift(real(ifft2(fft2(I) .* fft2(rot90(IR,2)))));

% ATT: is fftshift the right thing to do????

ccOut.CC = CC;




%%  FastCrossCorr: peak shift

% Center of cross correlation
% (seems to work for both even and odd image sizes)
cx = floor(sx/2);
cy = floor(sy/2);

[~, di] = max(CC(:));
d = [mod(di,sx)-cx,floor(di/sx)+1-cy];

% save pixel-resolution shift
ccOut.d = d;


%% FastCrossCorr: sub-pixel resolution

%  Refine location with sub-pixel resolution by fitting a polynomial
%  surface (paraboloid)
%
%  Do first a general fit, take over the the px0 and p0y fit results as
%  input for the fit of a paraboloid. This fit contains the center of the
%  curve as fit parameter and provides an estimate of the error in the peak
%  position. (Need to do general fit first to help finding the solution for
%  the paraboloid).

if strcmpi(ccIn.YNfit,'yes')
    
    % coordinates of small cropped region around peak
    nxf2 = ccIn.fit.sx2;
    nyf2 = ccIn.fit.sy2;
    xc = -nxf2:nxf2;
    yc = -nyf2:nyf2;
    [YC,XC] = meshgrid(yc,xc);
    
    % cross correlation cropped around maximum
    cxf = cx + d(1);
    cyf = cy + d(2);
    
    ic1 = cxf-nxf2;
    ic2 = cxf+nxf2;
    ic3 = cyf-nyf2;
    ic4 = cyf+nyf2;
        
    if ic1<1 || ic2>sx || ic3<1 || ic4>sy

        warning('[FastCrossCorr] warning: maximum too close to edge of image segment. No peak fitting performed.')
        
        % use pixel-resolution drift value
        ccOut.df1 = d;
        ccOut.df2 = d;
        
    else
        CCC = CC(ic1:ic2,ic3:ic4);
        
        [x,y,z] = prepareSurfaceData(XC,YC,CCC);
        
        % save for plotting cross correlation
        ccOut.XC = XC;
        ccOut.YC = YC;
        ccOut.CCC = CCC;
        
        
        % Fit 1
        
        if ccIn.fit.deg==4
            ft = fittype( 'poly44' );
        else
            ft = fittype( 'poly22' );
        end
        opts = fitoptions(ft);
        [fitresult1,gof1] = fit([x,y],z,ft,opts);
        cval1 = coeffvalues(fitresult1);
        
        % peak position fit 1 on fine grid
        xcf = -nxf2 : ccIn.fit.dpf : nxf2;
        ycf = -nyf2 : ccIn.fit.dpf : nyf2;
        [YF,XF] = meshgrid(xcf,ycf);
        xint = XF(:);
        yint = YF(:);
        
        qf = feval(fitresult1,[xint,yint]);
        [~,qi] = max(qf);
        df1 = [xint(qi),yint(qi)];
        
        % save fit and sub-pixel shift 1 (general polynomial)
        ccOut.fit1.result = fitresult1;
        ccOut.fit1.gof = gof1;
        ccOut.fit1.cval = cval1;
        ccOut.fit1.XF = XF;
        ccOut.fit1.YF = YF;
        ccOut.fit1.CCF = reshape(qf,size(XF));
        
        ccOut.df1 = d + df1;
        
        
        % Fit 2
        
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        
        if ccIn.fit.deg==4
            % have peak position as parameters a,b; the rest of the parameters
            % are polynomial coefficients
            % ft = fittype('c + d*(x-a)^2 + e*(y-b)^2 + f*(x-a)^4 + g*(y-b)^2', 'independent', {'x', 'y'}, 'dependent', 'z' );
            
            fstr = 'c + d*(cosd(h)*x-sind(h)*y-a)^2 + e*(sind(h)*x+cosd(h)*y-b)^2 + f*(cosd(h)*x-sind(h)*y-a)^4 + g*(sind(h)*x+cosd(h)*y-b)^2';
            ft = fittype(fstr, 'independent', {'x', 'y'}, 'dependent', 'z' );
            
            % polynomial coefficients as starting values + limits
            cv1 = [cval1(1),cval1(4),cval1(6),cval1(11),cval1(15)];
            lb1 = zeros([1,5]);
            ub1 = zeros([1,5]);
            for i=1:5
                if sign(cv1(i))>0
                    lb1(i) = 0;
                    ub1(i) = Inf;
                else
                    lb1(i) = -Inf;
                    ub1(i) = 0;
                end
            end
            
            opts.StartPoint = [df1(1),df1(2),cv1(1),cv1(2),cv1(3),cv1(4),cv1(5),0];
            opts.Lower = [-1,-1,lb1(1),lb1(2),lb1(3),lb1(4),lb1(5),-90];
            opts.Upper = [+1,+1,ub1(1),ub1(2),ub1(3),ub1(4),ub1(5),90];
            
        else
            % have peak position as parameters a,b; the rest of the parameters
            % are polynomial coefficients
            fstr = 'c + d*(cosd(f)*x-sind(f)*y-a)^2 + e*(sind(f)*x+cosd(f)*y-b)^2';
            ft = fittype(fstr, 'independent', {'x', 'y'}, 'dependent', 'z' );
            
            % polynomial coefficients as starting values + limits
            cv1 = [cval1(1),cval1(4),cval1(6)];
            lb1 = zeros([1,3]);
            ub1 = zeros([1,3]);
            for i=1:3
                if sign(cv1(i))>0
                    lb1(i) = 0;
                    ub1(i) = Inf;
                else
                    lb1(i) = -Inf;
                    ub1(i) = 0;
                end
            end
            
            opts.StartPoint = [df1(1),df1(2),cv1(1),cv1(2),cv1(3)];
            opts.Lower = [-1,-1,lb1(1),lb1(2),lb1(3)];
            opts.Upper = [+1,+1,ub1(1),ub1(2),ub1(3)];
            
        end
        
        
        [fitresult2,gof2] = fit([x,y],z,ft,opts);
        
        % shift value fit 2 + confidence values
        cval2 = coeffvalues(fitresult2);
        cconf = confint(fitresult2);
        
        % save data
        ccOut.fit2.result = fitresult2;
        ccOut.fit2.gof = gof2;
        ccOut.fit2.fstr = fstr;
        ccOut.fit2.cval = cval2;
        
        ccOut.df2 = d + cval2(1:2);
        
        ccOut.df2conf = zeros(2,2);
        ccOut.df2conf(:,1) = cconf(:,1) + d(1);
        ccOut.df2conf(:,2) = cconf(:,2) + d(2);
    end
end


%% FastCrossCorr: plot results

if strcmpi(ccIn.YNfit,'yes') && strcmpi(plotfit,'yes')
    
    figure('Name','Peak fit')
    set(gcf,'Position',[100,100,1200,400])
    
    % Make contour plot.
    subplot(1,3,1);
    if ccIn.plotfitnum==1
        plot( fitresult1,[x,y],z,'Style','Contour' );
    else
        plot( fitresult2,[x,y],z,'Style','Contour' );
    end
    xlabel x
    ylabel y
    grid on

    % Plot fit with data.
    subplot(1,3,2);
    if ccIn.plotfitnum==1
        plot( fitresult1,[x,y],z);
    else
        plot( fitresult2,[x,y],z);
    end
        xlabel x
    ylabel y
    zlabel cc
    grid on
    view(45,20);
    
    % Plot residuals.
    subplot(1,3,3);
    if ccIn.plotfitnum==1
        plot( fitresult1,[x,y],z, 'Style', 'Residual' );
    else
        plot( fitresult2,[x,y],z, 'Style', 'Residual' );
    end
        xlabel x
    ylabel y
    zlabel cc
    grid on
    view(45,20);
    
end



