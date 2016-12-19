function [yout,mse,deth,h1,h2] = PolarizationDemux( xin,mn,sps,polm, ...
    mu, ntaps, errid, iter, appLMS)
%POLARIZATIONDEMUX Polarization demultiplexing and channel equalization
%   algorithm includes constant modulus algorithm (CMA) and least mean
%   squared (LMS). Multi-stage configuration is applied.
%
%   [1] POLM indicates the type of transimtted signal, i.e., if the 
%       signal is not pol-muxed, then maximum ratio combination is applied.
%
%   [2] STAGE 1 is making one of two polars converge while the other one
%       remain the same. Stage 2 and 3 converge both polars.
%
%   [3] METHOD denotes ways to solve the singularity problem. Method 1
%       uses determination monitoring; method 2 utilizes the Jones matrix
%       relation. Method 0 turns off the de-singularity function
%
%   [4] The 'fast' version uses the .mex files while do not need to return
%       the value of tap coefficients while the .m files have to.
%
%   [5] The tap initialization can be various like [M*1.5]
%
%   Example:
%   
%   See also CMA_FILTER.c

%   Copyright2011 WANGDAWEI $16/3/2011$ 

%   Last Modified at PolyU

if size(xin,2)~=2
    error('POLARIZATIONDEMUX::CMAIN should have two columns')
end

x = DspAlg.Normalize(xin,mn) / (sqrt(mn)-1);

% make sure the tap number is odd
ntaps = ntaps + ~mod(ntaps,2);
nt = ntaps(1);

% taps initialization
halfnt = floor( nt/2);
if ~polm
    r = mean( x(:,1)./ x(:,2));
    M = rotpolar(1, r);
    x_mrc = x*M;
    h1 = zeros(nt,1);
    h1(halfnt+1) = 1;
else
    h1 = zeros(nt,2);
    h2 = zeros(nt,2);
    h1(halfnt+1,:) = [1 0];
    h2(halfnt+1,:) = [0 1];
end

% get constellation and apply periodic boundary condition
method = 0;
cstl = constellation(mn) / (sqrt(mn)-1);
if ~polm
    extx = [ x_mrc(end-halfnt+1:end,1); x_mrc(:,1); x_mrc(1:halfnt,1)];
else
    extx = [ x(end-halfnt+1:end,:); x; x(1:halfnt,:)];
end

% start LMS
if appLMS
    if polm
        for ii = 1:iter(1)
            [xx,mse,deth] = LMS_FILTER(extx,h1,h2,nt,mu(1),cstl,sps);
        end
    else
        for ii = 1:iter(1)
            [xx,mse,deth] = LMS_FILTER_sng(extx,h1,nt,mu(1),cstl,sps);
        end
    end
    yout = xx(  1:sps:end,:);
    mse  = mse( 1:sps:end  );
    deth = deth(1:sps:end  );
    return
end

% start all 3 stages
if polm
    for stage = 1:3
        if ntaps(stage)>nt
            % interp h1 h2
            h1 = interp1(h1,linspace(1,nt,ntaps(stage)),'linear');
            h2 = interp1(h2,linspace(1,nt,ntaps(stage)),'linear');
            % extend x
            halfnt = floor( ntaps(stage)/2 );
            extx = [ x(end-halfnt+1:end,:); x; x(1:halfnt,:) ];
            % update nt
            nt = ntaps(stage);
        end
        cm = get_radius(cstl,errid(stage),mn);
        for ii = 1:iter(stage)
            [xx,mse,deth] = CMA_FILTER(extx,h1,h2,nt,mu(stage),cm,sps,errid(stage),stage,method);
        end
    end
else
    for stage = 1:3
        if ntaps(stage)>nt
            % interp h1
            h1 = interp1(h1,linspace(1,nt,ntaps(stage)),'linear');
            % extend x
            halfnt = floor( ntaps(stage)/2);
            extx = [ x_mrc(end-halfnt+1:end,1); x_mrc; x_mrc(1:halfnt,1) ];
            % update nt
            nt = ntaps(stage);
        end
        cm = get_radius(cstl,errid(stage),mn);
        for ii = 1:iter(stage)
            [xx,mse,deth] = CMA_FILTER_sng(extx,h1,nt,mu(stage),cm,sps,errid(stage));
        end
    end
end

% format the output
yout = xx(  1:sps:end,:);
mse  = mse( 1:sps:end  );
deth = deth(1:sps:end  );


%% Radius Generation
function r = get_radius(c,id,mn)
if id==1 || id==7 || id==5 || id==6
    r = mean(abs(c).^4) / mean(abs(c).^2);
    return;
end
if id==2
    r = mean(real(c).^4) / mean(real(c).^2);
    return;
end
if id==3
    r = sqrt([1+1;1+9;9+9]/9);
    return;
end
if (id==4 && mn==16)
    r = [1+1;1+9;9+9]/9;
    return;
end
if (id==4 && mn==64)
    r = [1+1;1+9;1+25;9+9;9+25;9+49;25+25;25+49;49+49]/49;
    r = [1+1;1+9;1+25;9+25;9+49;25+49;49+49]/49;
    return;
end
if id==8
    r = sqrt([1+1;1+9;1+25;9+9;9+25;9+49;25+25;25+49;49+49]/49);
    return;
end


%% Maximum Ratio Combining (MRC)
function y = rotpolar( x, r)
if abs(r) < 0.5
    m = abs(r);
    delta = angle(r);
    alpha = m^2 / (m^2 + 1);
else
    m = abs( 1 / r);
    delta = -angle( 1 / r );
    alpha = 1 / ( m^2 + 1 );
end
M = [sqrt(alpha)*exp(-1j*delta) -sqrt(1-alpha)*exp(-1j*delta);...
     sqrt(1-alpha)               sqrt(alpha) ];
y = x * M;


%% Channel Parameter Monitoring
function test_h( h1, h2)
scatterplot(h1(:,1));
scatterplot(h1(:,2));
scatterplot(h2(:,1));
scatterplot(h2(:,2));
for ii=1:ntaps
    ev(:,ii) = eig([h1(ii,:);h2(ii,:)]);
end
scatterplot(ev(1,:));
scatterplot(ev(2,:));
figure; plot(angle(ev(1,:)));
figure; plot(angle(ev(2,:)));