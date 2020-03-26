% % conversion between SVY21 and Lat/Lon
% % From https://www.rubydoc.info/gems/SVY21/1.0.0/SVY21
% % further refer to https://www.linz.govt.nz/data/geodetic-services/coordinate-conversion/projection-conversions/transverse-mercator-transformation-formulae
% % see example below
% [lat, lon] = EN2LatLon([28001.642; 28001.642], [38744.572;38744.572])
% format long
% lat
% lon
% % in case degreeMinSec is preferred, use degrees2dms([lat, lon])
% degrees2dms([lat, lon])
function [lat, lon] = EN2LatLon(E, N)
RAD_RATIO = pi/180;
%ORIGIN_LATITUDE, Long 
lat0 = 1.366666*pi/180; 
lon0 = 103.833333*pi/180;
% FALSE_NORTHING , Easting
N0 = 38744.572;
E0 = 28001.642;
K0 = 1; % Central meridian scale factor.
A = 6378137.0; % Semi-major axis of reference ellipsoid.
F = 1.0 / 298.257223563; % Ellipsoidal flattening
K0 = 1; % Central meridian scale factor.
B = A*(1.0-F); %Semi-minor axis of reference ellipsoid.

N1 = N-N0;
E2 = (2.0 * F) - (F * F); 
E4 = E2*E2;
E6 = E4*E2;

A0 = 1.0 - (E2 / 4.0) - (3.0 * E4 / 64.0) - (5.0 * E6 / 256.0);
A2 = (3.0 / 8.0) * (E2 + (E4 / 4.0) + (15.0 * E6 / 128.0));
A4 = (15.0 / 256.0) * (E4 + (3.0 * E6 / 4.0));
A6 = 35.0 * E6 / 3072.0;
m = @(lat) A.*(A0.*lat-A2.*sin(2.*lat)+A4.*sin(4.*lat)-A6.*sin(6.*lat));
m0 = m(lat0);
m1 = m0+N1./K0;

n = (A-B)/(A+B);
n2 =n * n; 
n3 = n2 * n;
n4 = n2 * n2;
G = A * (1.0 - n) * (1.0 - n2) * (1.0 + (9.0 * n2 / 4.0) + (225.0 * n4 / 64.0)) * RAD_RATIO;
sigma = m1./G*RAD_RATIO;

phi1 = sigma+(3*n/2-27*n.^3/32).*sin(2*sigma)+(21/16*n.^2-55/32*n.^4).*sin(4*sigma)+...
    (151/96*n.^3).*sin(6.*sigma)+(1097/512.*n.^4).*sin(8*sigma);
rho1 = A*(1-E2)./(1-E2.*sin(phi1).^2).^1.5;
nu1 = A./sqrt(1-E2.*sin(phi1).^2);
psi = nu1/rho1;
t1 = tan(phi1);
E1 = E-E0;
x = E1./(K0.*nu1);

phiTerm1 = t1./(K0.*rho1).*(E1.*x/2);
phiTerm2 = t1./(K0.*rho1).*(E1.*x.^3/24).*(-4.*psi.^2+9.*psi.*(1-t1.^2)+12.*t1.^2);
phiTerm3 = t1./(K0.*rho1).*(E1.*x.^5/720).*(8*psi.^4.*(11-24.*t1.^2)-12.*psi.^3.*(21-71.*t1.^2)+...
    15.*psi.^2.*(15-98.*t1.^2+15.*t1.^4)+180.*psi.*(5.*t1.^2-3.*t1.^4)+360.*t1.^4);
phiTerm4 =  t1./(K0.*rho1).*(E1.*x.^7/40320).*(1385+3633.*51.^2+4095.*t1.^4+1575.*t1.^6);

lat = phi1 - phiTerm1 + phiTerm2 - phiTerm3 + phiTerm4;
lat = lat*180/pi;
lambTerm1 = x.*sec(phi1);
lambTerm2 = x.^3.*sec(phi1)/6.*(psi+2.*t1.^2);
lambTerm3 = x.^5.*sec(phi1)/120.*(-4.*psi.^3.*(1-6*t1.^2)+psi.^2*(9-68.*t1.^2)+72.*psi.*t1.^2+24.*t1.^4);
lambTerm4 = x.^7.*sec(phi1)/5040.*(61+662.*t1.^2+1320.*t1.^4+720.*t1.^6);
lon = lon0+ lambTerm1+ lambTerm2+ lambTerm3+ lambTerm4;
lon = lon*180/pi;
end
