% % conversion between SVY21 and Lat/Lon
% % From https://www.rubydoc.info/gems/SVY21/1.0.0/SVY21
% % further refer to https://www.linz.govt.nz/data/geodetic-services/coordinate-conversion/projection-conversions/transverse-mercator-transformation-formulae
% % see example below

% % in case degreeMinSec is input lat and long, use dms2degrees([lat, lon])
% dms2degrees([lat, lon])

% [E, N] = LatLon2EN([1.366666; 1.366666], [103.833333; 103.833333]); %input lat, long in degrees
% format long
% E
% N
function [E, N] = LatLon2EN(latData, lonData)
latData = latData*pi/180; 
lonData = lonData*pi/180; 
%ORIGIN_LATITUDE, Long 
lat0 = 1.366666*pi/180; 
lon0 = 103.833333*pi/180;
% FALSE_NORTHING , Easting
N0 = 38744.572;
E0 = 28001.642;

A = 6378137.0; % Semi-major axis of reference ellipsoid.
F = 1.0 / 298.257223563; % Ellipsoidal flattening
K0 = 1; % Central meridian scale factor.
B = A*(1.0-F); %Semi-minor axis of reference ellipsoid.

E2 = (2.0 * F) - (F * F); 
E4 = E2*E2;
E6 = E4*E2;

A0 = 1.0 - (E2 / 4.0) - (3.0 * E4 / 64.0) - (5.0 * E6 / 256.0);
A2 = (3.0 / 8.0) * (E2 + (E4 / 4.0) + (15.0 * E6 / 128.0));
A4 = (15.0 / 256.0) * (E4 + (3.0 * E6 / 4.0));
A6 = 35.0 * E6 / 3072.0;
m = @(lat) A.*(A0.*lat-A2.*sin(2.*lat)+A4.*sin(4.*lat)-A6.*sin(6.*lat));
m0 = m(lat0);
rho = A*(1-E2)./(1-E2.*sin(latData).^2).^1.5;
nu = A./sqrt(1-E2.*sin(latData).^2);
phi = nu./rho; 
t = tan(latData);
omega = lonData-lon0; 
NTerm1 = omega.^2/2.*nu.*sin(latData).*cos(latData);
NTerm2 = (1/24).*omega.^4.*nu.*sin(latData).*cos(latData).^3.*(4.*phi.^2+phi-t.^2);
NTerm3 = omega.^6/720.*nu.*sin(latData).*cos(latData).^5.*(8.*phi.^4.*(11-24.*t.^2)-...
    28*phi.^3.*(1-6.*t.^2)+phi.^2.*(1-32.*t.^2)-phi.*(2.*t.^2)+t.^4);
NTerm4 = omega.^8/40320.*nu.*sin(latData).*cos(latData).^7.*(1385-3111*t.^2+543.*t.^4-t.^6);
N = N0+K0*(m(latData)-m0+NTerm1+NTerm2+NTerm3+NTerm4);
ETerm1 = omega.^2/6.*cos(latData).^2.*(phi-t.^2);
ETerm2 = omega.^4/120.*cos(latData).^4.*(4*phi.^3.*(1-6*t.^2)+...
    phi.^2.*(1+8.*t.^2)-phi.*2.*t.^2+t.^4);
ETerm3 = omega.^6/5040.*cos(latData).^6.*(61-479.*t.^2+179.*t.^4-t.^6);
E = E0 + K0.*nu.*omega.*cos(latData).*(1+ETerm1+ETerm2+ETerm3);
end