% a conversion of easting and northing to ellipsoidal height for Singapore
% region. Coordinates based on SVY21, conversion formular based on below
% reference from a slide citing SLA's work.
% https://fig.net/resources/proceedings/2015/2015_07_vrfp_comm5/5A_Khoo_Singapore_Height_Datum.pdf
% For example: 
% Easting = 3.21542788e4;   
% Northing = 3.975825566e4;
% Elev = 11; % example ellipsoidal height of 11m from photogrammetry
% Hn = geoModN(Easting, Northing)
% HtSHD = Elev - Hn; % conversion to Singapore Height Datum
% for verification of HtSHD, please see 
% https://app1.sla.gov.sg/sirent/Services-SGeoid09.aspx


function Hn = geoModN(Easting, Northing)
Ee = (Easting - 4813.123)./43462.405;
Nn = (Northing - 25542.573)./24037.684;

Hn = 8.94184+2.08529.*Ee-0.16502.*Nn-0.661429.*Ee.^2-0.139884.*Nn.^2+...
    0.232462.*Ee.^6.*Nn;

end