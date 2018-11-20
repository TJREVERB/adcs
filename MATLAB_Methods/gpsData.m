% gpsData
% Preallocates data for GPS module

gps.alt = 408000;   % Altitude of satellite in meters (based on ISS altitude)
gps.lat = 0;        % Latitude of satellite in degrees (North latitude is positive, South latitude is negative)
gps.lon = 0;        % Longitude of satellite in degrees (East longitude is positive, West longitude is negative)
gps.time = decyear(2018,7,1); % Time in satellite orbit (change later to UT, current date is based on predicted launch date)