function d = juliandate(date)
%% Function description
%
% This function takes a date vector and outputs the Julian Calendar Date
%
% Inputs
% date	date vector with format [year,month,day,hour,minute,second] or [year,month,day]
%
% Outputs
% d	Julian date (days)
%
%%
yr = date(1);
mon = date(2);
day= date(3);

if length(date) == 6
    hr = date(4);
    min = date(5);
    sec = date(6);
    d = 367.0 * yr - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 ) + floor( 275 * mon / 9.0 ) + day + 1721013.5 + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;  
else
    d = 367.0 * yr - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 ) + floor( 275 * mon / 9.0 ) + day + 1721013.5;  
end
