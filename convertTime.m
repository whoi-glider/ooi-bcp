function [daytime] = convertTime(timeinsec)
% Converts time in 'seconds since 1900-01-01 0:0:0' to matlab date format
% where matlab format is days since January 0, 0000
   
secinday = 60*60*24;
[yr,mon,day,hr,mn,sec] = datevec(datestr(timeinsec./secinday));
daytime = datenum(yr + 1900,mon,day,hr,mn,sec);

end