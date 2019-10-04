function agecat = gAgeCat2num(gAgeCat)

if(strcmp('<1', strtrim(gAgeCat)))
    agecat = 0;
elseif(strcmp('1-9', strtrim(gAgeCat)))
    agecat = 1;
elseif(strcmp('10-19', strtrim(gAgeCat)))
    agecat = 2;
elseif(strcmp('20-29', strtrim(gAgeCat)))
    agecat = 3;
elseif(strcmp('30-39', strtrim(gAgeCat)))
    agecat = 4;
elseif(strcmp('40-49', strtrim(gAgeCat)))
    agecat = 5;
elseif(strcmp('50-59', strtrim(gAgeCat)))
    agecat = 6;
elseif(strcmp('60-69', strtrim(gAgeCat)))
    agecat = 7;
elseif(strcmp('70-79', strtrim(gAgeCat)))
    agecat = 8;
elseif(strcmp('80-101', strtrim(gAgeCat)))
    agecat = 9;  
end