function agecat = gAgeCat2char(gAgeCat)

if(gAgeCat == 0)
    agecat = '<1';
elseif(gAgeCat == 1)
    agecat = '1-9';
elseif(gAgeCat == 2)
    agecat = '10-19';
elseif(gAgeCat == 3)
    agecat = '20-29';
elseif(gAgeCat == 4)
    agecat = '30-39';
elseif(gAgeCat == 5)
    agecat = '40-49';
elseif(gAgeCat == 6)
    agecat = '50-59';
elseif(gAgeCat == 7)
    agecat = '60-69';
elseif(gAgeCat == 8)
    agecat = '70-79';
elseif(gAgeCat == 9)
    agecat = '80-101'; 
end