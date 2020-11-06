function [type] = motorclass(Impulse)
% lbf-sec

if Impulse > 9204.8
    type = 'I wish ';
elseif Impulse > 4602.4
    type = 'O';
elseif Impulse > 2301.20
    type = 'N';
elseif Impulse > 1150.6
    type = 'M';
elseif Impulse > 575.3
    type = 'L';
elseif Impulse > 287.65
    type = 'K';
elseif Impulse > 143.83
    type = 'J';
elseif Impulse > 71.92
    type = 'I';
elseif Impulse > 35.96
    type = 'H';
elseif Impulse > 17.92
    type = 'G';
elseif Impulse > 8.96
    type = 'F';
elseif Impulse > 4.48
    type = 'E';
elseif Impulse > 2.24
    type = 'D';
elseif Impulse > 1.12
    type = 'C';
elseif Impulse > .56
    type = 'B';
elseif Impulse > .29
    type = 'A';
end