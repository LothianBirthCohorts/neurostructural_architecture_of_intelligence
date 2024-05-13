% diverging winter-autumn colormap (white at 0)
function [cMap] = winter_autumn(n)

wintwhite = interp1([0;1], [0 1 0.5; 1 1 1], linspace(0,1,32));
wintwhite = wintwhite(2:32,:);

autwhite = interp1([0;1], [1 1 0; 1 1 1], linspace(0,1,32));
autwhite = autwhite(2:32,:);


wintall = [winter(128); wintwhite];
autall = [autumn(128); autwhite];
autall = flipud(autall);

cMap = [65/255 74/255 76/255; wintall; autall; 65/255 74/255 76/255];
wintwhite = interp1([0;1], [0 1 0.5; 1 1 1] , linspace(0,1,32));
wintwhite = wintwhite(2:32,:);
autwhite = interp1([0;1], [1 1 0; 1 1 1], linspace(0,1,32));
autwhite = autwhite(2:32,:);
wintall = [winter(128); wintwhite];
autall = [autumn(128); autwhite];
autall = flipud(autall);
darkred = interp1([0;1], [1 0 0; 74/255 4/255 4/255] , linspace(0,1,50));
darkblue = interp1([0;1], [0 0 128/255; 0 0 1], linspace(0,1,50));
cMap = [65/255 74/255 76/255; darkblue; wintall; autall; darkred; 65/255 74/255 76/255];

return