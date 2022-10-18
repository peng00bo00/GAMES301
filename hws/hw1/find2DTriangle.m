function [pjIdx, pkIdx, plIdx] = find2DTriangle(pjIdx, Theta)
%% find valid triangle on the plane
theta = Theta(pjIdx);
nRing = length(Theta);

if theta < pi
    thetaj = theta+pi;
else
    thetaj = theta-pi;
end

circle = circshift(1:nRing, -pjIdx);

%% check all the vertex on the plane
for k=1:nRing-1
    pkIdx = circle(k); plIdx = circle(k+1);
    thetak = Theta(pkIdx); thetal = Theta(plIdx);

    if isBetween(thetaj, thetak, thetal)
        return
    end
end

end

function ret = isBetween(thetaj, thetak, thetal)
%% a helper function to check if thetaj is between (thetak, thetal)
pj = [cos(thetaj), sin(thetaj)];
pk = [cos(thetak), sin(thetak)];
pl = [cos(thetal), sin(thetal)];

alpha = acos(dot(pj, pk));
beta  = acos(dot(pj, pl));
gamma = acos(dot(pk, pl));

if alpha < gamma && beta < gamma
    ret = 1;
else
    ret = 0;
end

end