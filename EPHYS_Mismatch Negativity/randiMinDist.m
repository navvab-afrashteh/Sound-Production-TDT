function X = randiMinDist(I1,I2,minDist,N)

x = randi([I1,I2], [5*(I2-I1),1]);
minDist = minDist+1;
% Initialize first point.
X = x(1);
% Try dropping down more points.
counter = 1;
k = 1;
while counter < N
    % Get a trial point.
    k = k+1;
    thisX = x(k);
    % See how far is is away from existing keeper points.
    dists = abs(thisX-X);
    if min(dists) >= minDist
        counter = counter + 1;
        X(counter) = thisX;
    end
end
X = X(:);