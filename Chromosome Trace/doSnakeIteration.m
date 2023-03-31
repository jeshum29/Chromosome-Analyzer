function snake = doSnakeIteration(snake, internalForceMatrix, forceImage, im, opts)


% Clamp contour to boundary
snake(:,1) = min(max(snake(:,1),1), size(forceImage,1));
snake(:,2) = min(max(snake(:,2),1), size(forceImage,2));
snake(:,3) = min(max(snake(:,3),1), size(forceImage,3));



% Get image force on the contour points
forceX = opts.imageWeight * interp3(squeeze(forceImage(:,:,:,1)), snake(:,2), snake(:,1), snake(:,3));
forceY = opts.imageWeight * interp3(squeeze(forceImage(:,:,:,2)), snake(:,2), snake(:,1), snake(:,3));
forceZ = opts.imageWeight * interp3(squeeze(forceImage(:,:,:,3)), snake(:,2), snake(:,1), snake(:,3));

forceX(isnan(forceX)) = 0;
forceY(isnan(forceY)) = 0;
forceZ(isnan(forceZ)) = 0;


% pull force at s = 0
ind     = 1;
tanVec  = snake(ind,:) - snake(ind+1,:);
tanVec  = tanVec ./ norm(tanVec);

intensityWeight = interp3(im, snake(ind,2), snake(ind,1), snake(ind,3)) / opts.imMean - opts.pullForceOffset;

pullForce = opts.pullForceWeight * intensityWeight * tanVec;

forceX(ind) = forceX(ind) + pullForce(1);
forceY(ind) = forceY(ind) + pullForce(2);
forceZ(ind) = forceZ(ind) + pullForce(3);

% pull force at s = L
ind     = size(snake,1);
tanVec  = snake(ind,:) - snake(ind-1,:);
tanVec  = tanVec ./ norm(tanVec);

intensityWeight = interp3(im, snake(ind,2), snake(ind,1), snake(ind,3)) / opts.imMean - opts.pullForceOffset;

pullForce = opts.pullForceWeight * intensityWeight * tanVec;

forceX(ind) = forceX(ind) + pullForce(1);
forceY(ind) = forceY(ind) + pullForce(2);
forceZ(ind) = forceZ(ind) + pullForce(3);


% update the positions
snake(:,1) = internalForceMatrix * (snake(:,1) + opts.dt * forceX);
snake(:,2) = internalForceMatrix * (snake(:,2) + opts.dt * forceY);
snake(:,3) = internalForceMatrix * (snake(:,3) + opts.dt * forceZ);


% Clamp contour to boundary
snake(:,1) = min(max(snake(:,1),1), size(forceImage,1));
snake(:,2) = min(max(snake(:,2),1), size(forceImage,2));
snake(:,3) = min(max(snake(:,3),1), size(forceImage,3));



