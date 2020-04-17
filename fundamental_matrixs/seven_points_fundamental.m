function [Fbest,best_inliers]=seven_points_fundamental(source_points,destination_points)
n=size(source_points,2);
max_iterations = inf;
iterations = 2000;
threshold = 3.0;
best_inliers = [];
confidence = 0.99;
for i = 1 : iterations
    indices = randperm(n, 7);
    x1 = source_points(:,indices);
    x2 = destination_points(:,indices);
    F = vgg_F_from_7pts_2img(x1,x2);
    if F == 0
        continue;
    end
    
    if size(F, 1) ~= 3
        continue;
    end    
    for fi = 1 : size(F, 3)
        inliers = [];
        Fi = F(:,:,fi);
        % Estimate the symmetric epipolar distance for each correspondences
        for j = 1 : n 
            l1 = destination_points(:,j)' * Fi;
            l2 = Fi * source_points(:,j);

            l1 = l1 / sqrt(l1(1)^2 + l1(2)^2);
            l2 = l2 / sqrt(l2(1)^2 + l2(2)^2);

            dist = abs(l1 * source_points(:,j)) + abs(destination_points(:,j)'* l2) * 0.5;

            if dist < threshold
                inliers = [inliers j];
            end
        end
        
        if length(best_inliers) < length(inliers)
            % Update inliers of the so-far-the-best model
            best_inliers = inliers;
            Fbest=Fi;
            % Update max iteration number
            max_iterations = log(1 - confidence) / log(1 - (length(best_inliers) / n)^5);
        end
    end
    
    if i > max_iterations
        break;
    end    
end
fprintf('Number of loaded points = %d\n', n);
fprintf('Number of found inliers = %d\n', length(best_inliers));
fprintf('Number of iterations = %d\n', i);