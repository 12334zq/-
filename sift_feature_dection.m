addpath(genpath('matchSIFT'));
addpath(genpath('data'));
addpath(genpath('fundamental_matrixs'));
addpath(genpath('homography'));
%∂¡Õº
image_i=J1;
image_j=J2;
% image_i=imread('KyotoA.png');
% image_j=imread('KyotoB.png');
%siftÃÿ’˜ºÏ≤‚
edge_thresh = 2.5;
[SIFTloc_i,SIFTdes_i] = vl_sift(single(rgb2gray(image_i)), 'edgethresh', edge_thresh) ;
[SIFTloc_j,SIFTdes_j] = vl_sift(single(rgb2gray(image_j)), 'edgethresh', edge_thresh) ;
%sift match
[matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j);
pts1=SIFTloc_i(1:2,matchPointsID_i)';%x,y
pts2=SIFTloc_j(1:2,matchPointsID_j)';
scale1=SIFTloc_i(3,matchPointsID_i)';
scale2=SIFTloc_j(3,matchPointsID_j)';
alpha1=SIFTloc_i(4,matchPointsID_i)';
alpha2=SIFTloc_j(4,matchPointsID_j)';
image = [image_i image_j];
figure,imshow(image)
hold on;
for i = 1 : length(pts1)
    color = rand(1,3);
    
    plot([pts1(i, 1), size(image_i, 2) + pts2(i, 1)], [pts1(i, 2) pts2(i, 2)], 'Color', color)
    
    scatter(pts1(i, 1), pts1(i, 2), 40, color, 'filled');
    scatter(size(image_i, 2) + pts2(i, 1), pts2(i, 2), 40, color,'filled');
end
%º∆À„µ•”¶æÿ’ÛDLT
ones1=ones(size(pts1,1),1);
M=[pts1,ones1,pts2,ones1,scale1,scale2,alpha1,alpha2];
H1 = normalizedDLT(pts1, pts2);
%%º∆À„≤–≤Ó
source_points=[pts1,ones1]';
destination_points=[pts2,ones1]';
pts2_t = H1 * source_points;
pts2_t = rdivide(pts2_t, pts2_t(3,:));
residuals = vecnorm(destination_points - pts2_t);
re_errors=mean(residuals);
%%2points
[H2,best_inliers,residuals2]= HomographyFromSIfT(M,image_i,image_j);
re_errors2=mean(residuals2);


%%ª˘¥°æÿ’Û
data=[pts1,pts2,alpha1-alpha2];
[F1,best_inliers1]=five_point_fundamental(data);
[F2,best_inliers2]=seven_points_fundamental(source_points,destination_points);
[F3, inliersIndex, status ]= estimateFundamentalMatrix(pts1, pts2);
%%º´œﬂ
  % Show the inliers in the first image.
 
  figure; 
  subplot(121); imshow(image_i); 
  title('Inliers and Epipolar Lines in First Image'); hold on;
  plot(pts1(inliersIndex,1),pts1(inliersIndex,2), 'go')
  % Compute the epipolar lines in the first image.
  epiLines = epipolarLine(F3', pts2(inliersIndex, :));
  % Compute the intersection points of the lines and the image border.
  pts = lineToBorderPoints(epiLines, size(image_i));

  % Show the epipolar lines in the first image
  line(pts(:, [1,3])', pts(:, [2,4])');

  % Show the inliers in the second image.
  
  subplot(122); imshow(image_j);
  title('Inliers and Epipolar Lines in Second Image'); hold on;
  plot(pts2(inliersIndex,1), pts2(inliersIndex,2), 'go')

  % Compute and show the epipolar lines in the second image.
  epiLines = epipolarLine(F3, pts1(inliersIndex, :));
  pts = lineToBorderPoints(epiLines, size(image_j));
  line(pts(:, [1,3])', pts(:, [2,4])');
%%
