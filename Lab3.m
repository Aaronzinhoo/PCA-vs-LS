load data_lab3
figure
imshow(im(:,:,[30,20,7]),[]) % for v2016+
len = 101^2;
img_vec = zeros(len,195);
%%
%1 PCA, Eigen Plots, PCA Scores
for i=1:195
    img_vec(:,i) = reshape(im(:,:,i),len,1);
end

img_vec = img_vec';
sample_mean = mean(img_vec,2);
size(img_vec)
Covar = cov(img_vec');
[eig_vec, eig_val] = eigs(Covar,194);
[seig_vec, seig_val] = eigs(Covar,1,'SM');
eig_vals = diag(eig_val);
eig_vals(end+1) = seig_val;
eig_vec(:,end+1) = seig_vec;

figure
plot(eig_vals(1:5,1))
title('Zoomed in Eigenvalue Plot')
xlabel('x')
ylabel('Eigen Value')
figure
plot(eig_vals(1:20,1))
title('Eigenvalue Plot')
xlabel('x')
ylabel('Eigen Value')
%%
% PCs are the col of XX' he cov matrix
% first few componets
PCA_score_first = zeros(101^2,8);
figure
for i=1:5
    subplot(2,3,i); 
    PCA_score_first(:,i) = img_vec'*eig_vec(:,i);
    imshow(mat2gray(reshape(PCA_score_first(:,i),101,101)))
    title([num2str(i), ' PCA Score'])
end

%Images of the first few PCA scores actually represent the image itself
%with segregation of the the different materials it looks like. All three main
%materials are differentiable within the first couple scores just like they
%are in the original image. The first couple PCAs are very similar to gray
%images of the actual picture but highlight different materials.
%After the first couple PCA scores though the images become less useful in their ability to
%discern materials. This implies that our data might lie in a low dimensional
%subspace. This follows from the eigenvalue plot which shows that the first
%two have the most varience while after that the value drops significantly. 

%%
% last few componets
figure
PCA_score_last = zeros(101^2,8);
for i=1:8
    subplot(2,4,i); 
    PCA_score_last(:,i) = img_vec'*eig_vec(:,186+i);
    imshow(mat2gray(reshape(PCA_score_last(:,i),101,101))) % for v2016+
    %imagesc(reshape(PCA_score_last(:,i),101,101))
    title([num2str(186+i), ' PCA Score'])
end
% The last few componets are uniterpretable since they almost explain no
% data at all. The images resemble white noise or static which seems corect
% since these points probably account for noise in the data rather than
% signals

%2. 
% The data leads me to believe the dimensionality of the image is really
% low around where the eigenvalues seem to drop in value. Since the first
% couple of eigenvalues accoutn for most of the varience it seems like the
% 2-D or 3-D subspace is where the data lies theoretically. Also only the
% first PCA scores account for a majority of the image's properties


%%
%3. Display of 2D and 3D scatterplots
figure
subplot(1,2,1);scatter(PCA_score_first(:,1), PCA_score_first(:,2), '.')
xlabel('1st PCA Score');
ylabel('2nd PCA Score');
title('2-D Scatter Plot of PCA Scores')
subplot(1,2,2);scatter3(PCA_score_first(:,1), PCA_score_first(:,2), PCA_score_first(:,3), '.')
xlabel('1st PCA Score');
ylabel('2nd PCA Score');
zlabel('3rd PCA Score');
title('3-D Scatter Plot of PCA Scores')

%3.
%Yes since the data can lie in the 2-D while explaining much of the data.
%The data seems to group into 2 different groups but these groups are
%already well defined on a 2-D subspace. There could be a third cluster
%which there should be but it is hard to tell. But considering the rest of
%the data it seems as though it will be only 2.

%%
%4. Display of 2D and 3D scatterplots
figure
for i=1:6
    % negative may or may not be needed.. just check graph
    subplot(2,3,i);plot((-1*eig_vec(:,i)))
    xlabel('x');
    ylabel('PC value');
    title([num2str(i), ' Evec'])
end
figure
for i=4:6
    subplot(1,3,i-3);plot(S(:,i-3))
    xlabel('x');
    ylabel('S value');
    title([num2str(i-3), ' S vector'])
end

%4
% In terms of viewing the eigenvectors we can see that after the 3rd PC the
% values oscillate around zero meaning they dont offer much in terms of
% varience. Although the values arent that similar the first three PCs 
% resemble the vectors of S but it seems S was able to capture the data much
% more accurately.
%%
%5. If eq. 3 holds then the theoretical subspace of the model will exist in
%3 dimensions since given the S matrix is a three dimsensional matrix that
%is a basis. 

%6. During an affine combination we restrict the linear combination of data
%to be one making the last term of the data dependent on the first two.
%This essentially mean that we reduce the space of our data to only a two
%dimensional subspace now.

%7. Once the positivity constraints are addded to the sum-to-one
%constraints we restrict the points still to a 2 dimensional subspace of
%only the first quadrant (since all values must be > 0). The dimension
%doesnt change from the previous answer just the area that we cover.

%8 Basically from this we can conclude that the true dimensionality of the 
%data provided is 2 since our data physically is trapped within two dimensions
% by the abundance with the third being decided by the first two componets


%9.
for i=1:3
    PCA_S1(:,i) = S(:,1)'*eig_vec(:,i); 
    PCA_S2(:,i) = S(:,2)'*eig_vec(:,i); 
    PCA_S3(:,i) = S(:,3)'*eig_vec(:,i); 
end
figure
hold on
scatter(PCA_S1(:,1), -PCA_S1(:,2),'filled', '^')
scatter(PCA_S2(:,1), -PCA_S2(:,2),'filled', '^')
scatter(PCA_S3(:,1), -PCA_S3(:,2),'filled', '^')
title('2-D Scatter Plot of PCAs')
xlabel('1st PC');
ylabel('2nd PC');
scatter(PCA_score_first(:,1), -PCA_score_first(:,2),'.')
title('2-D Scatter Plot of PCAs')
xlabel('1st PC');
ylabel('2nd PC');
hold off

figure
hold on
scatter3(PCA_S1(:,1), -PCA_S1(:,2), PCA_S1(:,3),'filled','^')
scatter3(PCA_S2(:,1), -PCA_S2(:,2), PCA_S2(:,3),'filled', '^')
scatter3(PCA_S3(:,1), -PCA_S3(:,2), PCA_S3(:,3) ,'filled', '^')
scatter3(PCA_score_first(:,1), -PCA_score_first(:,2), PCA_score_first(:,3), '.')
xlabel('1st PC');
ylabel('2nd PC');
zlabel('3rd PC');
title('3-D Scatter Plot of PCAs')
hold off

%9.
%The position of the makers lie on the 2-D space at the critical points of
%the plotted PCA scores when decomposed into the eigenvector space. 
%But when expanded into a 3D space the new points actually lie almost on a
%plane it seems. This means that S can approx be represented in 2D which
%follows what we previously thought from the questions before. That our
%data could be represented in 2D. Thus is it an appropriate basis to use.
%%
%10 get A without constraint
A = inv(S'*S)*S'*img_vec;
figure
for i=1:3
   subplot(1,3,i); imshow(reshape(A(i,:),101,101),[])
   title(sprintf('Material %d', i))
end
% The images correspond to the materials that each vector maps out in the
% original image. We see highlighted areas where the materials exist while
% the other materials are left a dull dark color meaning the abundance of the
%the material in this region is not apparent. Each image has different
% portions of the image that has color since they each correspond to one
% material. The edges are well defined for where one material is adjacent
% to another. This is similar to the PCA score images that we plotted earlier
% but each column of A corresponds to one material and its abundance in
% that pixel when you reshape it to an image. Obviously the last image maps
% water, but its interesting that the edges of the land stick out a little
% and that the water itself doesnt stick out as much as the other materials
% in their respective images.

%%
% 11 Now compute it with the constraint that A col must sum to one
% n = 10;
% timeForEachIteration = zeros(1, n);
% for k = 1:n
%     tic
%     Z = S'*S;
%     BLK = [Z ones(3,1);ones(1,3) 0];
%     BLK_inv = inv(BLK);
%     F = S'*img_vec;
%     G = ones(1,101^2);
%     FG = [F;G];
%     A_constr = BLK_inv(1:3,:)*FG;
%     L = BLK_inv(4,:)*FG;
%     timeForEachIteration(k) = toc;
% end
% figure
% hold on
% plot(timeForEachIteration)
tic
Z = S'*S;
BLK = [Z ones(3,1);ones(1,3) 0];
BLK_inv = inv(BLK);
F = S'*img_vec;
G = ones(1,101^2);
FG = [F;G];
A_constr = BLK_inv(1:3,:)*FG;
L = BLK_inv(4,:)*FG;
toc
figure
for i=1:3
   subplot(1,3,i);imshow(reshape(A_constr(i,:),101,101), [])
   title(sprintf('Material Constr %d', i))
   %imshow()
end

%11 for the new abundancies we get the first two materials similar to the
%first calulation of A unconstrained. The third material sticks out much
%more than the previous calculation. The edges are much more clear and
%defined than before where you could see the general outline of the water
%but not as strict. The edges that were on the land previously have also
%dissapered.

%%
%12 do both computations 
% recalc of 11 with reduced dim
% n = 10;
% timeForEachIteration = zeros(1, n);
% for k = 1:n
%     tic
%     S_hat = eig_vec(:,1:2)'*S;
%     Z = S_hat'*S_hat;
%     BLK = [Z ones(3,1);ones(1,3) 0];
%     BLK_inv = inv(BLK);
%     F = S_hat'*(PCA_score_first(:,1:2)');
%     G = ones(1,101^2);
%     FG = [F;G];
%     A_reduced = BLK_inv(1:3,:)*FG;
%     timeForEachIteration(k) = toc;
% end

% plot(timeForEachIteration)
% title('Reduced vs Regular Constraint')
% xlabel('Run')
% ylabel('Time (using tic toc)')
% hold off
tic
S_hat = eig_vec(:,1:2)'*S;
Z = S_hat'*S_hat;
BLK = [Z ones(3,1);ones(1,3) 0];
BLK_inv = inv(BLK);
F = S_hat'*(PCA_score_first(:,1:2)');
G = ones(1,101^2);
FG = [F;G];
A_reduced = BLK_inv(1:3,:)*FG;
toc
figure
for i=1:3 
   subplot(1,3,i);imshow((reshape(A_reduced(i,:),101,101)), [])
   title(sprintf('Constrained and Reduced %d', i))
   %imshow()
end

%The amount of calculations for 11 without dimension reduction is more
%heavy than what is needed for the dimension reduction version of 11. The
%problem that we face here normally is loss of data, but in this case the
%recreated images of A reduced calculation is actually very simialr to the
%original. Meaning that we get practically the same result, visually, with
%an almost .003 second decrease in computation time.

%%
%13 
RMSE=zeros(1,101^2);
L = 1/sqrt(195);
for i=1:101^2
   RMSE(1,i) = L*norm(img_vec(:,i) - S*A(:,i));
end
figure
subplot(1,2,1);imshow((reshape(RMSE(1,:),101,101)),[])
title('RMSE for 10')

for i=1:101^2
   RMSE_c(1,i) = L*norm(img_vec(:,i) - S*A_constr(:,i));
end
subplot(1,2,2);imshow((reshape(RMSE_c(1,:),101,101)),[])
title('RMSE for 11')

for i=1:101^2
   RMSE_r(1,i) = L*norm(img_vec(:,i) - S*A_reduced(:,i));
end
subplot(1,2,2);imshow((reshape(RMSE_c(1,:),101,101)),[])
title('RMSE for 11 reduced')
%%
image = reshape((S*A)',101,101,195);
image_c = reshape((S*A_constr)',101,101,195);
figure
imshow(image(:,:,[30,20,7]),[])
figure
imshow(image_c(:,:,[30,20,7]),[])

%%
mean(RMSE)
mean(RMSE_c)
mean(RMSE_r)

% The highest value of the RMSE for the constrained case could be
% attributed to the fact that constrained case makes assumptions on the
% data that shouldnt be made such as there existing only 3 materials and they 
%create a porportion that equal to one. In some cases this may weight some
%materials higher than they should be which may account for some of the
%error witnessed. Although the majority of the materials seem to be the ones 
%accounted for the edge between the river and land seem to have a lot of error. This
% could be because of other materials that exist here that the constrained
% version doesnt account well for. The best value is given by the constrained
% case also though. The mean of the RMSE for both though are pretty low but the image
% for the constrained case seems to dislpay where the materials exist a lot better. 
% which spectral unmixing is trying to accomplish. When you recreate the original 
% image with these S and As we see that the S*A unconstrained is much more
% accuarte than the constrained version, but the constrained version accoutns for the 
%edges of each material in a more definite way. The way the land is segregated
% could also come into play. Since the land is in blocks the porportions
% may not need to be enforced as much since the materials are all
% segregated.

