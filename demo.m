[filename,pathname]=uigetfile({'*.jpg'},'choose the picture');
str=[pathname, filename];

img = imread(str);
%img = rgb2gray(img);

img = im2double(img);
scale = 1.0;
img=imresize(img,scale);

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);
st = clock;
[boundary,label] = SCoW(R,G,B,0.5,144);
fprintf(' took %.5f second\n',etime(clock,st));

R(boundary>0)=255;
G(boundary>0)=0;
B(boundary>0)=0;

img(:,:,1) = R;
img(:,:,2) = G;
img(:,:,3) = B;

imshow(img,[]);
