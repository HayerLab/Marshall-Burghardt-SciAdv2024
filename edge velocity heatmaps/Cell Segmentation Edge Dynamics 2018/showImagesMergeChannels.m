function showImagesMergeChannels(r,g,b,adjustFlag)

if nargin<4
    adjustFlag=1;
end
if nargin<3
    b=zeros(size(r));
end
if nargin<2
    g=zeros(size(r));
end

im=0*repmat(r,[1 1 3]);

if adjustFlag>0
    im(:,:,1)=mat2gray(r, [prctile(vect(r),0.5) prctile(vect(r),99.5)]);
    im(:,:,2)=mat2gray(g, [prctile(vect(g),0.5) prctile(vect(g),99.5)]);
    im(:,:,3)=mat2gray(b, [prctile(vect(b),0.5) prctile(vect(b),99.5)]);
else
    im(:,:,1)=r;
    im(:,:,2)=g;
    im(:,:,3)=b;
end

imshow(im);