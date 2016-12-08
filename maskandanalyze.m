function[imagestats,maskedimage]=maskandanalyze(image)

doubleimage=im2double(image);
dilateimage=imdilate(doubleimage,strel('disk',15));
normimage=doubleimage./dilateimage;
maskimage=normimage>.3;
maskimage_2=imerode(maskimage,strel('disk',2));
maskimage_3=imdilate(maskimage_2,strel('disk',2));
imagestats=regionprops(maskimage_3,'Area','centroid','Eccentricity');
maskedimage=maskimage_3;
