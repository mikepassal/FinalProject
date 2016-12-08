function[IsWTarea,pvalarea,isWTdist,pvald,numofperox,]=EvaluatePeroxisome(image)

%% First Goal is to determine if the peroxisomes are WT or not. 
% Can look purely at size distribution and cluster phenomenon . First establish a baseline
% distribution for WT
isitproperformat=isa(image,'numeric');
isitmatrix=ismatrix(image);
msg='Input must be an 2 dimensional numeric matrix';
if isitproperformat=0|isitmatrix=0
    Error(msg)
end 
load('allimages.mat');
[WT5stats,WT5mask]=maskandanalyze(WT5);
[WT6stats,WT6mask]=maskandanalyze(WT6);
[WT7stats,WT7mask]=maskandanalyze(WT7);
WT5area=[WT5stats(:).Area];
WT6area=[WT6stats(:).Area];
WT7area=[WT7stats(:).Area];
allWTarea=horzcat(WT5area,WT6area,WT7area);
% Now we find the average area of our imput image
[imagestats,imagemask]=maskandanalyze(image);
imagearea=[imagestats(:).Area];
[isitwtarea,pvalarea]=ttest2(allWTarea,imagearea,'Vartype','unequal','Alpha',0.001); % If 1, not the same, if 0, the same
if isitwtarea=1
    IsWTarea='Sample has abnormally sized peroxisomes';
elseif isitwtarea=0
    IsWTarea='Sample has normal peroxisomes';
end
%% Now find average distance to nearest peroxisome for WT
WT5cent=extractcent(WT5stats);
[IDXWT5,DWT5]=knnsearch(WT5cent,WT5cent,'K',2);
WT5dist=DWT5(:,2);
WT6cent=extractcent(WT6stats);
[IDXWT6,DWT6]=knnsearch(WT6cent,WT6cent,'K',2);
WT6dist=DWT6(:,2);
WT7cent=extractcent(WT7stats);
[IDXWT7,DWT7]=knnsearch(WT7cent,WT7cent,'K',2);
WT7dist=DWT7(:,2);
allWTdist=vertcat(WT5dist,WT5dist,WT7dist);
imagecent=extractcent(imagestats);
[IDXimg,Dimg]=knnsearch(imagecent,imagecent,'K',2);
imagedist=Dimg(:,2);
[isitwtdist,pvald]=ttest2(allWTdist,imagedist,'Vartype','unequal','Alpha',0.01);
if isitwtdist=1
    IsWTdist='Sample has abnormal clustering';
elseif isitwtdist=0
    IsWTdist='Sample has normal clustering';
end
numofperox=