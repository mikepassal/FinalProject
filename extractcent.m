function fixedcent=extractcent(struct)
structcent=[struct.Centroid];
structcentroidlength=length(structcent);
h=1;
for zz=1:structcentroidlength
    if mod(zz,2)==0
        fixedcent(h,2)=structcent(1,zz);
        h=h+1;
    else
        fixedcent(h,1)=structcent(1,zz);
    end
end