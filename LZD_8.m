%Leading Zero Detector (8bits)
function nzeros=LZD_8(num)
nzeros=0;    
temp=8;
while temp
    if 1-num/2^(temp-1)>0
       temp=temp-1;
       nzeros=nzeros+1;
    else
       temp=0; 
    end 
end
end
