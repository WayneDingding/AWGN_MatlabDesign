%Leading Zero Detector (48bits)
%detect zeros per 8bits
%if 8bits are all zeros, jump into LZD_8 
function nzeros=LZD(num)
nzeros=zeros(size(num));
i=size(num,1);
while i 
    temp=48;
    numt=num(i);
    while temp
        if 1-numt/2^(temp-8)>0
           temp=temp-8;
           nzeros(i)=nzeros(i)+8;
        else
            nzeros(i)=nzeros(i)+LZD_8(numt/2^(temp-8));
            temp=0; 
        end 
    end
    i=i-1;
end

end

