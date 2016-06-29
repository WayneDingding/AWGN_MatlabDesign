function Tausworthe=taus(s0,s1,s2)

 b= floor(bitxor(s0*2^13,s0,'uint64')*2^-19);
 s0=bitand(bitxor((bitand(s0,4294967294,'uint32'))*2^12,b,'uint64'),2^32-1); 
 
 b= floor(bitxor(s1*2^2,s1,'uint64')*2^-25);
 s1=bitand(bitxor((bitand(s1,4294967288,'uint32'))*2^4,b,'uint64'),2^32-1); 
 
 b= floor(bitxor(s2*2^3,s2,'uint64')*2^-11);
 s2=bitand(bitxor((bitand(s2,4294967280,'uint32'))*2^17,b,'uint64'),2^32-1); 
 
 Tausworthe=bitxor(bitxor(s0,s1),s2);
end