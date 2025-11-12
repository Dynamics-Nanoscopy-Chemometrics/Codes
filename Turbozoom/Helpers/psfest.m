function psf=psfest(lng,sigma)

u=-lng:lng;
v=exp(-(u'/sigma).^ 2/2);
v=v/sum(v);
psf=zeros(lng,lng);
nv=length(v);

for i=1:lng
    
    l=max(i-lng,1);
    u=min(i+lng,lng);
    l2=max(l-i+lng+1,1);
    u2=min(u-i+lng+1,nv);
    psf(i,l:u)=v(l2:u2);
    
end