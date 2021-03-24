function quadraturePoints = getintergrate(nInt)
[xi, w] = GaussQuad(nInt,1);
lxi = length(xi);
quadraturePoints = zeros(lxi*lxi,3);
for i = 1:lxi
    for j = 1:lxi
        n = (i-1)*lxi+j;
        quadraturePoints(n,:)=[xi(i),xi(j),w(i)*w(j)];
    end
end
end