function distance=fr_segmentDistance(p1,p2,q1,q2)
% Exact minimum of the convex segment-pair quadratic on [0,1]^2.
% Check the interior stationary point and all four edge minima. Independent
% edge solutions avoid the sign/clamping errors in two-pass heuristics.
p1=reshape(p1,1,3);p2=reshape(p2,1,3);q1=reshape(q1,1,3);q2=reshape(q2,1,3);
u=p2-p1;v=q2-q1;w=p1-q1;
a=dot(u,u);b=dot(u,v);c=dot(v,v);d=dot(u,w);e=dot(v,w);
clip=@(x)max(0,min(1,x));
if a<realmin
 if c<realmin;distance=norm(w);else;distance=norm(w-clip(e/c)*v);end
 return
elseif c<realmin
 distance=norm(w+clip(-d/a)*u);return
end
st=[0,clip(e/c);1,clip((e+b)/c);clip(-d/a),0;clip((b-d)/a),1];
D=a*c-b*b;
if D>eps*max(a*c,realmin)
 s=(b*e-c*d)/D;t=(a*e-b*d)/D;
 if s>=0&&s<=1&&t>=0&&t<=1;st(end+1,:)=[s,t];end
end
distance=min(vecnorm(w+st(:,1).*u-st(:,2).*v,2,2));
end
