function [H,maxrow,maxcol] = ExpandH(Hc,Z);
% function [H,maxrow,maxcol] = ExpandH(Hc,Z);
%
% Expand compressed Hc matrix into full-size H matrix based on scalar size Z
% also return maximum number of 1's in a row and in a column as maxrow and maxcol

maxcol = max(sum(Hc >= 0));
maxrow = max(sum(Hc'>= 0));

H = zeros( Z*size(Hc) );
for row=1:size(Hc,1)
	for col=1:size(Hc,2)
		if (Hc(row,col)>-1)
			H((row-1)*Z+1:row*Z,(col-1)*Z+1:col*Z) = circshift(eye(Z),-Hc(row,col));
		end
	end
end
