function [Hc] = lift_matrix(orig_Hc,liftr,liftc)

%take original compressed matrix orig_Hc
%perform lifting where liftr and liftc are arrays with
%the row and column indeces in which the lifting matrix is 1
%
%return Hc, the lifted resulting compressed matrix

[r c] = size(orig_Hc);
oddr  = reshape([orig_Hc; -ones(r,c)],r,[]);
evenr = reshape([-ones(r,c); orig_Hc],r,[]);
Hc    = transpose( reshape( [transpose(oddr);transpose(evenr)],2*c,[]));

for i=1:length(liftr)
	row = liftr(i)*2-1;
	col = liftc(i)*2-1;
	tmp = Hc(row,col:col+1);
	Hc(row,col:col+1) = Hc(row+1,col:col+1);
	Hc(row+1,col:col+1) = tmp;
end
