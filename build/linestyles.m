linetype   = {'-.','-'};
markertype = {'x','*','+','o','s','d','v','^'};
colortype  = {'b','g','m','c','k'};

ls_tab = {};
for m=1:length(markertype)
	for c=1:length(colortype)
		ls_tab{end+1} = strcat(linetype{1},markertype{m},colortype{c});
		ls_tab{end+1} = strcat(linetype{2},markertype{m},colortype{c});
	end
end
ls_idx=1; %line-style index
