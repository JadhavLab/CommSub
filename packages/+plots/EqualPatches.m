

% this function creates equally sized patches of specified numbers

ax = gca;
y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
patch_color = ["blue", "red", "green"];

for n = 1:6
    x = [(n-1)*5+0.5 n*5+0.5 n*5+0.5 (n-1)*5+0.5];
    p = patch(x,y,patch_color(int16(ceil(n/2))));
    
    if mod(n,2)~=0
        p.FaceAlpha=0.3;
    else
        p.FaceAlpha=0.1;
        
    end
    legend('off'); 
end