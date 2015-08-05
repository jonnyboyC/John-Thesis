function symbols = dumb(derp)

symbols = cell(size(derp));

for i = 1:size(derp,1)
    for j = 1:size(derp,2)
        switch derp(i,j)
            case 0
                symbols{i,j} = '0';
            case 1
                symbols{i,j} = 'F_1';
            case 2
                symbols{i,j} = 'F_2';
            case 3
                symbols{i,j} = 'F_3';
            case 4
                symbols{i,j} = 'F_4';
            case 5
                symbols{i,j} = 'FB_3';
            case 6
                symbols{i,j} = 'FB_3';
            case 7
                symbols{i,j} = '0';
            case 8
                symbols{i,j} = 'C_2';
            case 9
                symbols{i,j} = 'C_4';
            case 10
                symbols{i,j} = 'BB_3';
            case 11
                symbols{i,j} = 'BB_4';
            case 12
                symbols{i,j} = 'B_1';
            case 13
                symbols{i,j} = 'B_2';
            case 14
                symbols{i,j} = 'B_3';
            case 15
                symbols{i,j} = 'B_4';
        end
    end
end