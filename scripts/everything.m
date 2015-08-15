clear direct

direct{1} = 'D:\thesis\PIVData\cavity\2005_09_14\M030f0000v000a';
direct{2} = 'D:\thesis\PIVData\cavity\2005_09_15\M030f0000v000a';
direct{3} = 'D:\thesis\PIVData\cavity\2005_09_15\M030f0000v000b';
direct{4} = 'D:\thesis\PIVData\cavity\2005_09_15\M030f0000v000c';
direct{5} = 'D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000a';
direct{6} = 'D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000b';
direct{7} = 'D:\thesis\PIVData\cavity\2005_09_19\M030f0000v000';
direct{8} = 'D:\thesis\PIVData\cavity\2005_09_20\M030f1830v400a';
direct{9} = 'D:\thesis\PIVData\cavity\2005_09_20\M030f1830v400b';
direct{10} = 'D:\thesis\PIVData\cavity\2005_09_20\M030f3000v300';
direct{11} = 'D:\thesis\PIVData\cavity\2005_10_14\M030f3000v300';
direct{12} = 'D:\thesis\PIVData\cavity\2005_10_18\M030f1610v250';
direct{13} = 'D:\thesis\PIVData\cavity\2005_10_18\M030f3920v250';
direct{14} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f0000v000b';
direct{15} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f2900v600b';
direct{16} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f3000v200';
direct{17} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f1610v170';
direct{18} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f2800v750';
direct{19} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f3000v150';
direct{20} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f3000v200';
direct{21} = 'D:\thesis\PIVData\cavity\2007_5_4\M030fM1BF4';
direct{22} = 'D:\thesis\PIVData\cavity\2007_5_4\M030fM1BF4comp';
direct{23} = 'D:\thesis\PIVData\cavity\2007_5_4\M030fWNcomp3';
% direct{24} = 'D:\thesis\PIVData\airfoil\Test';
% direct{25} = 'D:\thesis\PIVData\airfoil\Test_2';
% direct{26} = 'D:\thesis\PIVData\airfoil\Test_12';
% direct{27} = 'D:\thesis\PIVData\airfoil\Test_14';
% direct{28} = 'D:\thesis\PIVData\airfoil\Test_24';
% direct{29} = 'D:\thesis\PIVData\airfoil\Test_26';
% direct{30} = 'D:\thesis\PIVData\airfoil\Test_36';
% direct{31} = 'D:\thesis\PIVData\airfoil\Test_38';
% direct{32} = 'D:\thesis\PIVData\airfoil\Test_49';
% direct{33} = 'D:\thesis\PIVData\airfoil\Test_50';
% direct{34} = 'D:\thesis\PIVData\airfoil\Test_51';
% direct{35} = 'D:\thesis\PIVData\airfoil\Test_52';
% direct{36} = 'D:\thesis\PIVData\mixing layer\Forced';
% direct{37} = 'D:\thesis\PIVData\mixing layer\Baseline';
% direct{38} = 'D:\thesis\PIVData\DNS jet\streamwise';
% direct{39} = 'D:\thesis\PIVData\DNS jet\cross_stream1';
% direct{40} = 'D:\thesis\PIVData\DNS jet\cross_stream2';
% direct{41} = 'D:\thesis\PIVData\DNS jet\cross_stream3';
direct{42} = 'D:\thesis\PIVData\cavity\2005_09_14\M030f0000v000b';
direct{43} = 'D:\thesis\PIVData\cavity\2005_09_16\M030f1830v400';
direct{44} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f2900v600a';
% direct{45} = 'D:\thesis\PIVData\airfoil\Test_1';


big_R = [];
big_names = [];
list1 = {'fkm', 'fgmm', 'pkm', 'pgmm'};
list2 = {'tke', 'amp1', 'amp2', 'freq'};
list3 = {'mean', 'std', 'median'};
normal = {'8cluster', '10cluster', '12cluster'};

stuff = cell(41, 1);
for i = 1:length(direct)
    if ~exist([direct{i} filesep 'Scores'], 'dir')
        continue
    end
    if ~exist([direct{i} filesep 'Scores' filesep 'not steady'], 'dir')
        continue;
    end
    ext_direct = [direct{i} filesep 'Scores' filesep 'not steady'];
    
    folders = dir(ext_direct);
    for j = 3:length(folders)
        if ~folders(j).isdir
            continue;
        end
        if any(ismember(folders(j).name, normal))
            continue;
        end
        
        ext_xl_direct = [ext_direct filesep folders(j).name]; 
        files = dir(ext_xl_direct);
        [~, idx] = sort([files.datenum], 'descend');
        files = files(idx);
        
        for k = 1:length(files)
           if files(k).isdir
               continue;
           end
           vars = load([ext_xl_direct filesep files(k).name]);
           [R, names] = merge_struct2(vars.results_scores.relations, list1, list2, list3, 'R');
           big_R = [big_R, R];
           big_names = [big_names, names];
         end
    end
end
    
disp('you win!');

