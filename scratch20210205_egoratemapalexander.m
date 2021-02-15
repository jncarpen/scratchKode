load('D:\Data\Dataset\sample data\fromPablo\Archive\Example1_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example2_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example3_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example4_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example5_RawData.mat')

% arrange dataset
Jercog(1).P = Cell1_Trajectory;
Jercog(2).P = Cell2_Trajectory;
Jercog(3).P = Cell3_Trajectory;
Jercog(4).P = Cell4_Trajectory;
Jercog(5).P = Cell5_Trajectory;

Jercog(1).ST = Cell1_Spikes;
Jercog(2).ST = Cell2_Spikes;
Jercog(3).ST = Cell3_Spikes;
Jercog(4).ST = Cell4_Spikes;
Jercog(5).ST = Cell5_Spikes;

Jercog(1).ST = Cell1_Spikes(:,1);
Jercog(2).ST = Cell2_Spikes(:,1);
Jercog(3).ST = Cell3_Spikes(:,1);
Jercog(4).ST = Cell4_Spikes(:,1);
Jercog(5).ST = Cell5_Spikes(:,1);

for i = 1%:5
    root.x = Jercog(i).P(:,2);
    root.y = Jercog(i).P(:,3);
    root.md = get_MD(Jercog(i).P)';
    root.ts = Jercog(i).P(:,1);
    
    [train, ~] = binSpikes(Jercog(i).P(:,1), Jercog(i).ST);
    root.spike = train';
    
    [out] = EgocentricRatemap(root);
    plotEBC(root, out, i)
end



