
clear
clc

ftno = 5120;
cellno = 2595;
pert = 0.10;
%ftno = 2560;

%deleting the reclass folders
epath = './relabel/new_class/Early';
lpath = './relabel/new_class/Late';
vpath = './relabel/new_class/Viable';
eee = dir(epath);
lll = dir(lpath);
vvv = dir(vpath);

for fff = 1:length(eee)
    ffnam = eee(fff).name;
    delete(strcat(epath, '/' , ffnam))
end

for fff = 1:length(lll)
    ffnam = lll(fff).name;
    delete(strcat(lpath , '/' , ffnam))
end

for fff = 1:length(vvv)
    ffnam = vvv(fff).name;
    delete(strcat(vpath , '/' , ffnam))
end


supath = './relabel/supp';
nosupath = './relabel/nosupp';

%counting the support cells to preallocate for speed
celxls = dir(supath);
ssupclcnt = 0;
for fll = 3:height(celxls)
%
%
    supname = celxls(fll).name;
    suptable = table2cell(readtable(supname));
    clsiz=height(suptable);
    ssupclcnt = ssupclcnt + clsiz;
%
    ssupportcells( ssupclcnt-clsiz+1 : ssupclcnt ,: ) = suptable;
%
%
end


nocelxls = dir(nosupath);
nnosupclcnt=0;
for fll = 3:height(nocelxls)
%
%
    nosupname = nocelxls(fll).name;
    nosuptable = table2cell(readtable(nosupname));
    clsiz=height(nosuptable);
    nnosupclcnt = nnosupclcnt + clsiz;
%
    nnosupportcells( nnosupclcnt-clsiz+1 : nnosupclcnt ,: ) = nosuptable;
%
%
end


%counting cells in each support cell class
ssupcnt0=0;
ssupcnt1=0;
ssupcnt2=0;
for fll=1:ssupclcnt

    class=ssupportcells(fll,2);
    class=cell2mat(class);

    if class==0
        ssupcnt0=ssupcnt0+1;
    end

    if class==1
        ssupcnt1=ssupcnt1+1;
    end

    if class==2
        ssupcnt2=ssupcnt2+1;
    end

end

%counting cells in each nonsupport cell class
nnosupcnt0=0;
nnosupcnt1=0;
nnosupcnt2=0;
for fll=1:nnosupclcnt

    class=nnosupportcells(fll,2);
    class=cell2mat(class);

    if class==0
        nnosupcnt0=nnosupcnt0+1;
    end

    if class==1
        nnosupcnt1=nnosupcnt1+1;
    end

    if class==2
        nnosupcnt2=nnosupcnt2+1;
    end

end


%getting index of each support cell of each class
supinds0=zeros(ssupcnt0,1);
cnt0=0;
supinds1=zeros(ssupcnt1,1);
cnt1=0;
supinds2=zeros(ssupcnt2,1);
cnt2=0;
for fll=1:ssupclcnt

    class=ssupportcells(fll,2);
    class=cell2mat(class);

    if class==0
        cnt0=cnt0+1;
        supinds0(cnt0)=fll;
    end

    if class==1
        cnt1=cnt1+1;
        supinds1(cnt1)=fll;
    end

    if class==2
        cnt2=cnt2+1;
        supinds2(cnt2)=fll;
    end

end

%getting index of each nonsupport cell of each class
nnosupinds0=zeros(nnosupcnt0,1);
cnt0=0;
nnosupinds1=zeros(nnosupcnt1,1);
cnt1=0;
nnosupinds2=zeros(nnosupcnt2,1);
cnt2=0;
for fll=1:nnosupclcnt

    class=nnosupportcells(fll,2);
    class=cell2mat(class);

    if class==0
        cnt0=cnt0+1;
        nnosupinds0(cnt0)=fll;
    end

    if class==1
        cnt1=cnt1+1;
        nnosupinds1(cnt1)=fll;
    end

    if class==2
        cnt2=cnt2+1;
        nnosupinds2(cnt2)=fll;
    end

end

ssupclasses=cell2mat(ssupportcells(:,2));
ssupfeatures=cell2mat(ssupportcells(:,3:end));

ssuppaths0 = ssupportcells(supinds0);
ssupclasses0=ssupclasses(supinds0);
ssupfeatures0=ssupfeatures(supinds0,:);
ssupcells0=[ssuppaths0,num2cell(ssupclasses0),num2cell(ssupfeatures0)];

ssuppaths1 = ssupportcells(supinds1);
ssupclasses1=ssupclasses(supinds1);
ssupfeatures1=ssupfeatures(supinds1,:);
ssupcells1=[ssuppaths1,num2cell(ssupclasses1),num2cell(ssupfeatures1)];

ssuppaths2 = ssupportcells(supinds2);
ssupclasses2=ssupclasses(supinds2);
ssupfeatures2=ssupfeatures(supinds2,:);
ssupcells2=[ssuppaths2,num2cell(ssupclasses2),num2cell(ssupfeatures2)];

nnosupclasses=cell2mat(nnosupportcells(:,2));
nnosupfeatures=cell2mat(nnosupportcells(:,3:end));

nnonsuppaths0=nnosupportcells(nnosupinds0);
nnonsupclasses0=nnosupclasses(nnosupinds0);
nnonsupfeats0=nnosupfeatures(nnosupinds0,:);
nnonsupcells0=[nnonsuppaths0,num2cell(nnonsupclasses0),num2cell(nnonsupfeats0)];

nnonsuppaths1=nnosupportcells(nnosupinds1);
nnonsupclasses1=nnosupclasses(nnosupinds1);
nnonsupfeats1=nnosupfeatures(nnosupinds1,:);
nnonsupcells1=[nnonsuppaths1,num2cell(nnonsupclasses1),num2cell(nnonsupfeats1)];

nnonsuppaths2=nnosupportcells(nnosupinds2);
nnonsupclasses2=nnosupclasses(nnosupinds2);
nnonsupfeats2=nnosupfeatures(nnosupinds2,:);
nnonsupcells2=[nnonsuppaths2,num2cell(nnonsupclasses2),num2cell(nnonsupfeats2)];

cells0=[[ssuppaths0; nnonsuppaths0],[num2cell(ssupclasses0); num2cell(nnonsupclasses0)],[num2cell(ssupfeatures0); num2cell(nnonsupfeats0)]];
cells1=[[ssuppaths1; nnonsuppaths1],[num2cell(ssupclasses1); num2cell(nnonsupclasses1)],[num2cell(ssupfeatures1); num2cell(nnonsupfeats1)]];
cells2=[[ssuppaths2; nnonsuppaths2],[num2cell(ssupclasses2); num2cell(nnonsupclasses2)],[num2cell(ssupfeatures2); num2cell(nnonsupfeats2)]];
cellsclasses0=[ssupclasses0; nnonsupclasses0];
cellsclasses1=[ssupclasses1; nnonsupclasses1];
cellsclasses2=[ssupclasses2; nnonsupclasses2];

%if grouping in non support cells
%cellsfeat0=[ssupfeatures0;nnonsupfeats0];
cellsfeat0=[ssupfeatures0;zeros(size(nnonsupfeats0,1),ftno)];
%if grouping in non support cells
%cellsfeat1=[ssupfeatures1;nnonsupfeats1];
cellsfeat1=[ssupfeatures1;zeros(size(nnonsupfeats1,1),ftno)];
%if grouping in non support cells
%cellsfeat2=[ssupfeatures2;nnonsupfeats2];
cellsfeat2=[ssupfeatures2;zeros(size(nnonsupfeats2,1),ftno)];

cellcnt0=length(ssuppaths0)+length(nnonsuppaths0);
cellcnt1=length(ssuppaths1)+length(nnonsuppaths1);
cellcnt2=length(ssuppaths2)+length(nnonsuppaths2);


%first pearson correlation coefficient matrix for each class
%zscores of support cell groups
ssupzscore0 = zeros(cellcnt0, ftno);
for pcf = 1:cellcnt0
%
%
    afeat = cellsfeat0(pcf,: );
    zscor = (afeat - mean(afeat))/std(afeat);
    ssupzscore0(pcf,: ) = zscor;
%
%
end

ssupzscore1 = zeros(cellcnt1, ftno);
for pcf = 1:cellcnt1
%
%
    afeat = cellsfeat1(pcf,: );
    zscor = (afeat - mean(afeat))/std(afeat);
    ssupzscore1(pcf,: ) = zscor;
%
%
end

ssupzscore2 = zeros(cellcnt2, ftno);
for pcf = 1:cellcnt2
%
%
    afeat = cellsfeat2(pcf,: );
    zscor = (afeat - mean(afeat))/std(afeat);
    ssupzscore2(pcf,: ) = zscor;
%
%
end

corrmat0 = (1/(ftno-1))*(ssupzscore0*ssupzscore0');
corrmat1 = (1/(ftno-1))*(ssupzscore1*ssupzscore1');
corrmat2 = (1/(ftno-1))*(ssupzscore2*ssupzscore2');

%correlation density calculation
corrheav0 = heaviside(corrmat0 - 0.5);
corrheav1 = heaviside(corrmat1 - 0.5);
corrheav2 = heaviside(corrmat2 - 0.5);

%with non supports grouped in
%corrden0 = sum(corrheav0);
%corrden1 = sum(corrheav1);
%corrden2 = sum(corrheav2);

%without non supports grouped in
corrden0 = sum(corrheav0(1:ssupcnt0,1:ssupcnt0));
corrden1 = sum(corrheav1(1:ssupcnt1,1:ssupcnt1));
corrden2 = sum(corrheav2(1:ssupcnt2,1:ssupcnt2));

%enhanced correlation coefficient
%with non supports grouped in
%corcor0 = sum(corrden0.*corrmat0);
%corcor1 = sum(corrden1.*corrmat1);
%corcor2 = sum(corrden2.*corrmat2);

%without nonsupports grouped in
corcor0 = [sum(corrden0.*corrmat0(1:ssupcnt0,1:ssupcnt0)),zeros(1,nnosupcnt0)];
corcor1 = [sum(corrden1.*corrmat1(1:ssupcnt1,1:ssupcnt1)),zeros(1,nnosupcnt1)];
corcor2 = [sum(corrden2.*corrmat2(1:ssupcnt2,1:ssupcnt2)),zeros(1,nnosupcnt2)];

corcor0=corcor0/max(corcor0);
corcor1=corcor1/max(corcor1);
corcor2=corcor2/max(corcor2);


supportsterm0 = false(1, cellcnt0);
for crr = 1:cellcnt0
%terminal support cell groups
    celco0 = corcor0(crr);

    if celco0 >= pert
        supportsterm0(crr) = celco0;
    end
%
%
end

supportsterm1 = false(1, cellcnt1);
for crr = 1:cellcnt1
%terminal support cell groups
    celco1 = corcor1(crr);

    if celco1 >= pert
        supportsterm1(crr) = celco1;
    end
%
%
end

supportsterm2 = false(1, cellcnt2);
for crr = 1:cellcnt2
%terminal support cell groups
    celco2 = corcor2(crr);

    if celco2 >= pert
        supportsterm2(crr) = celco2;
    end
%
%
end

%organizing terminal support cell list in terms of paths, classes, features
suptermpaths0 = cells0(supportsterm0)';
suptermclasses0=cellsclasses0(supportsterm0);
suptermfeat0=cellsfeat0(supportsterm0,:);
suptermcells0=[suptermpaths0,num2cell(suptermclasses0),num2cell(suptermfeat0)];

suptermpaths1 = cells1(supportsterm1)';
suptermclasses1=cellsclasses1(supportsterm1);
suptermfeat1=cellsfeat1(supportsterm1,:);
suptermcells1=[suptermpaths1,num2cell(suptermclasses1),num2cell(suptermfeat1)];

suptermpaths2 = cells2(supportsterm2)';
suptermclasses2=cellsclasses2(supportsterm2);
suptermfeat2=cellsfeat2(supportsterm2,:);
suptermcells2=[suptermpaths2,num2cell(suptermclasses2),num2cell(suptermfeat2)];


%adding cells left out of terminal support cell list to non support cells
nonsup0=~supportsterm0;
nonsup1=~supportsterm1;
nonsup2=~supportsterm2;

nonsuppaths0=cells0(nonsup0)';
nonsupclasses0=cellsclasses0(nonsup0);
nonsupfeats0=cellsfeat0(nonsup0,:);

nonsuppaths1=cells1(nonsup1)';
nonsupclasses1=cellsclasses1(nonsup1);
nonsupfeats1=cellsfeat1(nonsup1,:);

nonsuppaths2=cells2(nonsup2)';
nonsupclasses2=cellsclasses2(nonsup2);
nonsupfeats2=cellsfeat2(nonsup2,:);

nosuppaths = [nonsuppaths0; nonsuppaths1; nonsuppaths2];
nosupclasses = [nonsupclasses0; nonsupclasses1; nonsupclasses2];
nosupfeatures = [nonsupfeats0; nonsupfeats1; nonsupfeats2];
nosupcells = [nosuppaths, num2cell(nosupclasses), num2cell(nosupfeatures)];
nosupcnt=length(nosuppaths);


%new zscores of terminal support cells
suptermzscore0 = zeros(length(suptermpaths0), ftno);
for pcf = 1:length(suptermpaths0)
%
%
    afeat = suptermfeat0(pcf,: );
    zscor = (afeat - mean(afeat))/std(afeat);
    suptermzscore0(pcf,: ) = zscor;
%
%
end

suptermzscore1 = zeros(length(suptermpaths1), ftno);
for pcf = 1:length(suptermpaths1)
%
%
    afeat = suptermfeat1(pcf,: );
    zscor = (afeat - mean(afeat))/std(afeat);
    suptermzscore1(pcf,: ) = zscor;
%
%
end

suptermzscore2 = zeros(length(suptermpaths2), ftno);
for pcf = 1:length(suptermpaths2)
%
%
    afeat = suptermfeat2(pcf,: );
    zscor = (afeat - mean(afeat))/std(afeat);
    suptermzscore2(pcf,: ) = zscor;
%
%
end
% 
% 
% %zscores of non support cells
% nosupzscore = zeros(nosupcnt, ftno);
% for pcf = 1:nosupcnt
% %
% %
%     afeat = nosupfeatures(pcf,: );
%     zscor = (afeat - mean(afeat))/std(afeat);
%     nosupzscore(pcf,: ) = zscor;
% %
% %
% end

%suptermzscore0=ssupzscore0(supportsterm0,:);
%suptermzscore1=ssupzscore1(supportsterm1,:);
%suptermzscore2=ssupzscore2(supportsterm2,:);

%nosupcorrmat0 = (1/(ftno-1))*(nosupzscore*suptermzscore0');
%nosupcorrmat1 = (1/(ftno-1))*(nosupzscore*suptermzscore1');
%nosupcorrmat2 = (1/(ftno-1))*(nosupzscore*suptermzscore2');

nosupcorrmat0 = (1/(ftno-1))*(nosupfeatures*suptermzscore0');
nosupcorrmat1 = (1/(ftno-1))*(nosupfeatures*suptermzscore1');
nosupcorrmat2 = (1/(ftno-1))*(nosupfeatures*suptermzscore2');

%correlation density between non support and terminal group calculation
nosupcorrheav0 = heaviside(nosupcorrmat0 - 0.5);
nosupcorrheav1 = heaviside(nosupcorrmat1 - 0.5);
nosupcorrheav2 = heaviside(nosupcorrmat2 - 0.5);

nosupcorrden0 = sum(nosupcorrheav0,2);
nosupcorrden1 = sum(nosupcorrheav1,2);
nosupcorrden2 = sum(nosupcorrheav2,2);

%enhanced correlation coefficient between non support of terminal support
%groups
nosupcorcor0 = nosupcorrden0.*sum(nosupcorrmat0,2);
nosupcorcor1 = nosupcorrden1.*sum(nosupcorrmat1,2);
nosupcorcor2 = nosupcorrden2.*sum(nosupcorrmat2,2);

%nosupcorcor0 = sum(nosupcorrmat0,2)/length(suptermpaths0);
%nosupcorcor1 = sum(nosupcorrmat1,2)/length(suptermpaths1);
%nosupcorcor2 = sum(nosupcorrmat2,2)/length(suptermpaths2);

cfold = 'C:\Users\alexp\Desktop\Images\SPcombined\';

for den=1:nosupcnt
%reclassifying nonsupport cells

    dens0=nosupcorcor0(den);
    dens1=nosupcorcor1(den);
    dens2=nosupcorcor2(den);

    celpathname=nosuppaths{den};
    [path, name, ext]=fileparts(celpathname);

    [maxdens,maxind]=max([dens0,dens1,dens2]);

    if maxind==1
        
        copyfile(strcat(cfold, name, ext), './relabel/new_class/Early')

    end

    if maxind==2
        
        copyfile(strcat(cfold, name, ext), './relabel/new_class/Late')

    end

    if maxind==3
        
        copyfile(strcat(cfold, name, ext), './relabel/new_class/Viable')

    end
end


for fl=1:length(suptermpaths0)
%adding terminal support cells to their original classes

    celpathname=suptermpaths0{fl};
    [path, name, ext]=fileparts(celpathname);
        
    copyfile(strcat(cfold, name, ext), './relabel/new_class/Early')

end

for fl=1:length(suptermpaths1)
%adding terminal support cells to their original classes

    celpathname=suptermpaths1{fl};
    [path, name, ext]=fileparts(celpathname);
        
    copyfile(strcat(cfold, name, ext), './relabel/new_class/Late')

end


for fl=1:length(suptermpaths2)
%adding terminal support cells to their original classes

    celpathname=suptermpaths2{fl};
    [path, name, ext]=fileparts(celpathname);
        
    copyfile(strcat(cfold, name, ext), './relabel/new_class/Viable')

end

