%clusters parameters (see sample csv) by local density (adapted from 'Clustering by fast search and find of density peaks', A. Rodriguez and A. Laio) 

clear
clc

delete('results_lcldens\*')

percent = 1.8;

%glcm_path = './results_lcldens/LabelB_bin.xlsx';
%glcm_path = './dataf/LabelB_bin.xlsx';
%glcm_path = './dataf/GLCMparameters_mean_90_24hr.xlsx';
glcm_path = './dataf/GLCMparameters_mean_90_48hr.xlsx';

glcm_table = readtable(glcm_path, 'VariableNamingRule', 'preserve');

%taking only the desired GLCM values (mean)
glcm_S = [glcm_table(:,1:15), glcm_table(:,17)];
glcm_P = [glcm_table(:,21:34), glcm_table(:,36)];

%converting table to matlab variable
new_glcm = [glcm_S, glcm_P];

glcm_min_arr = min(new_glcm{:,2:end});
glcm_min = [array2table(0), array2table(glcm_min_arr)];

%subtracting the minimum out of every column
for mincols = 2:width(glcm_min)
    
    new_glcm{:, mincols} = new_glcm{:, mincols} - glcm_min{:, mincols};
    
end

%finding the maximum value for each GLCM parameter
glcm_max_arr = max(new_glcm{:,2:end});
glcm_max = [array2table(0) , array2table(glcm_max_arr)];

%normalizing the GLCM parameter by maximum and minimums of each column
zero_glcm = new_glcm;
zero_glcm{:,2:end} = 0;
norm_glcm = zero_glcm;
for cols = 2:width(glcm_max)
    
    norm_col = new_glcm{:, cols} / glcm_max{:, cols};
    norm_glcm{:, cols} = norm_glcm{:, cols} + norm_col;
    
end

glcm_arr = norm_glcm{:,2:end};

prm_num = 30;

%counting the images and making an array to store the distances
img_cnt = length(glcm_arr);
diss = zeros(img_cnt^2, prm_num + 2);

for img = 1:img_cnt
    %looking at one image pair set of parameters
    params = glcm_arr(img, :);
    for img_otr = 1:img_cnt
        %looking at another image pair set of parameters
        params_otr = glcm_arr(img_otr, :);
    
        %taking the difference between single image parameter set in outer
        %loop and all the other image parameter sets
        param_dist = abs(params_otr - params);
        %param_avg = mean(param_dist);
        slice = img*img_cnt - img_cnt - img + img_otr;
        
        if img_otr > img
            diss(slice, 1) = diss(slice, 1) + img;
            diss(slice, 2) = diss(slice, 2) + img_otr;
            diss(slice, 3: prm_num + 2 ) = diss(slice, 3: prm_num + 2 ) + param_dist;
        end
    end
end

%performing vector addition on the parameters for each cell
DISS = diss(any(diss,2),:);
dissq = DISS(:,3:end).^2;
ds_fin = sqrt(sum(dissq,2));
xx = [DISS(:,1:2) ds_fin];

N=nchoosek(img_cnt,2);
% %cutoff distance
position = round(N*percent/100);
sda = sort(xx(:,3));
dc = sda(position);

ND=max(xx(:,2));
NL=max(xx(:,1));
if (NL>ND)
    ND=NL;
end
N=size(xx,1);
%prm_num = 1:(size(xx,2)-2);
for i=1:ND
    for j=1:ND
        dist(i,j)=0;
    end
end

for i=1:N
    ii=xx(i,1);
    jj=xx(i,2);
    dist(ii,jj) = xx(i,3);
    dist(jj,ii) = xx(i,3);
end

for i=1:ND
    rho(i)=0.;
end

for i=1:ND
    for j=1:ND
        disd(i,j)=0;
    end
end

% %old method for rho
% for i=1:ND-1
%    for j=i+1:ND
%        rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
%        rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
%    end
% end
% %old method for rho

%new method for rho
for i=1:ND
    for j=1:ND
        if dist(i,j) < dc
             disd(i,j) = 1;
        end
    end
end

rho = sum(disd)/2;
%new method for rho

maxd=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

for ii=2:ND
    delta(ordrho(ii))=maxd;
    for jj=1:ii-1
        if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
            delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
            nneigh(ordrho(ii))=ordrho(jj);
        end
    end
end

% Wrho = 0.08;
% Wdelta = 0.92;
% figmer = Wrho*rho + Wdelta*delta;

% figure(1)
% plot(rho(:),figmer(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')

delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

fid = fopen('DECISION_GRAPH', 'w');
for i=1:ND
    fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end


disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
%figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
figure(2)
% for i=1:ND
%     ind(i)=i;
%     gamma(i)=rho(i)*delta(i);
% end
%subplot(2,1,1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
%title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')


%subplot(2,1,1)
rect = drawrectangle;
rhomin=rect.Position(1);
deltamin=rect.Position(2);
NCLUST=0;
for i=1:ND
    cl(i)=-1;
end
for i=1:ND
    if ( (rho(i)>rhomin) && (delta(i)>deltamin))
        NCLUST=NCLUST+1;
        cl(i)=NCLUST;
        icl(NCLUST)=i;
    end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')

%assignation
for i=1:ND
    if (cl(ordrho(i))==-1)
        cl(ordrho(i))=cl(nneigh(ordrho(i)));
    end
end

for i=1:NCLUST
    nc=0;
    for j=1:ND
        if (cl(j)==i)
            nc=nc+1;
        end
    end
    fprintf('CLUSTER: %i CENTER: %i POINTS: %i \n', i,icl(i),nc);
end

%colormap hsv
%cmap=colormap;
cmap=[1 0 0; 0 0 1; 1 1 0; 0 1 0; 1 0 1; 0.929 0.694 0.125];
%for i=1:NCLUST %- for non-custom colormap
for ic=1:NCLUST
%    ic=int8((i*64.)/(NCLUST*1.)); %- for non custom colormap
    hold on
    %plot(rho(icl(ic)),delta(icl(ic)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
    %for non custom colormap
    plot(rho(icl(ic)),delta(icl(ic)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

% subplot(2,1,2)
% disp('Performing 2D nonclassical multidimensional scaling')
% Y1 = mdscale(dist, 2, 'criterion','sammon');
% title ('multidimensional scaling','FontSize',15.0)
% xlabel ('X')
% ylabel ('Y')

% for i=1:ND
%     A(i,1)=0.;
%     A(i,2)=0.;
% end
% for i=1:NCLUST
%     nn=0;
%     ic=int8((i*64.)/(NCLUST*1.));
%     for j=1:ND
%         if (cl(j)==i)
%             nn=nn+1;
%             A(nn,1)=Y1(j,1);
%             A(nn,2)=Y1(j,2);
%         end
%     end
%     hold on
%     plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
% end

glcm1 = glcm_arr(:,1)';
glcm2 = glcm_arr(:,2)';
glcm13 = glcm_arr(:,13)';
glcm15 = glcm_arr(:,15)';
 
idx3 = cl;

idx_1 = find(idx3==1);
idx_2 = find(idx3==2);
idx_3 = find(idx3==3);
idx_4 = find(idx3==4);
idx_5 = find(idx3==5);
idx_6 = find(idx3==6);
%idx_7 = find(idx3==7);

bin_list1 = glcm_table(idx_1,:);
bin_list2 = glcm_table(idx_2,:);
bin_list3 = glcm_table(idx_3,:);
bin_list4 = glcm_table(idx_4,:);
bin_list5 = glcm_table(idx_5,:);
bin_list6 = glcm_table(idx_6,:);
%bin_list7 = glcm_table(idx_7,:);

writetable(bin_list1, './results_lcldens/LabelA_bin.xlsx')
writetable(bin_list2, './results_lcldens/LabelB_bin.xlsx')
writetable(bin_list3, './results_lcldens/LabelC_bin.xlsx')
%writetable(bin_list4, './results_lcldens/LabelD_bin.xlsx')
%writetable(bin_list5,'./results_lcldens/LabelE_bin.xlsx')
%writetable(bin_list6,'./results_lcldens/LabelF_bin.xlsx')
%writetable(bin_list7,'./results_lcldens/LabelG_bin.xlsx')

%writetable(bin_list1, './results_lcldens/LabelAA_bin.xlsx')
%writetable(bin_list2, './results_lcldens/LabelAB_bin.xlsx')
%writetable(bin_list1, './results_lcldens/LabelBA_bin.xlsx')
%writetable(bin_list2, './results_lcldens/LabelBB_bin.xlsx')

for i=1:ND
    A1(i)=0.;
    A2(i)=0.;
    A13(i)=0.;
    A15(i)=0.;
end
for i=1:NCLUST
    nn=0;
    %ic=int8((i*64.)/(NCLUST*1.));
    for j=1:ND
       if (cl(j)==i)
            nn=nn+1;
            A1(nn)=glcm1(j);
            A2(nn)=glcm2(j);
            A13(nn)=glcm13(j);
            A15(nn)=glcm15(j);            
       end
    end
    figure(3)
    view(3)
    hold on
    %plot3(A13(1:nn),A1(1:nn),A15(1:nn),'o','MarkerSize', 2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:))
    %for non custom colormap
    plot3(A13(1:nn),A1(1:nn),A15(1:nn),'o','MarkerSize', 2,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:)) 
    alpha(0.1)
    grid on;
    %title('GLCM parameters')
    xlabel('s-CLS')
    ylabel('s-ASM')
    zlabel('s-MAP')
    figure(4)
    view(3)
    hold on
    %plot3(A2(1:nn),A15(1:nn),A13(1:nn),'o','MarkerSize', 2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:))
    %for non custom colormap
    plot3(A2(1:nn),A15(1:nn),A13(1:nn),'o','MarkerSize', 2,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:))
    alpha(0.1)
    grid on;
    %title('GLCM parameters')
    xlabel('s-CON')
    ylabel('s-MAP')
    zlabel('s-CLS')
end

rho1 = rho(icl(1))
rho2 = rho(icl(2))
%rho3 = rho(icl(3))
%rho4 = rho(icl(4))
%rho5 = rho(icl(5))
%rho6 = rho(icl(6))
%rho7 = rho(icl(end))

delta1 = delta(icl(1))
delta2 = delta(icl(2))
%delta3 = delta(icl(3))
%delta4 = delta(icl(4))
%delta5 = delta(icl(5))
%delta6 = delta(icl(6))
%delta7 = delta(icl(end))

A=length(find(idx3==1))
B=length(find(idx3==2))
C=length(find(idx3==3))
%D=length(find(idx3==4))
%E=length(find(idx3==5))
%F=length(find(idx3==6))
%G=length(find(idx3==7))

A/(A+B+C)
B/(A+B+C)
C/(A+B+C)
%D/(A+B+C+D+E+F)
%E/(A+B+C+D+E+F)
%F/(A+B+C+D+E+F)
%G/(A+B+C+D+E+F+G)




