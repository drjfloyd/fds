% McDermott
% 12-19-11
% flame_height2.m

close all
clear all

disp('flame_height2 ...')

outdir = '../../../out/Heskestad_Flame_Height/';
expdir = '../../../exp/Heskestad_Flame_Height/';

% list of line files
filename = {'Qs=p1_RI=05_line.csv',    'Qs=p1_RI=10_line.csv',    'Qs=p1_RI=20_line.csv';   ...
            'Qs=p2_RI=05_line.csv',    'Qs=p2_RI=10_line.csv',    'Qs=p2_RI=20_line.csv';   ...
            'Qs=p5_RI=05_line.csv',    'Qs=p5_RI=10_line.csv',    'Qs=p5_RI=20_line.csv';   ...
            'Qs=1_RI=05_line.csv',     'Qs=1_RI=10_line.csv',     'Qs=1_RI=20_line.csv';    ...
            'Qs=2_RI=05_line.csv',     'Qs=2_RI=10_line.csv',     'Qs=2_RI=20_line.csv';    ...
            'Qs=5_RI=05_line.csv',     'Qs=5_RI=10_line.csv',     'Qs=5_RI=20_line.csv';    ...
            'Qs=10_RI=05_line.csv',    'Qs=10_RI=10_line.csv',    'Qs=10_RI=20_line.csv';   ...
            'Qs=20_RI=05_line.csv',    'Qs=20_RI=10_line.csv',    'Qs=20_RI=20_line.csv';   ...
            'Qs=50_RI=05_line.csv',    'Qs=50_RI=10_line.csv',    'Qs=50_RI=20_line.csv';   ...
            'Qs=100_RI=05_line.csv',   'Qs=100_RI=10_line.csv',   'Qs=100_RI=20_line.csv';  ...
            'Qs=200_RI=05_line.csv',   'Qs=200_RI=10_line.csv',   'Qs=200_RI=20_line.csv';  ...
            'Qs=500_RI=05_line.csv',   'Qs=500_RI=10_line.csv',   'Qs=500_RI=20_line.csv';  ...
            'Qs=1000_RI=05_line.csv',  'Qs=1000_RI=10_line.csv',  'Qs=1000_RI=20_line.csv'; ...
            'Qs=2000_RI=05_line.csv',  'Qs=2000_RI=10_line.csv',  'Qs=2000_RI=20_line.csv'; ...
            'Qs=5000_RI=05_line.csv',  'Qs=5000_RI=10_line.csv',  'Qs=5000_RI=20_line.csv'; ...
            'Qs=10000_RI=05_line.csv', 'Qs=10000_RI=10_line.csv', 'Qs=10000_RI=20_line.csv'};

rho_inf = 1.2;
cp = 1;
T_inf = 293;
g = 9.81;
D = 1.13;
Qdot=[151 303 756 1513 3025 7564 15127 30255 75636 151273 302545 756363 1512725 3025450 7563625 15127250];

for i=1:16 % hrr loop
    for j=1:3 % resolution loop

        M = csvread([outdir,filename{i,j}],2,0);
        z = M(:,1); dz = z(2)-z(1);
        hrrpul = M(:,2);
        Qdot_line = sum(hrrpul)*dz;
        Qstar(i) = Qdot(i)/(rho_inf*cp*T_inf*sqrt(g)*D^(5/2));

        % determine flame height
        for n=1:length(z)
            hrr(n) = sum(hrrpul(1:n))*dz*Qdot(i)/Qdot_line; % cummulative heat release
        end

        k = find(hrr>(0.99)*Qdot(i),1);
        if (k>1)
            L_99(i,j) = z(k-1)+dz*((0.99)*Qdot(i)-hrr(k-1))/(hrr(k)-hrr(k-1));
        else
            L_99(i,j) = dz*(0.99)*Qdot(i)/hrr(k);
        end

        k = find(hrr>(0.95)*Qdot(i),1);
        if (k>1)
            L_95(i,j) = z(k-1)+dz*((0.95)*Qdot(i)-hrr(k-1))/(hrr(k)-hrr(k-1));
        else
            L_95(i,j) = dz*(0.95)*Qdot(i)/hrr(k);
        end

    end % resolution loop

end % hrr loop

fclose('all');

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M=importdata([expdir,'flame_lengths.csv'],',',1);
Steward             = M.data(:,5);
Becker_and_Liang    = M.data(:,6);
Cox_and_Chitty      = M.data(:,7);
Heskestad           = M.data(:,8);
Hasemi_and_Tokunaga = M.data(:,9);
Cetegen             = M.data(:,10);
Delichatsios        = M.data(:,11);

H(1)=loglog(Qstar,Steward,'-','Color',[.5 0 .5]); hold on
H(2)=loglog(Qstar,Becker_and_Liang,'g-');
H(3)=loglog(Qstar(2:16),Cox_and_Chitty(2:16),'c-');
H(4)=loglog(Qstar,Heskestad,'k-');
H(5)=loglog(Qstar(3:16),Hasemi_and_Tokunaga(3:16),'m-');
H(6)=loglog(Qstar,Cetegen,'-','Color',[.5 .5 0]);
H(7)=loglog(Qstar,Delichatsios,'-','Color',[0 .5 .5]);
H(8)=loglog(Qstar,max(L_99(:,1:3),[],2),'r--','Linewidth',2);
H(9)=loglog(Qstar,min(L_95(:,1:3),[],2),'b--','Linewidth',2);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

Dep_Title = '{\itQ}*';
Ind_Title = '{\itL}_f/{\itD}';
Plot_Title = 'Flame Height Variation';
Min_Ind = 0.05;
Max_Ind = 2e4;
Min_Dep = 0.1;
Max_Dep = 3e2;
Title_Position = [0.1 0.9];
X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
Y_Title_Position = 10^(log10(Min_Dep)+Title_Position(2)*(log10(Max_Dep)-log10(Min_Dep)));

xlabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([Min_Ind Max_Ind Min_Dep Max_Dep])
text(X_Title_Position,Y_Title_Position,...
    Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

legend_handle=legend(H,...
                    'Steward',...
                    'Becker & Liang',...
                    'Cox & Chitty',...
                    'Heskestad',...
                    'Hasemi & Tokunaga',...
                    'Cetegen',...
                    'Delichatsios',...
                    'Max FDS 99%',...
                    'Min FDS 95%',...
                    'Location','EastOutside');

set(legend_handle,'FontName',Font_Name,'Interpreter',Font_Interpreter)
set(legend_handle,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)
set(legend_handle,'Units',Paper_Units)
pos = get(legend_handle,'position');
set(legend_handle,'position',[Paper_Width pos(2:4)])

set(gca,'YTick',[1e-1 1e0 1e1 1e2 1e3])
set(gca,'XTick',[1e-1 1e0 1e1 1e2 1e3 1e4])

% add VerStr if file is available

git_file = [outdir,'Qs=10000_RI=05_git.txt'];
addverstr(gca,git_file,'loglog')

% print to pdf

Paper_Width=1.35*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
display('Printing plot Flame_Height2...')
print -dpdf '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Heskestad/Flame_Height2'

