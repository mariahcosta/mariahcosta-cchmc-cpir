function [scaledfids, sctype] = HFscale(traj,fids,sctype,Method_Params,FIDs,tmp_FIDMat)
% there are three different methods that can be used to scale the k-space
% data; the default is set to fifty
% case 'all' - scales each projection to the mean value of k-space 
% case 2 - 'HF' - scales k-space > than 50% radius to mean value 
% case 3 - 'first' - scales k-space so all pts within a breath are scaled
% based on the k0 of the first projection within each breath

%% MC Update - only getting 'all' method to work 
% MC 12/13/21

k0s = abs(tmp_FIDMat(1,:));
meank0 = mean(k0s);

% find number of projections that are left after deleting begin/end
% npro = sum(sum(abs(fids))>1);
npro = Method_Params.NPro;
traj1 = reshape(traj,[],npro,3);
% traj1 = traj1(:,:,1)
%MC Change 10/28/21
% traj1 = reshape(traj,503,[],3);
% traj1 = permute(traj,[2 3 1]); % permute to make indexing easier


switch sctype
    case 'all'
        disp('scaled all projections');
        scale_factor = meank0./abs(tmp_FIDMat(1,:)); % scale vector for each proj
        scale_factor = repmat(scale_factor,size(tmp_FIDMat,1),1); %scale matrix can be applied to all points
        scaledfids = tmp_FIDMat.*scale_factor;
        
        
    case 'HF'
        %% MC - 12/13
        disp('scaled outer 50% of projections');
        scale_factor = meank0./abs(tmp_FIDMat(1,:)); % scale vector for each proj
        % calculate distance from k0
        k0_dist = zeros(size(traj1)); % zero point for distance calculation

        % kspace is scaled to +/- 0.5 for Scott Robertson's code. We want to scale
        % k-space values more than 50 % of the way from k0 to the edge of k-space
        SRd = 0.5; % Scott Robertson's distance = radius of k-space
        dfromk0 = SRd/2; % half way from k0 to radius of k-space
        
        % distance matrix
        rad = sqrt((traj1(:,:,1)-k0_dist).^2+(traj1(:,:,2)-k0_dist).^2+(traj1(:,:,3)-k0_dist).^2);
        selmat = rad>dfromk0;
        scaledfids = tmp_FIDMat.*selmat(:,:,1).*scale_factor + tmp_FIDMat.*(1-selmat(:,:,1));
        
        
    case 'first'
%         disp('scaled the first projection of each breath');
%         % scale vector for first projection of each breath
%         scale_factor = meank0./abs(FIDs(1,1:Method_Params.ProjPerTrig:npro));
%         % generate scale matrix for scaling fids
%         % copy vector by number of projection per breath
%         scale_factor = repmat(scale_factor,Method_Params.ProjPerTrig,1);
%         % reshape the matrix into a vector
%         scale_factor = reshape(scale_factor,[],1)';
%         % make full scaling matrix
%         scale_factor = repmat(scale_factor,size(FIDs,1),1);
%         % due to repmat and reshape, the scale_factor may come out larger
%         % trim down scale factor to match the size of teh fids
%         scale_factor(:,size(FIDs,2)+1:end) = [];
%         % MC add 10/28/21
%         scale_factor = repmat(scale_factor,1,Method_Params.NSlices*Method_Params.Nbvalue);
%         scaledfids = FIDs.*scale_factor;


% MC attempt - not working 12/14/21
% current version is the same as 'all' method.. will fix later

% Only written for spiral diffusion sequence where Nslices is the inner 
% loop (all within a breath), then bvals, then Npro. This is because I am 
% worried the indexing won't work if all projections are captured in one 
% breath 

if strcmp(Method_Params.SequenceName, '<User:pn_spiral_NonIso>') ==1
    disp('scaled the first projection of each breath');
        % scale vector for first projection of each breath
        
        scale_factor = meank0./abs(tmp_FIDMat(1,1:Method_Params.ProjPerTrig:npro));
        % make full scaling matrix
        scale_factor = repmat(scale_factor,size(tmp_FIDMat,1),1);
        
        % due to repmat and reshape, the scale_factor may come out larger
        % trim down scale factor to match the size of teh fids
        scale_factor(:,size(tmp_FIDMat,2)+1:end) = [];
        
        % MC add 10/28/21
%         scale_factor = repmat(scale_factor,1,Method_Params.NSlices*Method_Params.Nbvalue);
        scaledfids = tmp_FIDMat.*scale_factor;
else
    scaledfids = tmp_FIDMat;
    disp('Error: unable to scale to first proj for this seq type - see HFscale_2D.m')
end


end
disp('done scaling...')





%% Extras: 

%% test a few things
% scale_factor = meank0./abs(fids(1,:)); % scale vector for each proj
% scale_factor = repmat(scale_factor,size(fids,1),1); %scale matrix can be applied to all points
% % calculate distance from k0
% % kspace is scaled to +/- 0.5 for Scott Robertson's code. We want to scale
% % k-space values more than 50 % of the way from k0 to the edge of k-space
% dfk0 = 0.5/2;
% traj1 = permute(traj,[2 3 1]); % permute to make indexing easier
% k0 = zeros(size(fids)); % zero point for distance calculation
% % distance matrix
% d = sqrt((traj1(:,:,1)-k0).^2+(traj1(:,:,2)-k0).^2+(traj1(:,:,3)-k0).^2);
% selmat = d>dfk0;
% scaledfids.*selmat.*scale_factor+fids.*(1-selmat);
% 
% % example that I think shows how the k-space can be scaled based on a
% % selection matrix
% scalefac = ones(size(A))*2;
% A = randi(10, 5,3); % orignial FIDs
% b = A>5; % selmat
% c = A.*b.*scalefac+A.*(1-b); % c is the scaled matrix

%% paper demonstration figure;
% HFtraj = zeros(size(traj1));
%     
% trajxHF(:,:,1) = traj1(:,:,1).*selmat;
% trajyHF(:,:,2) = traj1(:,:,2).*selmat;
% trajzHF(:,:,3) = traj1(:,:,3).*selmat;
% 
% % trajxLF = traj1(:,:,1).*1-selmat;
% % trajyLF = traj1(:,:,2).*1-selmat;
% % trajzLF = traj1(:,:,3).*1-selmat;
% 
% seltraj(seltraj==0) = [];
% HFtraj(:,:,kk) = temptraj.*selmat;
% % HFtraj(HFtraj == 0) = [];
% rpros = [10 78];
% % plot3(squeeze(Raw.Traj(1,:,rpros)),squeeze(Raw.Traj(2,:,rpros)),squeeze(Raw.Traj(3,:,rpros)),'o')
% plot3(squeeze(HFtraj(:,rpros,1)),squeeze(HFtraj(:,rpros,2)),squeeze(HFtraj(:,rpros,3)),'o')
% figure;
% a = [1:5000];
% plot3(trajx(roz(a),coz(a)),trajy(roz(a),coz(a)),trajz(roz(a),coz(a)),'o')
% 
% figure;
% hold on;
% % plot LF
% pro1 = randi(length(trajx),3000,1);
% plot3(trajxLF(1:8:143,pro1),trajyLF(1:8:143,pro1),trajzLF(1:8:143,pro1),'k.');
% set(gca,'XLim',[-0.5 0.5],'YLim',[-0.5 0.5],'ZLim',[-0.5 0.5])
% set(gca,'View',[0,0])
% 
% % plot HF
% pro = randi(length(trajx),5,1);
% plot3(trajxHF(144:end,pro),trajyHF(144:end,pro),trajzHF(144:end,pro),'o');

%% show results for a single projection and for k0
% show = 1;
% if show
%     figure('Name',sprintf('%s scaling',sctype),'NumberTitle','off');
%     ax0 = subplot(2,2,[1,3]);
%     title('single random proj');
%     hold on;
%     proj = randi(npro,1,1);
%     plot(abs(fids(:,proj)));
%     plot(abs(scaledfids(:,proj)));
%     hold off;
%     ax1 = subplot(2,2,2);
%     hold on;
%     title('input k0');
%     plot(abs(fids(1,:)),'Color',[0.6 0.6 0.6]);
%     plot(abs(fids(1,:)),'bo');
%     hold off;
%     ax2 = subplot(2,2,4);
%     title('scaled k0');
%     hold on;
%     plot(abs(scaledfids(1,:)),'Color',[0.6 0.6 0.6]);
%     plot(abs(scaledfids(1,:)),'ro');
%     hold off;
%     linkaxes([ax1 ax2],'x');
    
end
