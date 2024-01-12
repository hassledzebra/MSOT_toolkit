
%data = load('/Volumes/Samsung_T5/ZH 2-DG 2-LMP/Scan_6/output/msot_processed.mat');
data = load('E:\fat\Scan_1\output\unmix_processed.mat');
outputfolder = 'E:\fat\Scan_1\output';
img = data.unmix_img;
div = round(size(img,2)/4);
img = img(:,2*div+1:3*div,:,:);
w = data.w;

w_interest = 800;
ind = find(abs(w-w_interest) < 1e-5);

dce = squeeze(img(:,:,ind,:));
data = dce;
save(strcat(outputfolder,'m1_wat.mat'),'data')

data = load('E:/ZH 2-DG 2-LMP/Scan_6/output/unmix_processed.mat');
img = data.unmix_img;
img = img(:,1:div,:,:);

dce = squeeze(img(:,:,ind,:));
data = dce;
save(strcat(outputfolder,'m1_Hb'),'data')

data = load('E:/ZH 2-DG 2-LMP/Scan_6/output/unmix_processed.mat');
img = data.unmix_img;
img = img(:,div+1:2*div,:,:);


dce = squeeze(img(:,:,ind,:));
data = dce;
save(strcat(outputfolder,'m1_HbO2'),'data')

data = load('E:/ZH 2-DG 2-LMP/Scan_6/output/unmix_processed.mat');
img = data.unmix_img;
img = img(:,3*div+1:end,:,:);


dce = squeeze(img(:,:,ind,:));
data = dce;
save(strcat(outputfolder,'m1_BAT'),'data')
%%
figure
imagesc(img(:,:,4, 179),[0,0.01])
colormap hot
colorbar