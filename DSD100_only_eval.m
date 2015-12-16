% This script is intended to perform the evaluation of the separation
% quality on the Demixing Secrets Dataset 100 (DSD100).
% It was initially prepared for the task of "professionally-produced
% music recordings (MUS)" of the fifth community-based Signal Separation
% Evaluation Campaign (SiSEC 2015) (https://sisec.inria.fr/).
%
% This function should be used along with the DSD100. The file
% "DSD100_only_eval.m" should be placed in a root folder, along with the
% folder "DSD100" and the folder, called YOURFOLDER here, containing your
% results to evaluate.
%
% The folder YOURFOLDER should have exactly the same structure as the
% DSD100/Sources folder, the matching is case sensitive. In each directory,
% there should be the separated sources estimates whose quality is to be
% evaluated. There is the possibility of not including all sources, and
% also to include the "accompaniment" source, defined as the sum of all
% sources except vocals.
%
% The evaluation function should then be called simply as follows:
%   DSD100_only_eval.m
% The function loops over all the 100 songs of the MSD100 data set, and,
% for each song, for both the "Dev" and "Test" subsets, performs evaluation
% using the BSS Eval toolbox 3.0 (included in this function) and saves the
% results (i.e. SDR, ISR, SIR, and SAR) in the file "results.mat," including
% the song name. The function also saves the results for all the songs in a
% single file "resultX.mat" to the root folder, along with this function.

% We would like to thank Emmanuel Vincent for giving us the permission to
% use the BSS Eval toolbox 3.0;
%
% If you use this script, please reference the following paper
%
%@inproceedings{SiSEC2015,
%  TITLE = {{The 2015 Signal Separation Evaluation Campaign}},
%  AUTHOR = {N. Ono and Z. Rafii and D. Kitamura and N. Ito and A. Liutkus},
%  BOOKTITLE = {{International Conference on Latent Variable Analysis and Signal Separation  (LVA/ICA)}},
%  ADDRESS = {Liberec, France},
%  SERIES = {Latent Variable Analysis and Signal Separation},
%  VOLUME = {9237},
%  PAGES = {387-395},
%  YEAR = {2015},
%  MONTH = Aug,
%}
%
% A more adequate reference will soon be given at sisec.inria.fr
% please check there and provide the one provided.
%
% Original Author: Zafar Rafii, zafarrafii@gmail.com
% Last updated by A.Liutkus on December 7, 2015

function dsd100_only_eval

estimates_name = 'Estimates'; %here provide the name of the directory to eval

warning('off','all')
dataset_folder = fullfile(pwd,'.');
subsets_names = {'Dev','Test'};
sources_names = {'bass','drums','other','vocals','accompaniment'};
result = struct;
for subset_index = 1:2
    sources_folder = fullfile(dataset_folder,'Sources',subsets_names{subset_index});
    estimates_folder = fullfile(dataset_folder,estimates_name,subsets_names{subset_index});
    estimates_list = dir(estimates_folder);
    % Extract only those that are directories.
    estimates_list = estimates_list([estimates_list.isdir])
    estimates_names = {estimates_list.name}';
    estimates_names(ismember(estimates_names,{'.','..'})) = [];
    for song_index = 1:numel(estimates_names)
        disp([subsets_names{subset_index},' ',num2str(song_index),'/',num2str(50),' ',estimates_names(song_index)])

        sources_data = [];
        for source_index = 1:4
            source_file = fullfile(sources_folder,estimates_names(song_index),[sources_names{source_index},'.wav']);
            [source_data,source_sampling] = audioread(source_file);
            [sources_samples,sources_channels] = size(source_data);
            source_data = repmat(source_data,[1,3-sources_channels]);
            sources_data = cat(3,sources_data,source_data);
            clear source_data
        end
        sources_data = cat(3,sources_data,sum(sources_data(:,:,1:3),3));

        estimates_data = zeros(sources_samples,2,5);
        for estimate_index = 1:5
            estimate_file = fullfile(estimates_folder,estimates_names(song_index),[sources_names{estimate_index},'.wav']);
            if exist(estimate_file,'file') == 2
                estimate_data = audioread(estimate_file);
                [estimate_samples,estimate_channels] = size(estimate_data);
                estimate_data = repmat(estimate_data,[1,3-estimate_channels]);
                estimates_samples = min(size(estimates_data,1),estimate_samples);
                estimates_data = estimates_data(1:estimates_samples,:,:);
                estimates_data(:,:,estimate_index) = estimate_data(1:estimates_samples,:);
                clear estimate_data
            end
        end

        [SDR,ISR,SIR,SAR] = bss_eval(estimates_data,sources_data,30*source_sampling,15*source_sampling);
        clear estimates_data sources_data
        results = struct;
        results.name = estimates_names(song_index);
        for source_index = 1:5
            results.(sources_names{source_index}).sdr = SDR(source_index,:);
            results.(sources_names{source_index}).isr = ISR(source_index,:);
            results.(sources_names{source_index}).sir = SIR(source_index,:);
            results.(sources_names{source_index}).sar = SAR(source_index,:);
        end
        results_file = fullfile(estimates_folder,estimates_names(song_index),'results.mat');
        save(results_file,'results')
        result.(lower(subsets_names{subset_index}))(song_index).results  = results;
    end
end
result_file = fullfile(dataset_folder,[estimates_name,'.mat']);
save(result_file,'result')
warning('on','all')

function [SDR,ISR,SIR,SAR] = bss_eval(ie,i,win,ove)

[nsampl,~,nsrc] = size(ie);
nwin = floor((nsampl-win+1+ove)/ove);
SDR = zeros(nsrc,nwin);
ISR = zeros(nsrc,nwin);
SIR = zeros(nsrc,nwin);
SAR = zeros(nsrc,nwin);
for k = 1:nwin
    K = (k-1)*ove+1:(k-1)*ove+win;
    [SDR(:,k),ISR(:,k),SIR(:,k),SAR(:,k)] = bss_eval_images(ie(K,:,:),i(K,:,:));
end

function [SDR,ISR,SIR,SAR] = bss_eval_images(ie,i)

nsrc = size(ie,3);
SDR = zeros(nsrc,1);
ISR = zeros(nsrc,1);
SIR = zeros(nsrc,1);
SAR = zeros(nsrc,1);
for j = 1:nsrc
    [s_true,e_spat,e_interf,e_artif] = bss_decomp_mtifilt(ie(:,:,j),i,j,512);
    [SDR(j,1),ISR(j,1),SIR(j,1),SAR(j,1)] = bss_image_crit(s_true,e_spat,e_interf,e_artif);
end

function [s_true,e_spat,e_interf,e_artif] = bss_decomp_mtifilt(se,s,j,flen)

nchan = size(se,2);
s_true = [s(:,:,j);zeros(flen-1,nchan)];
e_spat = project(se,s(:,:,j),flen)-s_true;
e_interf = project(se,s,flen)-s_true-e_spat;
e_artif = [se;zeros(flen-1,nchan)]-s_true-e_spat-e_interf;

function sproj = project(se,s,flen)

[nsampl,nchan,nsrc] = size(s);
s = reshape(s,[nsampl,nchan*nsrc]);
s = [s;zeros(flen-1,nchan*nsrc)];
se = [se;zeros(flen-1,nchan)];
fftlen = 2^nextpow2(nsampl+flen-1);
sf = fft(s',fftlen,2);
sef = fft(se',fftlen,2);

G = zeros(nchan*nsrc*flen);
for k1 = 0:nchan*nsrc-1
    for k2 = 0:k1
        ssf = sf(k1+1,:).*conj(sf(k2+1,:));
        ssf = real(ifft(ssf));
        ss = toeplitz(ssf([1,fftlen:-1:fftlen-flen+2]),ssf(1:flen));
        G(k1*flen+1:k1*flen+flen,k2*flen+1:k2*flen+flen) = ss;
        G(k2*flen+1:k2*flen+flen,k1*flen+1:k1*flen+flen) = ss';
    end
end

D = zeros(nchan*nsrc*flen,nchan);
for k = 0:nchan*nsrc-1
    for i = 1:nchan
        ssef = sf(k+1,:).*conj(sef(i,:));
        ssef = real(ifft(ssef,[],2));
        D(k*flen+1:k*flen+flen,i) = ssef(:,[1,fftlen:-1:fftlen-flen+2])';
    end
end

C = G\D;
C = reshape(C,flen,nchan*nsrc,nchan);
sproj = zeros(nsampl+flen-1,nchan);
for k = 1:nchan*nsrc
    for i = 1:nchan
        sproj(:,i) = sproj(:,i)+fftfilt(C(:,k,i),s(:,k));
    end
end

function [SDR,ISR,SIR,SAR] = bss_image_crit(s_true,e_spat,e_interf,e_artif)

s_true = s_true(:);
e_spat = e_spat(:);
e_interf = e_interf(:);
e_artif = e_artif(:);
SDR = 10*log10(sum(s_true.^2)/sum((e_spat+e_interf+e_artif).^2));
ISR = 10*log10(sum(s_true.^2)/sum(e_spat.^2));
SIR = 10*log10(sum((s_true+e_spat).^2)/sum(e_interf.^2));
SAR = 10*log10(sum((s_true+e_spat+e_interf).^2)/sum(e_artif.^2));
