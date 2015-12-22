% This script is intended to perform the full evaluation of the separation
% quality of your method on the Demixing Secrets Dataset 100 (DSD100).
%
% It was initially prepared for the task of "professionally-produced
% music recordings (MUS)" of the fifth community-based Signal Separation
% Evaluation Campaign (SiSEC 2015) (https://sisec.inria.fr/).
%
% This function should be used along with the Demixing Secret Dataset 100
% (DSD100) for the purposed of music source separation.
%
% The separation function should be named "myfunction.m," placed in the
% root folder, and have the following syntax:
%   [bass, drums, other, vocals, accompaniment] = myfunction(mixture, fs)
% where "mixture" is a matrix of size [#samples, #channels] corresponding
% to the mixture, "fs" is the corresponding sampling frequency in Hz, and
% "bass," "drums," "other," "vocals," and "accompaniment" are matrices of
% same size as the mixture corresponding to the estimates, i.e., the bass,
% the drums, the other instruments, the vocals, and the full accompaniment
% (i.e., bass+drums+other), respectively. If one or more sources are not
% meant to be estimated, the function should return an empty matrix (i.e.,
% []). Any other parameter of the algorithm should be defined internally.
%
% The evaluation function should then be called simply as follows:
%  DSD100_separate_and_eval.m
%
% The function loops over all the 100 songs of the MSD100 data set, and,
% for each song, for both the "Dev" and "Test" subsets, performs source
% separation on the mixture "mixture.wav" from the folder "Mixtures" using
% the separation function "myfunction.m" and saves the estimates as
% "bass.wav," "drums.wav," "other.wav," and "vocals.wav" (if estimated) to
% the folder "Estimates." The function then measures the separation
% performance using the BSS Eval toolbox 3.0 (included in this function)
% and the sources from the folder "Sources," and saves the results (i.e.,
% SDR, ISR, SIR, and SAR) in the file "results.mat," including the song
% name and the processing time, along with the estimates to the folder
% "Estimates". The function also saves the results for all the songs in a
% single file "result.mat" to the root folder, along with this function.
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

function DSD100_separate_and_eval

method_name = 'MYMETHOD'; % change your method name here

warning('off','all')
dataset_folder = fullfile('./DSD100');

subsets_names = {'Dev','Test'};
sources_names = {'bass','drums','other','vocals','accompaniment'};
estimates_folder = fullfile(pwd,sprintf('Estimates%s',method_name));

mkdir(estimates_folder)
result = struct;
for i = 1:numel(subsets_names)
    mixtures_folder = fullfile(dataset_folder,'Mixtures',subsets_names{i});
    mixtures_dir = dir(mixtures_folder);
    % Extract only those that are directories.
    mixtures_dir = mixtures_dir([mixtures_dir.isdir]);
    mixtures_names = {mixtures_dir.name};
    mixtures_names(ismember(mixtures_names,{'.','..'})) = [];

    n = numel(mixtures_names);
    estimates_folder = fullfile(pwd,sprintf('Estimates%s',method_name),subsets_names{i});
    mkdir(estimates_folder)
    sources_folder = fullfile(dataset_folder,'Sources',subsets_names{i});
    for j = 1:n
        disp([subsets_names{i},': ',num2str(j),'/',num2str(n),' ',mixtures_names{j}])
        mixture_file = fullfile(mixtures_folder,mixtures_names{j},'mixture.wav');
        [mixture,fs] = wavread(mixture_file); %#ok<*DWVRD>

        [l,c] = size(mixture);
        tic

        [output{1},output{2},output{3},output{4},output{5}] = myfunction(mixture,fs);
        time = toc;
        clear mixture

        estimate_folder = fullfile(estimates_folder,mixtures_names{j});
        mkdir(estimate_folder)
        estimates = zeros(l,c,0);
        for k = 1:5
            estimate = output{1};
            output(1) = [];
            if ~isempty(estimate)
                estimate_file = fullfile(estimate_folder,[sources_names{k},'.wav']);
                wavwrite(estimate,fs,estimate_file) %#ok<*DWVWR>
                estimate = repmat(estimate,[1,3-size(estimate,2)]);

                estimate = estimate(1:l,1:c);
            else
                estimate = zeros(l,c);
            end
            estimates = cat(3,estimates,estimate);
            clear estimate
        end

        source_folder = fullfile(sources_folder,mixtures_names{j});
        sources = zeros(l,c,0);
        accompaniment = zeros(l,c);
        for k = 1:3
            source_file = fullfile(source_folder,[sources_names{k},'.wav']);
            source = wavread(source_file);
            source = repmat(source,[1,3-size(source,2)]);
            sources = cat(3,sources,source);
            accompaniment = accompaniment+source;
            clear source
        end
        source_file = fullfile(source_folder,[sources_names{4},'.wav']);
        source = wavread(source_file);
        source = repmat(source,[1,3-size(source,2)]);
        sources = cat(3,sources,source,accompaniment);
        clear source accompaniment

        [SDR,ISR,SIR,SAR] = bss_eval(estimates,sources,30*fs,15*fs);
        clear estimates sources
        results = struct;
        results.name = mixtures_names{j};
        results.time = time;
        for k = 1:5
            results.(sources_names{k}).sdr = SDR(k,:);
            results.(sources_names{k}).isr = ISR(k,:);
            results.(sources_names{k}).sir = SIR(k,:);
            results.(sources_names{k}).sar = SAR(k,:);
        end
        results_file = fullfile(estimate_folder,'results.mat');
        save(results_file,'results')
        result.(lower(subsets_names{i}))(j).results  = results;
    end
end
result_file = fullfile(estimates_folder,sprintf('result%s.mat',method_name));
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
