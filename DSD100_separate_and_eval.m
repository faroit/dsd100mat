% This script is intended to perform the full evaluation of the separation
% quality of your method on the Demixing Secrets Dataset 100 (DSD100).
%
% It is primarily prepared for the task of "professionally-produced
% music recordings (MUS)" of community-based Signal Separation
% Evaluation Campaign (SiSEC) (https://sisec.inria.fr/).
%
% This function should be used along with the Demixing Secret Dataset 100
% (DSD100) for the purpose of music source separation . The file
% "DSD100_separate_and_eval_parallel.m" should be placed in a root folder,
% along with the folder "DSD100" and the separation function to be evaluated.
%
% If you have the parallel toolbox, this code will automatically parallelize
% evaluation.
%
% The separation function should be named "myfunction.m," placed in the
% same root folder, and have the following syntax:
%   [bass, drums, other, vocals, accompaniment] = myfunction(mixture, fs)
% where "mixture" is a matrix of size [#samples, #channels] corresponding
% to the mixture stereo waveform, "fs" is the corresponding sampling
% frequency in Hz, and "bass," "drums," "other," "vocals," and "accompaniment"
% are matrices of same size as the mixture corresponding to the estimates,
% i.e., the bass, the drums, the other instruments, the vocals, and the full
% accompaniment (i.e., bass+drums+other), respectively. If one or more sources
% are not meant to be estimated, the function should return an empty matrix
% (i.e., []). Any other parameter of the algorithm should be defined internally.
%
% The evaluation function should then be called simply as follows:
%  DSD100_separate_and_eval.m
%
% The function loops over all the 100 songs of the DSD100 data set, and,
% for each song, for both the "Dev" and "Test" subsets, performs source
% separation on the mixture "mixture.wav" from the folder "Mixtures" using
% the separation function "myfunction.m" and saves the estimates as
% "bass.wav," "drums.wav," "other.wav,", "vocals.wav" and "accompaniment.wav"
% (if estimated) to the folder "EstimatesMETHOD" where METHOD is the name
% of your method.
%
% The function then measures the separation performance using the BSS Eval
% toolbox 3.0 (included in this function) and the sources from the folder
% "Sources," and saves the results (i.e., SDR, ISR, SIR, and SAR) in the file
% "results.mat," including the song name and the processing time, along with
% the estimates to the folder "Estimates". The function also saves the results
% for all the songs in a single file "resultMETHOD.mat" to the root folder, along
% with this function.
%  A first evaluation is performed for the 4 sources vocals/bass/drums
% and other, and a second is performed for accompaniment.
%
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
% Original Author: Zafar Rafii,
% Updated by A. Liutkus

function DSD100_separate_and_eval_parallel

% here, provide the root folder for the DSD100 dataset.
dataset_folder = 'DSD100';

% here provide the name of the directory where the results for your methods
% will be saved
base_estimates_directory = 'link/to/your/results/rootfolder';

% Your method name
% will save the estimates in base_estimates_directory/method_name/...
method_name='YOURMETHOD';

subsets_names = {'Dev','Test'};
sources_names = {'bass','drums','other','vocals','accompaniment'};
estimates_folder = fullfile(base_estimates_directory,sprintf('Estimates%s',method_name));

mkdir(estimates_folder)
result = struct;
warning('off','all')

%comment this line if you don't want to use parallel toolbox. This may
%happen for memory usage reasons. Also, modify this line if you don't want
%to use the maximum number of cores, which is the default behaviour.
parpool('local');

%loop over the subsets: dev and test
for i = 1:numel(subsets_names)
    mixtures_folder = fullfile(dataset_folder,'Mixtures',subsets_names{i});
    mixtures_dir = dir(mixtures_folder);
    n = numel(mixtures_dir)-2;
    estimates_folder = fullfile(pwd,sprintf('Estimates%s',method_name),subsets_names{i});
    mkdir(estimates_folder)
    sources_folder = fullfile(dataset_folder,'Sources',subsets_names{i});
    results_set=cell(50,1);

    parfor j = 1:n		%change to a standard for loop if you don't want parallelization
        %load mixture
        mixture_name = mixtures_dir(j+2).name;
        disp([subsets_names{i},': ',num2str(j),'/',num2str(n),' ',mixture_name])
        mixture_file = fullfile(mixtures_folder,mixture_name,'mixture.wav');
        [mixture,fs] = wavread(mixture_file);


        [l,c] = size(mixture);

        output = cell(4,1);

        %call the separation function myfunction.m, that should be designed
        %and written by the user. you can optionally change this name to
        %any other one here: this is the only place your separation
        %function is called.
        tic
        [output{1},output{2},output{3},output{4},output{5}] = myfunction(mixture,fs);
        time = toc;
        mixture = [];

        %store the estimates
        estimate_folder = fullfile(estimates_folder,mixture_name);
        mkdir(estimate_folder)
        estimates = zeros(l,c,0);
        for k = 1:5
            estimate = output{1};
            output(1) = [];
            if ~isempty(estimate)
                estimate_file = fullfile(estimate_folder,[sources_names{k},'.wav']);
                wavwrite(estimate,fs,estimate_file)
                estimate = repmat(estimate,[1,3-size(estimate,2)]);

                estimate = estimate(1:l,1:c);
            else
                estimate = zeros(l,c);
            end
            estimates = cat(3,estimates,estimate);
            estimate = [];
        end

        %now load the sources references
        source_folder = fullfile(sources_folder,mixture_name);
        sources = zeros(l,c,0);
        accompaniment = zeros(l,c);
        for k = 1:3
            source_file = fullfile(source_folder,[sources_names{k},'.wav']);
            source = wavread(source_file);
            source = repmat(source,[1,3-size(source,2)]);
            sources = cat(3,sources,source);
            accompaniment = accompaniment+source;
            source = [];
        end
        source_file = fullfile(source_folder,[sources_names{4},'.wav']);
        source = wavread(source_file);
        source = repmat(source,[1,3-size(source,2)]);
        sources = cat(3,sources,source,accompaniment);
        source = [];
        accompaniment= [];


        %estimate quality for accompaniment
        [SDRac,ISRac,SIRac,SARac] = bss_eval(estimates(:,:,4:5) ,...
            sources(:,:,4:5),...
            30*fs,15*fs);

        %estimate quality for sources
        [SDR,ISR,SIR,SAR] = bss_eval(estimates(:,:,1:4),sources(:,:,1:4),...
            30*fs,15*fs);

        %append the accompaniment quality to the sources quality
        SDR(end+1,:) = SDRac(2,:);
        ISR(end+1,:) = ISRac(2,:);
        SIR(end+1,:) = SIRac(2,:);
        SAR(end+1,:) = SARac(2,:);


        estimates = [];
        sources= [];

        results_set{j} = struct;
        results_set{j}.name = mixture_name;
        results_set{j}.time = time;
        for k = 1:5
            results_set{j}.(sources_names{k}).sdr = SDR(k,:);
            results_set{j}.(sources_names{k}).isr = ISR(k,:);
            results_set{j}.(sources_names{k}).sir = SIR(k,:);
            results_set{j}.(sources_names{k}).sar = SAR(k,:);
        end
    end
    % now  store the results in the subset folder
    for song_index=1:50
        mixture_name = mixtures_dir(song_index+2).name;
        results_file = fullfile(estimates_folder,mixture_name,'results.mat');
        results = results_set{song_index};
        save(results_file,'results')
        result.(lower(subsets_names{i}))(song_index).results =results_set{song_index};
    end
end
%save the globar results in the root folder
result_file = fullfile(dataset_folder,sprintf('result%s.mat',method_name));
save(result_file,'result')
warning('on','all')


%bss eval code
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
