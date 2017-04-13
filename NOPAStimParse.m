%Script for parsing and assigning stimuli for NOPA study
%
% script needs to take a list of "correctly paired" stimuli
% and a list of files from a directory of stimuli
% and generate a series of txt files that contain a subset
% of correctly paired, and incorrectly paired stimuli from the list of
% possible options. Each txt file represents an individual block, so the
% stimulus pairings need to be tracked from block to block such that over
% the total number of blocks each correct stimulus pair is tested 50 times,
% (50 * 15 = 750 correct stimulus presentations)

%first things first, what do we need...?
%...turns out we need a txt file for each block giving stimulus file names
%for EVERY trial in the block. This because each stimulus presentation will
%rely on a random color map being applied to a shape in an MS-associates
%cohort, wow, ok
% -----------------------------
%|          .txt file          |
%| trigs       stims           |
%|1 10 5 2 Stim1.bmp Stim1.wav |
%|2 8 14 8 Stim2.bmp Stim2.wav |
%|2 2 7 12 Stim3.bmp Stim3.wav |
%|etc...                       |
%|                             |
%|                             |
%|                             |
%-------------------------------
%just gonna start
clear
nstim = 15;
nreps = 50;

blocks_ncorrect = 5;
blocks_nwrong = 5;
blocks_nodd = 1;
blocks_nreps = 10;

corrpairs = {'VS01' '14';
    'VS02' '01';
    'VS03' '02';
    'VS04' '12';
    'VS05' '10';
    'VS06' '11';
    'VS07' '09';
    'VS08' '07';
    'VS09' '13';
    'VS10' '15';
    'VS11' '08';
    'VS12' '06';
    'VS13' '05';
    'VS14' '04';
    'VS15' '03';};

total_blocks = (nstim*nreps)/(blocks_ncorrect*blocks_nreps);

stimdir = uigetdir;
bmp_filelist = ls(fullfile(stimdir,'*bmp'));
wav_filelist = ls(fullfile(stimdir,'*wav'));

pat = '_+'; %This and pretty much all regexp calls will work only with the present NOPA filenaming structure (4/3/17)
for i = 1:length(bmp_filelist);  
    imfeats(i,:) = regexp(bmp_filelist(i,:),pat,'split');
end
im_forms = cell2mat(imfeats(:,1));
form_feats = unique(imfeats(:,1));
color_feats = unique(imfeats(:,2));

pat = '0(.*)\.'; % I kida love regexp
for i = 1:length(wav_filelist);
   temp = regexp(wav_filelist(i,:),pat,'tokens');
   soundfeats(i,:) = temp{:};
end
aud_feats = cell2mat(soundfeats);

blocksneeded = repmat(nreps/blocks_nreps,nstim,1);
prev_correct = zeros(blocks_ncorrect,1);
oddidx = zeros(total_blocks,1);
fudge = 5;
oddmode = 0; %0 = visual, 1 = auditory, this is just the starting point, will flip each block

for i = 1:total_blocks
    
    curr_correct = zeros(blocks_ncorrect,1);
    for j = 1:blocks_ncorrect
        candi = randi(15,1);
        while ~isempty(find(candi==prev_correct)) || ~isempty(find(candi==curr_correct)) || blocksneeded(candi)<fudge
            candi = randi(15,1);
        end
        curr_correct(j) = candi;
        blocksneeded(candi) = blocksneeded(candi)-1;
    end
    if mod(i,3) == 0
        fudge = fudge - 1;
    end
    correct_pairs{i,:} = corrpairs(curr_correct,:);
    prev_correct = curr_correct;
    
    for j = 1:blocks_nodd
        candi = randi(15,1);
        while ~isempty(find(candi==curr_correct)) || ~isempty(find(candi==oddidx))
            candi = randi(15,1);
        end
        oddidx(i) = candi;
    end
    
    seats_taken =[curr_correct' oddidx(i)];
    curr_wrong = 1:1:nstim;
    curr_wrong(ismember(curr_wrong,seats_taken)) = [];
    wrong_pairs{i,:} = corrpairs(curr_wrong,:);
    
    if ~oddmode;
        temp_free = [wrong_pairs{i}(:,2)' corrpairs(oddidx(i),2)];
        ix = randperm(length(temp_free),length(temp_free));
        dandi = temp_free(ix(1:end-1));
        
        while ~isempty(find(strcmp(dandi,wrong_pairs{i}(:,2)')))
            ix = randperm(length(temp_free),length(temp_free));
            dandi = temp_free(ix(1:end-1));
        end
        
        wrong_pairs{i}(:,2) = dandi';
        shuff = randperm(length(wrong_pairs{i}),blocks_nwrong);
        wrong_pairs{i} = wrong_pairs{i}(shuff,:);
        
        oddstim{i} = corrpairs(oddidx(i),1);
        image_locs = find(strcmp(oddstim{i},im_forms));
        theseones = randi(length(image_locs),1,blocks_nreps);
        oddimages = bmp_filelist(image_locs(theseones),:);
        theseones = randi(length(aud_feats),1,blocks_nreps);
        oddaudios = wav_filelist(theseones,:);
        thisone = randi(length(image_locs),1,1);
        oddtrainer = bmp_filelist(image_locs(thisone),:);
        oddmode = ~oddmode;
    else
        temp_free = [wrong_pairs{i}(:,1)' corrpairs(oddidx(i),1)];
        ix = randperm(length(temp_free),length(temp_free));
        dandi = temp_free(ix(1:end-1));
        
        while ~isempty(find(strcmp(dandi,wrong_pairs{i}(:,1)')))
            ix = randperm(length(temp_free),length(temp_free));
            dandi = temp_free(ix(1:end-1));
        end
        
        wrong_pairs{i}(:,1) = dandi';
        shuff = randperm(length(wrong_pairs{i}),blocks_nwrong);
        wrong_pairs{i} = wrong_pairs{i}(shuff,:);
        
        oddstim{i} = corrpairs(oddidx(i),2);
        theseones = randi(length(bmp_filelist),1,blocks_nreps);
        oddimages = bmp_filelist(theseones,:);
        oddaudios = repmat(wav_filelist(find(strcmp(oddstim{i},aud_feats)),:),blocks_nreps,1);
        oddtrainer = wav_filelist(find(strcmp(oddstim{i},aud_feats)),:);
        oddmode = ~oddmode;
    end
    fullset{i,:} = cat(1,correct_pairs{i},wrong_pairs{i});
    codeset = cat(1,repmat('1',length(correct_pairs{i}),1),repmat('2',length(wrong_pairs{i}),1));
   
    imageout = oddimages;
    audioout = oddaudios;
    code1 = repmat('3',blocks_nreps,1);
    code2 = [];
    code4 = [];
    for j = 1:blocks_nreps
        pat = '\_+';
        temp = regexp(oddimages(j,:),pat,'split');
        code2(j,1) = find(strcmp(temp{1},form_feats));
        pat = '0(.*)\.';
        temp = regexp(oddaudios(j,:),pat,'tokens');
        code4(j,1) = find(strcmp(temp{1},aud_feats));
    end
    
    for j = 1:length(fullset{i})
        image_locs = find(strcmp(fullset{i}(j,1),im_forms));
        theseones = randi(length(image_locs),1,blocks_nreps);
        imageout = [imageout; bmp_filelist(image_locs(theseones),:)]; 
        audioout = [audioout; repmat(wav_filelist(find(strcmp(fullset{i}(j,2),aud_feats)),:),blocks_nreps,1)];
        code1 = [code1; repmat(codeset(j,:),blocks_nreps,1)];
        code2 = [code2; repmat(find(strcmp(fullset{i}(j,1),form_feats)),blocks_nreps,1)];
        code4 = [code4; repmat(find(strcmp(fullset{i}(j,2),aud_feats)),blocks_nreps,1)];
    end
    code3 = [];
    for j = 1:length(imageout)
        pat = '\_+';
        temp = regexp(imageout(j,:),pat,'split');
        code3(j) = find(strcmp(temp{2},color_feats));
    end
    
    check = code1;
    while strcmp('3',check(1))
        check = code1;
        shuff = randperm(length(imageout),length(imageout));
        check = check(shuff);
    end
    
    code1 = code1(shuff); code2 = code2(shuff); code3 = code3(shuff); code4 = code4(shuff);
    imageout = imageout(shuff,:); audioout = audioout(shuff,:);
    outfiles(i,:) = {code1; code2; code3'; code4; imageout; audioout};
    
    fileID = fopen(['Block' num2str(i) '.txt'],'w');
    fprintf(fileID,'%d %s\r\n',~oddmode,oddtrainer);
    for j = 1:length(outfiles{1,1})
        fprintf(fileID,'%s %d %d %d %s %s\r\n',outfiles{i,1}(j),outfiles{i,2}(j),outfiles{i,3}(j),outfiles{i,4}(j),...
            outfiles{i,5}(j,:),outfiles{i,6}(j,:));
    end
    fclose(fileID);
end