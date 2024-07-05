tic

clc
clear all
close all

% dosya okuma icin
listing = dir('Kimia99_DB');
num_rec = zeros(9,11); 

%---------------------------------------------------
% Egitim asamasi
%---------------------------------------------------
for train_no = 1:11;
    
    train_no
    
    for class_indx = 1:9
        
        class_indx
        
        tmp_indx = (class_indx-1)*11 + train_no;
        tmp_indx = tmp_indx+2;
        filename = listing(tmp_indx).name;
        
        bw = im2bw(imread(sprintf('Kimia99_DB/%s',filename)),0); %binary
        
        tmp_bw = imresize(bw, [64 64]);
        bw1{class_indx} = zeros(64,64); % Her sýnýfýn eðitim imgesini saklýyor.
        bw1{class_indx} = tmp_bw;
        
        s = regionprops(bw1{class_indx},{...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation',...
            'PixelIdxList',...
            'PixelList',...
            'Area'});
        
        % Gürültüden dolayý birden fazla elips çýktýðýndan en büyük alana
        % sahip elips seçimi için kullanýlýyor.
        for num_ellipse = 1:size(s,1)
            size_ellipse(num_ellipse) = s(num_ellipse).Area;
        end
        [max_area,max_indx] = max(size_ellipse);
               
       
        %---------------------------------
        % Minör eksen köþeleri bulmak için
        %---------------------------------
        vertex(1,1) = s(max_indx,1).Centroid(1)-(s(max_indx,1).MinorAxisLength/2)*sind(s(max_indx,1).Orientation);
        vertex(1,2) = s(max_indx,1).Centroid(2)-(s(max_indx,1).MinorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        vertex(1,3) = s(max_indx,1).Centroid(1)+(s(max_indx,1).MinorAxisLength/2)*sind(s(max_indx,1).Orientation);
        vertex(1,4) = s(max_indx,1).Centroid(2)+(s(max_indx,1).MinorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        
        %---------------------------------
        % Majör eksen köþeleri bulmak için
        %---------------------------------
        vertex(1,5) = s(max_indx,1).Centroid(1)+(s(max_indx,1).MajorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        vertex(1,6) = s(max_indx,1).Centroid(2)-(s(max_indx,1).MajorAxisLength/2)*sind(s(max_indx,1).Orientation);
        vertex(1,7) = s(max_indx,1).Centroid(1)-(s(max_indx,1).MajorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        vertex(1,8) = s(max_indx,1).Centroid(2)+(s(max_indx,1).MajorAxisLength/2)*sind(s(max_indx,1).Orientation);
        vr{class_indx} = [vertex(1,1) vertex(1,2); vertex(1,5) vertex(1,6); vertex(1,3) vertex(1,4); vertex(1,7) vertex(1,8)];...
            clear s; clear size_ellipse;
    end
    
    %---------------------------------------------------
    % Test asamasi
    %---------------------------------------------------
    % Orjinal birinci resim için ayný iþlemlerin tekrarý
    %---------------------------------------------------
    
    tst_no = 0;
    
    for class_indx = 1:9
        
        for tst_indx = 1:11
            
            if tst_indx == train_no
                continue
            end
            
            % Ýmgeleri sýrasý ile çekme
            tst_no = tst_no+1;
            tmp_indx = (class_indx-1)*11 + tst_indx;
            tmp_indx = tmp_indx+2;
            filename = listing(tmp_indx).name;
            bw = im2bw(imread(sprintf('Kimia99_DB/%s',filename)),0);
            tmp_bw = imresize(bw, [64 64]);
            bw2{tst_no} = zeros(64,64);
            bw2{tst_no} = tmp_bw;
            
            s = regionprops(bw2{tst_no},{...
                'Centroid',...
                'MajorAxisLength',...
                'MinorAxisLength',...
                'Orientation',...
                'PixelIdxList',...
                'PixelList',...
                'Area'});
            
            % Ýmge içinde çýkan muhtemel gürültülere uydurulan elipsleri
            % ihmal etmek için
            size_ellipse = 0;
            for num_ellipse = 1:size(s,1)
                size_ellipse(num_ellipse) = s(num_ellipse).Area;
            end
            [max_area,max_indx] = max(size_ellipse);
            
            %---------------------------------
            % Minör eksen köþeleri bulmak için
            %---------------------------------
            orgvertex(1,1) = s(max_indx,1).Centroid(1)-(s(max_indx,1).MinorAxisLength/2)*sind(s(max_indx,1).Orientation);
            orgvertex(1,2) = s(max_indx,1).Centroid(2)-(s(max_indx,1).MinorAxisLength/2)*cosd(s(max_indx,1).Orientation);
            orgvertex(1,3) = s(max_indx,1).Centroid(1)+(s(max_indx,1).MinorAxisLength/2)*sind(s(max_indx,1).Orientation);
            orgvertex(1,4) = s(max_indx,1).Centroid(2)+(s(max_indx,1).MinorAxisLength/2)*cosd(s(max_indx,1).Orientation);
            %---------------------------------
            % Majör eksen köþeleri bulmak için
            %---------------------------------
            orgvertex(1,5) = s(max_indx,1).Centroid(1)+(s(max_indx,1).MajorAxisLength/2)*cosd(s(max_indx,1).Orientation);
            orgvertex(1,6) = s(max_indx,1).Centroid(2)-(s(max_indx,1).MajorAxisLength/2)*sind(s(max_indx,1).Orientation);
            orgvertex(1,7) = s(max_indx,1).Centroid(1)-(s(max_indx,1).MajorAxisLength/2)*cosd(s(max_indx,1).Orientation);
            orgvertex(1,8) = s(max_indx,1).Centroid(2)+(s(max_indx,1).MajorAxisLength/2)*sind(s(max_indx,1).Orientation);
            %             orgvertex = round(orgvertex);
            %orgvr: test
            orgvr{tst_no} = [orgvertex(1,1) orgvertex(1,2); orgvertex(1,5) orgvertex(1,6); orgvertex(1,3) orgvertex(1,4); orgvertex(1,7) orgvertex(1,8)];
            
            %------------------------------------
            % M Dönüþüm matrisinin hesaplanmasý
            %------------------------------------
            cntr_indx = zeros(9,4);
            
            for org_class_indx = 1:9
                imshow(bw2{tst_no},'InitialMagnification','fit')
                hold on
                
                cs_orgvr{1} = orgvr{tst_no};
                plot(orgvr{1,1}(:,1),orgvr{1,1}(:,2),'yo','Linewidth',5)
                
                % Test imgesinin vertexlerini dairesel olarak kaydýrýyoruz.
                for i = 2:4
                    cs_orgvr{1,i} = circshift(cs_orgvr{i-1},1);
                end
                % Eðitim imgesinin vertexlerini dairesel olarak kaydýrýyoruz.
                cs_vr{1} = vr{org_class_indx}(1:4,:);
                for i = 2:4
                    cs_vr{1,i} = circshift(cs_vr{i-1},1);
                end
                
                tmp_k = 0;
                for i = 1:1
                    coor_1 = cs_orgvr{1,i};
                    coor_im1 = coor_1(1:3,:)';
                    coor_im1(3,:) = 1;
                    % Dört farklý köþe sýralamasý deneniyor.
                    for j = 1:4
                        tmp_k = tmp_k+1;
                        coor_2 = cs_vr{1,j};
                        coor_im2 = coor_2(1:3,:)';
                        coor_im2(3,:) = 1;
                        
                        M = coor_im1 * pinv(coor_im2); %coor_im1: Test imgesinin köþe koordinatlarý
                        tst_coor2 = M*coor_im2; % coor_im2: Eðitim imgesinin köþe koordinatlarý
                        
                        org_s = regionprops(bw1{org_class_indx},{...
                            'Centroid',...
                            'MajorAxisLength',...
                            'MinorAxisLength',...
                            'Orientation',...
                            'PixelIdxList',...
                            'PixelList',...
                            'Area'...
                            'SubarrayIdx'});
                        
                        for num_ellipse = 1:size(org_s,1)
                            org_size_ellipse(num_ellipse) = org_s(num_ellipse).Area;
                        end
                        [org_max_area,org_max_indx] = max(org_size_ellipse);
                        org_image = bw1{org_class_indx};
                        
                        hold on
                        plot(vr{1,org_class_indx}(:,1),vr{1,org_class_indx}(:,2),'g+','Linewidth',5)
                        
                        tmp_coor = tst_coor2';
                        plot(tmp_coor(:,1),tmp_coor(:,2),'bp','Linewidth',2)
                        
                        PL = org_s(org_max_indx).PixelList; % pixel coordinates of the training image
                        PL2 = s(max_indx).PixelList; % pixel coordinates of the test image
                        
                        ln1 = length(PL); % the number of the white pixels in training image
                        ln2 = length(PL2); % the number of the white pixels in test image
                        PL_hom = [PL ones(ln1,1)];% adding ones
                        PL_hom_pro = round(10*M*PL_hom')'; % projection of the training image white pixels onto
                        % the test image's white pixels.
                        PL_pro = PL_hom_pro(:,1:2); % removing ones
                        PL_pro = unique(PL_pro,'rows'); % removing the repetitions
                        
                        PL_hom_pro_2 = round(M*PL_hom')'; % projection of the training image white pixels onto
                        % the test image's white pixels.
                        PL_pro_2 = PL_hom_pro_2(:,1:2); % removing ones
                        PL_pro_2 = unique(PL_pro_2,'rows'); % removing the repetitions
                        
                        I = find(PL_pro(:,1)>640);
                        PL_pro(I,1) = 640;
                        J = find(PL_pro(:,1)<10);
                        PL_pro(J,1) = 10;
                        I = find(PL_pro(:,2)>640);
                        PL_pro(I,2) = 640;
                        J = find(PL_pro(:,2)<10);
                        PL_pro(J,2) = 10;
                        
                        hold on
                        plot(PL_pro_2(:,1),PL_pro_2(:,2), 'r.' ) % ploting the projected points
                        ln_pro = length(PL_pro); % the number of white pixels after the projection
                        
                        % counting the number of overlapping white pixels
                        cntr = 0;
                        cntr_tmp = 0;
                        pixel_cntr = 0;
                        PL2_new = 10*PL2;
                        
                        for i=1:ln_pro
                            
                            for h = 1:ln2
                                
                                if ( (((PL2_new(h,1)-10) <= PL_pro(i,1)) & ((PL_pro(i,1) <= (PL2_new(h,1)+10)))) & ...
                                    (((PL2_new(h,2)-10) <= PL_pro(i,2)) & ((PL_pro(i,2) <= (PL2_new(h,2)+10)))) ) 
                                
                                cntr = cntr + 1;
                                cntr_tmp(cntr,1) = PL_pro(i,1); 
                                cntr_tmp(cntr,2) = PL_pro(i,2);                               
                                   
                               end
                           
                            end
                           
                        end
                        
                        cntr_remove = unique(cntr_tmp,'rows');                        
                        cntr_indx(org_class_indx,tmp_k) = length(cntr_remove)/ln1;
                        
                    end
                    
                end
                
                max_cntr = max(cntr_indx(org_class_indx,:));
                overlap_rate(org_class_indx,tst_no) = max_cntr;
                clear org_s; clear org_size_ellipse;
  
                figure   
                close all
            end
            
            tmp_rate = overlap_rate(:,tst_no);
            [max_val,max_cls_indx] = max(tmp_rate);
            rec_indx(class_indx,tst_indx) = max_cls_indx
            
            if max_cls_indx == class_indx
                num_rec(class_indx,train_no) = num_rec(class_indx,train_no) + 1
            end
            clear s, clear size_ellipse;
        end
    end  
    assigment_matrix{train_no} = rec_indx;
    
    elapsed_time = toc
end