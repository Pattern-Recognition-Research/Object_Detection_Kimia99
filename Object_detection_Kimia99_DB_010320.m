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
        
        bw = im2bw(imread(sprintf('Kimia99_DB/%s',filename)),0);
        
        tmp_bw = imresize(bw, [256 256]);
        bw1{class_indx} = zeros(256,256);
        bw1{class_indx} = tmp_bw;
        image_size = size(bw1{class_indx});
%             imshow(bw1{class_indx})
        
        s = regionprops(bw1{class_indx},{...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation',...
            'PixelIdxList',...
            'PixelList',...
            'Area'});
        
        for num_ellipse = 1:size(s,1)
            
            size_ellipse(num_ellipse) = s(num_ellipse).Area;
        end
        
        [max_area,max_indx] = max(size_ellipse);
        
        t = linspace(0,2*pi,100);
        coef = 1;
        hold on
        
        % for k = 1:length(s)
        a = s(max_indx).MajorAxisLength/2;
        b = s(max_indx).MinorAxisLength/2;
        Xc = s(max_indx).Centroid(1);
        Yc = s(max_indx).Centroid(2);
        phi = deg2rad(-s(max_indx).Orientation);
        
        a = a*coef;
        b = b*coef;
        
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi;
        hold on
        
        hammer = s(max_indx).PixelList;
        cnt = s(max_indx).Centroid;
        
        %---------------------------------
        % Minör eksen köþeleri bulmak için
        %---------------------------------
        vertex(1,1) = s(max_indx,1).Centroid(1)-(s(max_indx,1).MinorAxisLength/2)*sind(s(max_indx,1).Orientation);
        vertex(1,2) = s(max_indx,1).Centroid(2)-(s(max_indx,1).MinorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        % vertex(2,1) = s(2,1).Centroid(1)-(s(2,1).MinorAxisLength/2)*sind(s(2,1).Orientation);
        % vertex(2,2) = s(2,1).Centroid(2)-(s(2,1).MinorAxisLength/2)*cosd(s(2,1).Orientation);
        
        vertex(1,3) = s(max_indx,1).Centroid(1)+(s(max_indx,1).MinorAxisLength/2)*sind(s(max_indx,1).Orientation);
        vertex(1,4) = s(max_indx,1).Centroid(2)+(s(max_indx,1).MinorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        % vertex(2,3) = s(2,1).Centroid(1)+(s(2,1).MinorAxisLength/2)*sind(s(2,1).Orientation);
        % vertex(2,4) = s(2,1).Centroid(2)+(s(2,1).MinorAxisLength/2)*cosd(s(2,1).Orientation);
        
        %---------------------------------
        % Majör eksen köþeleri bulmak için
        %---------------------------------
        
        vertex(1,5) = s(max_indx,1).Centroid(1)+(s(max_indx,1).MajorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        vertex(1,6) = s(max_indx,1).Centroid(2)-(s(max_indx,1).MajorAxisLength/2)*sind(s(max_indx,1).Orientation);
        % vertex(2,5) = s(2,1).Centroid(1)+(s(2,1).MajorAxisLength/2)*cosd(s(2,1).Orientation);
        % vertex(2,6) = s(2,1).Centroid(2)-(s(2,1).MajorAxisLength/2)*sind(s(2,1).Orientation);
        
        vertex(1,7) = s(max_indx,1).Centroid(1)-(s(max_indx,1).MajorAxisLength/2)*cosd(s(max_indx,1).Orientation);
        vertex(1,8) = s(max_indx,1).Centroid(2)+(s(max_indx,1).MajorAxisLength/2)*sind(s(max_indx,1).Orientation);
        % vertex(2,7) = s(2,1).Centroid(1)-(s(2,1).MajorAxisLength/2)*cosd(s(2,1).Orientation);
        % vertex(2,8) = s(2,1).Centroid(2)+(s(2,1).MajorAxisLength/2)*sind(s(2,1).Orientation);
        
        vertex = round(vertex);
        
        vr{class_indx} = [vertex(1,1) vertex(1,2); vertex(1,5) vertex(1,6); vertex(1,3) vertex(1,4); vertex(1,7) vertex(1,8)];...
            %     vertex(2,1) vertex(2,2); vertex(2,5) vertex(2,6); vertex(2,3) vertex(2,4); vertex(2,7) vertex(2,8)];
        
        %     plot(vr{class_indx}(:,1),vr{class_indx}(:,2),'yo','Linewidth',2)
        
        %---------------------------------
        % Elipslerin alanýný bulmak için
        %---------------------------------
        % area(1,1) = pi*s(1,1).MinorAxisLength*s(1,1).MajorAxisLength/4;
        % area(2,1) = pi*s(2,1).MinorAxisLength*s(2,1).MajorAxisLength/4;
        
        clear s; clear size_ellipse;
        
    end % class_indx = 1:70
    
    %---------------------------------------------------
    % Test asamasi
    %---------------------------------------------------
    % Orjinal birinci resim için ayný iþlemlerin tekrarý
    %---------------------------------------------------
    
    tst_no = 0;
    
    for class_indx = 5:5 %1:9
        
        for tst_indx = 1:1 %1:11
            
            if tst_indx == train_no
                continue
            end
            
            tst_no = tst_no+1
            tmp_indx = (class_indx-1)*11 + tst_indx;
            tmp_indx = tmp_indx+2;
            filename = listing(tmp_indx).name;
            
            bw = im2bw(imread(sprintf('Kimia99_DB/%s',filename)),0);
            
            tmp_bw = imresize(bw, [256 256]);
            bw2{tst_no} = zeros(256,256);
            bw2{tst_no} = tmp_bw;
            
            %         bw2{tst_no} = imresize(bw2{tst_no}, [256 256]);
            
            image_size = size(bw2{tst_no});
            
                    imshow(bw)
            
            s = regionprops(bw2{tst_no},{...
                'Centroid',...
                'MajorAxisLength',...
                'MinorAxisLength',...
                'Orientation',...
                'PixelIdxList',...
                'PixelList',...
                'Area'});
            
            size_ellipse = 0;
            
            for num_ellipse = 1:size(s,1)
                
                size_ellipse(num_ellipse) = s(num_ellipse).Area;
                
            end
            
            [max_area,max_indx] = max(size_ellipse);
            
                   figure
                   imshow(bw2{tst_no},'InitialMagnification','fit')
            
            t = linspace(0,2*pi,100);
            coef = 1;
            hold on
            
            % for k = 1:length(s)
            a = s(max_indx).MajorAxisLength/2;
            b = s(max_indx).MinorAxisLength/2;
            Xc = s(max_indx).Centroid(1);
            Yc = s(max_indx).Centroid(2);
            phi = deg2rad(-s(max_indx).Orientation);
            
            a = a*coef;
            b = b*coef;
            
            x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
            y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
            
            plot(x,y,'r','Linewidth',5)
            
%             x_pts(k,:) = x;
%             y_pts(k,:) = y;
%             end
                    hold off
            
                    hold on
            
            hammer = s(max_indx).PixelList;
            cnt = s(max_indx).Centroid;
                    plot(hammer(:,1),hammer(:,2),'c.')
                    plot(cnt(1,1),cnt(1,2),'ro','Linewidth',2)
            
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
            orgvertex = round(orgvertex);
            
            orgvr{tst_no} = [orgvertex(1,1) orgvertex(1,2); orgvertex(1,5) orgvertex(1,6); orgvertex(1,3) orgvertex(1,4); orgvertex(1,7) orgvertex(1,8)];

            %---------------------------------
            % Elipslerin alanýný bulmak için
            %---------------------------------
            orgarea(1,1) = pi*s(max_indx,1).MinorAxisLength*s(max_indx,1).MajorAxisLength/4;
            
            %------------------------------------
            % M Dönüþüm matrisinin hesaplanmasý
            %------------------------------------
            
            cntr_indx = zeros(9,16);
            
            for org_class_indx = 6:6 %1:9
                
                cs_orgvr{1} = orgvr{tst_no};
                
                for i = 2:4
                    
                    cs_orgvr{1,i} = circshift(cs_orgvr{i-1},1);
                    
                end
                
                cs_vr{1} = vr{org_class_indx}(1:4,:);
                
                for i = 2:4
                    
                    cs_vr{1,i} = circshift(cs_vr{i-1},1);
                    
                end
                
                tmp_k = 0;
                
                for i = 1:1 %1:4
                    
                    coor_1 = cs_orgvr{1,i};
                    coor_im1 = coor_1(1:3,:)';
                    coor_im1(3,:) = 1;
                    
                    for j = 1:4
                        
                        tmp_k = tmp_k+1;
                        
                        coor_2 = cs_vr{1,j};
                        coor_im2 = coor_2(1:3,:)';
                        coor_im2(3,:) = 1;
                        
                        M = coor_im2 * inv(coor_im1);
                        
                        tst_coor2 = M*coor_im1;
                        
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
                        image_2 = bw2{tst_no};
                        
                        
                        figure
                        imshow(bw1{org_class_indx},'InitialMagnification','fit')
                        
                        t = linspace(0,2*pi,100);
                        coef = 1;
                        hold on

                        % for k = 1:length(s)
                        a = org_s(max_indx).MajorAxisLength/2;
                        b = org_s(max_indx).MinorAxisLength/2;
                        Xc = org_s(max_indx).Centroid(1);
                        Yc = org_s(max_indx).Centroid(2);
                        phi = deg2rad(-org_s(max_indx).Orientation);

                        a = a*coef;
                        b = b*coef;

                        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
                        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);

                        plot(x,y,'r','Linewidth',5)
                        
                        
                        PL = org_s(org_max_indx).PixelList; % pixel coordinates of the first image
                        PL2 = s(max_indx).PixelList; % pixel coordinates of the second image
                        
                        ln1 = length(PL); % the number of the white pixels in image 1
                        ln2 = length(PL2); % the number of the white pixels in image 2
                        
                        if ln1 > ln2 % It is used to determine the large image
                            
                            PL_hom = [PL ones(ln1,1)];% adding ones
                            PL_hom_pro = round(M*PL_hom')'; % projection of the second image white pixels onto
                            % the first image's white pixels.
                            
                            PL_pro = PL_hom_pro(:,1:2); % removing ones
                            
                            PL_pro = unique(PL_pro,'rows'); % removing the repetitions
                            
                           I = find(PL_pro(:,1)>256);
                            PL_pro(I,1) = 256;
                            J = find(PL_pro(:,1)<1);
                            PL_pro(I,1) = 1;
                            I = find(PL_pro(:,2)>256);
                            PL_pro(I,1) = 256;
                            J = find(PL_pro(:,2)<1);
                            PL_pro(I,1) = 1;
                            figure
                            imshow(bw2{tst_no},'InitialMagnification','fit')
                            hold on
                            
                            
                            plot(PL_pro(:,1),PL_pro(:,2), 'r.' ) % ploting the projected points
                            
                            ln_pro = length(PL_pro); % the number of white pixels after the projection
                            
                            % counting the number of overlapping white pixels
                            cntr = 0;
                            pixel_cntr = 0;
                            
                            for i=1:ln_pro
                                cntr = cntr + sum(PL2(:, 1) == PL_pro(i,1) & PL2(:, 2) == PL_pro(i,2));                    
                           end
                            
                            cntr_indx(org_class_indx,tmp_k) = cntr/ln_pro;
                            % overlap_rate(org_class_indx,tst_no) = cntr/size(PL,1);
                            
                        else
                            PL2_hom = [PL2 ones(ln2,1)];% adding ones
                            PL2_hom_pro = round(inv(M)*PL2_hom')'; % projection of the second image white pixels onto
                            % the first image's white pixels.
                            PL2_pro = PL2_hom_pro(:,1:2); % removing ones
                            PL2_pro = unique(PL2_pro,'rows'); % removing the repetitions

                            I = find(PL2_pro(:,1)>256);
                            PL2_pro(I,1) = 256;
                            J = find(PL2_pro(:,1)<1);
                            PL2_pro(I,1) = 1;
                            I = find(PL2_pro(:,2)>256);
                            PL2_pro(I,1) = 256;
                            J = find(PL2_pro(:,2)<1);
                            PL2_pro(I,1) = 1;
                            figure
                            imshow(bw2{tst_no},'InitialMagnification','fit')
                            hold on
                            plot(PL2_pro(:,1),PL2_pro(:,2), 'r.' ) % ploting the projected points
                            
                            ln2_pro = length(PL2_pro); % the number of white pixels after the projection
                            
                            % counting the number of overlapping white pixels
                            cntr = 0;
                            pixel_cntr = 0;
                            
                            for i=1:ln2_pro

                                cntr = cntr + sum(PL(:, 1) == PL2_pro(i,1) & PL(:, 2) == PL2_pro(i,2));
                            end
                            
                            cntr_indx(org_class_indx,tmp_k) = cntr/ln2_pro;

                            
                        end
                        
                    end
                    
                end
                
                max_cntr = max(cntr_indx(org_class_indx,:));
                overlap_rate(org_class_indx,tst_no) = max_cntr;
                
                clear org_s; clear org_size_ellipse;
            end
            
            tmp_rate = overlap_rate(:,tst_no);
            
            [max_val,max_cls_indx] = max(tmp_rate);
            rec_indx(class_indx,tst_indx) = max_cls_indx;
            
            if max_cls_indx == class_indx
                
                num_rec(class_indx,train_no) = num_rec(class_indx,train_no) + 1
            end
            
%             close all
            clear s, clear size_ellipse;
        end
        
    end
    
    assigment_matrix{train_no} = rec_indx;
    
    elapsed_time = toc
    
end