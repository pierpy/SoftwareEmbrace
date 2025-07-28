% Open emt files of BTS system and track the ball for the rally times.
% Output: in Trial are the time-points of the period of play (Trial{k} refers to the k-th file)    

[file,path] = uigetfile('.emt');

filelist = dir(path);
kfile=0;

for jjk = 3:length(filelist)    %% Loop over the emt files
    file = filelist(jjk).name;

    if strcmp(file(end-2:end),'emt') & ~strcmp(file(end-8:end-4),'_STOP')
        kfile = kfile+1;
        fname = strcat(path,file);
        ifp = fopen(fname,'r');

        for k=1:25
            fscanf(ifp,'%s',1)
        end

        k = 1;
        while ~feof(ifp)
            x = fscanf(ifp,'%d\t',1);
            x
            if ~isempty(x)
                frame(k) = x;
                Time(k) = fscanf(ifp,'%f\t',1);
                X(k) = fscanf(ifp,'%f\t',1);
                Y(k) = fscanf(ifp,'%f\t',1);
                Z(k) = fscanf(ifp,'%f\t',1);
                k = k+1;
            end
        end
        fclose('all')

        [pks,locs] = findpeaks(Y, 'MinPeakHeight',0.16);
     
        sold = 220;
        ydiff = diff(locs);

        ind = find(ydiff>sold);


        Rally = [];
        for k=2:length(ind)
            Rally{k} = locs(ind(k-1)+1):locs(ind(k));
            Rally{k}(find(Y(Rally{k})<0)) = [];
        end
        Rally{1} = locs(1):locs(ind(1));
        nn = length(Rally)+1;
        Rally{nn} = locs(ind(end)+1):locs(end);

        for k=length(Rally):-1:1
            if length(Rally{k})<=3
                Rally{k} = [];
            end
        end

        ind1 = [];
        for k=1:length(Rally)
            ind1 = [ind1 Rally{k}];
        end

        for k=1:length(Rally)
            if ~isempty(Rally{k})
                if Rally{k}(end)-Rally{k}(1) < 250
                    Rally{k}=[];
                end
            end
        end

        close all
        figure
        plot(Time,Y,'b')
        hold on
        for k=1:length(Rally)
            ind1 = Rally{k};
            plot(Time(ind1),Y(ind1),'r')

        end
        hold off
        clear X Y Z
        fn1 = strcat(fname(1:end-8),'START_STOP.emt');
        ifp = fopen(fn1,'r');

        for k=1:25
            fscanf(ifp,'%s',1)
        end

        k = 1;
        while ~feof(ifp)
            x = fscanf(ifp,'%d\t',1);
            x
            if ~isempty(x)
                Time1(k) = fscanf(ifp,'%f\t',1);
                Start(k) = fscanf(ifp,'%f\t',1);
                End(k) = fscanf(ifp,'%f\t',1);
                k = k+1;
            end
        end
        fclose('all')
        ind = find(diff(End)>1);
        Offset = find(abs((Time1(ind(1))-Time))<0.002);
        if length(Offset>1)
            Offset = Offset(1);
        end
        for k=1:length(Rally)
            if ~isempty(Rally{k})
                Rally{k} = Rally{k}-Offset;
            end
        end

        TimeLimit = [];
        for k=1:length(Rally)
            if ~isempty(Rally{k})
                TimeLimit = [TimeLimit; Rally{k}(1) Rally{k}(end)];
            end
            Trial{kfile}.TimeLimit = TimeLimit;
            Trial{kfile}.Cond = file(18:end-9);
        end
        clear TimeLimit Time Time1 Rally
    end
end
