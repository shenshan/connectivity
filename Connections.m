%{
connectivity.Connections (computed) # create a GUI to input connectivity
-> connectivity.Slice
-----

%}

classdef Connections < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = connectivity.Slice
    end
    
	methods(Access=protected)
        
		function makeTuples(self, key)
            
            % fetch all the cells in this slice
            cells = fetch(connectivity.Cell & key);
            
            % GUI to insert values
            cnames = cell(1,length(cells));
            rnames = cell(1,length(cells));
            for ii = 1:length(cells)
                cnames{ii} = num2str(cells(ii).cell_id);
                rnames{ii} = num2str(cells(ii).cell_id);
            end
            data = zeros(length(cells),length(cells));
            for ii = 1:length(cells)
                data(ii,ii) = nan;
            end
            f = figure('Position',[200,200,650,250]);
            t = uitable('Parent',f,'Data',data,'ColumnName',cnames,'RowName',rnames,'Position',[20,20,600,200],'ColumnEditable',true(1,length(cells)));
            fetch1(connectivity.Slice & key, 'slice_name')
            in = input('Please press Enter to continue:');
            data = get(t,'Data');
            close(f)
            
            % insert tuples to database
            cnt = 0;
            for ii = 1:length(cells)
                for jj = 1:length(cells)
                    if isnan(data(ii,jj))
                        continue
                    end
                    cnt = cnt+1;
                    tuple1 = key;
                    tuple2 = key;
                    tuple = key;
                    tuple.pair_id = cnt;
                    tuple.connected = data(ii,jj);
                    
                    tuple1.cell_id = cells(ii).cell_id;
                    tuple1.role = 'from';
                    tuple1.pair_id = cnt;
                    tuple2.cell_id = cells(jj).cell_id;
                    tuple2.role = 'to';
                    tuple2.pair_id = cnt;
                    
                    if data(ii,jj)==2
                        tuple1.role = 'EC';
                        tuple2.role = 'EC';
                        data(jj,ii) = NaN;
                    end
                    insert(connectivity.CellTestedPair,tuple)
                    insert(connectivity.ConnectMembership,tuple1)
                    insert(connectivity.ConnectMembership,tuple2)
                end
            end
            self.insert(key)       
        end
        
    end
    
    

end