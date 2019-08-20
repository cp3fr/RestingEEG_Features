function [EEG,isvalid] = fix_event_structure(EEG)

tbl = struct2table(EEG.event);
tbl.delay = diff([tbl.latency;EEG.pnts]);
tbl.select = ones(size(tbl,1),1);

ind = contains(tbl.type,'90');
idx = find(ind,1,'last');

if idx~=1
  tbl.select(1:idx-1)=0;
end

ind = contains(tbl.type,'20') & tbl.delay<8000;

tbl.select(ind)=0;

tbl = tbl(logical(tbl.select),1:end-2);

EEG.event = table2struct(tbl);

if sum(contains(tbl.type,'20'))==5 && sum(contains(tbl.type,'30'))==5
  isvalid = true;
else
  isvalid = false;
  tbl
end