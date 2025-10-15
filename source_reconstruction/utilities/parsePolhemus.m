function parsePolhemus(f_name)
%Enter the .hsp and .elp file name without extension as a string,
%E.g. parsehspPolhemus  ('001_gt_2009_06_26').
%Output is the shape.mat file saved in the current workspace.

pendStr1='.elp';
pendStr2='.hsp';

f_name1=[f_name,pendStr1];
f_name2=[f_name,pendStr2];

fid1 = fopen(f_name1);
C = fscanf(fid1,'%c');
fclose(fid1);

E = regexprep(C,'\n','xx'); 
E = regexprep(E,'\t','yy'); 

returnsi = strfind(E,'xx');
tabsi = strfind(E,'yy');
sensornamesi = strfind(E,'%N');
fiducialsstarti = strfind(E,'%F');
lastfidendi = strfind(E(fiducialsstarti(3):fiducialsstarti(length(fiducialsstarti))+100),'xx');
fiducialsendi = fiducialsstarti(1)+strfind(E(fiducialsstarti(1):fiducialsstarti(length(fiducialsstarti))+lastfidendi(1)),'xx');

NASION = E(fiducialsstarti(1)+4:fiducialsendi(1)-2);
NASION = regexprep(NASION,'yy','\t');
NASION = str2num(NASION);

LPA = E(fiducialsstarti(2)+4:fiducialsendi(2)-2);
LPA = regexprep(LPA,'yy','\t');
LPA = str2num(LPA);

RPA = E(fiducialsstarti(3)+4:fiducialsendi(3)-2);
RPA = regexprep(RPA,'yy','\t');
RPA = str2num(RPA);

LPAredstarti = strfind(E,'0-RED');
LPAredendi = strfind(E(LPAredstarti(1):LPAredstarti(length(LPAredstarti))+46),'xx');
LPAred = E(LPAredstarti(1)+11:LPAredstarti(1)+LPAredendi(2)-2);
LPAred = regexprep(LPAred,'yy','\t');
LPAred = str2num(LPAred);

RPAyelstarti = strfind(E,'1-YELLOW');
RPAyelendi = strfind(E(RPAyelstarti(1):RPAyelstarti(length(RPAyelstarti))+50),'xx');
% RPAyelendi = strfind(E(RPAyelstarti(1):RPAyelstarti(length(RPAyelstarti))+47),'xx');
RPAyel = E(RPAyelstarti(1)+14:RPAyelstarti(1)+RPAyelendi(2)-2);
RPAyel = regexprep(RPAyel,'yy','\t');
RPAyel = str2num(RPAyel);

PFbluestarti = strfind(E,'2-BLUE');
PFblueendi = strfind(E(PFbluestarti(1):PFbluestarti(length(PFbluestarti))+46),'xx');
PFblue = E(PFbluestarti(1)+11:PFbluestarti(1)+PFblueendi(2)-2);
PFblue = regexprep(PFblue,'yy','\t');
PFblue = str2num(PFblue);

LPFwhstarti = strfind(E,'3-WHITE');
LPFwhendi = strfind(E(LPFwhstarti(1):LPFwhstarti(length(LPFwhstarti))+46),'xx');
LPFwh = E(LPFwhstarti(1)+12:LPFwhstarti(1)+LPFwhendi(2)-2);
LPFwh = regexprep(LPFwh,'yy','\t');
LPFwh = str2num(LPFwh);

RPFblackstarti = strfind(E,'4-BLACK');
RPFblackendi = strfind(E(RPFblackstarti(1):end),'xx');
RPFblack = E(RPFblackstarti(1)+12:RPFblackstarti(1)+RPFblackendi(2)-2);
RPFblack = regexprep(RPFblack,'yy','\t');
RPFblack = str2num(RPFblack);

allfids = [NASION;LPA;RPA;LPAred;RPAyel;PFblue;LPFwh;RPFblack];
fidslabels = {'NASION';'LPA';'RPA';'LPAred';'RPAyel';'PFblue';'LPFwh';'RPFblack'};

fid2 = fopen(f_name2);
C = fscanf(fid2,'%c');
fclose(fid2);
E = regexprep(C,'\n','xx'); %replace returns with "xx"
E = regexprep(E,'\t','yy'); %replace tabs with "yy"
returnsi = strfind(E,'xx');
tabsi = strfind(E,'yy');

headshapestarti = strfind(E,'position of digitized points');
headshapestartii = strfind(E(headshapestarti(1):end),'xx');
headshape = E(headshapestarti(1)+headshapestartii(2)+2:end);
headshape = regexprep(headshape,'yy','\t');
headshape = regexprep(headshape,'xx',';');
headshape = str2num(headshape);

shape.pnt = headshape;
shape.fid.pnt = allfids;
shape.fid.label = fidslabels;

%convert to BESA style coordinates so can use the .pos file or sensor
%config from .con
shape.pnt = cat(2,fliplr(shape.pnt(:,1:2)),shape.pnt(:,3)).*1000;
shape.pnt = shape.pnt(1:length(shape.pnt)-15,:); % get rid of nose points may want to alter or comment this depending on your digitisation
%shape.pnt = shape.pnt*1000;
neg = shape.pnt(:,2)*-1;
shape.pnt(:,2) = neg;

shape.fid.pnt = cat(2,fliplr(shape.fid.pnt(:,1:2)),shape.fid.pnt(:,3)).*1000;
%shape.fid.pnt = shape.fid.pnt*1000;
neg2 = shape.fid.pnt(:,2)*-1;
shape.fid.pnt(:,2) = neg2;
shape.unit='mm';


new_name2 = ['shape.mat'];
save (new_name2,'shape');
clear all



