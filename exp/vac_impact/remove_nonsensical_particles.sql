delete from par where vac = 0 and (vac_rate != 0);
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
update par set seed = realization + 1;
