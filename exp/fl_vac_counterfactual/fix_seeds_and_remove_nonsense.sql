update par set seed = 1 + (realization % 100);
delete from par where state != 0 and vac == 0;
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
