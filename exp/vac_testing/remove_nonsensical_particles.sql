delete from par where vac = 0 and (vac_eff != 0.5 or vac_cov != 0.5);
delete from par where pvl = 1 and vac_day = 0;
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
update par set seed = realization;
