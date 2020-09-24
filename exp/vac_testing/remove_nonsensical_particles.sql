delete from par where vac = 0 and (vac_eff != 0.1 or vac_cov != 0.5);
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
