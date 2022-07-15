update par set seed = realization + 1;

-- delete particles with ring vax but no baseline vax
delete from par where vac == 0 and (active_vax_strat != 0 or dose_file != 1);
delete from par where vac == 1 and active_vax_strat == 0 and dose_file != 1;

delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
