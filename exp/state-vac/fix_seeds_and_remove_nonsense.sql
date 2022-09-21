update par set seed = realization + 1;

--select quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con, count(*) from par group by quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con;

-- delete particles with ring vax but no baseline vax
delete from par where pas_vac > 0 and act_vac > 0;
delete from par where pas_vac == 0 and state > 0;
delete from par where pas_vac != pas_alloc;
delete from par where act_vac != act_alloc;

delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);

vacuum;
analyze;
