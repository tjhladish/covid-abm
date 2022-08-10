update par set seed = realization + 1;

--select quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con, count(*) from par group by quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con;

-- delete particles with ring vax but no baseline vax
delete from par where inf_con == 2 and (pas_alloc < 3 and act_alloc < 3);
delete from par where (pas_vac == 0 and pas_alloc != 0) or (pas_vac != 0 and pas_alloc == 0);
delete from par where (act_vac == 0 and act_alloc != 0) or (act_vac != 0 and act_alloc == 0);
delete from par where (pas_alloc == 3 or act_alloc == 3) and ((pas_vac == 0 and act_vac == 0) or (pas_vac != 0 and act_vac != 0));
delete from par where pas_vac == 0 and (act_alloc == 1 or act_alloc == 2);
delete from par where (act_vac == 1 and act_alloc == 2) or (act_vac == 2 and act_alloc == 1);
delete from par where act_vac != 0 and pas_alloc == 2;

delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);

update job set status = 'P' where serial in (select serial from par where pas_alloc == 2 or act_alloc == 2);

vacuum;
analyze;
