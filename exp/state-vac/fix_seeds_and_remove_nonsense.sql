update par set seed = realization + 1;

-- delete extra particles
delete from par where pas_vac == 0 and act_vac > 0;
delete from par where act_vac > 0 and state != 0;
delete from par where pas_vac == 0 and act_vac == 0 and state != 0;

delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);

vacuum;
analyze;
