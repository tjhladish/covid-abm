DROP TABLE IF EXISTS infection_history;
CREATE TABLE infection_history(
    inf TEXT PRIMARY KEY,
    inf_place_id INTEGER,
    inf_by_id INTEGER,
    inf_owner_id INTEGER,
    infected_time INTEGER,
    infectious_start INTEGER,
    infectious_end INTEGER,
    symp_start INTEGER,
    symp_end INTEGER,
    sev_start INTEGER,
    sev_end INTEGER,
    hosp_time INTEGER,
    crit_start INTEGER,
    crit_end INTEGER,
    icu_time INTEGER,
    death_time INTEGER,
    strain TEXT,
    rel_infectiousness INTEGER,
    detected TEXT,
    num_sec_infs INTEGER,
    FOREIGN KEY(detected) REFERENCES infection_detection(inf_det)
);

DROP TABLE IF EXISTS secondary_infections;
CREATE TABLE secondary_infections(
    parent_inf TEXT,
    offspring_inf TEXT,
    FOREIGN KEY(parent_inf) REFERENCES infection_history(inf),
    FOREIGN KEY(offspring_inf) REFERENCES infection_history(inf)
);

DROP TABLE IF EXISTS infection_detection;
CREATE TABLE infection_detection(
    inf_det TEXT PRIMARY KEY,
    inf TEXT,
    detected_state TEXT,
    reported_time INTEGER,
    FOREIGN KEY(inf) REFERENCES infection_history(inf)
);

.mode csv
