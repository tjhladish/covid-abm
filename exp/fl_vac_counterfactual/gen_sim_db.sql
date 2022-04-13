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
    FOREIGN KEY(detected) REFERENCES infection_detection(inf_det),
    FOREIGN KEY(inf_owner_id) REFERENCES vaccination_history(p_id)
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

DROP TABLE IF EXISTS vaccination_history;
CREATE TABLE vaccination_history(
    p_id INTEGER,
    p_age_bin INTEGER,
    dose INTEGER,
    vax_sim_day INTEGER,
    vax_date TEXT,
    FOREIGN KEY(p_id) REFERENCES infection_history(inf_owner_id),
    FOREIGN KEY(p_age_bin) REFERENCES age_bins(bin_min),
    FOREIGN KEY(vax_sim_day) REFERENCES doses_available(sim_day),
    FOREIGN KEY(vax_date) REFERENCES doses_available(date_str)
);

DROP TABLE IF EXISTS age_bins;
CREATE TABLE age_bins(
    bin_min INTEGER,
    bin_pop INTEGER
);

DROP TABLE IF EXISTS doses_available;
CREATE TABLE doses_available (
    sim_day INTEGER,
    date_str TEXT,
    dose INTEGER,
    bin INTEGER,
    doses INTEGER,
    FOREIGN KEY(sim_day) REFERENCES vaccination_history(vax_sim_day),
    FOREIGN KEY(date_str) REFERENCES vaccination_history(vax_date)
);

.mode csv
