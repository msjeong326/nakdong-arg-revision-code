import numpy as np
import pandas as pd

PCS = [f"Comm_PC{i}" for i in range(1, 11)]

EXPOSURES = {
    "Comm_PC1": ("Microbial community", "community"),
    "MC_ppb": ("Microcystin", "toxin"),
    "bloom": ("Bloom (binary)", "toxin"),
    "Chl_a_mg_m3": ("Chlorophyll-a", "toxin"),
    "water_temp_C": ("Water temperature", "env"),
    "TN_mg_L": ("Total nitrogen", "env"),
    "TP_mg_L": ("Total phosphorus", "env"),
    "DO_mg_L_wq": ("Dissolved oxygen", "env"),
    "pH_wq": ("pH", "env"),
    "BOD_mg_L": ("BOD", "env"),
    "COD_mg_L": ("COD", "env"),
    "EC_uS_cm": ("Conductivity", "env"),
    "NH3N_mg_L": ("Ammonia-N", "env"),
    "NO3N_mg_L": ("Nitrate-N", "env"),
    "TOC_mg_L": ("TOC", "env"),
    "SS_mg_L": ("Suspended solids", "env"),
    "air_temp_avg_C": ("Air temperature", "met"),
    "humidity_avg_pct": ("Humidity", "met"),
    "precip_sum_mm": ("Precipitation", "met"),
    "wind_speed_avg_ms": ("Wind speed", "met"),
    "sunshine_avg_hr": ("Sunshine", "met"),
    "pressure_avg_hPa": ("Pressure", "met"),
    "discharge_mean_m3s": ("Discharge", "hydro"),
    "storage_rate_pct": ("Storage rate", "hydro"),
}

CORE_CONTROLS = [
    "water_temp_C",
    "Chl_a_mg_m3",
    "TN_mg_L",
    "TP_mg_L",
    "DO_mg_L_wq",
    "pH_wq",
]

GCM_DRIVERS = [
    ("Bloom", "bloom", "toxin"),
    ("Microcystin", "MC_ppb", "toxin"),
    ("Chlorophyll a", "Chl_a_mg_m3", "toxin"),
    ("Water temperature", "water_temp_C", "env"),
    ("TP", "TP_mg_L", "env"),
    ("TN", "TN_mg_L", "env"),
    ("pH", "pH_wq", "env"),
    ("DO", "DO_mg_L_wq", "env"),
    ("BOD", "BOD_mg_L", "env"),
    ("COD", "COD_mg_L", "env"),
    ("EC", "EC_uS_cm", "env"),
    ("NH3-N", "NH3N_mg_L", "env"),
    ("NO3-N", "NO3N_mg_L", "env"),
    ("TOC", "TOC_mg_L", "env"),
    ("SS", "SS_mg_L", "env"),
    ("Air temperature", "air_temp_avg_C", "met"),
    ("Humidity", "humidity_avg_pct", "met"),
    ("Precipitation", "precip_sum_mm", "met"),
    ("Wind speed", "wind_speed_avg_ms", "met"),
    ("Sunshine", "sunshine_avg_hr", "met"),
    ("Pressure", "pressure_avg_hPa", "met"),
    ("Discharge", "discharge_mean_m3s", "hydro"),
    ("Storage rate", "storage_rate_pct", "hydro"),
]

KNOCKOFF_ENV = [
    ("bloom", "bloom", "toxin"),
    ("microcystin", "MC_ppb", "toxin"),
    ("chlorophyll", "Chl_a_mg_m3", "toxin"),
    ("water temp", "water_temp_C", "env"),
    ("TP", "TP_mg_L", "env"),
    ("TN", "TN_mg_L", "env"),
    ("pH", "pH_wq", "env"),
    ("DO", "DO_mg_L_wq", "env"),
    ("BOD", "BOD_mg_L", "env"),
    ("COD", "COD_mg_L", "env"),
    ("EC", "EC_uS_cm", "env"),
    ("NH3N", "NH3N_mg_L", "env"),
    ("NO3N", "NO3N_mg_L", "env"),
    ("TOC", "TOC_mg_L", "env"),
    ("SS", "SS_mg_L", "env"),
    ("air temp", "air_temp_avg_C", "met"),
    ("humidity", "humidity_avg_pct", "met"),
    ("precip", "precip_sum_mm", "met"),
    ("wind", "wind_speed_avg_ms", "met"),
    ("sunshine", "sunshine_avg_hr", "met"),
    ("pressure", "pressure_avg_hPa", "met"),
    ("discharge", "discharge_mean_m3s", "hydro"),
    ("storage", "storage_rate_pct", "hydro"),
]

GCOMP_FEATURES = PCS + [
    "water_temp_C",
    "MC_ppb",
    "Chl_a_mg_m3",
    "TOC_mg_L",
    "NO3N_mg_L",
    "discharge_mean_m3s",
    "bloom",
    "TN_mg_L",
    "TP_mg_L",
    "DO_mg_L_wq",
    "pH_wq",
]

GCOMP_DRIVERS = [
    ("Microbial community", "Comm_PC1", "quantile"),
    ("Water temperature", "water_temp_C", "quantile"),
    ("TOC", "TOC_mg_L", "quantile"),
    ("Nitrate-N", "NO3N_mg_L", "quantile"),
    ("Microcystin", "MC_ppb", "quantile"),
    ("Chlorophyll-a", "Chl_a_mg_m3", "quantile"),
    ("Discharge", "discharge_mean_m3s", "quantile"),
    ("Bloom", "bloom", "binary"),
]

def season_design(df: pd.DataFrame, controls: list[str]) -> np.ndarray:
    season = pd.get_dummies(df["season"], drop_first=True, dtype=float)
    if controls:
        return np.column_stack([df[controls].to_numpy(dtype=float), season.to_numpy()])
    return season.to_numpy()
