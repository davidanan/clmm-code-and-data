# Data Dictionary

## General Conventions

- `Event_ID`, `Camera_ID`, `Site_ID`, `SeqID`, and `SeqSegID` are pseudonymous
  identifiers created for analysis and release.
- Event dates are: Event 1, 12 August 2023; Event 2, 20 July 2024; Event 3,
  26 July 2024; Event 4, 19 August 2024; and Event 5, 1 July 2025.
- Timestamps are local civil time in Dalian, China (China Standard Time, UTC+08:00).
- Observed and simulated inundation levels use the ordered integer scale 0--4.

## `data/paired_observed_simulated_inundation.csv`

This is the reduced model-input table used by both CLMM scripts in `code/`.
It contains 20,825 paired minute-level records from five events, 60 sites, and
164 continuous sequence segments.

| Column | Description | Unit or values |
| --- | --- | --- |
| `source_row_id` | Stable row identifier inherited from the source CLMM table. | Integer identifier |
| `Event_ID` | Rainfall-event identifier; Event 1 is the reference event in both CLMMs. | 1--5 |
| `Site_ID` | Monitoring-site identifier. | Integer identifier |
| `SeqID` | Observation-sequence identifier before temporal-gap splitting. | Character identifier |
| `SeqSegID` | Continuous sequence-segment identifier after temporal-gap splitting. | Character identifier |
| `Actual_Time` | Observation timestamp as retained in the model-input table. | Local time, `YYYY/M/D HH:MM` |
| `t` | Consecutive within-segment time index used by the CLMM workflow. | 1 min step |
| `Flood_Phase` | Observation-defined flood phase. | `rising`, `peak`, or `recession` |
| `Y_obs` | Expert-verified observed inundation level. | Ordered level 0--4 |
| `Y_sim` | Simulated inundation level after ordinal conversion. | Ordered level 0--4 |
| `Rain_Int` | Instantaneous rainfall intensity (I) assigned from the nearest gauge. | mm/min |
| `Rain_lag_30` | Antecedent 30 min rainfall accumulation (AR30). | mm |
| `DEM` | Local elevation (ELV). | m |
| `Slope` | Local slope (SLP). | % |
| `TPI` | Topographic position index (TPI). | m |
| `Distance_manhole` | Distance to the nearest manhole (DMH). | m |
| `Road_density` | Road length per unit area within a 200 m buffer (RD). | km/km² |
| `Diameter` | Diameter of the nearest drainage pipe (DIA). | m |
| `Imp` | Fraction of impervious area. | 0–1 |
| `sim_depth_m` | Simulated continuous inundation depth retained from `YsimDepth`. | m |

## `data/verified_video_inundation_levels.csv`

This is the release table of expert-verified observations. It contains 36,490
minute-level records from 222 event-specific camera records at 124 monitored sites.

| Column | Description | Unit or values |
| --- | --- | --- |
| `Event_ID` | Rainfall-event identifier. | 1--5 |
| `Actual_Time` | Observation timestamp formatted as ISO 8601 local time. | `YYYY-MM-DDTHH:MM:SS` |
| `Camera_ID` | Event-specific camera-record identifier (CID); the original field name is retained. | Integer identifier |
| `Site_ID` | Monitoring-site identifier shared across events where applicable. | Integer identifier |
| `Y_obs` | Expert-verified observed inundation level. | Ordered level 0--4 |

Retention and privacy rules:

- rows were retained if and only if `Y_obs` was not missing;
- no additional observation filtering was applied after the `Y_obs` screen;
- precise latitude and longitude fields were removed from the public release;
- raw surveillance images and videos are not included.

The derived data are licensed under [CC BY 4.0](../DATA_LICENSE.md).
