# Contribution: Aftershock Data Template

This directory contains a standardized CSV template for recording aftershock sequences, developed as part of the analysis of the 1979 UK earthquake.

## Purpose

The template facilitates systematic recording of aftershock parameters (time, location, depth, magnitude) in a format that can be easily ingested by seismological analysis tools such as EQcorrscan.

## File Format

`aftershock_data_template.csv` is a comma‑separated values file with the following columns:

| Column | Description | Format / Example |
|--------|-------------|------------------|
| `event_id` | Unique identifier for each aftershock | e.g., `AS001` |
| `date` | Date of the event | `YYYY‑MM‑DD` |
| `time` | Time of the event (UTC) | `HH:MM:SS` |
| `latitude` | Latitude in decimal degrees (north positive) | `55.1234` |
| `longitude` | Longitude in decimal degrees (east positive) | `‑2.5678` |
| `depth_km` | Depth in kilometers | `12.5` |
| `magnitude` | Magnitude value | `3.2` |
| `magnitude_type` | Magnitude scale used (`ML`, `Mw`, `Mb`, etc.) | `ML` |
| `location_quality` | Qualitative assessment of location uncertainty | `good` / `fair` / `poor` |
| `notes` | Free‑text remarks, source references, etc. | |

## Usage

1. Copy the template file and rename it for your specific aftershock sequence.
2. Fill in the rows with the observed aftershock data.
3. Use the CSV with Python scripts (e.g., pandas, ObsPy) to load and analyze the catalog.
4. The template is particularly suited for studying aftershock decay (Omori’s law), magnitude‑frequency distributions, and spatial clustering.

## Example

The file includes one example row (AS001) that illustrates the expected format. Replace the placeholder values with actual measurements.

## Related Work

This template was created for the analysis of the 1979 earthquake in northern England and southern Scotland. The same format can be applied to other historical or recent aftershock sequences.

## Contributing

If you have suggestions for additional columns or improvements, please open an issue or submit a pull request.