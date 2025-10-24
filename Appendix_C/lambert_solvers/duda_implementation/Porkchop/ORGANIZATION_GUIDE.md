# Porkchop Analysis - Clean Organization

## 📁 Directory Structure

### Main Files (Keep These)

**Core Modules:**
- `porkchop.py` - Core orbital mechanics calculator
- `quick_extended_porkchop.py` - Final porkchop plot generator (BEST VERSION)

**Final Results:**
- `quick_extended_total.png` - **MAIN RESULT: Total ΔV porkchop plot**
- `quick_extended_departure.png` - Departure ΔV component
- `quick_extended_arrival.png` - Arrival ΔV component

**Documentation:**
- `ORGANIZATION_GUIDE.md` - This file (organization guide)
- `README.md` - Project overview

### Files to Archive/Remove

**Development Files (can be moved to archive):**
- `simple_porkchop_test.py`
- `medium_porkchop.py` 
- `comprehensive_porkchop.py`
- `extended_porkchop.py`
- `porkchop_with_save.py`
- `porkchop_plot.py`
- `data_access.py`
- `detailed_analysis.py`
- `organize_files.py`

**Legacy Plots (can be removed):**
- `porkchop_total_deltav.png`
- `porkchop_departure_deltav.png`
- `porkchop_arrival_deltav.png`
- `extended_porkchop_total_deltav.png`
- `extended_porkchop_departure_deltav.png`
- `extended_porkchop_arrival_deltav.png`

**Generated Data:**
- `orbit_data.txt` - Keep for reference
- `__pycache__/` - Can be removed (auto-generated)

## 🚀 Quick Start Guide

### To Run the Analysis:
```bash
python quick_extended_porkchop.py
```

### To View Results:
- Open `quick_extended_total.png` for the main porkchop plot
- Shows optimal transfer: **4.020 km/s ΔV at 147.5 min departure, 96h flight time**

### To Understand the Orbits:
```bash
python porkchop.py
```

## 🎯 Key Results Summary

**Optimal LEO-to-Moon Transfer Solution:**
- **Departure Time**: 147.5 minutes (2.46 hours, ~1.67 LEO orbits)
- **Time of Flight**: 96.0 hours (4.0 days)
- **Total ΔV**: 4.020 km/s
  - Departure ΔV: 3.146 km/s (78.3% of total)
  - Arrival ΔV: 0.874 km/s (21.7% of total)

**Analysis Parameters:**
- LEO altitude: 200 km (88.5 minute period)
- Moon simulation: Extended to 99.7 hours (proper coverage)
- Grid resolution: 25 × 20 = 500 transfer cases analyzed
- Valid solutions: 471/500 cases

## 🧹 Cleanup Recommendations

### Keep These Files:
1. `porkchop.py` - Core calculator
2. `quick_extended_porkchop.py` - Best porkchop generator  
3. `quick_extended_*.png` - Final results (3 files)
4. `orbit_data.txt` - Generated data
5. `README.md` - Documentation
6. This `ORGANIZATION_GUIDE.md` - Organization guide

### Can Remove/Archive:
- All other `.py` files (development versions)
- All other `.png` files (legacy plots)
- `__pycache__/` directory
- `SIMULATION_UPDATE_SUMMARY.md` (if not needed)

This leaves you with **8-9 essential files** instead of 20+ files!