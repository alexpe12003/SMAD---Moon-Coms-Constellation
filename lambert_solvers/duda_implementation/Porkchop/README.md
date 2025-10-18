# 🚀 LEO-to-Moon Porkchop Transfer Analysis

**A comprehensive orbital mechanics analysis tool for calculating optimal transfer trajectories between Low Earth Orbit (LEO) and the Moon using Lambert's problem.**

## 📁 Project Structure

### 🎯 `final_results/`
**Your main results are here!**
- **`quick_extended_total.png`** - 🎨 **MAIN PORKCHOP PLOT** showing total ΔV requirements
- **`quick_extended_departure.png`** - Departure ΔV component analysis
- **`quick_extended_arrival.png`** - Arrival ΔV component analysis  
- **`orbit_data.txt`** - Generated orbital position data (LEO & Moon)

### ⚙️ `core_modules/`
**Essential code to run the analysis**
- **`porkchop.py`** - Core orbital mechanics calculator and orbit simulator
- **`quick_extended_porkchop.py`** - 🚀 **MAIN PORKCHOP GENERATOR** (run this!)

### 🔧 `development_files/`
**Development history and alternative implementations**
- Testing scripts and earlier versions
- Different resolution implementations
- Development utilities

### 📊 `legacy_plots/`
**Earlier version plots for comparison**
- Previous porkchop plots from development iterations
- Plots with limited Moon simulation coverage

### 📚 `documentation/`
**Project documentation and summaries**
- Technical documentation and development history

## 🚀 Quick Start

### 1. **View Results** (No coding required)
Open `final_results/quick_extended_total.png` to see the porkchop plot.

### 2. **Run New Analysis**
```bash
cd core_modules
python quick_extended_porkchop.py
```

### 3. **Understand the Orbits**
```bash
cd core_modules  
python porkchop.py
```

## 🎯 Key Results Summary

### **Optimal LEO-to-Moon Transfer Solution:**
- **🕐 Departure Time**: 147.5 minutes (2.46 hours into simulation)
- **⏱️ Time of Flight**: 96.0 hours (4.0 days)
- **🚀 Total ΔV**: 4.020 km/s
  - **🛫 Departure ΔV**: 3.146 km/s (78.3% of total)
  - **🛬 Arrival ΔV**: 0.874 km/s (21.7% of total)

### **Analysis Parameters:**
- **🌍 LEO Orbit**: 200 km altitude, 88.5 minute period
- **🌙 Moon Simulation**: Extended to 99.7 hours (proper coverage for all transfers)
- **📊 Grid Resolution**: 25 departure times × 20 time-of-flight values = 500 scenarios
- **✅ Valid Solutions**: 471/500 transfer cases found

### **Physical Interpretation:**
- **Front-loaded energy**: 78% of ΔV required at LEO departure
- **Efficient arrival**: Minimal insertion burn needed at Moon
- **4-day transfer**: Good balance of efficiency and flight time
- **Departure timing**: Optimal at ~1.67 LEO orbits into simulation

## 📊 How to Read the Porkchop Plot

The porkchop plot (`final_results/quick_extended_total.png`) shows:
- **X-axis**: Departure time from LEO (minutes)
- **Y-axis**: Time of flight (hours) 
- **Colors**: Total ΔV required (darker = lower ΔV = better)
- **White lines**: ΔV contours with values
- **Red star**: Optimal solution (minimum ΔV)
- **Cyan lines**: LEO orbit markers (every 88.5 minutes)

## 🛠️ Technical Details

### **Orbital Mechanics:**
- **LEO**: Circular orbit at 6,578 km radius (200 km altitude)
- **Moon**: Circular orbit at 384,400 km radius (average distance)
- **Coordinate System**: 2D coplanar orbits, Earth-centered
- **Transfer Method**: Lambert's problem solution for two-body dynamics

### **Simulation Coverage:**
- **LEO Data**: 1-minute intervals for one complete orbit (88.5 min)
- **Moon Data**: Extended simulation covering 99.7 hours (15.1% of lunar orbit)
- **Transfer Analysis**: Proper Moon position interpolation for all arrival times

### **Delta-V Calculation:**
- **Departure**: |V_transfer - V_circular_LEO|
- **Arrival**: |V_transfer - V_circular_Moon|
- **Total**: Departure ΔV + Arrival ΔV

## 🔬 Lambert Problem Implementation

Uses the Duda implementation of Lambert's problem to solve for transfer orbits:
- **Inputs**: Initial position (LEO), final position (Moon), time of flight
- **Outputs**: Required velocity vectors at departure and arrival
- **Method**: Universal variable formulation for robust convergence

## 📈 Mission Planning Insights

### **Optimal Strategy:**
- **Launch Window**: Depart at 147.5 minutes (specific LEO phase)
- **Transfer Type**: 4-day Hohmann-like transfer
- **Propulsion**: Most ΔV required at LEO departure (plan accordingly)
- **Arrival**: Minimal lunar insertion burn needed

### **Trade-offs:**
- **Shorter transfers** (1-2 days): Higher ΔV but faster arrival
- **Longer transfers** (5+ days): Slightly higher ΔV, more complex timing
- **Departure timing**: ±30 minutes can significantly affect ΔV

---

## 🤝 Usage Notes

**For mission planning**: Use the optimal solution parameters directly.

**For research**: Modify `core_modules/quick_extended_porkchop.py` to explore different:
- Departure windows
- Time of flight ranges  
- Orbital altitudes
- Grid resolutions

**For visualization**: All plots saved as high-resolution PNG files suitable for presentations.

---

*Generated from comprehensive Lambert problem analysis - October 2025*