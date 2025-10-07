# Lambda1 Calculation Failures: Analysis and Solution

## Problem Description

The lunar transfer trajectory calculations were failing for lambda1 values between 295° and 330° with a `ValueError: math domain error`.

## Root Cause Analysis

The error occurred in the `lunar_soi_calculations()` function at this line:
```python
epsilon2_rad = math.asin((vm/v2)*math.cos(lambda1_rad) - (v1_kms/v2)*math.cos(lambda1_rad + gamma1_rad - phi1_rad))
```

### Mathematical Issue

The `math.asin()` function requires its argument to be in the range [-1, 1]. For lambda1 values 295°-330°, the calculated argument exceeded 1.0:

| Lambda1 | asin argument | Status |
|---------|---------------|--------|
| 290°    | 0.979950     | ✓ Valid |
| 295°    | 1.023484     | ❌ Invalid |
| 300°    | 1.055818     | ❌ Invalid |
| 305°    | 1.076725     | ❌ Invalid |
| 310°    | 1.085833     | ❌ Invalid |
| 315°    | 1.083002     | ❌ Invalid |
| 320°    | 1.068212     | ❌ Invalid |
| 325°    | 1.041885     | ❌ Invalid |
| 330°    | 1.004472     | ❌ Invalid |
| 335°    | 0.996445     | ✓ Valid |

### Physical Interpretation

These lambda1 values correspond to **physically impossible trajectory geometries** where:
- The spacecraft's approach angle to the Moon's SOI
- Combined with the Moon's orbital motion
- Results in impossible relative velocity vectors

This suggests that for the given transfer parameters (R₀ = 1.05 DU, V₀ = 1.372 DU/TU), certain lunar encounter geometries cannot be achieved.

## Solution Implemented

### 1. Domain Validation
Added explicit validation before calling `math.asin()`:
```python
# Calculate the argument for arcsin and check for valid domain
asin_argument = (vm/v2)*math.cos(lambda1_rad) - (v1_kms/v2)*math.cos(lambda1_rad + gamma1_rad - phi1_rad)

# Check if the argument is within the valid domain [-1, 1] for arcsin
if not (-1.0 <= asin_argument <= 1.0):
    raise ValueError(f"Invalid flight path geometry for λ₁ = {lambda1_deg}°. "
                    f"arcsin argument = {asin_argument:.6f} is outside [-1, 1] range. "
                    f"This represents a physically impossible trajectory geometry.")
```

### 2. Improved Error Handling
Enhanced the parametric study to gracefully handle these cases:
```python
except ValueError as e:
    if "Invalid flight path geometry" in str(e):
        print(f"Skipping λ₁ = {lambda1}° (physically impossible geometry)")
    else:
        print(f"Warning: Calculation failed for lambda1 = {lambda1}°: {str(e)}")
    continue
```

### 3. Summary Statistics
Added reporting of skipped values:
```python
skipped_count = len(lambda1_values) - len(valid_lambda1_values)
if skipped_count > 0:
    print(f"Skipped {skipped_count} lambda1 values due to physically impossible geometries.")
```

## Results

The parametric study now runs successfully:
- **Total lambda1 values tested**: 73 (0° to 360° in 5° steps)
- **Valid calculations**: 65 data points
- **Skipped values**: 8 (λ₁ = 295°, 300°, 305°, 310°, 315°, 320°, 325°, 330°)
- **Success rate**: 89%

## Key Insights

1. **Physical Constraints**: Not all lambda1 values are achievable with given transfer parameters
2. **Orbital Mechanics**: The failure region (295°-330°) represents geometries where the Moon's orbital motion and spacecraft approach vector cannot be reconciled
3. **Mission Design**: These constraints should be considered when planning lunar transfer missions

## Validation

The solution maintains:
- ✅ Mathematical rigor (proper domain validation)
- ✅ Physical accuracy (excludes impossible geometries)  
- ✅ Computational robustness (graceful error handling)
- ✅ User feedback (clear reporting of limitations)

## Alternative Approaches Considered

1. **Clamp Method**: `asin_arg = max(-1.0, min(1.0, asin_arg))`
   - ❌ Would produce mathematically valid but physically meaningless results

2. **Skip Method**: Current implementation
   - ✅ Maintains physical integrity of calculations

3. **Parameter Adjustment**: Modify R₀ or V₀ to eliminate invalid regions
   - ⚠️ Would change mission requirements

The current solution (Skip Method) was chosen as it preserves the physical accuracy of the analysis while providing robust computational behavior.