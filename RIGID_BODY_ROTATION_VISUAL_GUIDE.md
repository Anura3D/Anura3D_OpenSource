# Visual Explanation of Rigid Body Rotation Fix

## Before the Fix

```
    Rotating Cylinder
         ω
         ↻
        ___
       /   \
      |  C  |  ← Centroid (C) moving with velocity v_trans
       \___/
         ↑
         |
        MP   ← Material Point

Contact velocity calculation:
v_contact = v_trans  ❌ WRONG - ignores rotation!
```

The old code only used the translational velocity `v_trans` of the rigid body centroid, completely ignoring the rotational component at the contact point.

## After the Fix

```
    Rotating Cylinder
         ω
         ↻
        ___
       /   \      v_rot = ω × r
      |  C  |     ↗
       \___/     
         ↑       
         |      ↗ v_total = v_trans + v_rot
        MP   

Where:
- r = position vector from centroid C to contact point MP
- ω = angular velocity (scalar in 2D)
- v_rot = ω × r = ω * (-ry, rx)  in 2D

Contact velocity calculation:
v_contact = v_trans + ω × r  ✓ CORRECT
```

## Mathematical Details

### 2D Cross Product
In 2D plane strain, the velocity at contact point is:

```
v_point = [vx_trans] + ω * [-ry]
          [vy_trans]       [ rx]

where:
rx = xp - xc  (x-component of radius vector)
ry = yp - yc  (y-component of radius vector)
```

### Example: Rotating Cylinder

Consider a cylinder rotating counter-clockwise (ω > 0) with no translation:

```
              N (North)
                ↑
                |
    W ←─────────●─────────→ E
                |
                ↓
              S (South)

For a material point at North (ry > 0, rx = 0):
v = ω * (-ry, 0) = velocity pointing WEST ←

For a material point at East (rx > 0, ry = 0):
v = ω * (0, rx) = velocity pointing NORTH ↑

For a material point at South (ry < 0, rx = 0):
v = ω * (-ry, 0) = velocity pointing EAST →

For a material point at West (rx < 0, ry = 0):
v = ω * (0, rx) = velocity pointing SOUTH ↓
```

This creates the circular motion pattern expected for rotation!

## Impact on Forces

### Relative Velocity
```
v_relative = v_particle - v_rigid_body_at_contact
           = v_particle - (v_trans + ω × r)
```

### Tangential Force (Friction)
The friction force depends on the relative sliding velocity in the tangential direction:

```
v_tangent = v_relative - (v_relative · n) * n

Before fix: Ignores ω × r component
After fix:  Correctly accounts for rotational velocity
```

### Normal Force (Damping)
The damping component of normal force depends on relative velocity:

```
F_damping = -c * (v_relative · n) * n

Before fix: Incorrect damping due to missing rotational component
After fix:  Correct damping based on true relative velocity
```

## Example Scenarios

### Scenario 1: Cylinder with MPs Inside (No Friction)
- **Before**: MPs would not experience correct centripetal acceleration
- **After**: MPs correctly interact with rotating wall

### Scenario 2: Rotating Drum (With Friction)
- **Before**: MPs sliding on rotating surface would experience incorrect friction
- **After**: Friction force correctly opposes relative tangential motion

### Scenario 3: Combined Translation + Rotation
- **Before**: Only translation component considered
- **After**: Full velocity field (translation + rotation) properly accounted for

## Code Location
File: `src/StructuralElements.for`
Subroutine: `LevelSetContact()`
Lines: ~995-1015 (velocity calculation)
       ~1088-1120 (force calculation)
