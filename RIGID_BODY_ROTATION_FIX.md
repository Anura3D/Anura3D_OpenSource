# Rigid Body Rotation Fix for Thin Rigid Bodies

## Problem Description

The thin rigid body contact code in `StructuralElements.for` was only considering the translational velocity of rigid bodies when calculating contact forces and tangential friction forces. This meant that for a rotating cylinder with material points inside, the rotational velocity component at the contact point was ignored, leading to incorrect contact behavior.

## Solution

The fix adds the rotational velocity component to the rigid body velocity calculation at each contact point. For a rigid body in 2D plane strain, the velocity at any point is given by:

```
v_point = v_translation + ω × r
```

Where:
- `v_translation` is the translational velocity of the rigid body centroid
- `ω` is the angular velocity (scalar in 2D)
- `r` is the position vector from the centroid to the contact point

In 2D, the cross product simplifies to:
```
v_point = [vx_trans, vy_trans] + ω * [-ry, rx]
```

## Changes Made

### 1. Added New Local Variables (Line ~957)
```fortran
real(REAL_TYPE) :: ..., RigidBodyVelAtPoint(2), RelativeVel(2), RadiusVec(2)
```

### 2. Calculate Rigid Body Velocity at Contact Point (After line ~988)
Added calculation that includes both translational and rotational components:
```fortran
! Calculate rigid body velocity at contact point including rotation
! v_point = v_translation + omega x r
! In 2D: v_point = v_translation + omega * (-ry, rx)
if (.not. ISAXISYMMETRIC) then
    RadiusVec(1) = Xp - DataStructureRigidBody(J)%Centroid(1)
    RadiusVec(2) = Yp - DataStructureRigidBody(J)%Centroid(2)
    RigidBodyVelAtPoint(1) = DataStructureRigidBody(J)%Velocity(1) + &
                             DataStructureRigidBody(J)%AngularVelocity * (-RadiusVec(2))
    RigidBodyVelAtPoint(2) = DataStructureRigidBody(J)%Velocity(2) + &
                             DataStructureRigidBody(J)%AngularVelocity * RadiusVec(1)
else
    ! For axisymmetric, only translational velocity (no rotation in this plane)
    RigidBodyVelAtPoint = DataStructureRigidBody(J)%Velocity
end if
```

### 3. Updated All Relative Velocity Calculations
Replaced all instances of:
```fortran
VelocityArray(I,:) - DataStructureRigidBody(J)%Velocity
```

With proper relative velocity calculation:
```fortran
RelativeVel = VelocityArray(I,:) - RigidBodyVelAtPoint
```

This affects:
- Crossing detection for time-stepping control
- Tangent vector calculation
- Normal force damping term
- Tangential force friction calculation
- All commented-out alternative friction models (for consistency)

### 4. Applied to Both Plane Strain and Axisymmetric Cases
The fix is applied in both the plane strain (`if (.not. ISAXISYMMETRIC)`) and axisymmetric (`if (ISAXISYMMETRIC)`) sections of the contact routine.

## Impact

This fix ensures that:
1. **Rotating cylinders** with material points inside will now properly account for the tangential velocity due to rotation
2. **Friction forces** will be calculated based on the correct relative velocity between MPs and the rotating surface
3. **Normal damping forces** will use the correct relative normal velocity
4. **Time-stepping control** will account for the full velocity field when checking for potential boundary crossings

## Test Recommendations

To validate the fix, test cases should include:
1. A rotating cylinder filled with material points
2. Material points sliding on a rotating disk
3. A spinning drum with internal MPs
4. Combined translation and rotation scenarios

Expected behavior: Material points should experience tangential friction forces that resist the relative sliding motion, including the component due to rigid body rotation.

## Related Files
- Modified: `src/StructuralElements.for` (subroutine `LevelSetContact`)
- Data structures: Uses existing `DataStructureRigidBody(J)%AngularVelocity` and `Centroid` members
