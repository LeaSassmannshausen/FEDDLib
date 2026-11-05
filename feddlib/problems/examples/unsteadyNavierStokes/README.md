# Artery inflow (`BC Type = parabolic_artery`)

This 3D boundary condition prescribes an elliptical parabola:

`u(x,t) = MaxVelocity * ramp(t) * max(0, 1 - (s/a)^2 - (q/b)^2) * direction`.

Here `s` and `q` are the inlet-plane coordinates relative to its centre. The
inflow direction is normalized and must point into the fluid. The principal
axis is projected onto the inlet plane and normalized; its cross product with
the direction defines the second axis. Equal radii give a circular profile.
`MaxVelocity` is the peak speed, not the mean speed. A positive `Inflow Ramp
Time` applies a half-cosine ramp from zero to the requested speed; zero disables
the ramp. All lengths use the mesh coordinate units.

Set the following entries under `Parameter`. These values fit the inlet of
`uniform_try_1_optimized_flags_1_2_3_10.mesh`:

```xml
<Parameter name="BC Type" type="string" value="parabolic_artery"/>
<Parameter name="MaxVelocity" type="double" value="1.5"/>
<Parameter name="Inflow Flag" type="int" value="2"/>
<Parameter name="Wall Flag" type="int" value="1"/>
<Parameter name="Inflow Center X" type="double" value="0.03632148"/>
<Parameter name="Inflow Center Y" type="double" value="0.08433928"/>
<Parameter name="Inflow Center Z" type="double" value="0.0"/>
<Parameter name="Inflow Direction X" type="double" value="0.0"/>
<Parameter name="Inflow Direction Y" type="double" value="0.0"/>
<Parameter name="Inflow Direction Z" type="double" value="1.0"/>
<Parameter name="Inflow Axis X" type="double" value="0.77287016"/>
<Parameter name="Inflow Axis Y" type="double" value="-0.63456419"/>
<Parameter name="Inflow Axis Z" type="double" value="0.0"/>
<Parameter name="Inflow Radius 1" type="double" value="1.40549310"/>
<Parameter name="Inflow Radius 2" type="double" value="1.23376006"/>
<Parameter name="Inflow Ramp Time" type="double" value="0.1"/>
```

The supplied case has inlet 2 at `z=0`, no-slip wall 1, and outlet 3 at `z=8`.
No Dirichlet condition is imposed at the outlet. The elliptical fit approximates
the slightly irregular inlet rim; nodes marked as wall remain exactly zero.

The build-directory case is saved as `parametersProblem_parabolic_artery.xml`,
so future builds that copy the default `parametersProblem.xml` do not remove it.
Run from `build/feddlib/problems/examples/unsteadyNavierStokes`:

```sh
mpirun -np 4 ./problems_unsteadyNavierStokes.exe --problemfile=parametersProblem_parabolic_artery.xml
```
