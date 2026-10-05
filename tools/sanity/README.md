# Sanity checks (cyl-3d ports)

Light scripts, run on a login node in minutes. Inputs cut from each flow's frozen DSMC
field with awk (one pass of ~100 MB, ~10 s per flow):

    P=~/orcd/scratch/cbgb-sparta-3d/cyl-ports/sccm$F/post/frozen.grid
    # on-axis cells (r < 0.6 mm): x n T
    awk 'NR>9 && ($3*$3+$4*$4)<3.6e-7 {print $2,$11,$15}' $P > $AXIS_DIR/axis_$F.txt
    # azimuthal average in 0.4 mm (x, r) bins, volume/density weighted: x r n u_x u_r T
    awk 'NR>9 { r=sqrt($3*$3+$4*$4); vol=($8-$5)*($9-$6)*($10-$7); ix=int(($2+0.0224)/0.0004); ir=int(r/0.0004); k=ix" "ir; ur=(r>0)?($3*$13+$4*$14)/r:0; W[k]+=vol; N[k]+=$11*vol; U[k]+=$11*$12*vol; R[k]+=$11*ur*vol; T[k]+=$11*$15*vol } END { for (k in W) if (N[k]>0) { split(k,a," "); printf "%.4f %.4f %.4e %.3f %.3f %.3f\n", -0.0224+(a[1]+0.5)*0.0004, (a[2]+0.5)*0.0004, N[k]/W[k], U[k]/N[k], R[k]/N[k], T[k]/N[k] } }' $P > $AXIS_DIR/field2d_$F.txt

Then (flows 10 15 20 30 40):

    AXIS_DIR=... python tools/sanity/range_check.py            # depth, collisions, dE/dx vs hand calc
    AXIS_DIR=... python tools/sanity/flight_mc.py              # start-point vs exact free-path sampling
    PYTHONPATH=tools AXIS_DIR=... python tools/sanity/more_checks.py   # checks 1-5

Results: results/sanity/.
